MODULE historyFile

  USE nrtype
  USE dataTypes,         ONLY: STRFLX
  USE var_lookup,        ONLY: ixRFLX, nVarsRFLX
  USE var_lookup,        ONLY: ixHFLX, nVarsHFLX
  USE globalData,        ONLY: idxSUM,idxIRF,idxKWT, &
                               idxKW,idxMC,idxDW
  USE globalData,        ONLY: meta_rflx
  USE globalData,        ONLY: meta_hflx
  USE globalData,        ONLY: pid, nNodes
  USE globalData,        ONLY: masterproc
  USE globalData,        ONLY: mpicom_route
  USE globalData,        ONLY: pio_netcdf_format
  USE globalData,        ONLY: pio_typename
  USE globalData,        ONLY: pio_numiotasks
  USE globalData,        ONLY: pio_rearranger
  USE globalData,        ONLY: pio_root
  USE globalData,        ONLY: pio_stride
  USE pio_utils
  USE perf_mod,          ONLY: t_startf,t_stopf ! timing start/stop (GPTL library)

  implicit none
  integer(i4b),parameter     :: recordDim=-999       ! record dimension Indicator

  public:: histFile

  type :: histFile
    private
    character(300)        :: fname               ! netCDF name
    integer(i4b)          :: iTime=0             ! time step in output netCDF
    logical(lgt)          :: fileStatus=.false.  ! flag to indicate history output netcdf is open
    logical(lgt)          :: gageOutput=.false.  ! flag to indicate this is at-gage-only output (== output subset of reaches)
    type(iosystem_desc_t) :: pioSystem           ! PIO system (this does not have to be initialize for each history file?)
    type(file_desc_t)     :: pioFileDesc         ! PIO data identifying the file
    type(io_desc_t)       :: ioDescRchFlux       ! PIO domain decomposition data for reach flux [nRch]
    type(io_desc_t)       :: ioDescHruFlux       ! PIO domain decomposition data for hru runoff [nHRU]

    CONTAINS

      generic,    public :: createNC => createNC_rch, createNC_rch_hru
      generic,    public :: set_compdof => set_compdof_rch, set_compdof_rch_hru
      procedure,  public :: openNC
      procedure,  public :: fileOpen
      generic,    public :: write_loc => write_loc_rch, write_loc_rch_hru
      procedure,  public :: write_flux
      procedure,  public :: closeNC
      procedure,  public :: cleanup
      procedure, private :: cleanup_rch
      procedure, private :: cleanup_hru
      procedure, private :: write_flux_hru
      procedure, private :: write_flux_rch
      procedure, private :: createNC_rch
      procedure, private :: createNC_rch_hru
      procedure, private :: set_compdof_rch
      procedure, private :: set_compdof_rch_hru
      procedure, private :: write_loc_rch
      procedure, private :: write_loc_rch_hru

  end type histFile

  INTERFACE histFile
    module procedure constructor
  END INTERFACE histFile

  CONTAINS

    ! -----------------------------------------------------
    ! set pio local-global mapping index for reach
    ! -----------------------------------------------------
    FUNCTION constructor(fname, pioSys, gageOutput) RESULT(instHistFile)
      implicit none
      type(histFile)                              :: instHistFile
      character(*),                    intent(in) :: fname
      logical(lgt),          optional, intent(in) :: gageOutput
      type(iosystem_desc_t), optional, intent(in) :: pioSys

      instHistFile%fname = fname

      if (present(gageOutput)) then
        instHistFile%gageOutput = gageOutput
      end if

      ! pio initialization - pioSystem
      if (present(pioSys)) then
          instHistFile%pioSystem=pioSys
      else
        pio_numiotasks = nNodes/pio_stride
        call pio_sys_init(pid, mpicom_route,          & ! input: MPI related parameters
                          pio_stride, pio_numiotasks, & ! input: PIO related parameters
                          pio_rearranger, pio_root,   & ! input: PIO related parameters
                          instHistFile%pioSystem)       ! output: PIO system descriptors
      end if

    END FUNCTION constructor

    ! -----------------------------------------------------
    ! set pio local-global mapping index for reach
    ! -----------------------------------------------------
    SUBROUTINE set_compdof_rch(this, compdof_rch, nRch_in)
      implicit none
      ! Argument variables
      class(histFile),          intent(inout)  :: this
      integer(i4b),allocatable, intent(in)     :: compdof_rch(:)   ! Mapping of the storage order for the computational decomposition
      integer(i4b),             intent(in)     :: nRch_in          ! total number of reaches

      ! initialze iodescRchFlux
      call pio_decomp(this%pioSystem,    & ! input: pio system descriptor
                      ncd_double,        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [nRch_in],         & ! input: dimension length == global array size
                      compdof_rch,       & ! input: local->global mapping
                      this%ioDescRchFlux)

    END SUBROUTINE

    ! -----------------------------------------------------
    ! set pio local-global mapping index for reach and hru
    ! -----------------------------------------------------
    SUBROUTINE set_compdof_rch_hru(this, compdof_rch, compdof_hru, nRch_in, nHRU_in)
      implicit none
      ! Argument variables
      class(histFile),          intent(inout)  :: this
      integer(i4b),allocatable, intent(in)     :: compdof_rch(:)   ! Mapping of the storage order for the computational decomposition
      integer(i4b),allocatable, intent(in)     :: compdof_hru(:)   ! Mapping of the storage order for the computational decomposition
      integer(i4b),             intent(in)     :: nHRU_in          ! total number of HRUs
      integer(i4b),             intent(in)     :: nRch_in          ! total number of reaches

      ! initialze iodescRchFlux and ioDescHruFlux
      call pio_decomp(this%pioSystem,    & ! input: pio system descriptor
                      ncd_double,        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [nRch_in],         & ! input: dimension length == global array size
                      compdof_rch,       & ! input: local->global mapping
                      this%ioDescRchFlux)

      call pio_decomp(this%pioSystem,    & ! input: pio system descriptor
                      ncd_double,        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [nHRU_in],         & ! input: dimension length == global array size
                      compdof_hru,       & ! input: local->global mapping
                      this%ioDescHruFlux)
    END SUBROUTINE

    ! ---------------------------------------------------------------
    ! Create netCDF and define dimension and variables - only reach
    ! ---------------------------------------------------------------
    SUBROUTINE createNC_rch(this, nRch_in, ierr, message)

      USE var_lookup, ONLY: ixQdims
      USE globalData, ONLY: meta_qDims
      USE public_var, ONLY: calendar          ! calendar name
      USE public_var, ONLY: time_units        ! time units (seconds, hours, or days)

      implicit none
      ! Argument variables
      class(histFile),          intent(inout)  :: this
      integer(i4b),             intent(in)     :: nRch_in          ! total number of reaches
      integer(i4b),             intent(out)    :: ierr             ! error code
      character(*),             intent(out)    :: message          ! error message
      ! local variables
      character(strLen)                        :: cmessage         ! error message of downwind routine
      integer(i4b)                             :: ixDim,iVar       ! dimension, and variable index
      integer(i4b)                             :: dim_array(20)    ! dimension id array (max. 20 dimensions
      integer(i4b)                             :: nDims            ! number of dimension

      ierr=0; message='createNC_rch/'

      ! 1. Create new netCDF - initialize pioFileDesc under pioSystem
      call createFile(this%pioSystem, trim(this%fname), pio_typename, pio_netcdf_format, this%pioFileDesc, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! 2. Define dimension
      call def_dim(this%pioFileDesc, trim(meta_qDims(ixQdims%time)%dimName), recordDim, meta_qDims(ixQdims%time)%dimId)
      call def_dim(this%pioFileDesc, trim(meta_qDims(ixQdims%seg)%dimName),  nRch_in,   meta_qDims(ixQdims%seg)%dimId)
      call def_dim(this%pioFileDesc, trim(meta_qDims(ixQdims%ens)%dimName),  1,         meta_qDims(ixQdims%ens)%dimId)

      ! 3. Define variables
      ! --- time variable
      call def_var(this%pioFileDesc,                           &                                        ! pio file descriptor
                  trim(meta_qDims(ixQdims%time)%dimName),      &                                        ! variable name
                  ncd_float,                                   &                                        ! variable type
                  ierr, cmessage,                              &                                        ! error handle
                  pioDimId=[meta_qDims(ixQdims%time)%dimId],   &                                        ! dimension array
                  vdesc=trim(meta_qDims(ixQdims%time)%dimName), vunit=trim(time_units), vcal=calendar)  ! optional attributes
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! --- reach ID variables
      call def_var(this%pioFileDesc, 'reachID', ncd_int, ierr, cmessage, pioDimId=[meta_qDims(ixQdims%seg)%dimId], vdesc='reach ID', vunit='-')
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! ---- flux variables
      do iVar=1,nVarsRFLX
        if (.not.meta_rflx(iVar)%varFile) cycle

        ! define dimension ID array
        nDims = size(meta_rflx(iVar)%varDim)
        do ixDim = 1, nDims
          dim_array(ixDim) = meta_qDims(meta_rflx(iVar)%varDim(ixDim))%dimId
        end do

        ! define variable
        call def_var(this%pioFileDesc,            &                 ! pio file descriptor
                     meta_rflx(iVar)%varName,     &                 ! variable name
                     ncd_float,                   &                 ! dimension array and type
                     ierr, cmessage,              &                 ! error handling
                     pioDimId=dim_array(1:nDims), &                 ! dimension id
                     vdesc=meta_rflx(iVar)%varDesc, vunit=meta_rflx(iVar)%varUnit)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end do

      ! 4. End definitions
      call end_def(this%pioFileDesc, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    END SUBROUTINE createNC_rch

    ! ---------------------------------------------------------------
    ! Create netCDF and define dimension and variables - hru & reach
    ! ---------------------------------------------------------------
    SUBROUTINE createNC_rch_hru(this, nRch_in, nHRU_in, ierr, message)

      USE var_lookup, ONLY: ixQdims
      USE globalData, ONLY: meta_qDims
      USE public_var, ONLY: calendar          ! calendar name
      USE public_var, ONLY: time_units        ! time units (seconds, hours, or days)

      implicit none
      ! Argument variables
      class(histFile),       intent(inout)  :: this
      integer(i4b),             intent(in)     :: nRch_in
      integer(i4b),             intent(in)     :: nHru_in
      integer(i4b),             intent(out)    :: ierr             ! error code
      character(*),             intent(out)    :: message          ! error message
      ! local variables
      character(strLen)                        :: cmessage         ! error message of downwind routine
      integer(i4b)                             :: ixDim,iVar       ! dimension, and variable index
      integer(i4b)                             :: dim_array(20)    ! dimension id array (max. 20 dimensions
      integer(i4b)                             :: nDims            ! number of dimension

      ierr=0; message='createNC_rch_hru/'

      ! 1. Create new netCDF
      call createFile(this%pioSystem, trim(this%fname), pio_typename, pio_netcdf_format, this%pioFileDesc, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! 2. Define dimension
      call def_dim(this%pioFileDesc, trim(meta_qDims(ixQdims%time)%dimName), recordDim, meta_qDims(ixQdims%time)%dimId)
      call def_dim(this%pioFileDesc, trim(meta_qDims(ixQdims%seg)%dimName),  nRch_in,   meta_qDims(ixQdims%seg)%dimId)
      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        call def_dim(this%pioFileDesc, trim(meta_qDims(ixQdims%hru)%dimName), nHru_in, meta_qDims(ixQdims%hru)%dimId)
      end if
      call def_dim(this%pioFileDesc, trim(meta_qDims(ixQdims%ens)%dimName),  1,         meta_qDims(ixQdims%ens)%dimId)

      ! 3. Define variables
      ! --- time variable
      call def_var(this%pioFileDesc,                           &                                        ! pio file descriptor
                  trim(meta_qDims(ixQdims%time)%dimName),      &                                        ! variable name
                  ncd_float,                                   &                                        ! variable type
                  ierr, cmessage,                              &                                        ! error handle
                  pioDimId=[meta_qDims(ixQdims%time)%dimId],   &                                        ! dimension array
                  vdesc=trim(meta_qDims(ixQdims%time)%dimName), vunit=trim(time_units), vcal=calendar)  ! optional attributes
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! --- hru ID variables
      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        call def_var(this%pioFileDesc, 'basinID', ncd_int, ierr, cmessage, pioDimId=[meta_qDims(ixQdims%hru)%dimId], vdesc='basin ID', vunit='-')
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

      ! --- reach ID variables
      call def_var(this%pioFileDesc, 'reachID', ncd_int, ierr, cmessage, pioDimId=[meta_qDims(ixQdims%seg)%dimId], vdesc='reach ID', vunit='-')
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! ---- hru flux variables
      do iVar=1,nVarsHFLX
        if (.not.meta_hflx(iVar)%varFile) cycle

        ! define dimension ID array
        nDims = size(meta_hflx(iVar)%varDim)
        do ixDim = 1, nDims
          dim_array(ixDim) = meta_qDims(meta_hflx(iVar)%varDim(ixDim))%dimId
        end do

        ! define variable
        call def_var(this%pioFileDesc,            &                 ! pio file descriptor
                     meta_hflx(iVar)%varName,     &                 ! variable name
                     ncd_float,                   &                 ! dimension array and type
                     ierr, cmessage,              &                 ! error handling
                     pioDimId=dim_array(1:nDims), &                 ! dimension id
                     vdesc=meta_hflx(iVar)%varDesc, vunit=meta_hflx(iVar)%varUnit)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end do

      ! ---- rch flux variables
      do iVar=1,nVarsRFLX
        if (.not.meta_rflx(iVar)%varFile) cycle

        ! define dimension ID array
        nDims = size(meta_rflx(iVar)%varDim)
        do ixDim = 1, nDims
          dim_array(ixDim) = meta_qDims(meta_rflx(iVar)%varDim(ixDim))%dimId
        end do

        ! define variable
        call def_var(this%pioFileDesc,            &                 ! pio file descriptor
                     meta_rflx(iVar)%varName,     &                 ! variable name
                     ncd_float,                   &                 ! dimension array and type
                     ierr, cmessage,              &                 ! error handling
                     pioDimId=dim_array(1:nDims), &                 ! dimension id
                     vdesc=meta_rflx(iVar)%varDesc, vunit=meta_rflx(iVar)%varUnit)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end do

      ! 4. End definitions
      call end_def(this%pioFileDesc, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    END SUBROUTINE createNC_rch_hru

    ! ---------------------------------
    ! open netCDF
    ! ---------------------------------
    SUBROUTINE openNC(this, ierr, message)
      implicit none
      ! Argument variables
      class(histFile), intent(inout)  :: this
      integer(i4b),    intent(out)    :: ierr             ! error code
      character(*),    intent(out)    :: message          ! error message
      ! local variables
      character(strLen)               :: cmessage         ! error message of downwind routine

      ierr=0; message='openNC/'

      call openFile(this%pioSystem, this%pioFileDesc, trim(this%fname), pio_typename, ncd_write, this%FileStatus, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call inq_dim_len(this%pioFileDesc, 'time', this%iTime)

    END SUBROUTINE openNC

    ! ---------------------------------
    ! Check if file is open or not
    ! ---------------------------------
    logical(lgt) FUNCTION fileOpen(this)
      implicit none
      class(histFile), intent(inout) :: this

      fileOpen = this%fileStatus
    END FUNCTION

    ! ---------------------------------
    ! close netCDF
    ! ---------------------------------
    SUBROUTINE closeNC(this)
      implicit none
      class(histFile), intent(inout) :: this

      if (this%fileStatus) then
        call closeFile(this%pioFileDesc, this%fileStatus)
      endif
    END SUBROUTINE closeNC

    ! ---------------------------------
    ! clean up decomposition
    ! ---------------------------------
    SUBROUTINE cleanup(this)
      implicit none
      class(histFile),           intent(inout) :: this

      if (.not.this%gageOutput) then
        if (meta_hflx(ixHFLX%basRunoff)%varFile) then
          call this%cleanup_hru()
        endif
      end if

      call this%cleanup_rch()

    END SUBROUTINE cleanup

    ! ---------------------------------
    ! Cleanup rch resources
    ! ---------------------------------
    SUBROUTINE cleanup_hru(this)
      implicit none
      class(histFile), intent(inout) :: this

      call freeDecomp(this%pioFileDesc, this%ioDescHruFlux)

    END SUBROUTINE cleanup_hru

    ! ---------------------------------
    ! Cleanup hru resources
    ! ---------------------------------
    SUBROUTINE cleanup_rch(this)
      implicit none
      class(histFile), intent(inout) :: this

      call freeDecomp(this%pioFileDesc, this%ioDescRchFlux)

    END SUBROUTINE cleanup_rch

    ! ---------------------------------
    ! writing location (hru and/or rch)
    ! ---------------------------------
    SUBROUTINE write_loc_rch_hru(this, reach_id, basin_id, ierr, message)

      implicit none
      ! Argument variables
      class(histFile),           intent(inout) :: this
      integer(i4b), allocatable, intent(in)    :: reach_id(:)   !
      integer(i4b), allocatable, intent(in)    :: basin_id(:)   !
      integer(i4b),              intent(out)   :: ierr          ! error code
      character(*),              intent(out)   :: message       ! error message
      ! local variables
      integer(i4b)                             :: nElem         !
      character(len=strLen)                    :: cmessage      ! error message of downwind routine

      ierr=0; message='write_loc_rch_hru/'

      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        nElem = size(basin_id)
        call write_netcdf(this%pioFileDesc, 'basinID', basin_id, [1], [nElem], ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

      nElem = size(reach_id)
      call write_netcdf(this%pioFileDesc, 'reachID', reach_id, [1], [nElem], ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    END SUBROUTINE write_loc_rch_hru

    SUBROUTINE write_loc_rch(this, reach_id, ierr, message)

      implicit none
      ! Argument variables
      class(histFile),           intent(inout) :: this
      integer(i4b), allocatable, intent(in)    :: reach_id(:)   !
      integer(i4b),              intent(out)   :: ierr          ! error code
      character(*),              intent(out)   :: message       ! error message
      ! local variables
      integer(i4b)                             :: nElem         !
      character(len=strLen)                    :: cmessage      ! error message of downwind routine

      ierr=0; message='write_loc_rch/'

      nElem = size(reach_id)
      call write_netcdf(this%pioFileDesc, 'reachID', reach_id, [1], [nElem], ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    END SUBROUTINE

    ! ---------------------------------
    ! writing flux (hru and/or rch)
    ! ---------------------------------
    SUBROUTINE write_flux(this, time, index_write, ierr, message)

      implicit none
      ! Argument variables
      class(histFile),           intent(inout) :: this
      real(dp),                  intent(in)    :: time
      integer(i4b), allocatable, intent(in)    :: index_write(:)
      integer(i4b),              intent(out)   :: ierr             ! error code
      character(*),              intent(out)   :: message          ! error message
      ! local variables
      character(len=strLen)                    :: cmessage         ! error message of downwind routine

      ierr=0; message='write_flux/'

      this%iTime = this%iTime + 1 ! this is only line to increment time step index

      ! write time -- note time is just carried across from the input
      call write_netcdf(this%pioFileDesc, 'time', [time], [this%iTime], [1], ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      if (.not.this%gageOutput) then
        if (meta_hflx(ixHFLX%basRunoff)%varFile) then
          call this%write_flux_hru(ierr, cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        endif
      end if

      call this%write_flux_rch(index_write, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call sync_file(this%pioFileDesc, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    END SUBROUTINE write_flux


    SUBROUTINE write_flux_hru(this, ierr, message)
      ! currently only basin Runoff to be output
      implicit none
      ! Argument variables
      class(histFile),           intent(inout) :: this
      integer(i4b),              intent(out)   :: ierr             ! error code
      character(*),              intent(out)   :: message          ! error message
      ! local variables
      real(dp),    allocatable                 :: basinRunoff(:)
      character(len=strLen)                    :: cmessage         ! error message of downwind routine

      ierr=0; message='write_flux_hru/'

      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        call get_proc_flux(ierr, cmessage, basinRunoff=basinRunoff)

        ! write the basin runoff at HRU (unit: the same as runoff input)
        call write_pnetcdf_recdim(this%pioFileDesc, 'basRunoff',basinRunoff, this%ioDescHruFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

    END SUBROUTINE write_flux_hru


    SUBROUTINE write_flux_rch(this, index_write, ierr, message)

      USE globalData, ONLY: RCHFLX_trib

      implicit none
      ! Argument variables
      class(histFile),           intent(inout) :: this
      integer(i4b), allocatable, intent(in)    :: index_write(:)
      integer(i4b),              intent(out)   :: ierr             ! error code
      character(*),              intent(out)   :: message          ! error message
      ! local variables
      real(dp),    allocatable                 :: array_temp(:)
      integer(i4b)                             :: ix
      integer(i4b)                             :: nRch_write
      character(len=strLen)                    :: cmessage         ! error message of downwind routine

      ierr=0; message='write_flux_rch/'

      nRch_write = size(index_write)
      if (nRch_write==1 .and. index_write(1)==-9999) then
        nRch_write=0
      end if

      if (nRch_write>0) then
        allocate(array_temp(nRch_write), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif
      else
        allocate(array_temp(1), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif
      end if

      ! write instataneous local runoff in each stream segment (m3/s)
      if (meta_rflx(ixRFLX%instRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = RCHFLX_trib(1,index_write)%BASIN_QI
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'instRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      ! write routed local runoff in each stream segment (m3/s)
      if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = RCHFLX_trib(1,index_write)%BASIN_QR(1)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'dlayRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile) then
        if (nRch_write>0) then
          do ix=1,nRch_write
            array_temp(ix) = RCHFLX_trib(1,index_write(ix))%ROUTE(idxSUM)%REACH_Q
          end do
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'sumUpstreamRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      call t_startf ('output/write_flux/kwt')
      if (meta_rflx(ixRFLX%KWTroutedRunoff)%varFile) then
        if (nRch_write>0) then
          do ix=1,nRch_write
            array_temp(ix) = RCHFLX_trib(1,index_write(ix))%ROUTE(idxKWT)%REACH_Q
          end do
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'KWTroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif
      call t_stopf ('output/write_flux/kwt')

      call t_startf ('output/write_flux/irf')
      if (meta_rflx(ixRFLX%IRFroutedRunoff)%varFile) then
        if (nRch_write>0) then
          do ix=1,nRch_write
            array_temp(ix) = RCHFLX_trib(1,index_write(ix))%ROUTE(idxIRF)%REACH_Q
          end do
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'IRFroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif
      call t_stopf ('output/write_flux/irf')

      if (meta_rflx(ixRFLX%KWroutedRunoff)%varFile) then
        if (nRch_write>0) then
          do ix=1,nRch_write
            array_temp(ix) = RCHFLX_trib(1,index_write(ix))%ROUTE(idxKW)%REACH_Q
          end do
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'KWroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%MCroutedRunoff)%varFile) then
        if (nRch_write>0) then
          do ix=1,nRch_write
            array_temp(ix) = RCHFLX_trib(1,index_write(ix))%ROUTE(idxMC)%REACH_Q
          end do
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'MCroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%DWroutedRunoff)%varFile) then
        if (nRch_write>0) then
          do ix=1,nRch_write
            array_temp(ix) = RCHFLX_trib(1,index_write(ix))%ROUTE(idxDW)%REACH_Q
          end do
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'DWroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%volume)%varFile) then
        if (nRch_write>0) then
          do ix=1,nRch_write
            array_temp(ix) = RCHFLX_trib(1,index_write(ix))%ROUTE(idxIRF)%REACH_VOL(1)
          end do
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'volume', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

    END SUBROUTINE write_flux_rch

    ! ---------------------------------
    ! non-type-bound subroutine - organize flux data array or structure into per processor
    ! ---------------------------------
    SUBROUTINE get_proc_flux(ierr, message, basinRunoff)

      USE globalData, ONLY: nHRU_mainstem       ! number of mainstem HRUs
      USE globalData, ONLY: basinRunoff_main    ! mainstem only HRU runoff
      USE globalData, ONLY: basinRunoff_trib    ! tributary only HRU runoff
      USE globalData, ONLY: hru_per_proc        ! number of hrus assigned to each proc (size = num of procs+1)

      implicit none
      ! Argument variables
      real(dp),    allocatable, optional, intent(out) :: basinRunoff(:)
      integer(i4b),                       intent(out) :: ierr             ! error code
      character(*),                       intent(out) :: message          ! error message
      ! local variables
      integer(i4b)                          :: nHRU_local
      character(strLen)                     :: cmessage         ! error message of downwind routine

      if (present(basinRunoff)) then
        if (masterproc) then
          nHRU_local = nHRU_mainstem + hru_per_proc(0)
          allocate(basinRunoff(nHRU_local), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [basinRunoff]'; return; endif
          if (nHRU_mainstem>0) basinRunoff(1:nHRU_mainstem) = basinRunoff_main(1:nHRU_mainstem)
          if (hru_per_proc(0)>0) basinRunoff(nHRU_mainstem+1:nHRU_local) = basinRunoff_trib(:)
        else
          nHRU_local = hru_per_proc(pid)
          allocate(basinRunoff(nHRU_local), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [basinRunoff]'; return; endif
          basinRunoff = basinRunoff_trib
        endif
      end if

    END SUBROUTINE get_proc_flux

END MODULE historyFile
