MODULE historyFile

  USE nrtype
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
  USE globalData,        ONLY: version
  USE globalData,        ONLY: gitBranch
  USE globalData,        ONLY: gitHash
  USE histVars_data,     ONLY: histVars
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
    type(iosystem_desc_t) :: pioSys              ! PIO system (this does not have to be initialize for each history file?)
    type(file_desc_t)     :: pioFileDesc         ! PIO data identifying the file
    type(io_desc_t)       :: ioDescRchFlux       ! PIO domain decomposition data for reach flux [nRch]
    type(io_desc_t)       :: ioDescHruFlux       ! PIO domain decomposition data for hru runoff [nHRU]

    CONTAINS

      generic,    public :: set_compdof => set_compdof_rch, set_compdof_rch_hru
      procedure,  public :: createNC
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
      instHistFile%fileStatus = .false.

      if (present(gageOutput)) then
        instHistFile%gageOutput = gageOutput
      end if

      ! pio initialization - pioSystem
      if (present(pioSys)) then
          instHistFile%pioSys=pioSys
      else
        pio_numiotasks = nNodes/pio_stride
        call pio_sys_init(pid, mpicom_route,          & ! input: MPI related parameters
                          pio_stride, pio_numiotasks, & ! input: PIO related parameters
                          pio_rearranger, pio_root,   & ! input: PIO related parameters
                          instHistFile%pioSys)       ! output: PIO system descriptors
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
      call pio_decomp(this%pioSys,       & ! input: pio system descriptor
                      ncd_float,         & ! input: data type (pio_int, pio_real, pio_double, pio_char)
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
      call pio_decomp(this%pioSys,       & ! input: pio system descriptor
                      ncd_float,         & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [nRch_in],         & ! input: dimension length == global array size
                      compdof_rch,       & ! input: local->global mapping
                      this%ioDescRchFlux)

      call pio_decomp(this%pioSys,       & ! input: pio system descriptor
                      ncd_float,         & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [nHRU_in],         & ! input: dimension length == global array size
                      compdof_hru,       & ! input: local->global mapping
                      this%ioDescHruFlux)
    END SUBROUTINE

    ! ---------------------------------------------------------------
    ! Create history netCDF and define dimension and variables
    ! ---------------------------------------------------------------
    SUBROUTINE createNC(this, ierr, message, nRch_in, nHRU_in)

      USE var_lookup, ONLY: ixQdims
      USE globalData, ONLY: meta_qDims
      USE public_var, ONLY: calendar          ! calendar name
      USE public_var, ONLY: time_units        ! time units (seconds, hours, or days)

      implicit none
      ! Argument variables
      class(histFile),          intent(inout)  :: this
      integer(i4b),             intent(out)    :: ierr             ! error code
      character(*),             intent(out)    :: message          ! error message
      integer(i4b), optional,   intent(in)     :: nRch_in
      integer(i4b), optional,   intent(in)     :: nHRU_in
      ! local variables
      character(strLen)                        :: cmessage         ! error message of downwind routine
      integer(i4b)                             :: nRch_local
      integer(i4b)                             :: nHRU_local
      integer(i4b)                             :: ixDim,iVar       ! dimension, and variable index
      integer(i4b)                             :: dim_array(20)    ! dimension id array (max. 20 dimensions
      integer(i4b)                             :: nDims            ! number of dimension

      ierr=0; message='createNC/'

      nRch_local=0
      if (present(nRch_in)) then
        nRch_local=nRch_in
      end if

      nHRU_local=0
      if (present(nHRU_in)) then
        nHRU_local=nHRU_in
      end if

      ! gauge only history file does not include hru fluxes
      if ( this%gageOutput) then
        nHRU_local=0
      end if

      ! 1. Create new netCDF
      call createFile(this%pioSys, this%fname, pio_typename, pio_netcdf_format, this%pioFileDesc, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! 2. Define dimension
      call def_dim(this%pioFileDesc, meta_qDims(ixQdims%time)%dimName,   recordDim, meta_qDims(ixQdims%time)%dimId)
      call def_dim(this%pioFileDesc, meta_qDims(ixQdims%tbound)%dimName, 2,         meta_qDims(ixQdims%tbound)%dimId)
      call def_dim(this%pioFileDesc, meta_qDims(ixQdims%ens)%dimName,    1,         meta_qDims(ixQdims%ens)%dimId)
      if (nRch_local>0) then
        call def_dim(this%pioFileDesc, meta_qDims(ixQdims%seg)%dimName, nRch_in, meta_qDims(ixQdims%seg)%dimId)
      end if
      if (nHRU_local>0) then ! gauge only history file does not include hru fluxes
        call def_dim(this%pioFileDesc, meta_qDims(ixQdims%hru)%dimName, nHru_in, meta_qDims(ixQdims%hru)%dimId)
      end if

      ! 3. Define variables
      ! --- time variable
      call def_var(this%pioFileDesc,                           &                           ! pio file descriptor
                  meta_qDims(ixQdims%time)%dimName,            &                           ! variable name
                  ncd_double,                                  &                           ! variable type
                  ierr, cmessage,                              &                           ! error handle
                  pioDimId=[meta_qDims(ixQdims%time)%dimId],   &                           ! dimension array
                  vdesc=meta_qDims(ixQdims%time)%dimName, vunit=time_units, vcal=calendar) ! optional attributes
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! --- time bound variable
      call def_var(this%pioFileDesc,                           &                                  ! pio file descriptor
                  'time_bounds',                               &                                  ! variable name
                  ncd_double,                                  &                                  ! variable type
                  ierr, cmessage,                              &                                  ! error handle
                  pioDimId=[meta_qDims(ixQdims%tbound)%dimId, meta_qDims(ixQdims%time)%dimId], &  ! dimension array
                  vdesc='time interval endpoints', vunit=time_units, vcal=calendar)               ! optional attributes
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! ---- basinID and hru flux variables
      if (nHRU_local>0) then
        call def_var(this%pioFileDesc, 'basinID', ncd_int, ierr, cmessage, &
                     pioDimId=[meta_qDims(ixQdims%hru)%dimId], vdesc='basin ID', vunit='-')
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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
                       meta_hflx(iVar)%varType,     &                 ! variable type
                       ierr, cmessage,              &                 ! error handling
                       pioDimId=dim_array(1:nDims), &                 ! dimension id
                       vdesc=meta_hflx(iVar)%varDesc, vunit=meta_hflx(iVar)%varUnit)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        end do
      end if

      ! --- reach ID and rch flux variables
      if (nRch_local>0) then
        call def_var(this%pioFileDesc, 'reachID', ncd_int, ierr, cmessage, &
                     pioDimId=[meta_qDims(ixQdims%seg)%dimId], vdesc='reach ID', vunit='-')
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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
                       meta_rflx(iVar)%varType,     &                 ! variable type
                       ierr, cmessage,              &                 ! error handling
                       pioDimId=dim_array(1:nDims), &                 ! dimension id
                       vdesc=meta_rflx(iVar)%varDesc, vunit=meta_rflx(iVar)%varUnit)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        end do
      end if

      call put_attr(this%pioFileDesc, ncd_global, 'mizuRoute-version', version ,ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call put_attr(this%pioFileDesc, ncd_global, 'gitBranch', gitBranch ,ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call put_attr(this%pioFileDesc, ncd_global, 'gitHash', gitHash ,ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! 4. End definitions
      call end_def(this%pioFileDesc, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    END SUBROUTINE createNC

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

      call openFile(this%pioSys, this%pioFileDesc, trim(this%fname), pio_typename, ncd_write, this%FileStatus, ierr, cmessage)
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

      if (this%fileOpen() ) then
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

      !call freeDecomp(this%pioFileDesc, this%ioDescHruFlux)

    END SUBROUTINE cleanup_hru

    ! ---------------------------------
    ! Cleanup hru resources
    ! ---------------------------------
    SUBROUTINE cleanup_rch(this)
      implicit none
      class(histFile), intent(inout) :: this

      !call freeDecomp(this%pioFileDesc, this%ioDescRchFlux)

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
    SUBROUTINE write_flux(this, hVars_local, index_write, ierr, message)

      implicit none
      ! Argument variables
      class(histFile),           intent(inout) :: this
      type(histVars),            intent(in)    :: hVars_local
      integer(i4b), allocatable, intent(in)    :: index_write(:)
      integer(i4b),              intent(out)   :: ierr             ! error code
      character(*),              intent(out)   :: message          ! error message
      ! local variables
      character(len=strLen)                    :: cmessage         ! error message of downwind routine

      ierr=0; message='write_flux/'

      this%iTime = this%iTime + 1 ! this is only line to increment time step index

      ! write time -- note time is just carried across from the input
      call write_netcdf(this%pioFileDesc, 'time', [hVars_local%timeVar(1)], [this%iTime], [1], ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call write_netcdf(this%pioFileDesc, 'time_bounds', hVars_local%timeVar, [1,this%iTime], [2,1], ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      if (.not.this%gageOutput) then
        if (meta_hflx(ixHFLX%basRunoff)%varFile) then
          call this%write_flux_hru(hVars_local, ierr, cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        endif
      end if

      call this%write_flux_rch(hVars_local, index_write, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call sync_file(this%pioFileDesc, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    END SUBROUTINE write_flux


    SUBROUTINE write_flux_hru(this, hVars_local, ierr, message)

      ! Write HRU flux variables
      ! currently only basin hru Runoff [m/s]

      implicit none
      ! Argument variables
      class(histFile),           intent(inout) :: this
      type(histVars),            intent(in)    :: hVars_local
      integer(i4b),              intent(out)   :: ierr             ! error code
      character(*),              intent(out)   :: message          ! error message
      ! local variables
      real(sp),    allocatable                 :: array_temp(:)
      character(len=strLen)                    :: cmessage         ! error message of downwind routine

      ierr=0; message='write_flux_hru/'

      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        allocate(array_temp(hVars_local%nHru), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif
        array_temp = hVars_local%basRunoff
        call write_pnetcdf_recdim(this%pioFileDesc, 'basRunoff', array_temp, this%ioDescHruFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

    END SUBROUTINE write_flux_hru


    SUBROUTINE write_flux_rch(this, hVars_local, index_write, ierr, message)

      implicit none
      ! Argument variables
      class(histFile),           intent(inout) :: this
      type(histVars),            intent(in)    :: hVars_local
      integer(i4b), allocatable, intent(in)    :: index_write(:)
      integer(i4b),              intent(out)   :: ierr             ! error code
      character(*),              intent(out)   :: message          ! error message
      ! local variables
      real(sp),    allocatable                 :: array_temp(:)
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
          array_temp(1:nRch_write) = hVars_local%instRunoff(index_write)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'instRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      ! write routed local runoff in each stream segment (m3/s)
      if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = hVars_local%dlayRunoff(index_write)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'dlayRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxSUM)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'sumUpstreamRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%KWTroutedRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxKWT)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'KWTroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      call t_startf ('output/write_flux/irf')
      if (meta_rflx(ixRFLX%IRFroutedRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxIRF)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'IRFroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif
      call t_stopf ('output/write_flux/irf')

      if (meta_rflx(ixRFLX%KWroutedRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxKW)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'KWroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%MCroutedRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxMC)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'MCroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%DWroutedRunoff)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxDW)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'DWroutedRunoff', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      if (meta_rflx(ixRFLX%IRFvolume)%varFile) then
        if (nRch_write>0) then
          array_temp(1:nRch_write) = hVars_local%volume(index_write, idxIRF)
        end if
        call write_pnetcdf_recdim(this%pioFileDesc, 'IRFvolume', array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

    END SUBROUTINE write_flux_rch

END MODULE historyFile
