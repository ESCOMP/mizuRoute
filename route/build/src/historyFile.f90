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
  USE globalData,        ONLY: pioSystem
  USE globalData,        ONLY: pio_netcdf_format
  USE globalData,        ONLY: pio_typename
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
    character(FileStrLen)          :: fname               ! netCDF name
    integer(i4b)                   :: iTime=0             ! time step in output netCDF
    logical(lgt)                   :: fileStatus=.false.  ! flag to indicate history output netcdf is open
    logical(lgt)                   :: gageOutput=.false.  ! flag to indicate this is at-gage-only output (== output subset of reaches)
    type(file_desc_t)              :: pioFileDesc         ! PIO data identifying the file
    type(io_desc_t),      pointer  :: ioDescRchFlux       ! PIO domain decomposition data for reach flux [nRch]
    type(io_desc_t),      pointer  :: ioDescHruFlux       ! PIO domain decomposition data for hru runoff [nHRU]

    CONTAINS

      generic,    public :: set_compdof => set_compdof_rch, set_compdof_rch_hru
      procedure,  public :: createNC
      procedure,  public :: openNC
      procedure,  public :: fileOpen
      generic,    public :: write_loc => write_loc_rch, write_loc_rch_hru
      procedure,  public :: write_flux
      procedure,  public :: closeNC
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
    FUNCTION constructor(fname, gageOutput) RESULT(instHistFile)
      implicit none
      type(histFile)                              :: instHistFile
      character(*),                    intent(in) :: fname
      logical(lgt),          optional, intent(in) :: gageOutput

      instHistFile%fname = fname
      instHistFile%fileStatus = .false.

      if (present(gageOutput)) then
        instHistFile%gageOutput = gageOutput
      end if

    END FUNCTION constructor

    ! -----------------------------------------------------
    ! set pio local-global mapping index for reach
    ! -----------------------------------------------------
    SUBROUTINE set_compdof_rch(this, ioDescRchFlux)
      implicit none
      ! Argument variables
      class(histFile),          intent(inout)  :: this
      type(io_desc_t), target,  intent(in)     :: ioDescRchFlux

      this%ioDescRchFlux => ioDescRchFlux

    END SUBROUTINE

    ! -----------------------------------------------------
    ! set pio local-global mapping index for reach and hru
    ! -----------------------------------------------------
    SUBROUTINE set_compdof_rch_hru(this, ioDescRchFlux, ioDescHruFlux)
      implicit none
      ! Argument variables
      class(histFile),          intent(inout)  :: this
      type(io_desc_t), target,  intent(in)     :: ioDescRchFlux
      type(io_desc_t), target,  intent(in)     :: ioDescHruFlux

      this%ioDescRchFlux => ioDescRchFlux
      this%ioDescHruFlux => ioDescHruFlux

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
      call createFile(pioSystem, this%fname, pio_typename, pio_netcdf_format, this%pioFileDesc, ierr, cmessage)
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

      call openFile(pioSystem, this%pioFileDesc, trim(this%fname), pio_typename, ncd_write, this%FileStatus, ierr, cmessage)
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
      integer(i4b)                             :: iVar             ! variable index
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

      do iVar=1,nVarsRFLX

        if (.not.meta_rflx(iVar)%varFile) cycle

        if (nRch_write>0) then
          select case(iVar)
            case(ixRFLX%instRunoff)
              array_temp(1:nRch_write) = hVars_local%instRunoff(index_write)
            case(ixRFLX%dlayRunoff)
              array_temp(1:nRch_write) = hVars_local%dlayRunoff(index_write)
            case(ixRFLX%sumUpstreamRunoff)
              array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxSUM)
            case(ixRFLX%KWTroutedRunoff)
              array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxKWT)
            case(ixRFLX%IRFroutedRunoff)
              array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxIRF)
            case(ixRFLX%KWroutedRunoff)
              array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxKW)
            case(ixRFLX%MCroutedRunoff)
              array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxMC)
            case(ixRFLX%DWroutedRunoff)
              array_temp(1:nRch_write) = hVars_local%discharge(index_write, idxDW)
            case(ixRFLX%KWTvolume)
              array_temp(1:nRch_write) = hVars_local%volume(index_write, idxKWT)
            case(ixRFLX%IRFvolume)
              array_temp(1:nRch_write) = hVars_local%volume(index_write, idxIRF)
            case(ixRFLX%KWvolume)
              array_temp(1:nRch_write) = hVars_local%volume(index_write, idxKW)
            case(ixRFLX%MCvolume)
              array_temp(1:nRch_write) = hVars_local%volume(index_write, idxMC)
            case(ixRFLX%DWvolume)
              array_temp(1:nRch_write) = hVars_local%volume(index_write, idxDW)
            case default
              write(message,'(2A,1X,G0,1X,1A)') trim(message), 'Invalid RFLX variable index:',iVar,'. Check ixRFLX in var_lookup.f90'
              ierr=81; return
          end select
        end if

        call write_pnetcdf_recdim(this%pioFileDesc, meta_rflx(iVar)%varName, array_temp, this%ioDescRchFlux, this%iTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      end do

    END SUBROUTINE write_flux_rch

END MODULE historyFile
