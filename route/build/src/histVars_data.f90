MODULE histVars_data

 ! history file variable data

  USE nrtype
  USE dataTypes,         ONLY: STRFLX
  USE var_lookup,        ONLY: ixRFLX, nVarsRFLX
  USE var_lookup,        ONLY: ixHFLX, nVarsHFLX
  USE public_var,        ONLY: accumRunoff                ! accum runoff ID = 0
  USE public_var,        ONLY: impulseResponseFunc        ! IRF routing ID = 1
  USE public_var,        ONLY: kinematicWaveTracking      ! KWT routing ID = 2
  USE public_var,        ONLY: kinematicWave              ! KW routing ID = 3
  USE public_var,        ONLY: muskingumCunge             ! MC routing ID = 4
  USE public_var,        ONLY: diffusiveWave              ! DW routing ID = 5
  USE globalData,        ONLY: nRoutes
  USE globalData,        ONLY: routeMethods
  USE globalData,        ONLY: meta_rflx, meta_hflx
  USE globalData,        ONLY: idxSUM,idxIRF,idxKWT, &
                               idxKW,idxMC,idxDW
  USE globalData,        ONLY: nRch, nRch_mainstem, nRch_trib
  USE globalData,        ONLY: nHRU, nHRU_mainstem, nHRU_trib
  USE globalData,        ONLY: nTribOutlet
  USE globalData,        ONLY: hru_per_proc        ! number of hrus assigned to each proc (size = num of procs
  USE globalData,        ONLY: rch_per_proc        ! number of reaches assigned to each proc (size = num of procs)
  ! pio stuff
  USE globalData,        ONLY: pid, nNodes
  USE globalData,        ONLY: masterproc
  USE globalData,        ONLY: mpicom_route
  USE globalData,        ONLY: pio_typename
  USE globalData,        ONLY: pioSystem
  USE nr_utils,          ONLY: arth
  USE ncio_utils,        ONLY: get_nc
  USE pio_utils,         ONLY: file_desc_t, io_desc_t
  USE pio_utils,         ONLY: ncd_float, ncd_double, ncd_nowrite
  USE pio_utils,         ONLY: pio_decomp
  USE pio_utils,         ONLY: read_dist_array
  USE pio_utils,         ONLY: openFile
  USE pio_utils,         ONLY: closeFile

  implicit none

  public:: histVars

  type :: histVars
    ! --------
    ! output varialbes: basRunoff, instRunoff, dlayRunoff, discharge, elev, volume
    ! if new variables need to be written, need to add here AND each procedure
    ! --------
    integer(i4b)             :: nt                 ! number of time steps over variables are aggregated
    integer(i4b)             :: nHru               ! number of
    integer(i4b)             :: nRch               ! number of
    real(dp)                 :: timeVar            ! time variable [unit] at nt = 1
    real(dp), allocatable    :: basRunoff(:)       ! catchment average runoff [m/s] [nHru]
    real(dp), allocatable    :: instRunoff(:)      ! instantaneous lateral inflow into each river/lake [m3/s]  [nRch]
    real(dp), allocatable    :: dlayRunoff(:)      ! lateral inflow into river or lake [m3/s] for each reach [nRch]
    real(dp), allocatable    :: discharge(:,:)     ! river/lake discharge [m3/s] for each reach/lake and routing method [nRch,nMethod]
    real(dp), allocatable    :: elev(:,:)          ! river/lake surface water elevation [m] for each reach/lake and routing method [nRch,nMethod]
    real(dp), allocatable    :: volume(:,:)        ! river/lake volume [m3] for each reach/lake and routing method [nRch,nMethod]

    CONTAINS

      procedure,  public :: read_restart  ! initialize from restart file
      procedure,  public :: aggregate  ! Accumulating output variables
      procedure,  public :: finalize   ! Compute aggregated values (currently only mean) to be written in netCDF
      procedure,  public :: refresh    ! Reset output arrays to zero
      procedure,  public :: clean      ! Deallocate all the output array memories

  end type histVars

  INTERFACE histVars
    module procedure constructor
  END INTERFACE histVars

  CONTAINS

    ! -----------------------------------------------------
    ! constructor - instantiate history output data structure
    ! -----------------------------------------------------
    FUNCTION constructor(nHru_local, nRch_local, ierr, message) RESULT(instHistVar)

      implicit none
      ! argument variables
      integer(i4b),   intent(in)      :: nHru_local          ! number of hru in each proc
      integer(i4b),   intent(in)      :: nRch_local          ! number of Rch in each proc
      integer(i4b),   intent(out)     :: ierr                ! error code
      character(*),   intent(out)     :: message             ! error message
      type(histVars)                  :: instHistVar
      ! local variables
      character(strLen)               :: cmessage            ! error message from subroutine

      ierr=0; message='initHistVars/'

      instHistVar%nt   = 0
      instHistVar%nHru = nHRU_local
      instHistVar%nRch = nRch_local

      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        allocate(instHistVar%basRunoff(nHRU_local), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%basRunoff]'; return; endif
        instHistVar%basRunoff(1:nHRU_local) = 0._dp
      end if

      if (meta_rflx(ixRFLX%instRunoff)%varFile) then
        allocate(instHistVar%instRunoff(nRch_local), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%instRunoff]'; return; endif
        instHistVar%instRunoff(1:nRch_local) = 0._dp
      end if

      if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
        allocate(instHistVar%dlayRunoff(nRch_local), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%instRunoff]'; return; endif
        instHistVar%dlayRunoff(1:nRch_local) = 0._dp
      end if

      if (nRoutes>0) then ! this should be number of methods that ouput
        allocate(instHistVar%discharge(nRch_local, nRoutes), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%discharge]'; return; endif
        instHistVar%discharge(1:nRch_local, 1:nRoutes) = 0._dp

        allocate(instHistVar%volume(nRch_local, nRoutes), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%volume]'; return; endif
        instHistVar%volume(1:nRch_local, 1:nRoutes) = 0._dp
      end if

    END FUNCTION constructor

    ! ---------------------------------------------------------------
    ! accumulate data
    ! ---------------------------------------------------------------
    SUBROUTINE aggregate(this,             & ! inout:
                         timeVar_local,    & ! input: time variables current
                         basRunoff_local,  & ! input:
                         RCHFLX_local,     & ! input:
                         ierr, message)      ! output: error handling

      implicit none
      ! argument variables
      class(histVars),   intent(inout)   :: this
      real(dp),          intent(in)      :: timeVar_local
      real(dp),          intent(in)      :: basRunoff_local(:)
      type(STRFLX),      intent(in)      :: RCHFLX_local(:,:)
      integer(i4b),      intent(out)     :: ierr                ! error code
      character(*),      intent(out)     :: message             ! error message
      ! local variables
      integer(i4b)                       :: nHRU_input
      integer(i4b)                       :: nRch_input
      integer(i4b)                       :: ix, iRoute          ! loop indices
      integer(i4b)                       :: idxMethod           ! temporal method index

      ierr=0; message='aggregate/'

      ! -- increment number of sample
      this%nt = this%nt + 1

      if (this%nt == 1) then
        this%timeVar = timeVar_local
      end if

      ! -- array size checks - input data vs history output buffer
      ! hru and reach size in input data
      nHru_input = size(basRunoff_local)
      nRch_input = size(RCHFLX_local(1,:))
      if (nHru_input/=this%nHru) then
        write(message,'(2A,G0,A,G0)') trim(message),'history buffer hru size:',this%nHru,'/= input data hru size:',nHru_input
        ierr=81; return
      end if
      if (nRch_input/=this%nRch) then
        write(message,'(2A,G0,A,G0)') trim(message),'history buffer reach size:',this%nRch,'/= input data reach size:',nRch_input
        ierr=81; return
      end if

      ! ---- aggregate
      ! 1. basin runoff
      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        this%basRunoff(1:this%nHRU) = this%basRunoff(1:this%nHRU) + basRunoff_local(1:this%nHru)
      end if

      ! 2. instantaneous runoff into reach
      if (meta_rflx(ixRFLX%instRunoff)%varFile) then
        this%instRunoff(1:this%nRch) = this%instRunoff(1:this%nRch) + RCHFLX_local(1,1:this%nRch)%BASIN_QI
      end if

      ! 3. delayed runoff into reach
      if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
        this%dlayRunoff(1:this%nRch) = this%dlayRunoff(1:this%nRch) + RCHFLX_local(1,1:this%nRch)%BASIN_QR(1)
      end if

      ! 4. discharge and volume
      do iRoute = 1, nRoutes
        select case(routeMethods(iRoute))
          case(accumRunoff);           idxMethod=idxSUM
          case(impulseResponseFunc);   idxMethod=idxIRF
          case(kinematicWaveTracking); idxMethod=idxKWT
          case(kinematicWave);         idxMethod=idxKW
          case(muskingumCunge);        idxMethod=idxMC
          case(diffusiveWave);         idxMethod=idxDW
          case default
            write(message,'(2A,X,G0,X,A)') trim(message), 'routing method index:',routeMethods(iRoute), 'must be 0-5'
            ierr=81; return
        end select

        do ix=1,this%nRch
          this%discharge(ix,iRoute) = this%discharge(ix,iRoute) + RCHFLX_local(1,ix)%ROUTE(idxMethod)%REACH_Q
          this%volume(ix,iRoute)    = this%volume(ix,iRoute) + RCHFLX_local(1,ix)%ROUTE(idxMethod)%REACH_VOL(1)
        end do
      end do

    END SUBROUTINE aggregate

    ! ---------------------------------------------------------------
    ! finalize output variable (compute mean)
    ! ---------------------------------------------------------------
    SUBROUTINE finalize(this)

      implicit none
      ! Argument variables
      class(histVars), intent(inout)  :: this

      ! ---- aggregate
      ! 1. basin runoff
      if (allocated(this%basRunoff)) then
        this%basRunoff = this%basRunoff/real(this%nt, kind=dp)
      end if

      ! 2. instantaneous runoff into reach
      if (allocated(this%instRunoff)) then
        this%instRunoff = this%instRunoff/real(this%nt, kind=dp)
      end if

      ! 3. delayed runoff into reach
      if (allocated(this%dlayRunoff)) then
        this%dlayRunoff = this%dlayRunoff/real(this%nt, kind=dp)
      end if

      ! 4. discharge
      if (allocated(this%discharge)) then
        this%discharge = this%discharge/real(this%nt, kind=dp)
      end if

      ! 5. volume
      if (allocated(this%volume)) then
        this%volume    = this%volume/real(this%nt, kind=dp)
      end if

    END SUBROUTINE finalize

    ! ---------------------------------------------------------------
    ! re-initialzae intantiated data structure
    ! ---------------------------------------------------------------
    SUBROUTINE refresh(this)

      implicit none
      ! Argument variables
      class(histVars), intent(inout)  :: this

      this%nt = 0

      if (allocated(this%basRunoff))  this%basRunoff = 0._dp
      if (allocated(this%instRunoff)) this%instRunoff = 0._dp
      if (allocated(this%dlayRunoff)) this%dlayRunoff = 0._dp
      if (allocated(this%discharge))  this%discharge = 0._dp
      if (allocated(this%volume))     this%volume = 0._dp

    END SUBROUTINE refresh

    ! ---------------------------------------------------------------
    ! Release memory (deallocate array)
    ! ---------------------------------------------------------------
    SUBROUTINE clean(this)

      implicit none
      ! Argument variables
      class(histVars), intent(inout)  :: this

      if (allocated(this%basRunoff))  deallocate(this%basRunoff)
      if (allocated(this%instRunoff)) deallocate(this%instRunoff)
      if (allocated(this%dlayRunoff)) deallocate(this%dlayRunoff)
      if (allocated(this%discharge))  deallocate(this%discharge)
      if (allocated(this%volume))     deallocate(this%volume)

    END SUBROUTINE clean

    ! ---------------------------------------------------------------
    ! instantiate histVars by reading restart file
    ! ---------------------------------------------------------------
    SUBROUTINE read_restart(this, restart_name, ierr, message)

      implicit none
      ! Argument variables
      class(histVars)                    :: this
      character(*),         intent(in)   :: restart_name          ! restart file name
      integer(i4b),         intent(out)  :: ierr                  ! error code
      character(len=strLen),intent(out)  :: message               ! error message
      ! local variable
      character(len=strLen)              :: cmessage              ! error message from subroutines
      real(dp), allocatable              :: array_tmp(:)          ! temp array
      integer(i4b)                       :: ixRoute, ix1, ix2     ! loop index
      integer(i4b)                       :: ixFlow, ixVol         ! temporal method index
      logical(lgt)                       :: FileStatus            ! file open or close
      integer(i4b), allocatable          :: compdof_rch(:)        ! global reach index for PIO reading
      integer(i4b), allocatable          :: compdof_hru(:)        ! globa hru index for PIO reading
      type(file_desc_t)                  :: pioFileDesc           ! pio file handle
      type(io_desc_t)                    :: ioDescRchHist         ! PIO domain decomposition data for reach flux [nRch]
      type(io_desc_t)                    :: ioDescHruHist         ! PIO domain decomposition data for hru runoff [nHRU]

      ierr=0; message='read_restart/'

      call openFile(pioSystem, pioFileDesc, trim(restart_name), pio_typename, ncd_nowrite, FileStatus, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call get_nc(trim(restart_name), 'history_time', this%timeVar, 1, ierr, cmessage) ! get restart datetime
      call get_nc(trim(restart_name), 'nt', this%nt, 1, ierr, cmessage)                ! get number of aggretation time steps

      ! get number of HRUs and Reaches for each core
      if (masterproc) then
        this%nHru =nHru_mainstem+nHru_trib
        this%nRch =nRch_mainstem+nRch_trib+nTribOutlet
        allocate(array_tmp(nRch_mainstem+nRch_trib))
      else
        this%nHru = nHru_trib
        this%nRch = nRch_trib
        allocate(array_tmp(nRch_trib))
      end if

      ! global reach and hru index
      if (masterproc) then
        allocate(compdof_rch(sum(rch_per_proc(-1:pid))))
        allocate(compdof_hru(sum(hru_per_proc(-1:pid))))
      else
        allocate(compdof_rch(rch_per_proc(pid)))
        allocate(compdof_hru(hru_per_proc(pid)))
      end if

      ! For reach flux/volume
      if (masterproc) then
        ix1 = 1
      else
        ix1 = sum(rch_per_proc(-1:pid-1))+1
      endif
      ix2 = sum(rch_per_proc(-1:pid))
      compdof_rch = arth(ix1, 1, ix2-ix1+1)

      ! For HRU flux/volume
      if (masterproc) then
        ix1 = 1
      else
        ix1 = sum(hru_per_proc(-1:pid-1))+1
      endif
      ix2 = sum(hru_per_proc(-1:pid))
      compdof_hru = arth(ix1, 1, ix2-ix1+1)

      ! initialze ioDesc
      call pio_decomp(pioSystem,         & ! input: pio system descriptor
                      ncd_double,        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [nRch],            & ! input: dimension length == global array size
                      compdof_rch,       & ! input: local->global mapping
                      ioDescRchHist)

      call pio_decomp(pioSystem,         & ! input: pio system descriptor
                      ncd_double,        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [nHRU],            & ! input: dimension length == global array size
                      compdof_hru,       & ! input: local->global mapping
                      ioDescHruHist)

      ! reading basin runoff (hru dimension)
      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        allocate(this%basRunoff(this%nHru), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%basRunoff]'; return; endif

        call read_dist_array(pioFileDesc, meta_hflx(ixHFLX%basRunoff)%varName, this%basRunoff, ioDescHruHist, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

      ! reading instantaneous runoff (rch dimension)
      if (meta_rflx(ixRFLX%instRunoff)%varFile) then
        allocate(this%instRunoff(this%nRch), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%instRunoff]'; return; endif

        call read_dist_array(pioFileDesc, meta_rflx(ixRFLX%instRunoff)%varName, array_tmp, ioDescRchHist, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
        if (masterproc) then
          this%instRunoff(1:nRch_mainstem) = array_tmp(1:nRch_mainstem)
          this%instRunoff(nRch_mainstem+nTribOutlet+1:this%nRch) = array_tmp(nRch_mainstem+1:nRch_mainstem+nRch_trib)
        else
          this%instRunoff = array_tmp
        end if
      end if

      ! reading overland routed runoff (rch dimension)
      if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
        allocate(this%dlayRunoff(this%nRch), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%dlayRunoff]'; return; endif

        call read_dist_array(pioFileDesc, meta_rflx(ixRFLX%dlayRunoff)%varName, array_tmp, ioDescRchHist, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
        if (masterproc) then
          this%dlayRunoff(1:nRch_mainstem) = array_tmp(1:nRch_mainstem)
          this%dlayRunoff(nRch_mainstem+nTribOutlet+1:this%nRch) = array_tmp(nRch_mainstem+1:nRch_mainstem+nRch_trib)
        else
          this%dlayRunoff = array_tmp
        end if
      end if

      ! reading reach routed runoff and volume
      if (nRoutes>0) then
        allocate(this%discharge(this%nRch, nRoutes), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%discharge]'; return; endif
        allocate(this%volume(this%nRch, nRoutes), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%discharge]'; return; endif

        this%discharge = 0._dp
        this%volume    = 0._dp

        do ixRoute=1,nRoutes
          select case(routeMethods(ixRoute))
            case(accumRunoff)
              ixFlow=ixRFLX%sumUpstreamRunoff
            case(impulseResponseFunc)
              ixFlow=ixRFLX%IRFroutedRunoff
              ixVol=ixRFLX%IRFvolume
            case(kinematicWaveTracking)
              ixFlow=ixRFLX%KWTroutedRunoff
              ixVol=ixRFLX%KWTvolume
            case(kinematicWave)
              ixFlow=ixRFLX%KWroutedRunoff
              ixVol=ixRFLX%KWvolume
            case(muskingumCunge)
              ixFlow=ixRFLX%MCroutedRunoff
              ixVol=ixRFLX%MCvolume
            case(diffusiveWave)
              ixFlow=ixRFLX%DWroutedRunoff
              ixVol=ixRFLX%DWvolume
            case default
              write(message,'(2A,X,G0,X,A)') trim(message), 'routing method index:',routeMethods(ixRoute), 'must be 0-5'
              ierr=81; return
          end select

          if (meta_rflx(ixFlow)%varFile) then
            call read_dist_array(pioFileDesc, meta_rflx(ixFlow)%varName, array_tmp, ioDescRchHist, ierr, cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
          end if
          ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
          if (masterproc) then
            this%discharge(1:nRch_mainstem, ixRoute) = array_tmp(1:nRch_mainstem)
            this%discharge(nRch_mainstem+nTribOutlet+1:this%nRch, ixRoute) = array_tmp(nRch_mainstem+1:nRch_mainstem+nRch_trib)
          else
            this%discharge(:,ixRoute) = array_tmp
          end if

          if (routeMethods(ixRoute)==accumRunoff) cycle  ! accumuRunoff has only discharge

          if (meta_rflx(ixVol)%varFile) then
            call read_dist_array(pioFileDesc, meta_rflx(ixVol)%varName, this%volume(1:nRch,ixRoute), ioDescRchHist, ierr, cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
          end if
          ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
          if (masterproc) then
            this%volume(1:nRch_mainstem, ixRoute) = array_tmp(1:nRch_mainstem)
            this%volume(nRch_mainstem+nTribOutlet+1:this%nRch,:) = this%volume(nRch_mainstem+1:nRch_mainstem+nRch_trib,:)
          else
            this%volume(:,ixRoute) = array_tmp
          end if
        end do
      end if

      call closeFile(pioFileDesc, FileStatus)

    END SUBROUTINE read_restart

END MODULE histVars_data
