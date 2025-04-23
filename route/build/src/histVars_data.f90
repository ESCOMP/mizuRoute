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
  USE public_var,        ONLY: outputInflow               ! logical for outputting upstream inflow in history file
  USE public_var,        ONLY: floodplain                 ! logical for outputting upstream inflow in history file
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
  USE globalData,        ONLY: pioSystem
  USE globalData,        ONLY: ioDesc_hru_double
  USE globalData,        ONLY: ioDesc_hist_rch_double
  USE public_var,        ONLY: pio_typename
  USE ncio_utils,        ONLY: get_nc
  USE pio_utils,         ONLY: file_desc_t
  USE pio_utils,         ONLY: ncd_nowrite
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
    integer(i4b)             :: nHru               ! number of HRUs
    integer(i4b)             :: nRch               ! number of Reaches
    real(dp)                 :: timeVar(2)         ! time variables [time unit] at end points of current aggregaed step
    real(dp), allocatable    :: basRunoff(:)       ! catchment average runoff [m/s] [nHru]
    real(dp), allocatable    :: instRunoff(:)      ! instantaneous lateral inflow into each river/lake [m3/s]  [nRch]
    real(dp), allocatable    :: dlayRunoff(:)      ! lateral inflow into river or lake [m3/s] for each reach [nRch]
    real(dp), allocatable    :: discharge(:,:)     ! river/lake discharge [m3/s] for each reach/lake and routing method [nRch,nMethod]
    real(dp), allocatable    :: inflow(:,:)        ! inflow from upstream rivers/lakes [m3/s] for each reach/lake and routing method [nRch,nMethod]
    real(dp), allocatable    :: waterHeight(:,:)   ! river/lake surface water elevation [m] for each reach/lake and routing method [nRch,nMethod]
    real(dp), allocatable    :: floodVolume(:,:)   ! river/lake volume [m3] for each reach/lake and routing method [nRch,nMethod]
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
        allocate(instHistVar%basRunoff(nHRU_local), source=0._dp, stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%basRunoff]'; return; endif
      end if

      if (meta_rflx(ixRFLX%instRunoff)%varFile) then
        allocate(instHistVar%instRunoff(nRch_local), source=0._dp, stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%instRunoff]'; return; endif
      end if

      if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
        allocate(instHistVar%dlayRunoff(nRch_local), source=0._dp, stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%instRunoff]'; return; endif
      end if

      if (nRoutes>0) then ! this should be number of methods that ouput
        allocate(instHistVar%discharge(nRch_local, nRoutes), source=0._dp, stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%discharge]'; return; endif

        allocate(instHistVar%volume(nRch_local, nRoutes), source=0._dp, stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%volume]'; return; endif

        if (floodplain) then
          allocate(instHistVar%floodVolume(nRch_local, nRoutes), source=0._dp, stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%floodVolume]'; return; endif

          allocate(instHistVar%waterHeight(nRch_local, nRoutes), source=0._dp, stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%floodVolume]'; return; endif
        end if

        if (outputInflow) then
          allocate(instHistVar%inflow(nRch_local, nRoutes), source=0._dp, stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%inflow]'; return; endif
        end if
      end if

    END FUNCTION constructor

    ! ---------------------------------------------------------------
    ! accumulate data
    ! ---------------------------------------------------------------
    SUBROUTINE aggregate(this,             & ! inout:
                         timeVar_local,    & ! input: time variables [time unit] at endpoints of current simulation step
                         basRunoff_local,  & ! input: HRU average runoff depth [m/s]
                         RCHFLX_local,     & ! input: Reach flux data structure
                         ierr, message)      ! output: error handling

      implicit none
      ! argument variables
      class(histVars),   intent(inout)   :: this
      real(dp),          intent(in)      :: timeVar_local(2)    ! time variables [time unit] at endpoints of current simulation step
      real(dp),          intent(in)      :: basRunoff_local(:)  ! HRU average runoff depth [m/s]
      type(STRFLX),      intent(in)      :: RCHFLX_local(:,:)   ! Reach flux data structure
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
        this%timeVar(1) = timeVar_local(1)
        this%timeVar(2) = timeVar_local(2)
      else
        this%timeVar(2) = timeVar_local(2)
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
            write(message,'(2A,1X,G0,1X,A)') trim(message), 'routing method index:',routeMethods(iRoute), 'must be 0-5'
            ierr=81; return
        end select

        do ix=1,this%nRch
          this%discharge(ix,iRoute) = this%discharge(ix,iRoute) + RCHFLX_local(1,ix)%ROUTE(idxMethod)%REACH_Q
          this%volume(ix,iRoute)    = this%volume(ix,iRoute) + RCHFLX_local(1,ix)%ROUTE(idxMethod)%REACH_VOL(1)
          if (floodplain) then
            this%floodVolume(ix,iRoute)    = this%floodVolume(ix,iRoute) + RCHFLX_local(1,ix)%ROUTE(idxMethod)%FLOOD_VOL(1)
            this%waterHeight(ix,iRoute)    = this%waterHeight(ix,iRoute) + RCHFLX_local(1,ix)%ROUTE(idxMethod)%REACH_ELE
          end if
          if (outputInflow) then
            this%inflow(ix,iRoute)  = this%inflow(ix,iRoute) + RCHFLX_local(1,ix)%ROUTE(idxMethod)%REACH_INFLOW
          end if
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
        this%volume = this%volume/real(this%nt, kind=dp)
      end if

      ! 6. volume
      if (allocated(this%waterHeight)) then
        this%waterHeight = this%waterHeight/real(this%nt, kind=dp)
      end if

      ! 7. volume
      if (allocated(this%floodVolume)) then
        this%floodVolume = this%floodVolume/real(this%nt, kind=dp)
      end if

      ! 8. inflow
      if (allocated(this%inflow)) then
        this%inflow = this%inflow/real(this%nt, kind=dp)
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

      if (allocated(this%basRunoff))   this%basRunoff = 0._dp
      if (allocated(this%instRunoff))  this%instRunoff = 0._dp
      if (allocated(this%dlayRunoff))  this%dlayRunoff = 0._dp
      if (allocated(this%discharge))   this%discharge = 0._dp
      if (allocated(this%volume))      this%volume = 0._dp
      if (allocated(this%waterHeight)) this%waterHeight = 0._dp
      if (allocated(this%floodVolume)) this%floodVolume = 0._dp
      if (allocated(this%inflow))      this%inflow = 0._dp

    END SUBROUTINE refresh

    ! ---------------------------------------------------------------
    ! Release memory (deallocate array)
    ! ---------------------------------------------------------------
    SUBROUTINE clean(this)

      implicit none
      ! Argument variables
      class(histVars), intent(inout)  :: this

      if (allocated(this%basRunoff))   deallocate(this%basRunoff)
      if (allocated(this%instRunoff))  deallocate(this%instRunoff)
      if (allocated(this%dlayRunoff))  deallocate(this%dlayRunoff)
      if (allocated(this%discharge))   deallocate(this%discharge)
      if (allocated(this%volume))      deallocate(this%volume)
      if (allocated(this%waterHeight)) deallocate(this%waterHeight)
      if (allocated(this%floodVolume)) deallocate(this%floodVolume)
      if (allocated(this%inflow))      deallocate(this%inflow)

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
      integer(i4b)                       :: ixRoute               ! loop index
      integer(i4b)                       :: ixFlow, ixVol         ! temporal discharge, volume variable indices
      integer(i4b)                       :: ixWaterH              ! temporal water height variable index
      integer(i4b)                       :: ixFloodV              ! temporal flood volume variable index
      integer(i4b)                       :: ixInflow              ! temporal inflow variable index
      logical(lgt)                       :: FileStatus            ! file open or close
      type(file_desc_t)                  :: pioFileDesc           ! pio file handle

      ierr=0; message='read_restart/'

      call openFile(pioSystem, pioFileDesc, trim(restart_name), pio_typename, ncd_nowrite, FileStatus, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call get_nc(trim(restart_name), 'history_time', this%timeVar, 1, 2, ierr, cmessage) ! get time variables at time step endpoints
      call get_nc(trim(restart_name), 'nt', this%nt, 1, ierr, cmessage)                   ! get the number of sampled time steps

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

      ! reading basin runoff (hru dimension)
      if (meta_hflx(ixHFLX%basRunoff)%varFile) then
        allocate(this%basRunoff(this%nHru), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%basRunoff]'; return; endif

        call read_dist_array(pioFileDesc, meta_hflx(ixHFLX%basRunoff)%varName, this%basRunoff, ioDesc_hru_double, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

      ! reading instantaneous runoff (rch dimension)
      if (meta_rflx(ixRFLX%instRunoff)%varFile) then
        allocate(this%instRunoff(this%nRch), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%instRunoff]'; return; endif

        call read_dist_array(pioFileDesc, meta_rflx(ixRFLX%instRunoff)%varName, array_tmp, ioDesc_hist_rch_double, ierr, cmessage)
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

        call read_dist_array(pioFileDesc, meta_rflx(ixRFLX%dlayRunoff)%varName, array_tmp, ioDesc_hist_rch_double, ierr, cmessage)
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
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%volume]'; return; endif

        this%discharge = 0._dp
        this%volume    = 0._dp

        if (floodplain) then
          allocate(this%floodVolume(this%nRch, nRoutes), source=0.0_dp, stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%floodVolume]'; return; endif

          allocate(this%waterHeight(this%nRch, nRoutes), source=0.0_dp, stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%waterHeight]'; return; endif
        end if

        if (outputInflow) then
          allocate(this%inflow(this%nRch, nRoutes), source=0.0_dp, stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [hVars%volume]'; return; endif
        end if

        do ixRoute=1,nRoutes
          select case(routeMethods(ixRoute))
            case(accumRunoff)
              ixFlow=ixRFLX%sumUpstreamRunoff
            case(impulseResponseFunc)
              ixFlow=ixRFLX%IRFroutedRunoff
              ixVol=ixRFLX%IRFvolume
              ixWaterH=ixRFLX%IRFheight
              ixFloodV=ixRFLX%IRFfloodVolume
              ixInflow=ixRFLX%IRFinflow
            case(kinematicWaveTracking)
              ixFlow=ixRFLX%KWTroutedRunoff
              ixVol=ixRFLX%KWTvolume
              ixWaterH=ixRFLX%KWTheight
              ixFloodV=ixRFLX%KWTfloodVolume
              ixInflow=ixRFLX%KWTinflow
            case(kinematicWave)
              ixFlow=ixRFLX%KWroutedRunoff
              ixVol=ixRFLX%KWvolume
              ixWaterH=ixRFLX%KWheight
              ixFloodV=ixRFLX%KWfloodVolume
              ixInflow=ixRFLX%KWinflow
            case(muskingumCunge)
              ixFlow=ixRFLX%MCroutedRunoff
              ixVol=ixRFLX%MCvolume
              ixWaterH=ixRFLX%MCheight
              ixFloodV=ixRFLX%MCfloodVolume
              ixInflow=ixRFLX%MCinflow
            case(diffusiveWave)
              ixFlow=ixRFLX%DWroutedRunoff
              ixVol=ixRFLX%DWvolume
              ixWaterH=ixRFLX%DWheight
              ixFloodV=ixRFLX%DWfloodVolume
              ixInflow=ixRFLX%DWinflow
            case default
              write(message,'(2A,1X,G0,1X,A)') trim(message), 'routing method index:',routeMethods(ixRoute), 'must be 0-5'
              ierr=81; return
          end select

          if (meta_rflx(ixFlow)%varFile) then
            call read_dist_array(pioFileDesc, meta_rflx(ixFlow)%varName, array_tmp, ioDesc_hist_rch_double, ierr, cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
            ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
            if (masterproc) then
              this%discharge(1:nRch_mainstem, ixRoute) = array_tmp(1:nRch_mainstem)
              this%discharge(nRch_mainstem+nTribOutlet+1:this%nRch, ixRoute) = array_tmp(nRch_mainstem+1:nRch_mainstem+nRch_trib)
            else
              this%discharge(:,ixRoute) = array_tmp
            end if
          end if

          if (routeMethods(ixRoute)==accumRunoff) cycle  ! accumuRunoff has only discharge

          if (meta_rflx(ixVol)%varFile) then
            call read_dist_array(pioFileDesc, meta_rflx(ixVol)%varName, array_tmp, ioDesc_hist_rch_double, ierr, cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
            ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
            if (masterproc) then
              this%volume(1:nRch_mainstem, ixRoute) = array_tmp(1:nRch_mainstem)
              this%volume(nRch_mainstem+nTribOutlet+1:this%nRch, ixRoute) = array_tmp(nRch_mainstem+1:nRch_mainstem+nRch_trib)
            else
              this%volume(:,ixRoute) = array_tmp
            end if
          end if

          if (meta_rflx(ixFloodV)%varFile) then
            call read_dist_array(pioFileDesc, meta_rflx(ixFloodV)%varName, array_tmp, ioDesc_hist_rch_double, ierr, cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
            ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
            if (masterproc) then
              this%floodVolume(1:nRch_mainstem, ixRoute) = array_tmp(1:nRch_mainstem)
              this%floodVolume(nRch_mainstem+nTribOutlet+1:this%nRch, ixRoute) = array_tmp(nRch_mainstem+1:nRch_mainstem+nRch_trib)
            else
              this%floodVolume(:,ixRoute) = array_tmp
            end if
          end if

          if (meta_rflx(ixWaterH)%varFile) then
            call read_dist_array(pioFileDesc, meta_rflx(ixWaterH)%varName, array_tmp, ioDesc_hist_rch_double, ierr, cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
            ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
            if (masterproc) then
              this%waterHeight(1:nRch_mainstem, ixRoute) = array_tmp(1:nRch_mainstem)
              this%waterHeight(nRch_mainstem+nTribOutlet+1:this%nRch, ixRoute) = array_tmp(nRch_mainstem+1:nRch_mainstem+nRch_trib)
            else
              this%waterHeight(:,ixRoute) = array_tmp
            end if
          end if

          if (meta_rflx(ixInflow)%varFile) then
            call read_dist_array(pioFileDesc, meta_rflx(ixInflow)%varName, array_tmp, ioDesc_hist_rch_double, ierr, cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
            ! need to shift tributary part in main core to after halo reaches (nTribOutlet)
            if (masterproc) then
              this%inflow(1:nRch_mainstem, ixRoute) = array_tmp(1:nRch_mainstem)
              this%inflow(nRch_mainstem+nTribOutlet+1:this%nRch, ixRoute) = array_tmp(nRch_mainstem+1:nRch_mainstem+nRch_trib)
            else
              this%inflow(:,ixRoute) = array_tmp
            end if
          end if
        end do ! end of ixRoute loop
      end if ! end of nRoute>0 if-statement

      call closeFile(pioFileDesc, FileStatus)

    END SUBROUTINE read_restart

END MODULE histVars_data
