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
    real(sp), allocatable    :: basRunoff(:)       ! catchment average runoff [m/s] [nHru]
    real(sp), allocatable    :: instRunoff(:)      ! instantaneous lateral inflow into each river/lake [m3/s]  [nRch]
    real(sp), allocatable    :: dlayRunoff(:)      ! lateral inflow into river or lake [m3/s] for each reach [nRch]
    real(sp), allocatable    :: discharge(:,:)     ! river/lake discharge [m3/s] for each reach/lake and routing method [nRch,nMethod]
    real(sp), allocatable    :: elev(:,:)          ! river/lake surface water elevation [m] for each reach/lake and routing method [nRch,nMethod]
    real(sp), allocatable    :: volume(:,:)        ! river/lake volume [m3] for each reach/lake and routing method [nRch,nMethod]

    CONTAINS

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
        instHistVar%basRunoff(1:nHRU_local) = 0._sp
      end if

      if (meta_rflx(ixRFLX%instRunoff)%varFile) then
        allocate(instHistVar%instRunoff(nRch_local), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%instRunoff]'; return; endif
        instHistVar%instRunoff(1:nRch_local) = 0._sp
      end if

      if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
        allocate(instHistVar%dlayRunoff(nRch_local), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%instRunoff]'; return; endif
        instHistVar%dlayRunoff(1:nRch_local) = 0._sp
      end if

      if (nRoutes>0) then ! this should be number of methods that ouput
        allocate(instHistVar%discharge(nRch_local, nRoutes), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%discharge]'; return; endif
        instHistVar%discharge(1:nRch_local, 1:nRoutes) = 0._sp

        allocate(instHistVar%volume(nRch_local, nRoutes), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [instHistVar%volume]'; return; endif
        instHistVar%volume(1:nRch_local, 1:nRoutes) = 0._sp
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
        write(message,'(2A,I,A,I)') trim(message),'history buffer hru size:',this%nHru,'/= input data hru size:',nHru_input
        ierr=81; return
      end if
      if (nRch_input/=this%nRch) then
        write(message,'(2A,I,A,I)') trim(message),'history buffer reach size:',this%nRch,'/= input data reach size:',nRch_input
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
            write(message,'(2A,X,I,X,A)') trim(message), 'routing method index:',routeMethods(iRoute), 'must be 0-5'
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
        this%basRunoff = this%basRunoff/real(this%nt, kind=sp)
      end if

      ! 2. instantaneous runoff into reach
      if (allocated(this%instRunoff)) then
        this%instRunoff = this%instRunoff/real(this%nt, kind=sp)
      end if

      ! 3. delayed runoff into reach
      if (allocated(this%dlayRunoff)) then
        this%dlayRunoff = this%dlayRunoff/real(this%nt, kind=sp)
      end if

      ! 4. discharge
      if (allocated(this%discharge)) then
        this%discharge = this%discharge/real(this%nt, kind=sp)
      end if

      ! 5. volume
      if (allocated(this%volume)) then
        this%volume    = this%volume/real(this%nt, kind=sp)
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

      if (allocated(this%basRunoff))  this%basRunoff = 0._sp
      if (allocated(this%instRunoff)) this%instRunoff = 0._sp
      if (allocated(this%dlayRunoff)) this%dlayRunoff = 0._sp
      if (allocated(this%discharge))  this%discharge = 0._sp
      if (allocated(this%volume))     this%volume = 0._sp

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

END MODULE histVars_data
