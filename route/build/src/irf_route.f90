MODULE irf_route_module

USE nrtype
! data type
USE dataTypes, ONLY: STRFLX           ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO          ! Network topology
USE dataTypes, ONLY: RCHPRP           ! Reach parameter
USE dataTypes, ONLY: irfRCH           ! irf specific state data structure
 USE dataTypes,ONLY: subbasin_omp     ! mainstem+tributary data structures
! global parameters
USE public_var, ONLY: iulog           ! i/o logical unit number
USE public_var, ONLY: realMissing     ! missing value for real number
USE public_var, ONLY: integerMissing  ! missing value for integer number
USE public_var, ONLY: dt              ! routing time step duration [sec]
USE public_var, ONLY: qmodOption      ! qmod option (use 1==direct insertion)
USE public_var, ONLY: ntsQmodStop     ! number of time steps for which direct insertion is performed
USE public_var, ONLY: QerrTrend       ! temporal discharge error trend: 1->constant,2->linear, 3->logistic
USE globalData, ONLY: nThreads        ! number of threads used for openMP
USE globalData, ONLY: idxIRF          ! index of IRF method
! subroutines: general
USE model_finalize, ONLY : handle_err

implicit none

private
public::irf_route

CONTAINS

 ! *********************************************************************
 ! subroutine: perform network UH routing
 ! *********************************************************************
 SUBROUTINE irf_route(iEns,          &  ! input: index of runoff ensemble to be processed
                      river_basin,   &  ! input: river basin information (mainstem, tributary outlet etc.)
                      ixDesire,      &  ! input: reachID to be checked by on-screen pringing
                      NETOPO_in,     &  ! input: reach topology data structure
                      RPARAM_in,     &  ! input: reach parameter data structure
                      RCHFLX_out,    &  ! inout: reach flux data structure
                      ierr, message, &  ! output: error control
                      ixSubRch)         ! optional input: subset of reach indices to be processed

 implicit none
 ! argument variables
 integer(i4b),       intent(in)                  :: iEns                ! runoff ensemble to be routed
 type(subbasin_omp), intent(in),    allocatable  :: river_basin(:)      ! river basin information (mainstem, tributary outlet etc.)
 integer(i4b),       intent(in)                  :: ixDesire            ! index of the reach for verbose output ! Output
 type(RCHTOPO),      intent(in),    allocatable  :: NETOPO_in(:)        ! River Network topology
 type(RCHPRP),       intent(in),    allocatable  :: RPARAM_in(:)        ! River reach parameter
 TYPE(STRFLX),       intent(inout)               :: RCHFLX_out(:,:)     ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b),       intent(out)                 :: ierr                ! error code
 character(*),       intent(out)                 :: message             ! error message
 integer(i4b),       intent(in),    optional     :: ixSubRch(:)         ! subset of reach indices to be processed
 ! Local variables
 character(len=strLen)                           :: cmessage            ! error message from subroutine
 logical(lgt),                      allocatable  :: doRoute(:)          ! logical to indicate which reaches are processed
 integer(i4b)                                    :: nOrder              ! number of stream order
 integer(i4b)                                    :: nTrib               ! number of tributary basins
 integer(i4b)                                    :: nSeg                ! number of reaches in the network
 integer(i4b)                                    :: iSeg, jSeg          ! loop indices - reach
 integer(i4b)                                    :: iTrib               ! loop indices - branch
 integer(i4b)                                    :: ix                  ! loop indices stream order
 ! variables needed for timing
 !integer(i4b)                                    :: omp_get_thread_num
 !integer(i4b), allocatable                       :: ixThread(:)         ! thread id
 !integer*8,    allocatable                       :: openMPend(:)        ! time for the start of the parallelization section
 !integer*8,    allocatable                       :: timeTribStart(:)    ! time Tributaries start
 !real(dp),     allocatable                       :: timeTrib(:)         ! time spent on each Tributary

 ierr=0; message='irf_route/'

 ! number of reach check
 if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
  ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
 endif

 nSeg = size(RCHFLX_out(iens,:))

 allocate(doRoute(nSeg), stat=ierr)

 if (present(ixSubRch))then
  doRoute(:)=.false.
  doRoute(ixSubRch) = .true. ! only subset of reaches are on
 else
  doRoute(:)=.true. ! every reach is on
 endif

 nOrder = size(river_basin)


 do ix = 1,nOrder

   nTrib=size(river_basin(ix)%branch)

!  allocate(ixThread(nTrib), openMPend(nTrib), timeTrib(nTrib), timeTribStart(nTrib), stat=ierr)
!  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to allocate space for Trib timing'; return; endif
!  timeTrib(:) = realMissing
!  ixThread(:) = integerMissing

  ! 1. Route tributary reaches (parallel)
!$OMP PARALLEL DO schedule(dynamic,1)                   &
!$OMP          private(jSeg, iSeg)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(doRoute)                          & ! data array shared
!$OMP          shared(NETOPO_in)                        & ! data structure shared
!$OMP          shared(RPARAM_in)                        & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(ix, iEns, ixDesire)               & ! indices shared
!$OMP          firstprivate(nTrib)
!!$OMP          shared(openMPend, nThreads)              & ! timing variables shared
!!$OMP          shared(timeTribStart)                    & ! timing variables shared
!!$OMP          shared(timeTrib)                         & ! timing variables shared
!!$OMP          shared(ixThread)                         & ! thread id array shared
   trib:do iTrib = 1,nTrib
!!$    ixThread(iTrib) = omp_get_thread_num()
!    call system_clock(timeTribStart(iTrib))
     seg:do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
       jSeg = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
       if (.not. doRoute(jSeg)) cycle
       call irf_rch(iEns, jSeg, ixDesire, NETOPO_IN, RPARAM_in, RCHFLX_out, ierr, cmessage)
       if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
     end do seg
!    call system_clock(openMPend(iTrib))
!    timeTrib(iTrib) = real(openMPend(iTrib)-timeTribStart(iTrib), kind(dp))
   end do trib
!$OMP END PARALLEL DO

!  write(*,'(a)') 'iTrib nSeg ixThread nThreads StartTime EndTime'
!  do iTrib=1,nTrib
!    write(*,'(4(i5,1x),2(I20,1x))') iTrib, river_basin(iOut)%branch(iTrib)%nRch, ixThread(iTrib), nThreads, timeTribStart(iTrib), openMPend(iTrib)
!  enddo
!  deallocate(ixThread, openMPend, timeTrib, timeTribStart, stat=ierr)
!  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to deallocate space for Trib timing'; return; endif

 end do ! basin loop

 end subroutine irf_route

 ! *********************************************************************
 ! subroutine: perform one segment route UH routing
 ! *********************************************************************
 SUBROUTINE irf_rch(iEns,         & ! input: index of runoff ensemble to be processed
                    segIndex,     & ! input: index of runoff ensemble to be processed
                    ixDesire,     & ! input: reachID to be checked by on-screen pringing
                    NETOPO_in,    & ! input: reach topology data structure
                    RPARAM_in,    & ! input: reach parameter data structure
                    RCHFLX_out,   & ! inout: reach flux data structure
                    ierr, message)  ! output: error control

 implicit none
 ! Argument variables
 integer(i4b), intent(in)                 :: iEns           ! runoff ensemble to be routed
 integer(i4b), intent(in)                 :: segIndex       ! segment where routing is performed
 integer(i4b), intent(in)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)   ! River Network topology
 type(RCHPRP), intent(in),    allocatable :: RPARAM_in(:)   ! River reach parameter
 type(STRFLX), intent(inout)              :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables
 real(dp)                                 :: q_upstream     ! total discharge at top of the reach being processed
 real(dp)                                 :: WB_error       ! water balance error [m3/s]
 real(dp)                                 :: Qcorrect       ! Discharge correction (when qmodOption=1) [m3/s]
 real(dp)                                 :: k              ! the logistic decay rate or steepness of the curve.
 real(dp)                                 :: y0             !
 real(dp)                                 :: x0             !
 integer(i4b)                             :: nUps           ! number of upstream segment
 integer(i4b)                             :: iUps           ! upstream reach index
 integer(i4b)                             :: iRch_ups       ! index of upstream reach in NETOPO
 integer(i4b)                             :: ntdh           ! number of time steps in IRF
 integer(i4b)                             :: itdh           ! loop index for unit hydrograph
 character(len=strLen)                    :: fmt1           ! format string
 character(len=strLen)                    :: cmessage       ! error message from subroutine
 integer(i4b),parameter                   :: const=1        ! error reduction type: step function
 integer(i4b),parameter                   :: linear=2       ! error reduction type: linear function
 integer(i4b),parameter                   :: logistic=3     ! error reduction type: logistic function
 integer(i4b),parameter                   :: exponential=4  ! error reduction type: exponential function

 ierr=0; message='irf_rch/'

 ! initialize future discharge array at first time
  if (.not.allocated(RCHFLX_out(iens,segIndex)%QFUTURE_IRF))then

   ntdh = size(NETOPO_in(segIndex)%UH)

   allocate(RCHFLX_out(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX_out(iens,segIndex)%QFUTURE_IRF'; return; endif

   RCHFLX_out(iens,segIndex)%QFUTURE_IRF(:) = 0._dp

  end if

  ! get discharge coming from upstream
  nUps = size(NETOPO_in(segIndex)%UREACHI)
  q_upstream = 0.0_dp
  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach

      if (qmodOption==1) then
        if (RCHFLX_out(iens,iRch_ups)%QOBS>0._dp) then ! there is observation
          RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%Qerror = RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%REACH_Q - RCHFLX_out(iens,iRch_ups)%QOBS ! compute error
        end if
        if (RCHFLX_out(iens,iRch_ups)%Qelapsed > ntsQmodStop) then
          RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%Qerror=0._dp
        end if
        if (RCHFLX_out(iens,iRch_ups)%Qelapsed <= ntsQmodStop) then
          select case(QerrTrend)
            case(const)
              Qcorrect = RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%Qerror
            case(linear)
              Qcorrect = RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%Qerror*(1._dp - real(RCHFLX_out(iens,iRch_ups)%Qelapsed,dp)/real(ntsQmodStop, dp))
            case(logistic)
              x0 =0.25; y0 =0.90
              k = log(1._dp/y0-1._dp)/(ntsQmodStop/2._dp-ntsQmodStop*x0)
              Qcorrect = RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%Qerror/(1._dp + exp(-k*(1._dp*RCHFLX_out(iens,iRch_ups)%Qelapsed-ntsQmodStop/2._dp)))
            case(exponential)
              if (RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%Qerror/=0._dp) then
                k = log(0.1_dp/abs(RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%Qerror))/(1._dp*ntsQmodStop)
                Qcorrect = RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%Qerror*exp(k*RCHFLX_out(iens,iRch_ups)%Qelapsed)
              else
                Qcorrect = 0._dp
              end if
            case default; message=trim(message)//'discharge error trend model must be 1(const),2(liear), or 3(logistic)'; ierr=81; return
          end select
        else
          Qcorrect=0._dp
        end if
        RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%REACH_Q = max(RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%REACH_Q-Qcorrect, 0._dp)
      end if

      q_upstream = q_upstream + RCHFLX_out(iens,iRch_ups)%ROUTE(idxIRF)%REACH_Q
    end do
  endif

  ! perform UH convolution
  call conv_upsbas_qr(NETOPO_in(segIndex)%UH,         &  ! input: reach unit hydrograph
                      q_upstream,                     &  ! input: total discharge at top of the reach being processed
                      RCHFLX_out(iens,segIndex),      &  ! inout: updated fluxes at reach
                      RCHFLX_out(iens,segIndex)%TAKE, &  ! input: abstraction(-)/injection(+) [m3/s]
                      RPARAM_in(segIndex)%MINFLOW,    &  ! input: minimum environmental flow [m3/s]
                      ierr, message)                     ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! check
  if(segIndex==ixDesire)then
    ntdh = size(NETOPO_in(segIndex)%UH)
    WB_error = q_upstream - RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q + RCHFLX_out(iens,segIndex)%BASIN_QR(1) &
               - (RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(0))/dt
    write(fmt1,'(A,I5,A)') '(A, 1X',ntdh,'(1X,F20.7))'
    write(*,'(2a)') new_line('a'),'** Check Impulse Response Function routing **'
    write(*,'(a,x,I10,x,I10)') ' Reach index & ID       =', segIndex, NETOPO_in(segIndex)%REACHID
    write(*,fmt1)              ' Unit-Hydrograph        =', (NETOPO_in(segIndex)%UH(itdh), itdh=1,ntdh)
    write(*,'(a)')             ' * inflow [m3/s]: q_upstream'
    write(*,'(a)')             '   local flow [m3/s]: RCHFLX_out(iens,segIndex)%BASIN_QR(1)'
    write(*,'(a)')             '   outflow [m3/s]: RCHFLX_out%ROUTE%REACH_Q'
    write(*,'(a)')             '   volume [m3] at end of previous [0] and current [1] time step: REACH_VOL'
    write(*,'(a)')             '   water balance error [m3/s]: WB error'
    write(*,'(a,x,G15.4)')     ' inflow     =', q_upstream
    write(*,'(a,x,G15.4)')     ' local flow =', RCHFLX_out(iens,segIndex)%BASIN_QR(1)
    write(*,'(a,x,G15.4)')     ' outflow    =', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q - RCHFLX_out(iens,segIndex)%BASIN_QR(1)
    write(*,'(a,x,G20.4)')     ' volume[0]  =', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(0)
    write(*,'(a,x,G20.4)')     ' volume[1]  =', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1)
    write(*,'(a,x,G15.4)')     ' WB error   =', WB_error
  endif

  if (RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) < 0) then
    write(iulog,'(A,X,G12.5,X,A,X,I9)') ' ---- NEGATIVE VOLUME [m3]= ', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1), 'at ', NETOPO_in(segIndex)%REACHID
!    RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) = 0._dp
  end if
 END SUBROUTINE irf_rch


 ! *********************************************************************
 ! subroutine: Compute delayed runoff from the upstream segments
 ! *********************************************************************
 SUBROUTINE conv_upsbas_qr(reach_uh,   &    ! input: reach unit hydrograph
                           q_upstream, &    ! input:
                           rflux,      &    ! input: input flux at reach
                           Qtake,      &    ! input: abstraction(-)/injection(+) [m3/s]
                           Qmin,       &    ! input: minimum environmental flow [m3/s]
                           ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Details: Convolute runoff volume of upstream at one reach at one time step
 ! ----------------------------------------------------------------------------------------

 implicit none
 ! Argument variables
 real(dp),     intent(in)               :: reach_uh(:)  ! reach unit hydrograph
 real(dp),     intent(in)               :: q_upstream   ! total discharge at top of the reach being processed
 real(dp),     intent(in)               :: Qtake        ! abstraction(-)/injection(+) [m3/s]
 real(dp),     intent(in)               :: Qmin         ! minimum environmental flow [m3/s]
 type(STRFLX), intent(inout)            :: rflux        ! current Reach fluxes
 integer(i4b), intent(out)              :: ierr         ! error code
 character(*), intent(out)              :: message      ! error message
 ! Local variables
 real(dp)                               :: QupMod       ! modified total discharge at top of the reach being processed
 real(dp)                               :: Qabs         ! maximum allowable water abstraction rate [m3/s]
 real(dp)                               :: Qmod         ! abstraction rate to be taken from outlet discharge [m3/s]
 integer(i4b)                           :: nTDH         ! number of UH data (i.e., number of future time step
 integer(i4b)                           :: iTDH         ! index of UH data

 ierr=0; message='conv_upsbas_qr/'

 ! Q injection, add at top of reach
 QupMod = q_upstream
 if (Qtake>0) then
   QupMod = QupMod+ Qtake
 end if

 ! place a fraction of runoff in future time steps
 nTDH = size(reach_uh) ! identify the number of future time steps of UH for a given segment
 do iTDH=1,nTDH
   rflux%QFUTURE_IRF(iTDH) = rflux%QFUTURE_IRF(iTDH) &
                             + reach_uh(iTDH)*QupMod
 enddo

 ! compute volume in reach
 rflux%ROUTE(idxIRF)%REACH_VOL(0) = rflux%ROUTE(idxIRF)%REACH_VOL(1)
 ! For very low flow condition, outflow - inflow > current storage, so limit outflow and adjust rflux%QFUTURE_IRF(1)
 rflux%QFUTURE_IRF(1) = min(rflux%ROUTE(idxIRF)%REACH_VOL(0)/dt + QupMod*0.999, rflux%QFUTURE_IRF(1))
 rflux%ROUTE(idxIRF)%REACH_VOL(1) = rflux%ROUTE(idxIRF)%REACH_VOL(0) + (QupMod - rflux%QFUTURE_IRF(1))*dt

 ! Add local routed flow at the bottom of reach
 rflux%ROUTE(idxIRF)%REACH_Q = rflux%QFUTURE_IRF(1) + rflux%BASIN_QR(1)

 ! Q abstraction
 ! Compute actual abstraction (Qabs) m3/s - values should be negative
 ! Compute abstraction (Qmod) m3 taken from outlet discharge (REACH_Q)
 ! Compute REACH_Q subtracted from Qmod abstraction
 ! Compute REACH_VOL subtracted from total abstraction minus abstraction from outlet discharge
 if (Qtake<0) then
   Qabs               = max(-(rflux%ROUTE(idxIRF)%REACH_VOL(1)/dt+rflux%ROUTE(idxIRF)%REACH_Q), Qtake)
   Qmod               = min(rflux%ROUTE(idxIRF)%REACH_VOL(1) + Qabs*dt, 0._dp)
   rflux%ROUTE(idxIRF)%REACH_Q      = max(rflux%ROUTE(idxIRF)%REACH_Q + Qmod/dt, Qmin)
   rflux%ROUTE(idxIRF)%REACH_VOL(1) = rflux%ROUTE(idxIRF)%REACH_VOL(1) + (Qabs*dt - Qmod)
 end if

 ! move array back   use eoshift
 rflux%QFUTURE_IRF=eoshift(rflux%QFUTURE_IRF,shift=1)

 rflux%QFUTURE_IRF(ntdh) = 0._dp

 END SUBROUTINE conv_upsbas_qr

END MODULE irf_route_module
