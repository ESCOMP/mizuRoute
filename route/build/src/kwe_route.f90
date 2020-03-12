MODULE kwe_route_module

! WORKING in PROGRESS

!numeric type
USE nrtype
! data types
USE dataTypes, ONLY: STRFLX            ! fluxes in each reach
USE dataTypes, ONLY: STRSTA            ! state in each reach
USE dataTypes, ONLY: RCHTOPO           ! Network topology
USE dataTypes, ONLY: RCHPRP            ! Reach parameter
USE dataTypes, ONLY: EKWRCH            ! Eulerian KW state
! global data
USE public_var, ONLY: iulog             ! i/o logical unit number
USE public_var, ONLY: verySmall         ! a very small value
USE public_var, ONLY: runoffMin         ! minimum runoff
USE public_var, ONLY: realMissing       ! missing value for real number
USE public_var, ONLY: integerMissing    ! missing value for integer number

! privary
implicit none
private

public::kwe_route

contains

 ! *********************************************************************
 ! subroutine: route kinematic waves with Euler solution through the river network
 ! *********************************************************************
 SUBROUTINE kwe_route(iens,                 & ! input: ensemble index
                      river_basin,          & ! input: river basin information (mainstem, tributary outlet etc.)
                      T0,T1,                & ! input: start and end of the time step
                      ixDesire,             & ! input: reachID to be checked by on-screen pringing
                      NETOPO_in,            & ! input: reach topology data structure
                      RPARAM_in,            & ! input: reach parameter data structure
                      RCHSTA_out,           & ! inout: reach state data structure
                      RCHFLX_out,           & ! inout: reach flux data structure
                      ierr,message,         & ! output: error control
                      ixSubRch)               ! optional input: subset of reach indices to be processed

   USE dataTypes, ONLY: subbasin_omp          ! mainstem+tributary data strucuture
   implicit none
   ! Input
   integer(i4b),       intent(in)                 :: iEns                 ! ensemble member
   type(subbasin_omp), intent(in),    allocatable :: river_basin(:)       ! river basin information (mainstem, tributary outlet etc.)
   real(dp),           intent(in)                 :: T0,T1                ! start and end of the time step (seconds)
   integer(i4b),       intent(in)                 :: ixDesire             ! index of the reach for verbose output
   type(RCHTOPO),      intent(in),    allocatable :: NETOPO_in(:)         ! River Network topology
   type(RCHPRP),       intent(in),    allocatable :: RPARAM_in(:)         ! River reach parameter
   ! inout
   type(STRSTA),       intent(inout), allocatable :: RCHSTA_out(:,:)      ! reach state data
   type(STRFLX),       intent(inout), allocatable :: RCHFLX_out(:,:)      ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   ! output variables
   integer(i4b),       intent(out)                :: ierr                 ! error code
   character(*),       intent(out)                :: message              ! error message
   ! input (optional)
   integer(i4b),       intent(in), optional       :: ixSubRch(:)          ! subset of reach indices to be processed
   ! local variables
   character(len=strLen)                          :: cmessage             ! error message for downwind routine
   logical(lgt),                      allocatable :: doRoute(:)           ! logical to indicate which reaches are processed
   integer(i4b)                                   :: LAKEFLAG=0           ! >0 if processing lakes
   integer(i4b)                                   :: nOrder               ! number of stream order
   integer(i4b)                                   :: nTrib                ! number of tributary basins
   integer(i4b)                                   :: nSeg                 ! number of reaches in the network
   integer(i4b)                                   :: iSeg, jSeg           ! loop indices - reach
   integer(i4b)                                   :: iTrib                ! loop indices - branch
   integer(i4b)                                   :: ix                   ! loop indices stream order
   ! variables needed for timing
   integer*8                                      :: cr                   ! rate
   integer*8                                      :: startTime,endTime    ! date/time for the start and end of the initialization
   real(dp)                                       :: elapsedTime          ! elapsed time for the process

   ! initialize error control
   ierr=0; message='kwe_route/'
   call system_clock(count_rate=cr)

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

   call system_clock(startTime)

   do ix = 1, nOrder

     nTrib=size(river_basin(ix)%branch)

     ! 1. Route tributary reaches (parallel)
!$OMP PARALLEL DO schedule(dynamic,1)                   & ! chunk size of 1
!$OMP          private(jSeg, iSeg)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(T0,T1)                            & ! private for a given thread
!$OMP          shared(LAKEFLAG)                         & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(doRoute)                          & ! data array shared
!$OMP          shared(NETOPO_in)                        & ! data structure shared
!$OMP          shared(RPARAM_in)                        & ! data structure shared
!$OMP          shared(RCHSTA_out)                       & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(ix, iEns, ixDesire)               & ! indices shared
!$OMP          firstprivate(nTrib)
     trib:do iTrib = 1,nTrib
       seg:do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
         jSeg  = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
         if (.not. doRoute(jSeg)) cycle
         ! route kinematic waves through the river network
         call kwe_rch(iEns,jSeg,           & ! input: array indices
                      ixDesire,            & ! input: index of the desired reach
                      T0,T1,               & ! input: start and end of the time step
                      LAKEFLAG,            & ! input: flag if lakes are to be processed
                      NETOPO_in,           & ! input: reach topology data structure
                      RPARAM_in,           & ! input: reach parameter data structure
                      RCHSTA_out,          & ! inout: reach state data structure
                      RCHFLX_out,          & ! inout: reach flux data structure
                      ierr,cmessage)         ! output: error control
         !if (ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
       end do  seg
     end do trib
!$OMP END PARALLEL DO

   end do

   call system_clock(endTime)
   elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
!   write(*,"(A,1PG15.7,A)") '  elapsed-time [routing/kwe] = ', elapsedTime, ' s'

 END SUBROUTINE kwe_route

 ! *********************************************************************
 ! subroutine: perform one segment route KW routing
 ! *********************************************************************
 SUBROUTINE kwe_rch(iEns, segIndex, & ! input: index of runoff ensemble to be processed
                    ixDesire,       & ! input: reachID to be checked by on-screen pringing
                    T0,T1,          & ! input: start and end of the time step
                    LAKEFLAG,       & ! input: flag if lakes are to be processed
                    NETOPO_in,      & ! input: reach topology data structure
                    RPARAM_in,      & ! input: reach parameter data structure
                    RCHSTA_out,     & ! inout: reach state data structure
                    RCHFLX_out,     & ! inout: reach flux data structure
                    ierr, message)    ! output: error control

 implicit none

 ! Input
 integer(i4b),  intent(in)                 :: iEns              ! runoff ensemble to be routed
 integer(i4b),  intent(in)                 :: segIndex          ! segment where routing is performed
 integer(i4b),  intent(in)                 :: ixDesire          ! index of the reach for verbose output
 real(dp),      intent(in)                 :: T0,T1             ! start and end of the time step (seconds)
 integer(i4b),  intent(in)                 :: LAKEFLAG          ! >0 if processing lakes
 type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
 type(RCHPRP),  intent(in),    allocatable :: RPARAM_in(:)      ! River reach parameter
 ! inout
 type(STRSTA),  intent(inout), allocatable :: RCHSTA_out(:,:)   ! reach state data
 type(STRFLX),  intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! Output
 integer(i4b),  intent(out)                :: ierr              ! error code
 character(*),  intent(out)                :: message           ! error message
 ! Local variables to
 type(STRFLX), allocatable                 :: uprflux(:)        ! upstream Reach fluxes
 integer(i4b)                              :: nUps              ! number of upstream segment
 integer(i4b)                              :: iUps              ! upstream reach index
 integer(i4b)                              :: iRch_ups          ! index of upstream reach in NETOPO
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ! initialize error control
 ierr=0; message='kwe_rch/'


 ! identify number of upstream segments of the reach being processed
 nUps = size(NETOPO_in(segIndex)%UREACHI)

 allocate(uprflux(nUps), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': uprflux'; return; endif

 if (nUps>0) then
   do iUps = 1,nUps
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     uprflux(iUps) = RCHFLX_out(iens,iRch_ups)
   end do
 endif

 ! perform river network UH routing
 call kw_euler(RPARAM_in(segIndex),         &    ! input: parameter at segIndex reach
               T0,T1,                       &    ! input: start and end of the time step
               uprflux,                     &    ! input: upstream reach fluxes
               RCHSTA_out(iens,segIndex),   &    ! inout:
               RCHFLX_out(iens,segIndex),   &    ! inout: updated fluxes at reach
               ierr, message)                    ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! Check True since now this reach now routed
 RCHFLX_out(iEns,segIndex)%isRoute=.True.

 ! check
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
  write(iulog,*) 'RCHFLX_out(iens,segIndex)%BASIN_QR(1),RCHFLX_out(iens,segIndex)%REACH_Q = ', &
                  RCHFLX_out(iens,segIndex)%BASIN_QR(1),RCHFLX_out(iens,segIndex)%REACH_Q
 endif

 END SUBROUTINE kwe_rch


 ! *********************************************************************
 ! subroutine: route kinematic waves at one segment
 ! *********************************************************************
 SUBROUTINE kw_euler(rch_param,     & ! input: river parameter data structure
                     T0,T1,         & ! input: start and end of the time step
                     uprflux,       & ! input: upstream reach fluxes
                     rstate,        & ! inout: reach state at a reach
                     rflux,         & ! inout: reach flux at a reach
                     ierr,message)
 ! ----------------------------------------------------------------------------------------

 ! https://www.fs.fed.us/rm/pubs_other/rmrs_2010_wang_l001.pdf
 ! STILL WORKING IN PROGRESS.

 IMPLICIT NONE
 ! Input
 type(RCHPRP), intent(in)                 :: rch_param    ! River reach parameter
 real(dp),     intent(in)                 :: T0,T1        ! start and end of the time step (seconds)
 type(STRFLX), intent(in), allocatable    :: uprflux(:)   ! upstream Reach fluxes
 ! Input/Output
 type(STRSTA), intent(inout)              :: rstate       ! curent reach states
 type(STRFLX), intent(inout)              :: rflux        ! current Reach fluxes
 ! Output
 integer(i4b), intent(out)                :: ierr         ! error code
 character(*), intent(out)                :: message      ! error message
 ! LOCAL VAIRABLES
 type(EKWRCH)                             :: ekw_tmp      ! temporal eulerian kw state
 real(DP)                                 :: beta         ! constant, 5/3
 real(DP)                                 :: alpha        ! sqrt(slope)(/mannings N* width)
 real(dp)                                 :: theta        ! Courant number [-]
 real(dp)                                 :: dT           ! interval of time step [sec]
 real(dp)                                 :: dX           ! length of segment [m]
 real(dp)                                 :: Abar         !
 real(dp)                                 :: Qbar         !
 integer(i4b)                             :: iUps,ix      ! loop index
 integer(i4b)                             :: nUps         ! number of all upstream segment
 real(dp), allocatable                    :: WC(:)        ! wave celerity
 character(len=strLen)                    :: cmessage     ! error message from subroutine

 ! initialize error control
 ierr=0; message='kw_euler/'

 nUps = size(uprflux)

 if (nUps>0) then

   ! (time:0:1, loc:0:1) 0-previous time step/inlet, 1-current time step/outlet.
   do ix = 0,1
     ekw_tmp%A(0,ix) = rstate%EKW_ROUTE%A(1,ix)
     ekw_tmp%Q(0,ix) = rstate%EKW_ROUTE%Q(1,ix)
   end do

   ekw_tmp%A(1,1) = 0._dp
   ekw_tmp%Q(1,1) = 0._dp

   ! Get the reach parameters
   ! A = (Q/alpha)**(1/beta)
   ! Q = alpha*A**beta
   alpha  = sqrt(rch_param%R_SLOPE)/(rch_param%R_MAN_N*rch_param%R_WIDTH**(2._dp/3._dp))
   beta = 5._dp/3._dp
   dX  = rch_param%RLENGTH
   dT = T1-T0

   ! compute flow rate and flow area at upstream end at current time step
   ekw_tmp%Q(1,0) = 0.0_dp
   do iUps = 1,nUps
     ekw_tmp%Q(1,0) = ekw_tmp%Q(1,0) + uprflux(iUps)%REACH_Q
   end do

   ekw_tmp%A(1,0) = (ekw_tmp%Q(1,0)/alpha)**(1/beta)

   ! Compute Courant number == ration of celarity to dx/dt
   Qbar = 0.5*(ekw_tmp%Q(0,1)+ekw_tmp%Q(1,0))
   Abar = 0.5*(ekw_tmp%A(0,1)+ekw_tmp%A(1,0))
   if (Abar < 0._dp .and. Abar >0._dp) then
     Qbar = runoffMin
     Abar = (Qbar/alpha)**(1/beta)
   endif
   theta = beta*(dT/dX)*(Qbar)/(Abar)

   ! ----------
   ! solve flow rate and flow area at downstream end at current time step
   ! ----------
   if (theta >= 1) then
     ekw_tmp%Q(1,1) = ekw_tmp%Q(1,0) - dX/dT*(ekw_tmp%A(1,0)-ekw_tmp%A(0,0)) + rflux%BASIN_QR(1)
     ekw_tmp%A(1,1) = (ekw_tmp%Q(1,1)/alpha)**(1/beta)
   else
     ekw_tmp%A(1,1) = ekw_tmp%A(1,0) + dT/dX*(ekw_tmp%Q(0,0) - ekw_tmp%Q(0,1))
     ekw_tmp%Q(1,1) = alpha*ekw_tmp%A(1,1)**beta + rflux%BASIN_QR(1)
     ekw_tmp%A(1,1) = (ekw_tmp%Q(1,1)/alpha)**(1/beta)
   end if

 else ! if head-water

   do ix = 0,1
     ekw_tmp%A(0,ix) = rstate%EKW_ROUTE%A(1,ix)
     ekw_tmp%Q(0,ix) = rstate%EKW_ROUTE%Q(1,ix)
   end do

   ekw_tmp%A(1,0) = 0._dp
   ekw_tmp%Q(1,0) = 0._dp

   ekw_tmp%Q(1,1) = rflux%BASIN_QR(1)
   ekw_tmp%A(1,1) = (ekw_tmp%Q(1,1)/alpha)**(1/beta)
 endif

 ! add catchment flow
 rflux%REACH_Q = ekw_tmp%Q(1,1)

 ! update state
 rstate%EKW_ROUTE%Q = ekw_tmp%Q
 rstate%EKW_ROUTE%A = ekw_tmp%A

 END SUBROUTINE kw_euler


END MODULE kwe_route_module
