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
 ! Need to handle empty chanel (at first time step and the headwate basin where no channel routing is performed)

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
 integer(i4b)                             :: NN           ! number of reach sub-segments
 real(DP)                                 :: ALFA         ! constant, 5/3
 real(DP)                                 :: K            ! sqrt(slope)/mannings N
 real(DP)                                 :: XMX          ! length of the stream segment
 real(dp)                                 :: dT           ! interval of time step [sec]
 real(dp)                                 :: dx           ! time duration [sec]
 real(dp)                                 :: q_upstream   ! total upstream discharge [m3/s]
 integer(i4b)                             :: iUps,ix      ! loop index
 integer(i4b)                             :: nUps         ! number of all upstream segment
 real(dp), allocatable                    :: qin(:)       ! interval of time step [sec]
 real(dp), allocatable                    :: WC(:)        ! wave celerity
 character(len=strLen)                    :: cmessage     ! error message from subroutine

 ! initialize error control
 ierr=0; message='kw_euler/'


 ! deallocate current reach sub q array
 if (allocated(rstate%EKW_ROUTE%qsub)) then

  NN = size(rstate%EKW_ROUTE%qsub,2)
  allocate(ekw_tmp%qsub(0:1,NN), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': '; return; endif
  ekw_tmp%qsub = rstate%EKW_ROUTE%qsub

  deallocate(rstate%EKW_ROUTE%qsub)

 else

  NN = 1
  allocate(ekw_tmp%qsub(0:1,NN))
  ekw_tmp%qsub(:,:) = 0._dp

 endif

 ! identify number of upstream segments of the reach being processed
 nUps = size(uprflux)

 ! Get the reach parameters
 ALFA = 5._dp/3._dp        ! should this be initialized here or in a parameter file?
 K    = sqrt(rch_param%R_SLOPE)/rch_param%R_MAN_N
 XMX  = rch_param%RLENGTH

 ! Identify the number of points to route
 dT = T1-T0
 dx = XMX/NN

 allocate(WC(NN), qin(NN), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': '; return; endif

 ! compute wave celerity for all flow points (array operation)
 WC(1:NN) = ALFA*K**(1._dp/ALFA)*ekw_tmp%qsub(0,1:NN)**((ALFA-1._dp)/ALFA)

 qin(1:NN) = 0._dp
 qin(NN)   = rflux%BASIN_QR(1)

 ! Find out total q at top of a segment
 q_upstream = 0.0_dp
 if(nUps>0)then
   do iUps = 1,nUps
     q_upstream = q_upstream + uprflux(iUps)%REACH_Q
   end do
 endif

 subrch:do ix = 1, NN
   if (ix ==1) then
     ekw_tmp%qsub(1,ix) = dT/dX*q_upstream + 1._dp/WC(ix)*ekw_tmp%qsub(0,ix) + dT*qin(ix)
   else
     ekw_tmp%qsub(1,ix) = dT/dX*ekw_tmp%qsub(1,ix-1)+1._dp/WC(ix)*ekw_tmp%qsub(0,ix) + dT*qin(ix)
   endif
   ekw_tmp%qsub(1,ix) = ekw_tmp%qsub(1,ix)/(dT/dX+1._dp/WC(ix))
 end do subrch

 rflux%REACH_Q = ekw_tmp%qsub(1,NN)


 allocate(rstate%EKW_ROUTE%qsub(0:1,NN), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': rstate%EKW_ROUTE%qsub'; return; endif

 rstate%EKW_ROUTE%qsub = ekw_tmp%qsub

 END SUBROUTINE kw_euler


END MODULE kwe_route_module
