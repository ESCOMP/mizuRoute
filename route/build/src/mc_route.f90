MODULE mc_route_module

! muskingum-cunge routing

USE nrtype
! data types
USE dataTypes,   ONLY: STRFLX            ! fluxes in each reach
USE dataTypes,   ONLY: STRSTA            ! state in each reach
USE dataTypes,   ONLY: RCHTOPO           ! Network topology
USE dataTypes,   ONLY: RCHPRP            ! Reach parameter
USE dataTypes,   ONLY: mcRCH             ! MC specific state data structure 
USE dataTypes,   ONLY: subbasin_omp      ! mainstem+tributary data strucuture
! global data
USE public_var,  ONLY: iulog             ! i/o logical unit number
USE public_var,  ONLY: realMissing       ! missing value for real number
USE public_var,  ONLY: integerMissing    ! missing value for integer number
USE globalData,  ONLY: idxMC             ! index of IRF method 
! subroutines: general
USE model_finalize, ONLY : handle_err

! privary
implicit none
private

public::mc_route

contains

 ! *********************************************************************
 ! subroutine: perform muskingum-cunge routing through the river network
 ! *********************************************************************
 SUBROUTINE mc_route(iens,                 & ! input: ensemble index
                     river_basin,          & ! input: river basin information (mainstem, tributary outlet etc.)
                     T0,T1,                & ! input: start and end of the time step
                     ixDesire,             & ! input: reachID to be checked by on-screen pringing
                     NETOPO_in,            & ! input: reach topology data structure
                     RPARAM_in,            & ! input: reach parameter data structure
                     RCHSTA_out,           & ! inout: reach state data structure
                     RCHFLX_out,           & ! inout: reach flux data structure
                     ierr,message,         & ! output: error control
                     ixSubRch)               ! optional input: subset of reach indices to be processed

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

   ierr=0; message='mc_route/'

   ! number of reach check
   if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
    ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
   endif

   nSeg = size(NETOPO_in)

   allocate(doRoute(nSeg), stat=ierr)

   if (present(ixSubRch))then
    doRoute(:)=.false.
    doRoute(ixSubRch) = .true. ! only subset of reaches are on
   else
    doRoute(:)=.true. ! every reach is on
   endif

   nOrder = size(river_basin)

   do ix = 1, nOrder

     nTrib=size(river_basin(ix)%branch)

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
         call mc_rch(iEns,jSeg,           & ! input: array indices
                     ixDesire,            & ! input: index of the desired reach
                     T0,T1,               & ! input: start and end of the time step
                     LAKEFLAG,            & ! input: flag if lakes are to be processed
                     NETOPO_in,           & ! input: reach topology data structure
                     RPARAM_in,           & ! input: reach parameter data structure
                     RCHSTA_out,          & ! inout: reach state data structure
                     RCHFLX_out,          & ! inout: reach flux data structure
                     ierr,cmessage)         ! output: error control
         if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
       end do  seg
     end do trib
!$OMP END PARALLEL DO

   end do

 END SUBROUTINE mc_route

 ! *********************************************************************
 ! subroutine: perform muskingum-cunge routing for one segment
 ! *********************************************************************
 SUBROUTINE mc_rch(iEns, segIndex, & ! input: index of runoff ensemble to be processed
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
 logical(lgt)                              :: doCheck           ! check details of variables
 logical(lgt)                              :: isHW              ! headwater basin?
 integer(i4b)                              :: nUps              ! number of upstream segment
 integer(i4b)                              :: iUps              ! upstream reach index
 integer(i4b)                              :: iRch_ups          ! index of upstream reach in NETOPO
 real(dp)                                  :: q_upstream        ! total discharge at top of the reach being processed
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ierr=0; message='mc_rch/'

 doCheck = .false.
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   doCheck = .true.
 end if

 ! get discharge coming from upstream
 nUps = size(NETOPO_in(segIndex)%UREACHI)
 isHW = .true.
 q_upstream = 0.0_dp
 if (nUps>0) then
   isHW = .false.
   do iUps = 1,nUps
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxMC)%REACH_Q
   end do
 endif

 if(doCheck)then
   write(iulog,'(A)') 'CHECK muskingum-cunge routing'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,X,I6,X,G12.5)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps),RCHFLX_out(iens, iRch_ups)%ROUTE(idxMC)%REACH_Q
     enddo
   end if
   write(iulog,'(A,X,G12.5)') ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1)=',RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
 endif

 ! solve muskingum-cunge alogorithm
 call muskingum_cunge(RPARAM_in(segIndex),                & ! input: parameter at segIndex reach
                      T0,T1,                              & ! input: start and end of the time step
                      q_upstream,                         & ! input: total discharge at top of the reach being processed
                      isHW,                               & ! input: is this headwater basin?
                      RCHSTA_out(iens,segIndex)%MC_ROUTE, & ! inout:
                      RCHFLX_out(iens,segIndex),          & ! inout: updated fluxes at reach
                      doCheck,                            & ! input: reach index to be examined
                      ierr, cmessage)                       ! output: error control
 if(ierr/=0)then
    write(message, '(A,X,I10,X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage)
    return
 endif

 if(doCheck)then
  write(iulog,'(A,X,G12.5)') ' RCHFLX_out(iens,segIndex)%REACH_Q=', RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_Q
 endif

 END SUBROUTINE mc_rch


 ! *********************************************************************
 ! subroutine: solve muskingum equation
 ! *********************************************************************
 SUBROUTINE muskingum_cunge(rch_param,     & ! input: river parameter data structure
                            T0,T1,         & ! input: start and end of the time step
                            q_upstream,    & ! input: discharge from upstream
                            isHW,          & ! input: is this headwater basin?
                            rstate,        & ! inout: reach state at a reach
                            rflux,         & ! inout: reach flux at a reach
                            doCheck,       & ! input: reach index to be examined
                            ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Perform muskingum-cunge routing
 !
 !
 !
 !
 ! *
 ! *
 ! *
 !
 ! state array:
 ! (time:0:1, loc:0:1) 0-previous time step/inlet, 1-current time step/outlet.
 ! Q or A(1,2,3,4): 1: (t=0,x=0), 2: (t=0,x=1), 3: (t=1,x=0), 4: (t=1,x=1)

 IMPLICIT NONE
 ! Input
 type(RCHPRP), intent(in)                 :: rch_param    ! River reach parameter
 real(dp),     intent(in)                 :: T0,T1        ! start and end of the time step (seconds)
 real(dp),     intent(in)                 :: q_upstream   ! total discharge at top of the reach being processed
 logical(lgt), intent(in)                 :: isHW         ! is this headwater basin?
 logical(lgt), intent(in)                 :: doCheck      ! reach index to be examined
 ! Input/Output
 type(mcRCH),  intent(inout)              :: rstate       ! curent reach states
 type(STRFLX), intent(inout)              :: rflux        ! current Reach fluxes
 ! Output
 integer(i4b), intent(out)                :: ierr         ! error code
 character(*), intent(out)                :: message      ! error message
 ! LOCAL VAIRABLES
 real(dp)                                 :: alpha        ! sqrt(slope)(/mannings N* width)
 real(dp)                                 :: beta         ! constant, 5/3
 real(dp)                                 :: theta        ! dT/dX
 real(dp)                                 :: X            ! X-factor in descreterized kinematic wave
 real(dp)                                 :: dt           ! interval of time step [sec]
 real(dp)                                 :: dx           ! length of segment [m]
 real(dp)                                 :: Q(0:1,0:1)   ! discharge at computational molecule
 real(dp)                                 :: Qbar         ! 3-point average discharge [m3/s]
 real(dp)                                 :: Abar         ! 3-point average flow area [m2]
 real(dp)                                 :: Vbar         ! 3-point average velocity [m/s]
 real(dp)                                 :: Ybar         ! 3-point average flow depth [m]
 real(dp)                                 :: B            ! flow top width [m]
 real(dp)                                 :: ck           ! kinematic wave celerity [m/s]
 real(dp)                                 :: Cn           ! Courant number [-]
 real(dp)                                 :: dTsub        ! time inteval for sut time-step [sec]
 real(dp)                                 :: C0,C1,C2     ! muskingum parameters
 real(dp), allocatable                    :: QoutLocal(:) ! out discharge [m3/s] at sub time step
 real(dp), allocatable                    :: QinLocal(:)  ! in discharge [m3/s] at sub time step
 integer(i4b)                             :: ix           ! loop index
 integer(i4b)                             :: ntSub        ! number of sub time-step
 character(len=strLen)                    :: cmessage     ! error message from subroutine
 real(dp), parameter                      :: Y = 0.5      ! muskingum parameter Y (this is fixed)

 ierr=0; message='muskingum-cunge/'

 Q(0,0) = rstate%molecule%Q(1) ! inflow at previous time step (t-1)
 Q(0,1) = rstate%molecule%Q(2) ! outflow at previous time step (t-1)
 Q(1,1) = realMissing

 if (.not. isHW) then

   ! Get the reach parameters
   ! A = (Q/alpha)**(1/beta)
   ! Q = alpha*A**beta
   alpha = sqrt(rch_param%R_SLOPE)/(rch_param%R_MAN_N*rch_param%R_WIDTH**(2._dp/3._dp))
   beta  = 5._dp/3._dp
   dx = rch_param%RLENGTH
   dt = T1-T0
   theta = dt/dx     ! [s/m]

   ! compute total flow rate and flow area at upstream end at current time step
   Q(1,0) = q_upstream

   if (doCheck) then
     write(iulog,'(4(A,X,G12.5))') ' length [m] =',rch_param%RLENGTH,'slope [-] =',rch_param%R_SLOPE,'channel width [m] =',rch_param%R_WIDTH,'manning coef =',rch_param%R_MAN_N
     write(iulog,'(3(A,X,G12.5))') ' Qin(t-1) Q(0,0)=',Q(0,0),' Qin(t) Q(0,1)=',Q(0,1),' Qout(t-1) Q(1,0)=',Q(1,0)
   end if

   ! first, using 3-point average in computational molecule, check Cournat number is less than 1, otherwise subcycle within one time step
   Qbar = (Q(0,0)+Q(1,0)+Q(0,1))/3.0  ! average discharge [m3/s]
   Abar = (Qbar/alpha)**(1/beta)      ! average flow area [m2] (from manning equation)
   Vbar = Qbar/Abar                   ! average velocity [m/s]
   ck   = beta*Vbar                   ! kinematic wave celerity [m/s]
   Cn   = ck*theta                    ! Courant number [-]

   ! time-step adjustment so Courant number is less than 1
   ntSub = 1
   dTsub = dt
   if (Cn>1.0_dp) then
     ntSub = ceiling(dt/dx*cK)
     dTsub = dt/ntSub
   end if
   if (doCheck) then
     write(iulog,'(A,X,I3,A,X,G12.5)') ' No. sub timestep=',nTsub,' sub time-step [sec]=',dTsub
   end if

   allocate(QoutLocal(0:ntSub), QinLocal(0:ntSub), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   QoutLocal(:) = realMissing
   QoutLocal(0) = Q(0,1)        ! outfloe last time step
   QinLocal(0)  = Q(0,0)        ! inflow at last time step
   QinLocal(1:ntSub)  = Q(1,0)  ! infllow at sub-time step in current time step

   ! solve outflow at each sub time step
   do ix = 1, nTsub

     Qbar = (QinLocal(ix)+QinLocal(ix-1)+QoutLocal(ix-1))/3.0 ! 3 point average discharge [m3/s]
     Abar = (Qbar/alpha)**(1/beta)                            ! flow area [m2] (manning equation)
     Ybar = Abar/rch_param%R_WIDTH                            ! flow depth [m] (rectangular channel)
     B = rch_param%R_WIDTH                                    ! top width at water level [m] (rectangular channel)
     ck = beta*(Qbar/Abar)                                    ! kinematic wave celerity [m/s]

     X = 0.5*(1.0 - Qbar/(B*rch_param%R_SLOPE*ck*dX))         ! X factor for descreterized kinematic wave equation
     Cn = ck*dTsub/dx                                         ! Courant number [-]

     C0 = (-X+Cn*(1-Y))/(1-X+Cn*(1-Y))
     C1 = (X+Cn*Y)/(1-X+Cn*(1-Y))
     C2 = (1-X-Cn*Y)/(1-X+Cn*(1-Y))

     QoutLocal(ix) = C0* QinLocal(ix)+ C1* QinLocal(ix-1)+ C2* QoutLocal(ix-1)
     QoutLocal(ix) = max(0.0, QoutLocal(ix))

     if (isnan(QoutLocal(ix))) then
       ierr=10; message=trim(message)//'QoutLocal is Nan; activate vodose for this segment for diagnosis';return
     end if

     if (doCheck) then
       write(iulog,'(A,I3,X,A,G12.5,X,A,G12.5)') '   sub time-step= ',ix,'Courant number= ',Cn, 'Q= ',QoutLocal(ix)
     end if
   end do

   Q(1,1) = sum(QoutLocal(1:nTsub))/real(nTsub,kind=dp)

 else ! if head-water

   Q(1,0) = 0._dp
   Q(1,1) = 0._dp

   if (doCheck) then
     write(iulog,'(A)')            ' This is headwater '
   endif

 endif

 ! compute volume
 rflux%ROUTE(idxMC)%REACH_VOL(0) = rflux%ROUTE(idxMC)%REACH_VOL(1)
 rflux%ROUTE(idxMC)%REACH_VOL(1) = rflux%ROUTE(idxMC)%REACH_VOL(0) + (Q(1,0)-Q(1,1))*dT
 rflux%ROUTE(idxMC)%REACH_VOL(1) = max(rflux%ROUTE(idxMC)%REACH_VOL(1), 0._dp)

 ! add catchment flow
 rflux%ROUTE(idxMC)%REACH_Q = Q(1,1)+rflux%BASIN_QR(1)

 if (doCheck) then
   write(iulog,'(A,X,G12.5)') ' Qout(t)=',Q(1,1)
 endif

 ! save inflow (index 1) and outflow (index 2) at current time step
 rstate%molecule%Q(1) = Q(1,0)
 rstate%molecule%Q(2) = Q(1,1)

 END SUBROUTINE muskingum_cunge


END MODULE mc_route_module
