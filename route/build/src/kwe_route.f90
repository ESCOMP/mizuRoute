MODULE kwe_route_module

USE nrtype
! data types
USE dataTypes, ONLY: STRFLX            ! fluxes in each reach
USE dataTypes, ONLY: STRSTA            ! state in each reach
USE dataTypes, ONLY: RCHTOPO           ! Network topology
USE dataTypes, ONLY: RCHPRP            ! Reach parameter
USE dataTypes, ONLY: EKWRCH            ! Eulerian KW state
! global data
USE public_var, ONLY: iulog             ! i/o logical unit number
USE public_var, ONLY: realMissing       ! missing value for real number
USE public_var, ONLY: integerMissing    ! missing value for integer number
! subroutines: general
USE perf_mod,  ONLY: t_startf,t_stopf   ! timing start/stop

! privary
implicit none
private

public::kwe_route

real(dp), parameter  :: critFactor=0.01

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
   USE model_utils, ONLY: handle_err

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

   ! initialize error control
   ierr=0; message='kwe_route/'

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

   call t_startf('route/kwe')

   ! route kinematic waves through the river network
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
         call kwe_rch(iEns,jSeg,           & ! input: array indices
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

   call t_stopf('route/kwe')

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
 logical(lgt)                              :: doCheck           ! check details of variables
 logical(lgt)                              :: isHW              ! headwater basin?
 integer(i4b)                              :: nUps              ! number of upstream segment
 integer(i4b)                              :: iUps              ! upstream reach index
 integer(i4b)                              :: iRch_ups          ! index of upstream reach in NETOPO
 real(dp)                                  :: q_upstream        ! total discharge at top of the reach being processed
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ! initialize error control
 ierr=0; message='kwe_rch/'

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
     q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%REACH_Q
   end do
 endif

 if(doCheck)then
   write(iulog,'(A)') 'CHECK Euler Kinematic wave'
   if (nUps>0) then
     do iUps = 1,nUps
       write(iulog,'(A,X,I6,X,G12.5)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps),RCHFLX_out(iens, NETOPO_in(segIndex)%UREACHK(iUps))%REACH_Q
     enddo
   end if
   write(iulog,'(A,X,G12.5)') ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1)=',RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
 endif

 ! perform river network KWE routing
 call kw_euler(RPARAM_in(segIndex),         &    ! input: parameter at segIndex reach
               T0,T1,                       &    ! input: start and end of the time step
               q_upstream,                  &    ! input: total discharge at top of the reach being processed
               isHW,                        &    ! input: is this headwater basin?
               RCHSTA_out(iens,segIndex),   &    ! inout:
               RCHFLX_out(iens,segIndex),   &    ! inout: updated fluxes at reach
               doCheck,                     &    ! input: reach index to be examined
               ierr, message)                    ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! Check True since now this reach now routed
 RCHFLX_out(iEns,segIndex)%isRoute=.True.

 if(doCheck)then
  write(iulog,'(A,X,G12.5)') ' RCHFLX_out(iens,segIndex)%REACH_Q=', RCHFLX_out(iens,segIndex)%REACH_Q
 endif

 END SUBROUTINE kwe_rch


 ! *********************************************************************
 ! subroutine: route kinematic waves at one segment
 ! *********************************************************************
 SUBROUTINE kw_euler(rch_param,     & ! input: river parameter data structure
                     T0,T1,         & ! input: start and end of the time step
                     q_upstream,    & ! input: discharge from upstream
                     isHW,          & ! input: is this headwater basin?
                     rstate,        & ! inout: reach state at a reach
                     rflux,         & ! inout: reach flux at a reach
                     doCheck,       & ! input: reach index to be examined
                     ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Kinematic wave equation is solved based on conservative form the equation
 !
 ! Method: Li, R.‐M., Simons, D. B., and Stevens, M. A. (1975), Nonlinear kinematic wave approximation for water routing,
 !         Water Resour. Res., 11( 2), 245– 252, doi:10.1029/WR011i002p00245
 !
 ! * Use analytical solution (eq 29 in paper) for Q using the 2nd order Talor series of nonlinear kinematic equation: theta*Q + alpha*Q^beta = omega
 ! * Use initial guess using explicit Euler solution of kinematic approximation equation to start iterative computation of Q
 ! * iterative Q computation till LHS ~= RHS
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
 type(STRSTA), intent(inout)              :: rstate       ! curent reach states
 type(STRFLX), intent(inout)              :: rflux        ! current Reach fluxes
 ! Output
 integer(i4b), intent(out)                :: ierr         ! error code
 character(*), intent(out)                :: message      ! error message
 ! LOCAL VAIRABLES
 real(dp)                                 :: alpha        ! sqrt(slope)(/mannings N* width)
 real(dp)                                 :: beta         ! constant, 5/3
 real(dp)                                 :: alpha1       ! sqrt(slope)(/mannings N* width)
 real(dp)                                 :: beta1        ! constant, 5/3
 real(dp)                                 :: theta        ! dT/dX
 real(dp)                                 :: omega        ! right-hand side of kw finite difference
 real(dp)                                 :: f0,f1,f2     ! values of function f, 1st and 2nd derivatives at solution
 real(dp)                                 :: X            !
 real(dp)                                 :: dT           ! interval of time step [sec]
 real(dp)                                 :: dX           ! length of segment [m]
 real(dp)                                 :: A(0:1,0:1)   !
 real(dp)                                 :: Q(0:1,0:1)   !
 real(dp)                                 :: Qtrial(2)    ! trial solution of kw equation
 real(dp)                                 :: Abar         !
 real(dp)                                 :: Qbar         !
 real(dp)                                 :: absErr(2)    ! absolute error of nonliear equation solution
 real(dp)                                 :: f0eval(2)    !
 integer(i4b)                             :: imin         ! index at minimum value
 integer(i4b)                             :: ix           ! loop index
 character(len=strLen)                    :: cmessage     ! error message from subroutine

 ierr=0; message='kw_euler/'

 !  current time and inlet  3 (1,0) -> previous time and inlet  1 (0,0)
 !  current time and outlet 4 (1,1) -> previous time and outlet 2 (0,1)
 do ix = 0,1
   A(0,ix) = rstate%EKW_ROUTE%A(ix+3)
   Q(0,ix) = rstate%EKW_ROUTE%Q(ix+3)
 end do

 if (.not. isHW) then

   A(1,1) = realMissing
   Q(1,1) = realMissing

   ! Get the reach parameters
   ! A = (Q/alpha)**(1/beta)
   ! Q = alpha*A**beta
   alpha = sqrt(rch_param%R_SLOPE)/(rch_param%R_MAN_N*rch_param%R_WIDTH**(2._dp/3._dp))
   beta  = 5._dp/3._dp
   beta1  = 1._dp/beta
   alpha1 = (1.0/alpha)**beta1
   dX = rch_param%RLENGTH
   dT = T1-T0
   theta = dT/dX

   ! compute total flow rate and flow area at upstream end at current time step
   Q(1,0) = q_upstream
   A(1,0) = (Q(1,0)/alpha)**(1/beta)

   if (doCheck) then
     write(iulog,'(3(A,X,G12.5))') ' R_SLOPE=',rch_param%R_SLOPE,' R_WIDTH=',rch_param%R_WIDTH,' R_MANN=',rch_param%R_MAN_N
     write(iulog,'(3(A,X,G12.5))') ' Q(0,0)=',Q(0,0),' Q(0,1)=',Q(0,1),' Q(1,0)=',Q(1,0)
     write(iulog,'(3(A,X,G12.5))') ' A(0,0)=',A(0,0),' A(0,1)=',A(0,1),' A(1,0)=',A(1,0)
   end if

   ! ----------
   ! solve flow rate and flow area at downstream end at current time step
   ! ----------
   ! initial guess
   Qbar   = (Q(0,1) + Q(1,0))/2._dp
   Q(1,1) = (theta*Q(1,0) + alpha1*beta1*Qbar**(beta1-1)*Q(0,1))/(theta + alpha1*beta1*Qbar**(beta1-1))

   omega = theta*Q(1,0)+alpha1*Q(0,1)**(beta1)

   f0eval(1) = theta*Q(1,1) + alpha1*Q(1,1)**beta1
   absErr(1) = abs(f0eval(1)-omega)

   if ( abs(Q(1,1)-0.0_dp) < epsilon(Q(1,1))) then
     Q(1,1) = omega/(theta+alpha1)
   else if (absErr(1) > critFactor*omega) then
     ! iterative solution
     do
       f0 = theta*Q(1,1) + alpha1*Q(1,1)**beta1
       f1 = theta + alpha1*beta1*Q(1,1)**(beta1-1)     ! 1st derivative of f w.r.t. Q
       f2 = alpha1*beta1*(beta1-1)*Q(1,1)**(beta1-2)   ! 2nd derivative of f w.r.t. Q

       X = (f1/f2)**2._dp - 2._dp*(f0-omega)/f2
       if (X<0) X=0._dp

       ! two solutions
       Qtrial(1) = abs(Q(1,1) - f1/f2 + sqrt(X))
       Qtrial(2) = abs(Q(1,1) - f1/f2 - sqrt(X))

       f0eval = theta*Qtrial + alpha1*Qtrial**beta1
       absErr = abs(f0eval-omega)
       imin   = minloc(absErr,DIM=1)
       Q(1,1) = Qtrial(imin)

       if (absErr(imin) < critFactor*omega) exit
     end do
   endif

   A(1,1) = (Q(1,1)/alpha)**(1/beta)

   if (doCheck) then
     write(iulog,'(2(A,X,G12.5))') ' A(1,1)=',A(1,1),' Q(1,1)=',Q(1,1)
   end if

 else ! if head-water

   A(1,0) = 0._dp
   Q(1,0) = 0._dp

   Q(1,1) = 0._dp
   A(1,1) = 0._dp

   if (doCheck) then
     write(iulog,'(A)')            ' This is headwater '
     write(iulog,'(2(A,X,G12.5))') ' A(1,1)=',A(1,1),' Q(1,1)=',Q(1,1)
   endif

 endif

 ! add catchment flow
 rflux%REACH_Q = Q(1,1)+rflux%BASIN_QR(1)

 ! update state
 do ix = 0,1
   rstate%EKW_ROUTE%Q(ix+1) = Q(0,ix)
   rstate%EKW_ROUTE%Q(ix+3) = Q(1,ix)

   rstate%EKW_ROUTE%A(ix+1) = A(0,ix)
   rstate%EKW_ROUTE%A(ix+3) = A(1,ix)
 end do

 END SUBROUTINE kw_euler


END MODULE kwe_route_module
