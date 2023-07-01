MODULE kw_route_module

! kinematic wave routing

USE nrtype
USE dataTypes,     ONLY: STRFLX          ! fluxes in each reach
USE dataTypes,     ONLY: STRSTA          ! state in each reach
USE dataTypes,     ONLY: RCHTOPO         ! Network topology
USE dataTypes,     ONLY: RCHPRP          ! Reach parameter
USE dataTypes,     ONLY: kwRCH           ! kw specific state data structure
USE public_var,    ONLY: iulog           ! i/o logical unit number
USE public_var,    ONLY: realMissing     ! missing value for real number
USE public_var,    ONLY: integerMissing  ! missing value for integer number
USE public_var,    ONLY: qmodOption      ! qmod option (use 1==direct insertion)
USE globalData,    ONLY: idxKW           ! routing method index for kinematic wwave
USE water_balance, ONLY: comp_reach_wb   ! compute water balance error
USE base_route,    ONLY: base_route_rch  ! base (abstract) reach routing method class

implicit none

private
public::kwe_route_rch

real(dp), parameter  :: critFactor=0.01

type, extends(base_route_rch) :: kwe_route_rch
 CONTAINS
   procedure, pass :: route => kw_rch
end type kwe_route_rch

CONTAINS

 ! *********************************************************************
 ! subroutine: perform one segment route KW routing
 ! *********************************************************************
 SUBROUTINE kw_rch(this,           &
                   iEns, segIndex, & ! input: index of runoff ensemble to be processed
                   ixDesire,       & ! input: reachID to be checked by on-screen pringing
                   T0,T1,          & ! input: start and end of the time step
                   NETOPO_in,      & ! input: reach topology data structure
                   RPARAM_in,      & ! input: reach parameter data structure
                   RCHSTA_out,     & ! inout: reach state data structure
                   RCHFLX_out,     & ! inout: reach flux data structure
                   ierr, message)    ! output: error control
 implicit none
 ! Argument variables
 class(kwe_route_rch)                      :: this
 integer(i4b),  intent(in)                 :: iEns              ! runoff ensemble to be routed
 integer(i4b),  intent(in)                 :: segIndex          ! segment where routing is performed
 integer(i4b),  intent(in)                 :: ixDesire          ! index of the reach for verbose output
 real(dp),      intent(in)                 :: T0,T1             ! start and end of the time step (seconds)
 type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
 type(RCHPRP),  intent(inout), allocatable :: RPARAM_in(:)      ! River reach parameter
 type(STRSTA),  intent(inout)              :: RCHSTA_out(:,:)   ! reach state data
 type(STRFLX),  intent(inout)              :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b),  intent(out)                :: ierr              ! error code
 character(*),  intent(out)                :: message           ! error message
 ! Local variables
 logical(lgt)                              :: verbose           ! check details of variables
 logical(lgt)                              :: isHW              ! headwater basin?
 integer(i4b)                              :: nUps              ! number of upstream segment
 integer(i4b)                              :: iUps              ! upstream reach index
 integer(i4b)                              :: iRch_ups          ! index of upstream reach in NETOPO
 real(dp)                                  :: q_upstream        ! total discharge at top of the reach being processed
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ierr=0; message='kw_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   verbose = .true.
 end if

 ! get discharge coming from upstream
 nUps = size(NETOPO_in(segIndex)%UREACHI)
 isHW = .true.
 q_upstream = 0.0_dp
 if (nUps>0) then
   isHW = .false.
   do iUps = 1,nUps
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     if (qmodOption==1 .and. RCHFLX_out(iens,iRch_ups)%Qobs/=realMissing) then
       RCHFLX_out(iens, iRch_ups)%ROUTE(idxKW)%REACH_Q = RCHFLX_out(iens,iRch_ups)%Qobs
     end if
     q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxKW)%REACH_Q
   end do
 endif

 if(verbose)then
   write(iulog,'(2A)') new_line('a'), '** CHECK Kinematic wave routing **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,X,I12,X,G12.5)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps),RCHFLX_out(iens, iRch_ups)%ROUTE(idxKW)%REACH_Q
     enddo
   end if
   write(iulog,'(A,X,G15.4)') ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1)=',RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
 endif

 ! perform river network KW routing
 call kinematic_wave(RPARAM_in(segIndex),                     & ! input: parameter at segIndex reach
                     T0,T1,                                   & ! input: start and end of the time step
                     q_upstream,                              & ! input: total discharge at top of the reach being processed
                     RCHFLX_out(iens,segIndex)%REACH_WM_FLUX, & ! input: abstraction(-)/injection(+) [m3/s]
                     RPARAM_in(segIndex)%MINFLOW,             & ! input: minimum environmental flow [m3/s]
                     isHW,                                    & ! input: is this headwater basin?
                     RCHSTA_out(iens,segIndex)%KW_ROUTE,      & ! inout:
                     RCHFLX_out(iens,segIndex),               & ! inout: updated fluxes at reach
                     verbose,                                 & ! input: reach index to be examined
                     ierr, cmessage)                            ! output: error control
 if(ierr/=0)then
   write(message, '(A,X,I12,X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage); return
 endif

 if(verbose)then
   write(iulog,'(A,X,G15.4)') ' RCHFLX_out(iens,segIndex)%REACH_Q=', RCHFLX_out(iens,segIndex)%ROUTE(idxKW)%REACH_Q
 endif

 if (RCHFLX_out(iens,segIndex)%ROUTE(idxKW)%REACH_VOL(1) < 0) then
   write(iulog,'(A,X,G12.5,X,A,X,I9)') ' ---- NEGATIVE VOLUME = ', RCHFLX_out(iens,segIndex)%ROUTE(idxKW)%REACH_VOL(1), 'at ', NETOPO_in(segIndex)%REACHID
 end if

 call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxKW, q_upstream, RCHFLX_out(iens,segIndex), verbose, lakeFlag=.false.)

 END SUBROUTINE kw_rch


 ! *********************************************************************
 ! subroutine: route kinematic waves at one segment
 ! *********************************************************************
 SUBROUTINE kinematic_wave(rch_param,     & ! input: river parameter data structure
                           T0,T1,         & ! input: start and end of the time step
                           q_upstream,    & ! input: discharge from upstream
                           Qtake,         & ! input: abstraction(-)/injection(+) [m3/s]
                           Qmin,          & ! input: minimum environmental flow [m3/s]
                           isHW,          & ! input: is this headwater basin?
                           rstate,        & ! inout: reach state at a reach
                           rflux,         & ! inout: reach flux at a reach
                           verbose,       & ! input: reach index to be examined
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

 implicit none
 ! argument variables
 type(RCHPRP), intent(in)                 :: rch_param    ! River reach parameter
 real(dp),     intent(in)                 :: T0,T1        ! start and end of the time step (seconds)
 real(dp),     intent(in)                 :: q_upstream   ! total discharge at top of the reach being processed
 real(dp),     intent(in)                 :: Qtake        ! abstraction(-)/injection(+) [m3/s]
 real(dp),     intent(in)                 :: Qmin         ! minimum environmental flow [m3/s]
 logical(lgt), intent(in)                 :: isHW         ! is this headwater basin?
 type(kwRCH),  intent(inout)              :: rstate       ! curent reach states
 type(STRFLX), intent(inout)              :: rflux        ! current Reach fluxes
 logical(lgt), intent(in)                 :: verbose      ! reach index to be examined
 integer(i4b), intent(out)                :: ierr         ! error code
 character(*), intent(out)                :: message      ! error message
 ! Local variables
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
 real(dp)                                 :: Q(0:1,0:1)   !
 real(dp)                                 :: Qtrial(2)    ! trial solution of kw equation
 real(dp)                                 :: Qbar         !
 real(dp)                                 :: absErr(2)    ! absolute error of nonliear equation solution
 real(dp)                                 :: f0eval(2)    !
 real(dp)                                 :: QupMod       ! modified total discharge at top of the reach being processed
 real(dp)                                 :: Qabs         ! maximum allowable water abstraction rate [m3/s]
 real(dp)                                 :: Qmod         ! abstraction rate to be taken from outlet discharge [m3/s]
 integer(i4b)                             :: imin         ! index at minimum value

 ierr=0; message='kinematic_wave/'

 Q(0,0) = rstate%molecule%Q(1) ! previous time and inlet  1 (0,0)
 Q(0,1) = rstate%molecule%Q(2) ! previous time and outlet 2 (0,1)
 dt = T1-T0

 ! Q injection, add at top of reach
 QupMod = q_upstream
 if (Qtake>0) then
   QupMod = QupMod+ Qtake
 end if

 if (.not. isHW) then

   Q(1,1) = realMissing ! current time and outlet 4 (1,1)

   ! Get the reach parameters
   ! A = (Q/alpha)**(1/beta)
   ! Q = alpha*A**beta
   alpha = sqrt(rch_param%R_SLOPE)/(rch_param%R_MAN_N*rch_param%R_WIDTH**(2._dp/3._dp))
   beta  = 5._dp/3._dp
   beta1  = 1._dp/beta
   alpha1 = (1.0/alpha)**beta1
   dX = rch_param%RLENGTH
   theta = dt/dX

   ! compute total flow rate and flow area at upstream end at current time step
   Q(1,0) = QupMod

   if (verbose) then
     write(iulog,'(A,X,G12.5)') ' length [m]        =',rch_param%RLENGTH
     write(iulog,'(A,X,G12.5)') ' slope [-]         =',rch_param%R_SLOPE
     write(iulog,'(A,X,G12.5)') ' channel width [m] =',rch_param%R_WIDTH
     write(iulog,'(A,X,G12.5)') ' manning coef. [-] =',rch_param%R_MAN_N
     write(iulog,'(A)')         ' Initial 3 point discharge [m3/s]: '
     write(iulog,'(3(A,X,G12.5))') ' Q(0,0)=',Q(0,0),' Q(0,1)=',Q(0,1),' Q(1,0)=',Q(1,0)
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

 else ! if head-water

   Q(1,0) = 0._dp
   Q(1,1) = 0._dp

   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif

 endif

 ! compute volume
 rflux%ROUTE(idxKW)%REACH_VOL(0) = rflux%ROUTE(idxKW)%REACH_VOL(1)
 ! For very low flow condition, outflow - inflow > current storage, so limit outflow and adjust Q(1,1)
 Q(1,1) = min(rflux%ROUTE(idxKW)%REACH_VOL(0)/dt + Q(1,0)*0.999, Q(1,1))
 rflux%ROUTE(idxKW)%REACH_VOL(1) = rflux%ROUTE(idxKW)%REACH_VOL(0) + (Q(1,0)-Q(1,1))*dt

 ! add catchment flow
 rflux%ROUTE(idxKW)%REACH_Q = Q(1,1)+rflux%BASIN_QR(1)

 if (verbose) then
   write(iulog,'(1(A,X,G15.4))') ' Q(1,1)=',Q(1,1)
 end if

 ! Q abstraction
 ! Compute actual abstraction (Qabs) m3/s - values should be negative
 ! Compute abstraction (Qmod) m3 taken from outlet discharge (REACH_Q)
 ! Compute REACH_Q subtracted from Qmod abstraction
 ! Compute REACH_VOL subtracted from total abstraction minus abstraction from outlet discharge
 if (Qtake<0) then
   Qabs = max(-(rflux%ROUTE(idxKW)%REACH_VOL(1)/dt+rflux%ROUTE(idxKW)%REACH_Q-Qmin), Qtake)
   Qmod = min(rflux%ROUTE(idxKW)%REACH_VOL(1) + Qabs*dt, 0._dp) ! Qtake taken from outflow portion, Qmod <=0

   rflux%ROUTE(idxKW)%REACH_Q      = rflux%ROUTE(idxKW)%REACH_Q + Qmod/dt
   rflux%ROUTE(idxKW)%REACH_VOL(1) = rflux%ROUTE(idxKW)%REACH_VOL(1) + (Qabs*dt - Qmod)

   ! modify computational molecule state (Q)
   Q(1,1) = Q(1,1) - max(abs(Qmod/dt)-rflux%BASIN_QR(1), 0._dp)
 end if

 ! update state
 rstate%molecule%Q(1) = Q(1,0)
 rstate%molecule%Q(2) = Q(1,1)

 END SUBROUTINE kinematic_wave

END MODULE kw_route_module
