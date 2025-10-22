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
USE public_var,    ONLY: desireId        ! ID or reach where detailed reach state is print in log
USE public_var,    ONLY: dt              ! simulation time step [sec]
USE public_var,    ONLY: is_flux_wm      ! logical water management components fluxes should be read
USE public_var,    ONLY: qmodOption      ! qmod option (use 1==direct insertion)
USE public_var,    ONLY: hw_drain_point  ! headwater catchment pour point (top_reach==1 or bottom_reach==2)
USE public_var,    ONLY: min_length_route! minimum reach length for routing to be performed.
USE globalData,    ONLY: idxKW           ! routing method index for kinematic wwave
USE water_balance, ONLY: comp_reach_wb   ! compute water balance error
USE base_route,    ONLY: base_route_rch  ! base (abstract) reach routing method class
USE hydraulic,     ONLY: flow_depth
USE hydraulic,     ONLY: water_height
USE hydraulic,     ONLY: Pwet
USE data_assimilation, ONLY: direct_insertion ! qmod option (use 1==direct insertion)

implicit none

private
public::kwe_route_rch

real(dp), parameter  :: critFactor=0.01
integer(i4b), parameter :: top_reach=1
integer(i4b), parameter :: bottom_reach=2

type, extends(base_route_rch) :: kwe_route_rch
 CONTAINS
   procedure, pass :: route => kw_rch
end type kwe_route_rch

CONTAINS

 ! *********************************************************************
 ! subroutine: perform one segment route KW routing
 ! *********************************************************************
 SUBROUTINE kw_rch(this,           & ! kwe_route_rch object to bound this procedure
                   segIndex,       & ! input: index of runoff reach to be processed
                   T0,T1,          & ! input: start and end of the time step
                   NETOPO_in,      & ! input: reach topology data structure
                   RPARAM_in,      & ! input: reach parameter data structure
                   RCHSTA_out,     & ! inout: reach state data structure
                   RCHFLX_out,     & ! inout: reach flux data structure
                   ierr, message)    ! output: error control
 implicit none
 ! Argument variables
 class(kwe_route_rch)                      :: this              ! kwe_route_rch object to bound this procedure
 integer(i4b),  intent(in)                 :: segIndex          ! segment where routing is performed
 real(dp),      intent(in)                 :: T0,T1             ! start and end of the time step (seconds)
 type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
 type(RCHPRP),  intent(inout), allocatable :: RPARAM_in(:)      ! River reach parameter
 type(STRSTA),  intent(inout)              :: RCHSTA_out(:)     ! reach state data
 type(STRFLX),  intent(inout)              :: RCHFLX_out(:)     ! Reach fluxes (space [reaches]) for decomposed domains
 integer(i4b),  intent(out)                :: ierr              ! error code
 character(*),  intent(out)                :: message           ! error message
 ! Local variables
 logical(lgt)                              :: verbose           ! check details of variables
 logical(lgt)                              :: isHW              ! headwater basin?
 integer(i4b)                              :: nUps              ! number of upstream segment
 integer(i4b)                              :: iUps              ! upstream reach index
 integer(i4b)                              :: iRch_ups          ! index of upstream reach in NETOPO
 real(dp)                                  :: Qlat              ! lateral flow into channel [m3/s]
 real(dp)                                  :: Qabs              ! maximum allowable water abstraction rate [m3/s]
 real(dp)                                  :: q_upstream        ! total discharge at top of the reach [m3/s]
 real(dp)                                  :: q_upstream_mod    ! total discharge at top of the reach after water abstraction [m3/s]
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ierr=0; message='kw_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHID == desireId) verbose = .true.

 ! get discharge coming from upstream
 nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
 isHW = .true.
 q_upstream = 0.0_dp

 Qabs = RCHFLX_out(segIndex)%REACH_WM_FLUX ! initial water abstraction (positive) or injection (negative)
 RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_WM_FLUX_actual = RCHFLX_out(segIndex)%REACH_WM_FLUX ! initialize actual water abstraction

 ! update volume at previous time step
 RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(0) = RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(1)

 if (nUps>0) then
   isHW = .false.
   do iUps = 1,nUps
     if (.not. NETOPO_in(segIndex)%goodBas(iUps)) cycle ! skip upstream reach which does not any flow due to zero total contributory areas
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     q_upstream = q_upstream + RCHFLX_out(iRch_ups)%ROUTE(idxKW)%REACH_Q
   end do
   q_upstream_mod  = q_upstream
   Qlat = RCHFLX_out(segIndex)%BASIN_QR(1)
 else ! headwater
   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif
   if (hw_drain_point==top_reach) then ! lateral flow is poured in a reach at the top
     q_upstream = q_upstream + RCHFLX_out(segIndex)%BASIN_QR(1)
     q_upstream_mod = q_upstream
     Qlat = 0._dp
   else if (hw_drain_point==bottom_reach) then ! lateral flow is poured in a reach at the top
     q_upstream_mod = q_upstream
     Qlat = RCHFLX_out(segIndex)%BASIN_QR(1)
   end if
 endif

 RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_INFLOW = q_upstream ! total inflow from the upstream reaches

 ! Water management - water injection or abstraction (irrigation or industrial/domestic water usage)
 ! For water abstraction, water is extracted from the following priorities:
 ! 1. existing storage(REACH_VOL(0), 2. upstream inflow , 3 lateral flow (BASIN_QR)
 if((RCHFLX_out(segIndex)%REACH_WM_FLUX /= realMissing).and.(is_flux_wm)) then
   if (Qabs > 0) then ! positive == abstraction
     if (RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(1)/dt > Qabs) then ! take out all abstraction from strorage
       RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(1) = RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(1) - Qabs*dt
     else ! if inital abstraction is greater than volume
       Qabs = Qabs - RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(1)/dt ! get residual Qabs after extracting from strorage
       RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(1) = 0._dp ! voluem gets 0
       if (q_upstream > Qabs) then ! then take out all residual abstraction from upstream inflow
         q_upstream_mod = q_upstream - Qabs
       else ! if residual abstraction is still greater than lateral flow
         Qabs = Qabs - q_upstream ! get residual abstraction after extracting upstream inflow and storage.
         q_upstream_mod = 0._dp ! upstream inflow gets 0 (all is gone to abstracted flow).
         if (Qlat > Qabs) then ! then take residual abstraction out from lateral flow
           Qlat = Qlat - Qabs
         else ! if residual abstraction is greater than upstream inflow
           Qabs = Qabs - Qlat ! take out residual abstraction from lateral flow
           Qlat = 0._dp ! lateral flow gets 0 (all are gone to abstraction)
           RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_WM_FLUX_actual = RCHFLX_out(segIndex)%REACH_WM_FLUX - Qabs
         end if
       end if
     end if
   else ! negative == injection
     Qlat = Qlat - Qabs
   endif
 endif

 if(verbose)then
   write(iulog,'(2A)') new_line('a'), '** CHECK Kinematic wave routing **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,1X,I12,1X,G12.5)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps), &
             RCHFLX_out(iRch_ups)%ROUTE(idxKW)%REACH_Q
     enddo
   end if
   write(iulog,'(A,1X,G15.4)') ' RCHFLX_out(segIndex)%BASIN_QR(1)=',RCHFLX_out(segIndex)%BASIN_QR(1)
 endif

 ! perform river network KW routing
 call kinematic_wave(RPARAM_in(segIndex),                     & ! input: parameter at segIndex reach
                     T0,T1,                                   & ! input: start and end of the time step
                     q_upstream_mod,                          & ! input: total discharge at top of the reach being processed
                     Qlat,                                    & ! input: lateral flow [m3/s]
                     isHW,                                    & ! input: is this headwater basin?
                     RCHSTA_out(segIndex)%KW_ROUTE,           & ! inout:
                     RCHFLX_out(segIndex),                    & ! inout: updated fluxes at reach
                     verbose,                                 & ! input: reach index to be examined
                     ierr, cmessage)                            ! output: error control
 if(ierr/=0)then
   write(message, '(A,1X,I12,1X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage); return
 endif

 if(verbose)then
   write(iulog,'(A,1X,G15.4)') ' RCHFLX_out(segIndex)%REACH_Q=', RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_Q
 endif

 if (RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(1) < 0) then
   write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE VOLUME = ', RCHFLX_out(segIndex)%ROUTE(idxKW)%REACH_VOL(1), &
         'at ', NETOPO_in(segIndex)%REACHID
 end if

 if (qmodOption==1) then
   call direct_insertion(segIndex,       & ! input: reach index
                         idxKW,          & ! input: routing method id for diffusive wave routing
                         NETOPO_in,      & ! input: reach topology data structure
                         RCHSTA_out,     & ! inout: reach state data structure
                         RCHFLX_out,     & ! inout: reach fluxes datq structure
                         ierr, cmessage)   ! output: error control
   if(ierr/=0)then
     write(message,'(A,X,I12,X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage); return
   endif
 end if

 if (qmodOption==0) then ! check reach water balance only if data assimilation is off
   call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxKW, q_upstream, Qlat, RCHFLX_out(segIndex), verbose, lakeFlag=.false.)
 end if

 END SUBROUTINE kw_rch


 ! *********************************************************************
 ! subroutine: route kinematic waves at one segment
 ! *********************************************************************
 SUBROUTINE kinematic_wave(rch_param,     & ! input: river parameter data structure
                           T0,T1,         & ! input: start and end of the time step
                           q_upstream,    & ! input: discharge from upstream
                           Qlat,          & ! input: lateral discharge into chaneel [m3/s]
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
 real(dp),     intent(in)                 :: Qlat         ! lateral discharge into chaneel [m3/s]
 logical(lgt), intent(in)                 :: isHW         ! is this headwater basin?
 type(kwRCH),  intent(inout)              :: rstate       ! curent reach states
 type(STRFLX), intent(inout)              :: rflux        ! current Reach fluxes
 logical(lgt), intent(in)                 :: verbose      ! reach index to be examined
 integer(i4b), intent(out)                :: ierr         ! error code
 character(*), intent(out)                :: message      ! error message
 ! Local variables
 real(dp)                                 :: depth        ! flow depth [m]
 real(dp)                                 :: p            ! wetness perimeter [m]
 real(dp)                                 :: alpha        ! sqrt(slope)(/mannings N* width)
 real(dp)                                 :: beta         ! constant, 5/3
 real(dp)                                 :: alpha1       ! sqrt(slope)(/mannings N* width)
 real(dp)                                 :: beta1        ! constant, 5/3
 real(dp)                                 :: theta        ! dT/dX
 real(dp)                                 :: omega        ! right-hand side of kw finite difference
 real(dp)                                 :: f0,f1,f2     ! values of function f, 1st and 2nd derivatives at solution
 real(dp)                                 :: X            !
 real(dp)                                 :: dX           ! length of segment [m]
 real(dp)                                 :: Q(0:1,0:1)   !
 real(dp)                                 :: Qtrial(2)    ! trial solution of kw equation
 real(dp)                                 :: Qbar         !
 real(dp)                                 :: absErr(2)    ! absolute error of nonliear equation solution
 real(dp)                                 :: f0eval(2)    !
 integer(i4b)                             :: imin         ! index at minimum value

 ierr=0; message='kinematic_wave/'

 Q(0,0) = rstate%molecule%Q(1) ! previous time and inlet  1 (0,0)
 Q(0,1) = rstate%molecule%Q(2) ! previous time and outlet 2 (0,1)

 associate(S         => rch_param%R_SLOPE,    & ! channel slope
           n         => rch_param%R_MAN_N,    & ! manning n
           bt        => rch_param%R_WIDTH,    & ! channel bottom width
           bankDepth => rch_param%R_DEPTH,    & ! bankfull depth
           zc        => rch_param%SIDE_SLOPE, & ! channel side slope
           zf        => rch_param%FLDP_SLOPE, & ! floodplain slope
           bankVol   => rch_param%R_STORAGE,  & ! bankful volume
           L         => rch_param%RLENGTH)      ! channel length

 if (.not. isHW .or. hw_drain_point==top_reach) then

   if (rch_param%RLENGTH > min_length_route) then
   ! compute total flow rate and flow area at upstream end at current time step
   Q(1,0) = q_upstream
   Q(1,1) = realMissing ! current time and outlet 4 (1,1)

   ! Get the reach parameters
   ! A = (Q/alpha)**(1/beta)
   ! Q = alpha*A**beta
   Qbar   = (Q(0,1) + Q(1,0))/2._dp
   depth = flow_depth(abs(Qbar), bt, zc, S, n, zf=zf, bankDepth=bankDepth) ! compute flow depth as normal depth (a function of flow)
   p = Pwet(depth, bt, zc, zf, bankDepth=bankDepth)
   alpha = sqrt(S)/(n*p**(2._dp/3._dp))
   beta  = 5._dp/3._dp
   beta1  = 1._dp/beta
   alpha1 = (1.0/alpha)**beta1
   theta = dt/L

   if (verbose) then
     write(iulog,'(A,1X,G12.5)') ' length [m]        =',rch_param%RLENGTH
     write(iulog,'(A,1X,G12.5)') ' slope [-]         =',rch_param%R_SLOPE
     write(iulog,'(A,1X,G12.5)') ' channel width [m] =',rch_param%R_WIDTH
     write(iulog,'(A,1X,G12.5)') ' manning coef. [-] =',rch_param%R_MAN_N
     write(iulog,'(A)')          ' Initial 3 point discharge [m3/s]: '
     write(iulog,'(3(A,1X,G12.5))') ' Q(0,0)=',Q(0,0),' Q(0,1)=',Q(0,1),' Q(1,0)=',Q(1,0)
   end if

   ! ----------
   ! solve flow rate and flow area at downstream end at current time step
   ! ----------
   ! initial guess
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
   else ! length < min_length_route: length is short enough to just pass upstream to downstream
     Q(1,0) = q_upstream
     Q(1,1) = q_upstream
   end if
 else ! if head-water and pour runnof to the bottom of reach

   Q(1,0) = 0._dp
   Q(1,1) = 0._dp

   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif

 endif

 if (rch_param%RLENGTH > min_length_route) then
 ! For very low flow condition, outflow - inflow > current storage, so limit outflow and adjust Q(1,1)
 Q(1,1) = min(rflux%ROUTE(idxKW)%REACH_VOL(1)/dt + Q(1,0)*0.999, Q(1,1))
 rflux%ROUTE(idxKW)%REACH_VOL(1) = rflux%ROUTE(idxKW)%REACH_VOL(1) + (Q(1,0)-Q(1,1))*dt
 end if

 ! if reach volume exceeds flood threshold volume, excess water is flooded volume.
 if (rflux%ROUTE(idxKW)%REACH_VOL(1) > bankVol) then
   rflux%ROUTE(idxKW)%FLOOD_VOL(1) = rflux%ROUTE(idxKW)%REACH_VOL(1) - bankVol  ! floodplain volume == overflow volume
 else
   rflux%ROUTE(idxKW)%FLOOD_VOL(1) = 0._dp
 end if
 ! compute surface water height [m]
 rflux%ROUTE(idxKW)%REACH_ELE = water_height(rflux%ROUTE(idxKW)%REACH_VOL(1)/L, bt, zc, zf=zf, bankDepth=bankDepth)

 ! add catchment flow
 rflux%ROUTE(idxKW)%REACH_Q = Q(1,1)+Qlat

 if (verbose) then
   write(iulog,'(1(A,1X,G15.4))') ' Q(1,1)=',Q(1,1)
 end if

 ! update state
 rstate%molecule%Q(1) = Q(1,0)
 rstate%molecule%Q(2) = Q(1,1)

 end associate

 END SUBROUTINE kinematic_wave

END MODULE kw_route_module
