MODULE dfw_route_module

! diffusive wave routing

USE nrtype
USE dataTypes,     ONLY: STRFLX          ! fluxes in each reach
USE dataTypes,     ONLY: STRSTA          ! state in each reach
USE dataTypes,     ONLY: RCHTOPO         ! Network topology
USE dataTypes,     ONLY: RCHPRP          ! Reach parameter
USE dataTypes,     ONLY: dwRCH           ! dw specific state data structure
USE public_var,    ONLY: iulog           ! i/o logical unit number
USE public_var,    ONLY: realMissing     ! missing value for real number
USE public_var,    ONLY: integerMissing  ! missing value for integer number
USE public_var,    ONLY: dt              ! simulation time step [sec]
USE public_var,    ONLY: is_flux_wm      ! logical water management components fluxes should be read
USE public_var,    ONLY: qmodOption      ! qmod option (use 1==direct insertion)
USE public_var,    ONLY: hw_drain_point  ! headwater catchment pour point (top_reach==1 or bottom_reach==2)
USE globalData,    ONLY: idxDW           ! routing method index for diffusive wave
USE water_balance, ONLY: comp_reach_wb   ! compute water balance error
USE base_route,    ONLY: base_route_rch  ! base (abstract) reach routing method class

implicit none

private
public::dfw_route_rch

integer(i4b), parameter :: top_reach=1
integer(i4b), parameter :: bottom_reach=2

type, extends(base_route_rch) :: dfw_route_rch
 CONTAINS
   procedure, pass :: route => dfw_rch
end type dfw_route_rch

CONTAINS

 ! *********************************************************************
 ! subroutine: perform diffusive wave routing for one segment
 ! *********************************************************************
 SUBROUTINE dfw_rch(this,           & ! dfw_route_rch object to bound this procedure
                    iens, segIndex, & ! input: index of runoff ensemble to be processed
                    ixDesire,       & ! input: reachID to be checked by on-screen pringing
                    T0,T1,          & ! input: start and end of the time step
                    NETOPO_in,      & ! input: reach topology data structure
                    RPARAM_in,      & ! input: reach parameter data structure
                    RCHSTA_out,     & ! inout: reach state data structure
                    RCHFLX_out,     & ! inout: reach flux data structure
                    ierr, message)    ! output: error control

 implicit none
 ! Argument variables
 class(dfw_route_rch)                      :: this              ! dfw_route_rch object to bound this procedure
 integer(i4b),  intent(in)                 :: iens              ! runoff ensemble to be routed
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
 real(dp)                                  :: Qlat              ! lateral flow into channel [m3/s]
 real(dp)                                  :: Qabs              ! maximum allowable water abstraction rate [m3/s]
 real(dp)                                  :: q_upstream        ! total discharge at top of the reach [m3/s]
 real(dp)                                  :: q_upstream_mod    ! total discharge at top of the reach after water abstraction [m3/s]
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ierr=0; message='dfw_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   verbose = .true.
 end if

 ! get discharge coming from upstream
 nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
 isHW = .true.
 q_upstream = 0.0_dp

 Qabs = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initial water abstraction (positive) or injection (negative)
 RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initialize actual water abstraction

 ! update volume at previous time step
 RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(0) = RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(1)

 if (nUps>0) then ! this hru is not headwater
   isHW = .false.
   do iUps = 1,nUps
     if (.not. NETOPO_in(segIndex)%goodBas(iUps)) cycle ! skip upstream reach which does not any flow due to zero total contributory areas
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     if (qmodOption==1 .and. RCHFLX_out(iens,iRch_ups)%Qobs/=realMissing) then
       RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q = RCHFLX_out(iens,iRch_ups)%Qobs
     end if
     q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q
   end do
   q_upstream_mod  = q_upstream
   Qlat = RCHFLX_out(iens,segIndex)%BASIN_QR(1)
 else ! headwater
   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif
   if (hw_drain_point==top_reach) then ! lateral flow is poured in a reach at the top
     q_upstream = q_upstream + RCHFLX_out(iens,segIndex)%BASIN_QR(1)
     q_upstream_mod = q_upstream
     Qlat = 0._dp
   else if (hw_drain_point==bottom_reach) then ! lateral flow is poured in a reach at the top
     q_upstream_mod = q_upstream
     Qlat = RCHFLX_out(iens,segIndex)%BASIN_QR(1)
   end if
 end if

 ! Water management - water injection or abstraction (irrigation or industrial/domestic water usage)
 ! For water abstraction, water is extracted from the following priorities:
 ! 1. existing storage(REACH_VOL(0), 2. upstream inflow , 3 lateral flow (BASIN_QR)
 if((RCHFLX_out(iens,segIndex)%REACH_WM_FLUX /= realMissing).and.(is_flux_wm)) then
   if (Qabs > 0) then ! positive == abstraction
     if (RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(1)/dt > Qabs) then ! take out all abstraction from strorage
       RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(1) - Qabs*dt
     else ! if inital abstraction is greater than volume
       Qabs = Qabs - RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(1)/dt ! get residual Qabs after extracting from strorage
       RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(1) = 0._dp ! voluem gets 0
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
           RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX - Qabs
         end if
       end if
     end if
   else ! negative == injection
     Qlat = Qlat - Qabs
   endif
 endif

 if(verbose)then
   write(iulog,'(2A)') new_line('a'), '** CHECK diffusive wave routing **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,1X,I12,1X,G15.4)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps), &
             RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q
     enddo
   end if
   write(iulog,'(A,1X,G15.4)') ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1)=',RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
 endif

 ! solve diffusive wave equation
 call diffusive_wave(RPARAM_in(segIndex),                     &  ! input: parameter at segIndex reach
                     T0,T1,                                   &  ! input: start and end of the time step
                     q_upstream_mod,                          &  ! input: total discharge at top of the reach being processed
                     Qlat,                                    &  ! input: lateral flow [m3/s]
                     isHW,                                    &  ! input: is this headwater basin?
                     RCHSTA_out(iens,segIndex)%DW_ROUTE,      &  ! inout:
                     RCHFLX_out(iens,segIndex),               &  ! inout: updated fluxes at reach
                     verbose,                                 &  ! input: reach index to be examined
                     ierr, cmessage)                             ! output: error control
 if(ierr/=0)then
    write(message, '(A,1X,I12,1X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage)
    return
 endif

 if(verbose)then
   write(iulog,'(A,1X,G15.4)') ' RCHFLX_out(iens,segIndex)%REACH_Q=', RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_Q
 endif

 if (RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(1) < 0) then
   write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE VOLUME = ', RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_VOL(1), &
         'at ', NETOPO_in(segIndex)%REACHID
 end if

 call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxDW, q_upstream, Qlat, RCHFLX_out(iens,segIndex), verbose, lakeFlag=.false.)

 END SUBROUTINE dfw_rch


 ! *********************************************************************
 ! subroutine: solve diffuisve wave equation
 ! *********************************************************************
 SUBROUTINE diffusive_wave(rch_param,     & ! input: river parameter data structure
                           T0,T1,         & ! input: start and end of the time step
                           q_upstream,    & ! input: discharge from upstream
                           Qlat,          & ! input: lateral discharge into chaneel [m3/s]
                           isHW,          & ! input: is this headwater basin?
                           rstate,        & ! inout: reach state at a reach
                           rflux,         & ! inout: reach flux at a reach
                           verbose,       & ! input: reach index to be examined
                           ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Solve linearlized diffusive wave equation per reach and time step.
 !  dQ/dt + ck*dQ/dx = dk*d2Q/dx2  - a)
 !
 !  ck (celerity) and dk (diffusivity) are computed with previous inflow and outflow and current inflow
 !
 !  1) dQ/dt   = (Q(t,x) - Q(t-1,x))/dt
 !  2) dQ/dx   = [(1-wck)(Q(t-1,x+1)-Q(t-1,x-1)) + wck*(Q(t,x+1)-Q(t,x-1))]/2dx
 !  3) d2Q/dx2 = [(1-wdk)(Q(t-1,x+1)-2Q(t-1,x)+Q(t-1,x-1)) + wdk*(Q(t,x+1)-2Q(t,x)+Q(t,x-1))]/2dx
 !
 !  upstream B.C:   Dirchlet BC with inflow at current time-step,t, from upstream basin
 !  downstream B.C: Neumann BC with prescribed Q gradient (Sbc)
 !                  dQ/dx|x=N = Sbc ->  4) Q(t,N)-Q(t,N-1)) = Sbc*dx
 !  Another downstream B.C option is absorbing boundary condition
 !                  dQ/dt|x=N + ck*dQ/dx|x=N = 0
 !
 !  Inserting 1), 2), 3) and 4) into a) and moving Q(t,*) terms to left-hand side and Q(t-1,*) terms to the righ-hand side
 !  results in tridiagonal matrix equation A*Q = b
 !  where A is [N x N] matrix, Q is [N x 1] vector to be solved (next time step Q) and b is [N x 1] vector
 !  N (nMolecule is used in the code) is the number of internal nodes including upstream and downstream boundaries
 !
 !  Since A is a tridiagonal matrix, the code stores only three diagnoal elements - upper, diagonal, and lower
 !  solving the matrix equation use thomas algorithm
 !
 ! ----------------------------------------------------------------------------------------
 USE globalData, ONLY : nMolecule   ! number of internal nodes for finite difference (including upstream and downstream boundaries)
 implicit none
 ! Argument variables
 type(RCHPRP), intent(in)        :: rch_param      ! River reach parameter
 real(dp),     intent(in)        :: T0,T1          ! start and end of the time step (seconds)
 real(dp),     intent(in)        :: q_upstream     ! total discharge at top of the reach being processed
 real(dp),     intent(in)        :: Qlat           ! lateral discharge into chaneel [m3/s]
 logical(lgt), intent(in)        :: isHW           ! is this headwater basin?
 type(dwRch),  intent(inout)     :: rstate         ! curent reach states
 type(STRFLX), intent(inout)     :: rflux          ! current Reach fluxes
 logical(lgt), intent(in)        :: verbose        ! reach index to be examined
 integer(i4b), intent(out)       :: ierr           ! error code
 character(*), intent(out)       :: message        ! error message
 ! Local variables
 real(dp)                        :: alpha          ! sqrt(slope)(/mannings N* width)
 real(dp)                        :: beta           ! constant, 5/3
 real(dp)                        :: Cd             ! Fourier number
 real(dp)                        :: Ca             ! Courant number
 real(dp)                        :: dx             ! length of segment [m]
 real(dp)                        :: Qbar           ! 3-point average discharge [m3/s]
 real(dp)                        :: Abar           ! 3-point average flow area [m2]
 real(dp)                        :: Vbar           ! 3-point average velocity [m/s]
 real(dp)                        :: ck             ! kinematic wave celerity [m/s]
 real(dp)                        :: dk             ! diffusivity [m2/s]
 real(dp)                        :: Sbc            ! neumann BC slope
 real(dp), allocatable           :: diagonal(:,:)  ! diagonal part of matrix
 real(dp), allocatable           :: b(:)           ! right-hand side of the matrix equation
 real(dp), allocatable           :: Qlocal(:,:)    ! sub-reach & sub-time step discharge at previous and current time step [m3/s]
 real(dp), allocatable           :: Qsolved(:)     ! solved Q at sub-reach at current time step [m3/s]
 real(dp), allocatable           :: Qprev(:)       ! sub-reach discharge at previous time step [m3/s]
 real(dp)                        :: dTsub          ! time inteval for sub time-step [sec]
 real(dp)                        :: wck            ! weight for advection
 real(dp)                        :: wdk            ! weight for diffusion
 real(dp)                        :: pcntReduc      ! flow profile adjustment based on storage [-]
 integer(i4b)                    :: ix,it          ! loop index
 integer(i4b)                    :: ntSub          ! number of sub time-step
 integer(i4b)                    :: Nx             ! number of internal reach segments
 integer(i4b)                    :: downstreamBC   ! method of B.C condition - absorbing or Neumann
 character(len=strLen)           :: fmt1           ! format string
 character(len=strLen)           :: cmessage       ! error message from subroutine
 ! Local parameters
 integer(i4b), parameter         :: absorbingBC=1
 integer(i4b), parameter         :: neumannBC=2    ! flux derivative w.r.t. distance at downstream boundary

 ierr=0; message='diffusive_wave/'

 ! hard-coded parameters
 downstreamBC = neumannBC

 ntSub = 1  ! number of sub-time step
 wck = 1.0  ! weight in advection term
 wdk = 1.0  ! weight in diffusion term 0.0-> fully explicit, 0.5-> Crank-Nicolson, 1.0 -> fully implicit

 Nx = nMolecule%DW_ROUTE - 1  ! Nx: number of internal reach segments

 if (.not. isHW .or. hw_drain_point==top_reach) then

   allocate(Qprev(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   allocate(b(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   allocate(diagonal(nMolecule%DW_ROUTE,3), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! initialize previous time step flow
   Qprev(1:nMolecule%DW_ROUTE) = rstate%molecule%Q     ! flow state at previous time step

   ! Get the reach parameters
   ! A = (Q/alpha)**(1/beta)
   ! Q = alpha*A**beta
   alpha = sqrt(rch_param%R_SLOPE)/(rch_param%R_MAN_N*rch_param%R_WIDTH**(2._dp/3._dp))
   beta  = 5._dp/3._dp
   dx = rch_param%RLENGTH/(Nx-1) ! one extra sub-segment beyond outlet

   if (verbose) then
     write(iulog,'(A,1X,G12.5)') ' length [m]        =',rch_param%RLENGTH
     write(iulog,'(A,1X,G12.5)') ' slope [-]         =',rch_param%R_SLOPE
     write(iulog,'(A,1X,G12.5)') ' channel width [m] =',rch_param%R_WIDTH
     write(iulog,'(A,1X,G12.5)') ' manning coef [-]  =',rch_param%R_MAN_N
   end if

   ! time-step adjustment so Courant number is less than 1
   dTsub = dt/ntSub

   if (verbose) then
     write(iulog,'(A,1X,I3,A,1X,G12.5)') ' No. sub timestep=',nTsub,' sub time-step [sec]=',dTsub
   end if

   allocate(Qlocal(1:nMolecule%DW_ROUTE, 0:1), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   allocate(Qsolved(1:nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   Qlocal(:,0:1) = realMissing
   Qlocal(1:nMolecule%DW_ROUTE, 0) = Qprev ! previous time step
   Qlocal(1,1)  = q_upstream     ! infllow at sub-time step in current time step

   do it = 1, nTsub

     Qbar = (Qlocal(1,1)+Qlocal(1,0)+Qlocal(nMolecule%DW_ROUTE-1,0))/3.0 ! 3 point average discharge [m3/s]
     Abar = (abs(Qbar)/alpha)**(1/beta)                            ! flow area [m2] (manning equation)
     Vbar = 0._dp
     if (Abar>0._dp) Vbar = Qbar/Abar                     ! average velocity [m/s]
     ck   = beta*Vbar                                     ! kinematic wave celerity [m/s]
     dk   = Qbar/(2*rch_param%R_WIDTH*rch_param%R_SLOPE)  ! diffusivity [m2/s]

     Cd = dk*dTsub/dx**2
     Ca = ck*dTsub/dx

     ! create a matrix - current time step
     ! populate tridiagonal elements
     ! diagonal
     diagonal(1,2)             = 1._dp
     diagonal(2:nMolecule%DW_ROUTE-1,2) = 2._dp + 4*wdk*Cd
     if (downstreamBC == absorbingBC) then
       diagonal(nMolecule%DW_ROUTE,2)     = 1._dp + wck*Ca
     else if (downstreamBC == neumannBC) then
       diagonal(nMolecule%DW_ROUTE,2)     = 1._dp
     end if

     ! upper
     diagonal(:,1)           = 0._dp
     diagonal(3:nMolecule%DW_ROUTE,1) = wck*Ca - 2._dp*wdk*Cd

     ! lower
     diagonal(:,3)             = 0._dp
     diagonal(1:nMolecule%DW_ROUTE-2,3) = -wck*Ca - 2._dp*wdk*Cd
     if (downstreamBC == absorbingBC) then
       diagonal(nMolecule%DW_ROUTE-1,3)   = -wck*Ca
     else if (downstreamBC == neumannBC) then
       diagonal(nMolecule%DW_ROUTE-1,3)     = -1._dp
     end if

     ! populate left-hand side
     ! upstream boundary condition
     b(1)             = Qlocal(1,1)
     ! downstream boundary condition
     if (downstreamBC == absorbingBC) then
       b(nMolecule%DW_ROUTE) = (1._dp-(1._dp-wck)*Ca)*Qlocal(nMolecule%DW_ROUTE,0) + (1-wck)*Ca*Qlocal(nMolecule%DW_ROUTE-1,0)
     else if (downstreamBC == neumannBC) then
       Sbc = (Qlocal(nMolecule%DW_ROUTE,0)-Qlocal(nMolecule%DW_ROUTE-1,0))
       b(nMolecule%DW_ROUTE)     = Sbc
     end if
     ! internal node points
     b(2:nMolecule%DW_ROUTE-1) = ((1._dp-wck)*Ca+2._dp*(1._dp-wdk))*Cd*Qlocal(1:nMolecule%DW_ROUTE-2,0)  &
                      + (2._dp-4._dp*(1._dp-wdk)*Cd)*Qlocal(2:nMolecule%DW_ROUTE-1,0)           &
                      - ((1._dp-wck)*Ca - (1._dp-wdk)*Cd)*Qlocal(3:nMolecule%DW_ROUTE,0)

     ! solve matrix equation - get updated Qlocal
     call TDMA(nMolecule%DW_ROUTE, diagonal, b, Qsolved)

     if (verbose) then
       write(fmt1,'(A,I5,A)') '(A,1X',nMolecule%DW_ROUTE,'(1X,G15.4))'
       write(iulog,fmt1) ' Q sub_reqch=', (Qsolved(ix), ix=1,nMolecule%DW_ROUTE)
     end if

     Qlocal(:,1) = Qsolved
     Qlocal(:,0) = Qlocal(:,1)
   end do

   ! For very low flow condition, outflow - inflow may exceed current storage, so limit outflow and adjust flow profile
   if (abs(Qlocal(nMolecule%DW_ROUTE-1,1))>0._dp) then
     pcntReduc = min((rflux%ROUTE(idxDW)%REACH_VOL(1)/dt + Qlocal(1,1) *0.999)/Qlocal(nMolecule%DW_ROUTE-1,1), 1._dp)
     Qlocal(2:nMolecule%DW_ROUTE,1) = Qlocal(2:nMolecule%DW_ROUTE,1)*pcntReduc
   end if

   rflux%ROUTE(idxDW)%REACH_VOL(1) = rflux%ROUTE(idxDW)%REACH_VOL(1) + (Qlocal(1,1) - Qlocal(nMolecule%DW_ROUTE-1,1))*dt

   ! store final outflow in data structure
   rflux%ROUTE(idxDW)%REACH_Q = Qlocal(nMolecule%DW_ROUTE-1,1) + Qlat

   ! update state
   rstate%molecule%Q = Qlocal(:,1)

   if (verbose) then
     write(fmt1,'(A,I5,A)') '(A,1X',nMolecule%DW_ROUTE,'(1X,G15.4))'
     write(iulog,'(A,1X,G12.5)') ' rflux%REACH_Q= ', rflux%ROUTE(idxDW)%REACH_Q
     write(iulog,fmt1) ' Qprev(1:nMolecule)= ', Qprev(1:nMolecule%DW_ROUTE)
     write(iulog,'(A,5(1X,G12.5))') ' Qbar, Abar, Vbar, ck, dk= ',Qbar, Abar, Vbar, ck, dk
     write(iulog,'(A,2(1X,G12.5))') ' Cd, Ca= ', Cd, Ca
     write(iulog,fmt1) ' diagonal(:,1)= ', diagonal(:,1)
     write(iulog,fmt1) ' diagonal(:,2)= ', diagonal(:,2)
     write(iulog,fmt1) ' diagonal(:,3)= ', diagonal(:,3)
     write(iulog,fmt1) ' b= ', b(1:nMolecule%DW_ROUTE)
   end if

 else ! if head-water and pour runnof to the bottom of reach

   rflux%ROUTE(idxDW)%REACH_Q = Qlat

   rflux%ROUTE(idxDW)%REACH_VOL(0) = 0._dp
   rflux%ROUTE(idxDW)%REACH_VOL(1) = 0._dp

   rstate%molecule%Q(1:nMolecule%DW_ROUTE) = 0._dp
   rstate%molecule%Q(nMolecule%DW_ROUTE)   = rflux%ROUTE(idxDW)%REACH_Q

   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif

 endif

 END SUBROUTINE diffusive_wave

 SUBROUTINE TDMA(NX,MAT,b,T)
 ! Solve tridiagonal matrix system of linear equation
 ! NX is the number of unknowns (gridblocks minus boundaries)
 ! Solve system of linear equations, A*T = b where A is tridiagonal matrix
 ! MAT = NX x 3 array holding tri-diagonal portion of A
 ! MAT(NX,1) - uppder diagonals for matrix A
 ! MAT(NX,2) - diagonals for matrix A
 ! MAT(NX,3) - lower diagonals for matrix A
 ! b(NX) - vector of the right hand coefficients b
 ! T(NX) - The solution matrix
 !
 ! example, A
 ! d u 0 0 0
 ! l d u 0 0
 ! 0 l d u 0
 ! 0 0 l d u
 ! 0 0 0 l d
 !
 ! MAT(:,1) = [0, u, u, u, u]
 ! MAT(:,2) = [d, d, d, d, d]
 ! MAT(:,3) = [l, l, l, l, 0]

   implicit none
   ! Argument variables
   integer(i4b),  intent(in)     :: NX     ! number of unknown (= number of matrix size, grid point minus two end points)
   real(dp),      intent(in)     :: MAT(NX,3)
   real(dp),      intent(in)     :: b(NX)
   real(dp),      intent(inout)  :: T(NX)
   ! Local variables
   integer(i4b)                  :: ix
   real(dp)                      :: U(NX)
   real(dp)                      :: D(NX)
   real(dp)                      :: L(NX)
   real(dp)                      :: b1(NX)
   real(dp)                      :: coef

   U(1:NX) = MAT(1:NX,1)
   D(1:NX) = MAT(1:NX,2)
   L(1:NX) = MAT(1:NX,3)
   b1(1:NX) = b(1:NX)
   do ix = 2, NX
     coef  = L(ix-1)/D(ix-1)
     D(ix) = D(ix)-coef*U(ix)
     b1(ix) = b1(ix)-coef*b1(ix-1)
   end do

   T(NX) = b1(NX)/D(NX) ! Starts the countdown of answers
   do ix = NX-1, 1, -1
       T(ix) = (b1(ix) - U(ix+1)*T(ix+1))/D(ix)
   end do

 END SUBROUTINE TDMA

END MODULE dfw_route_module
