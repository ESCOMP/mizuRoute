MODULE dfw_route_module

USE nrtype
! data types
USE dataTypes,   ONLY: STRFLX            ! fluxes in each reach
USE dataTypes,   ONLY: STRSTA            ! state in each reach
USE dataTypes,   ONLY: RCHTOPO           ! Network topology
USE dataTypes,   ONLY: RCHPRP            ! Reach parameter
USE dataTypes,   ONLY: dwRCH             ! dw specific state data structure
USE dataTypes,   ONLY: subbasin_omp      ! mainstem+tributary data strucuture
! global data
USE public_var,  ONLY: iulog             ! i/o logical unit number
USE public_var,  ONLY: realMissing       ! missing value for real number
USE public_var,  ONLY: integerMissing    ! missing value for integer number
USE public_var,  ONLY: qmodOption        ! qmod option (use 1==direct insertion)
USE globalData,  ONLY: idxDW
! subroutines: general
USE model_finalize, ONLY : handle_err

implicit none

private
public::dfw_route

CONTAINS

 ! *********************************************************************
 ! subroutine: perform diffusive wave routing through the river network
 ! *********************************************************************
 SUBROUTINE dfw_route(iens,                 & ! input: ensemble index
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
   ! Argument variables
   integer(i4b),       intent(in)                 :: iEns                 ! ensemble member
   type(subbasin_omp), intent(in),    allocatable :: river_basin(:)       ! river basin information (mainstem, tributary outlet etc.)
   real(dp),           intent(in)                 :: T0,T1                ! start and end of the time step (seconds)
   integer(i4b),       intent(in)                 :: ixDesire             ! index of the reach for verbose output
   type(RCHTOPO),      intent(in),    allocatable :: NETOPO_in(:)         ! River Network topology
   type(RCHPRP),       intent(in),    allocatable :: RPARAM_in(:)         ! River reach parameter
   type(STRSTA),       intent(inout)              :: RCHSTA_out(:,:)      ! reach state data
   type(STRFLX),       intent(inout)              :: RCHFLX_out(:,:)      ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   integer(i4b),       intent(out)                :: ierr                 ! error code
   character(*),       intent(out)                :: message              ! error message
   integer(i4b),       intent(in),    optional    :: ixSubRch(:)          ! subset of reach indices to be processed
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

   ierr=0; message='dfw_route/'

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
         call dfw_rch(iEns,jSeg,           & ! input: array indices
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

 END SUBROUTINE dfw_route

 ! *********************************************************************
 ! subroutine: perform diffusive wave routing for one segment
 ! *********************************************************************
 SUBROUTINE dfw_rch(iEns, segIndex, & ! input: index of runoff ensemble to be processed
                    ixDesire,       & ! input: reachID to be checked by on-screen pringing
                    T0,T1,          & ! input: start and end of the time step
                    LAKEFLAG,       & ! input: flag if lakes are to be processed
                    NETOPO_in,      & ! input: reach topology data structure
                    RPARAM_in,      & ! input: reach parameter data structure
                    RCHSTA_out,     & ! inout: reach state data structure
                    RCHFLX_out,     & ! inout: reach flux data structure
                    ierr, message)    ! output: error control

 implicit none
 ! Argument variables
 integer(i4b),  intent(in)                 :: iEns              ! runoff ensemble to be routed
 integer(i4b),  intent(in)                 :: segIndex          ! segment where routing is performed
 integer(i4b),  intent(in)                 :: ixDesire          ! index of the reach for verbose output
 real(dp),      intent(in)                 :: T0,T1             ! start and end of the time step (seconds)
 integer(i4b),  intent(in)                 :: LAKEFLAG          ! >0 if processing lakes
 type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
 type(RCHPRP),  intent(in),    allocatable :: RPARAM_in(:)      ! River reach parameter
 type(STRSTA),  intent(inout)              :: RCHSTA_out(:,:)   ! reach state data
 type(STRFLX),  intent(inout)              :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b),  intent(out)                :: ierr              ! error code
 character(*),  intent(out)                :: message           ! error message
 ! Local variables
 logical(lgt)                              :: doCheck           ! check details of variables
 logical(lgt)                              :: isHW              ! headwater basin?
 integer(i4b)                              :: nUps              ! number of upstream segment
 integer(i4b)                              :: iUps              ! upstream reach index
 integer(i4b)                              :: iRch_ups          ! index of upstream reach in NETOPO
 real(dp)                                  :: q_upstream        ! total discharge at top of the reach being processed
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ierr=0; message='dfw_rch/'

 doCheck = .false.
 if(segIndex == ixDesire)then
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
     if (qmodOption==1) then
       if (RCHFLX_out(iens,iRch_ups)%QOBS>0._dp) then
         RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%Qerror = RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q - RCHFLX_out(iens,iRch_ups)%QOBS ! compute error
         RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q = RCHFLX_out(iens,iRch_ups)%QOBS
       else
         RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q = max(RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q-RCHFLX_out(iens,iRch_ups)%ROUTE(idxDW)%Qerror, 0.0001)
       end if
     end if
     q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q
   end do
 endif

 if(doCheck)then
   write(iulog,'(2A)') new_line('a'), '** CHECK diffusive wave routing **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,X,I12,X,G15.4)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps),RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q
     enddo
   end if
   write(iulog,'(A,X,G15.4)') ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1)=',RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
 endif

 ! solve diffusive wave equation
 call diffusive_wave(RPARAM_in(segIndex),                &  ! input: parameter at segIndex reach
                     T0,T1,                              &  ! input: start and end of the time step
                     q_upstream,                         &  ! input: total discharge at top of the reach being processed
                     RCHFLX_out(iens,segIndex)%TAKE,     &  ! input: abstraction(-)/injection(+) [m3/s]
                     RPARAM_in(segIndex)%MINFLOW,        &  ! input: minimum environmental flow [m3/s]
                     isHW,                               &  ! input: is this headwater basin?
                     RCHSTA_out(iens,segIndex)%DW_ROUTE, &  ! inout:
                     RCHFLX_out(iens,segIndex),          &  ! inout: updated fluxes at reach
                     doCheck,                            &  ! input: reach index to be examined
                     ierr, cmessage)                        ! output: error control
 if(ierr/=0)then
    write(message, '(A,X,I12,X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage)
    return
 endif

 if(doCheck)then
   write(iulog,'(A,X,G15.4)') ' RCHFLX_out(iens,segIndex)%REACH_Q=', RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_Q
 endif

 END SUBROUTINE dfw_rch


 ! *********************************************************************
 ! subroutine: solve diffuisve wave equation
 ! *********************************************************************
 SUBROUTINE diffusive_wave(rch_param,     & ! input: river parameter data structure
                           T0,T1,         & ! input: start and end of the time step
                           q_upstream,    & ! input: discharge from upstream
                           Qtake,         & ! input: abstraction(-)/injection(+) [m3/s]
                           Qmin,          & ! input: minimum environmental flow [m3/s]
                           isHW,          & ! input: is this headwater basin?
                           rstate,        & ! inout: reach state at a reach
                           rflux,         & ! inout: reach flux at a reach
                           doCheck,       & ! input: reach index to be examined
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
 real(dp),     intent(in)        :: Qtake          ! abstraction(-)/injection(+) [m3/s]
 real(dp),     intent(in)        :: Qmin           ! minimum environmental flow [m3/s]
 logical(lgt), intent(in)        :: isHW           ! is this headwater basin?
 type(dwRch),  intent(inout)     :: rstate         ! curent reach states
 type(STRFLX), intent(inout)     :: rflux          ! current Reach fluxes
 logical(lgt), intent(in)        :: doCheck        ! reach index to be examined
 integer(i4b), intent(out)       :: ierr           ! error code
 character(*), intent(out)       :: message        ! error message
 ! Local variables
 real(dp)                        :: alpha          ! sqrt(slope)(/mannings N* width)
 real(dp)                        :: beta           ! constant, 5/3
 real(dp)                        :: Cd             ! Fourier number
 real(dp)                        :: Ca             ! Courant number
 real(dp)                        :: dt             ! interval of time step [sec]
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
 real(dp)                        :: QupMod         ! modified total discharge at top of the reach being processed
 real(dp)                        :: Qabs           ! maximum allowable water abstraction rate [m3/s]
 real(dp)                        :: Qmod           ! abstraction rate to be taken from outlet discharge [m3/s]
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
 dt = T1-T0

 Nx = nMolecule%DW_ROUTE - 1  ! Nx: number of internal reach segments

 ! Q injection, add at top of reach
 QupMod = q_upstream
 if (Qtake>0) then
   QupMod = QupMod+ Qtake
 end if

 if (.not. isHW) then

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
   dx = rch_param%RLENGTH/Nx

   if (doCheck) then
     write(iulog,'(A,X,G12.5)') ' length [m]        =',rch_param%RLENGTH
     write(iulog,'(A,X,G12.5)') ' slope [-]         =',rch_param%R_SLOPE
     write(iulog,'(A,X,G12.5)') ' channel width [m] =',rch_param%R_WIDTH
     write(iulog,'(A,X,G12.5)') ' manning coef [-]  =',rch_param%R_MAN_N
   end if

   ! time-step adjustment so Courant number is less than 1
   dTsub = dt/ntSub

   if (doCheck) then
     write(iulog,'(A,X,I3,A,X,G12.5)') ' No. sub timestep=',nTsub,' sub time-step [sec]=',dTsub
   end if

   allocate(Qlocal(1:nMolecule%DW_ROUTE, 0:1), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   allocate(Qsolved(1:nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   Qlocal(:,0:1) = realMissing
   Qlocal(1:nMolecule%DW_ROUTE, 0) = Qprev ! previous time step
   Qlocal(1,1)  = QupMod     ! infllow at sub-time step in current time step

   do it = 1, nTsub

     Qbar = (Qlocal(1,1)+Qlocal(1,0)+Qlocal(nMolecule%DW_ROUTE,0))/3.0 ! 3 point average discharge [m3/s]
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

     if (doCheck) then
       write(fmt1,'(A,I5,A)') '(A,1X',nMolecule%DW_ROUTE,'(1X,G15.4))'
       write(iulog,fmt1) ' Q sub_reqch=', (Qlocal(ix,1), ix=1,nMolecule%DW_ROUTE)
     end if

     Qlocal(:,1) = Qsolved
     Qlocal(:,0) = Qlocal(:,1)
   end do

   ! store final outflow in data structure
   rflux%ROUTE(idxDW)%REACH_Q = Qlocal(nMolecule%DW_ROUTE,1) + rflux%BASIN_QR(1)

   if (doCheck) then
     write(fmt1,'(A,I5,A)') '(A,1X',nMolecule%DW_ROUTE,'(1X,G15.4))'
     write(iulog,'(A,X,G12.5)') 'rflux%REACH_Q= ', rflux%ROUTE(idxDW)%REACH_Q
     write(iulog,fmt1) 'Qprev(1:nMolecule)= ', Qprev(1:nMolecule%DW_ROUTE)
     write(iulog,'(A,5(1X,G12.5))') 'Qbar, Abar, Vbar, ck, dk= ',Qbar, Abar, Vbar, ck, dk
     write(iulog,'(A,2(1X,G12.5))') 'Cd, Ca= ', Cd, Ca
     write(iulog,fmt1) 'diagonal(:,1)= ', diagonal(:,1)
     write(iulog,fmt1) 'diagonal(:,2)= ', diagonal(:,2)
     write(iulog,fmt1) 'diagonal(:,3)= ', diagonal(:,3)
     write(iulog,fmt1) 'b= ', b(1:nMolecule%DW_ROUTE)
   end if

   ! compute volume
   rflux%ROUTE(idxDW)%REACH_VOL(0) = rflux%ROUTE(idxDW)%REACH_VOL(1)
   rflux%ROUTE(idxDW)%REACH_VOL(1) = rflux%ROUTE(idxDW)%REACH_VOL(0) + (Qlocal(1,1) - Qlocal(nMolecule%DW_ROUTE,1))*dt

   ! update state
   rstate%molecule%Q = Qlocal(:,1)

 else ! if head-water

   rflux%ROUTE(idxDW)%REACH_Q = rflux%BASIN_QR(1)

   rflux%ROUTE(idxDW)%REACH_VOL(0) = 0._dp
   rflux%ROUTE(idxDW)%REACH_VOL(1) = 0._dp

   rstate%molecule%Q(1:nMolecule%DW_ROUTE) = 0._dp
   rstate%molecule%Q(nMolecule%DW_ROUTE)   = rflux%ROUTE(idxDW)%REACH_Q

   if (doCheck) then
     write(iulog,'(A)')            ' This is headwater '
   endif

 endif

 ! Q abstraction
 ! Compute actual abstraction (Qabs) m3/s - values should be negative
 ! Compute abstraction (Qmod) m3 taken from outlet discharge (REACH_Q)
 ! Compute REACH_Q subtracted from Qmod abstraction
 ! Compute REACH_VOL subtracted from total abstraction minus abstraction from outlet discharge
 if (Qtake<0) then
   Qabs = max(-(rflux%ROUTE(idxDW)%REACH_VOL(1)/dt+rflux%ROUTE(idxDW)%REACH_Q-Qmin), Qtake)
   Qabs = min(Qabs, 0._dp)  ! this should be <=0
   Qmod = min(rflux%ROUTE(idxDW)%REACH_VOL(1) + Qabs*dt, 0._dp) ! this should be <=0

   rflux%ROUTE(idxDW)%REACH_Q      = rflux%ROUTE(idxDW)%REACH_Q + Qmod/dt
   rflux%ROUTE(idxDW)%REACH_VOL(1) = rflux%ROUTE(idxDW)%REACH_VOL(1) + (Qabs*dt - Qmod)

   ! modify computational molecule state (Q)
   rstate%molecule%Q(nMolecule%DW_ROUTE) = Qlocal(nMolecule%DW_ROUTE,1) - max(abs(Qmod/dt)-rflux%BASIN_QR(1), 0._dp)
 end if

 if (doCheck) then
   write(iulog,'(A,X,G12.5)') ' Qout(t) =', rflux%ROUTE(idxDW)%REACH_Q
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
