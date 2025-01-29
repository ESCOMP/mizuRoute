MODULE tracer_module

! solve 1-D advection(-diffusion) equation from upstream to downstream

USE nrtype
USE dataTypes,     ONLY: STRFLX                ! fluxes in each reach
USE dataTypes,     ONLY: SUBRCH                ! states at sub-reach
USE dataTypes,     ONLY: STRSTA                ! state in each reach
USE dataTypes,     ONLY: RCHTOPO               ! Network topology
USE dataTypes,     ONLY: RCHPRP                ! Reach parameter
USE public_var,    ONLY: iulog                 ! i/o logical unit number
USE public_var,    ONLY: realMissing           ! missing value for real number
USE public_var,    ONLY: integerMissing        ! missing value for integer number
USE public_var,    ONLY: dt                    ! simulation time step [sec]
USE public_var,    ONLY: hw_drain_point        ! headwater catchment pour point (top_reach==1 or bottom_reach==2)
USE public_var,    ONLY: min_length_route      ! minimum reach length for routing to be performed.
USE public_var,    ONLY: kinematicWave
USE public_var,    ONLY: muskingumCunge
USE public_var,    ONLY: diffusiveWave
USE globalData,    ONLY: idxKW,idxMC,idxDW
USE water_balance, ONLY: comp_reach_wb         ! compute water balance error
USE hydraulic,     ONLY: flow_depth
USE hydraulic,     ONLY: flow_area

implicit none

private
public::constituent_rch

integer(i4b)            :: route_method=5    ! routing method used for water routing
integer(i4b), parameter :: top_reach=1
integer(i4b), parameter :: bottom_reach=2

CONTAINS

 ! *********************************************************************
 ! subroutine: perform diffusive wave routing for one segment
 ! *********************************************************************
 SUBROUTINE constituent_rch(iens, segIndex,   & ! input: index of runoff ensemble to be processed
                            ixDesire,         & ! input: reachID to be checked by on-screen pringing
                            NETOPO_in,        & ! input: reach topology data structure
                            RPARAM_in,        & ! input: reach parameter data structure
                            RCHSTA_out,       & ! inout: reach state data structure
                            RCHFLX_out,       & ! inout: reach flux data structure
                            ierr, message)      ! output: error control

 implicit none
 ! Argument variables
 integer(i4b),     intent(in)                 :: iens              ! runoff ensemble to be routed
 integer(i4b),     intent(in)                 :: segIndex          ! segment where routing is performed
 integer(i4b),     intent(in)                 :: ixDesire          ! index of the reach for verbose output
 type(RCHTOPO),    intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
 type(RCHPRP),     intent(inout), allocatable :: RPARAM_in(:)      ! River reach parameter
 type(STRSTA),     intent(inout)              :: RCHSTA_out(:,:)   ! reach state data
 type(STRFLX),     intent(inout)              :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b),     intent(out)                :: ierr              ! error code
 character(*),     intent(out)                :: message           ! error message
 ! Local variables
 logical(lgt)                                 :: verbose           ! check details of variables
 logical(lgt)                                 :: isHW              ! headwater basin?
 integer(i4b)                                 :: idxRoute          ! routing method index
 integer(i4b)                                 :: nUps              ! number of upstream segment
 integer(i4b)                                 :: iUps              ! upstream reach index
 integer(i4b)                                 :: iRch_ups          ! index of upstream reach in NETOPO
 type(SUBRCH)                                 :: subrch_state      ! curent states at sub-reaches
 real(dp)                                     :: Clat              ! lateral flow into channel [m3/s]
 real(dp)                                     :: Cupstream         ! total discharge at top of the reach [m3/s]
 character(len=strLen)                        :: cmessage          ! error message from subroutine

 ierr=0; message='constituent_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   verbose = .true.
 end if

 select case(route_method)
   case(muskingumCunge);  idxRoute = idxMC; subrch_state=RCHSTA_out(iens,segIndex)%MC_ROUTE%molecule
   case(kinematicWave);   idxRoute = idxKW; subrch_state=RCHSTA_out(iens,segIndex)%KW_ROUTE%molecule
   case(diffusiveWave);   idxRoute = idxDW; subrch_state=RCHSTA_out(iens,segIndex)%DW_ROUTE%molecule
   case default
     message=trim(message)//'route_method used for tracer needs to be 3(kinematic wave), 4(muskingum cunge), or 5(diffusive wave)'; ierr=81; return
 end select

 ! get mass flux from upstream
 nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
 isHW = .true.
 Cupstream = 0.0_dp

 ! update volume at previous time step
 RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%reach_mass(0) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%reach_mass(1)

 if (nUps>0) then ! this hru is not headwater
   isHW = .false.
   do iUps = 1,nUps
     if (.not. NETOPO_in(segIndex)%goodBas(iUps)) cycle ! skip upstream reach which does not any flow due to zero total contributory areas
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     Cupstream = Cupstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxRoute)%reach_cc
   end do
   Clat = RCHFLX_out(iens,segIndex)%basin_c
 else ! headwater
   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif
   if (hw_drain_point==top_reach) then ! lateral flow is poured in a reach at the top
     Cupstream = Cupstream + RCHFLX_out(iens,segIndex)%basin_c
     Clat = 0._dp
   else if (hw_drain_point==bottom_reach) then ! lateral flow is poured in a reach at the top
     Clat = RCHFLX_out(iens,segIndex)%basin_c
   end if
 end if

 if(verbose)then
   write(iulog,'(2A)') new_line('a'), '** CHECK tracer **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,1X,I12,1X,G15.4)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps), &
             RCHFLX_out(iens, iRch_ups)%ROUTE(idxRoute)%reach_cc
     enddo
   end if
   write(iulog,'(A,1X,G15.4)') ' RCHFLX_out(iEns,segIndex)%basin_c=',RCHFLX_out(iEns,segIndex)%basin_c
 endif

 ! solve diffusive wave equation
 call comp_mass_flux(RPARAM_in(segIndex),                     &  ! input: parameter at segIndex reach
                     Cupstream,                               &  ! input: total discharge at top of the reach being processed
                     Clat,                                    &  ! input: lateral flow [m3/s]
                     isHW,                                    &  ! input: is this headwater basin?
                     subrch_state,                            &  ! inout:
                     RCHFLX_out(iens,segIndex),               &  ! inout: updated fluxes at reach
                     route_method,                            &  ! input: routing method used
                     verbose,                                 &  ! input: reach index to be examined
                     ierr, cmessage)                             ! output: error control
 if(ierr/=0)then
    write(message, '(A,1X,I12,1X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage)
    return
 endif

 if(verbose)then
   write(iulog,'(A,1X,G15.4)') ' RCHFLX_out(iens,segIndex)%reach_cc=', RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%reach_cc
 endif

 if (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%reach_mass(1) < 0) then
   write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE MASS = ', RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%reach_mass(1), &
         'at ', NETOPO_in(segIndex)%REACHID
 end if

 END SUBROUTINE constituent_rch


 ! *********************************************************************
 ! subroutine: solve diffuisve wave equation
 ! *********************************************************************
 SUBROUTINE comp_mass_flux(rch_param,      & ! input: river parameter data structure
                           Cupstream,      & ! input: discharge from upstream
                           Clat,           & ! input: lateral discharge into chaneel [m3/s]
                           isHW,           & ! input: is this headwater basin?
                           subrch_state,   & ! inout: state at a sub-reach
                           rflux,          & ! inout: reach flux at a reach
                           route_method,   & ! input: routining method used
                           verbose,        & ! input: reach index to be examined
                           ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Solve linearlized diffusive wave equation per reach and time step.
 !  dQ/dt + ck*dQ/dx = dk*d2Q/dx2  - a)
 !
 !  ck (celerity) and dk (diffusivity) are computed with previous inflow and outflow and current inflow
 !
 ! ----------------------------------------------------------------------------------------
 USE globalData, ONLY : nMolecule   ! number of internal nodes for finite difference (including upstream and downstream boundaries)
 USE advection_diffusion, ONLY: solve_ade
 implicit none
 ! Argument variables
 type(RCHPRP), intent(in)        :: rch_param      ! River reach parameter
 real(dp),     intent(in)        :: Cupstream      ! total discharge at top of the reach being processed
 real(dp),     intent(in)        :: Clat           ! lateral discharge into chaneel [m3/s]
 logical(lgt), intent(in)        :: isHW           ! is this headwater basin?
 type(SUBRCH), intent(inout)     :: subrch_state   ! curent states at sub-reaches
 type(STRFLX), intent(inout)     :: rflux          ! current reach fluxes
 integer(i4b), intent(in)        :: route_method   ! routing method used for water routing
 logical(lgt), intent(in)        :: verbose        ! reach index to be examined
 integer(i4b), intent(out)       :: ierr           ! error code
 character(*), intent(out)       :: message        ! error message
 ! Local variables
 real(dp)                        :: Ybar           ! mean reach flow depth over sub-reach points [m]
 real(dp)                        :: Qbar           ! mean reach discharge over sub-reach points [m3/s]
 real(dp)                        :: Abar           ! mean reach flow area over sub-reach points [m2]
 real(dp)                        :: Vbar           ! cross-sectional mean water velocity [m/s]
 real(dp)                        :: dk             ! diffusivity [m2/s]
 real(dp), allocatable           :: Clocal(:,:)    ! sub-reach & sub-time step discharge at previous and current time step [m3/s]
 real(dp), allocatable           :: Cprev(:)       ! sub-reach discharge at previous time step [m3/s]
 real(dp)                        :: dTsub          ! time inteval for sub time-step [sec]
 real(dp)                        :: pcntReduc      ! flow profile adjustment based on storage [-]
 integer(i4b)                    :: idxRoute       ! routing method index
 integer(i4b)                    :: nMole          ! number of sub-nodes
 integer(i4b)                    :: it             ! loop index
 integer(i4b)                    :: ntSub          ! number of sub time-step
 character(len=strLen)           :: cmessage       ! error message from subroutine

 ierr=0; message='comp_mass_flux/'

 ntSub = 1  ! number of sub-time step

 select case(route_method)
   case(muskingumCunge);  idxRoute = idxMC; nMole=nMolecule%MC_ROUTE
   case(kinematicWave);   idxRoute = idxKW; nMole=nMolecule%KW_ROUTE
   case(diffusiveWave);   idxRoute = idxDW; nMole=nMolecule%DW_ROUTE
   case default
     message=trim(message)//'route_method used for tracer needs to be 3(kinematic wave), 4(muskingum cunge), or 5(diffusive wave)'; ierr=81; return
 end select

 associate(S         => rch_param%R_SLOPE,    & ! channel slope
           n         => rch_param%R_MAN_N,    & ! manning n
           bt        => rch_param%R_WIDTH,    & ! channel bottom width
           bankDepth => rch_param%R_DEPTH,    & ! bankfull depth
           zc        => rch_param%SIDE_SLOPE, & ! channel side slope
           zf        => rch_param%FLDP_SLOPE, & ! floodplain slope
           bankVol   => rch_param%R_STORAGE,  & ! bankful volume
           L         => rch_param%RLENGTH)      ! channel length

 if (.not. isHW .or. hw_drain_point==top_reach) then

   allocate(Cprev(nMole), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! initialize previous time step flow
   Cprev(1:nMole) = subrch_state%CC     ! flow state at previous time step

   if (verbose) then
     write(iulog,'(A,1X,G12.5)') ' length [m] =',L
   end if

   ! time-step adjustment so Courant number is less than 1
   dTsub = dt/ntSub

   if (verbose) then
     write(iulog,'(A,1X,I3,A,1X,G12.5)') ' No. sub timestep=',nTsub,' sub time-step [sec]=',dTsub
   end if

   allocate(Clocal(1:nMole, 0:1), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   do it = 1, nTsub
     ! use constant flow velocity
     Qbar = sum(subrch_state%Q)/nMole
     Ybar = flow_depth(Qbar, bt, zc, S, n, zf=zf, bankDepth=bankDepth)
     Abar = flow_area(abs(Ybar), bt, zc, zf=zf, bankDepth=bankDepth)
     Vbar =  Qbar/Abar
     dk = 0._dp   ! dk=0 -> no diffusion, only advection
     call solve_ade(L,                  & ! input: river parameter data structure
                    nMole,              & ! input: number of sub-segments
                    dTsub,              & ! input: time_step [sec]
                    Cupstream,          & ! input: quantity from upstream [unit of quantity]
                    Vbar,               & ! input: reach mean flow velocity [m/s]
                    dk,                 & ! input: diffusivity [m2/s]
                    Clat,               & ! input: lateral quantity into chaneel [unit of quantity]
                    Cprev,              & ! input: quantity at previous time step [unit of quantity]
                    Clocal,             & ! inout: quantity soloved at current time step [unit of quantity]
                    verbose,            & ! input: reach index to be examined
                    ierr,message)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end do

   ! For very low flow condition, outflow - inflow may exceed current storage, so limit outflow and adjust flow profile
   if (abs(Clocal(nMole-1,1))>0._dp) then
     pcntReduc = min((rflux%ROUTE(idxRoute)%reach_mass(1)/dt + Clocal(1,1) *0.999)/Clocal(nMole-1,1), 1._dp)
     Clocal(2:nMole,1) = Clocal(2:nMole,1)*pcntReduc
   end if

   rflux%ROUTE(idxRoute)%reach_mass(1) = rflux%ROUTE(idxRoute)%reach_mass(1) + (Cupstream - Clocal(nMole-1,1))*dt

   ! store final outflow in data structure
   rflux%ROUTE(idxRoute)%reach_cc = Clocal(nMole-1,1) + Clat

   ! update state
   subrch_state%CC = Clocal(:,1)

 else ! if head-water and pour runnof to the bottom of reach

   rflux%ROUTE(idxRoute)%reach_cc = Clat

   rflux%ROUTE(idxRoute)%reach_mass(0) = 0._dp
   rflux%ROUTE(idxRoute)%reach_mass(1) = 0._dp

   subrch_state%CC(1:nMole) = 0._dp
   subrch_state%CC(nMole)   = rflux%ROUTE(idxRoute)%reach_cc

   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif
 endif

 end associate

 END SUBROUTINE comp_mass_flux

END MODULE tracer_module
