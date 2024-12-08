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
USE public_var,    ONLY: min_length_route! minimum reach length for routing to be performed.
USE globalData,    ONLY: idxDW           ! routing method index for diffusive wave
USE water_balance, ONLY: comp_reach_wb   ! compute water balance error
USE base_route,    ONLY: base_route_rch  ! base (abstract) reach routing method class
USE hydraulic,     ONLY: flow_depth
USE hydraulic,     ONLY: water_height
USE hydraulic,     ONLY: celerity
USE hydraulic,     ONLY: diffusivity

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
 real(dp)                                  :: Qupstream        ! total discharge at top of the reach [m3/s]
 real(dp)                                  :: Qupstream_mod    ! total discharge at top of the reach after water abstraction [m3/s]
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ierr=0; message='dfw_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   verbose = .true.
 end if

 ! get discharge coming from upstream
 nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
 isHW = .true.
 Qupstream = 0.0_dp

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
     Qupstream = Qupstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxDW)%REACH_Q
   end do
   Qupstream_mod  = Qupstream
   Qlat = RCHFLX_out(iens,segIndex)%BASIN_QR(1)
 else ! headwater
   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif
   if (hw_drain_point==top_reach) then ! lateral flow is poured in a reach at the top
     Qupstream = Qupstream + RCHFLX_out(iens,segIndex)%BASIN_QR(1)
     Qupstream_mod = Qupstream
     Qlat = 0._dp
   else if (hw_drain_point==bottom_reach) then ! lateral flow is poured in a reach at the top
     Qupstream_mod = Qupstream
     Qlat = RCHFLX_out(iens,segIndex)%BASIN_QR(1)
   end if
 end if

 RCHFLX_out(iens,segIndex)%ROUTE(idxDW)%REACH_INFLOW = Qupstream ! total inflow from the upstream reaches

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
       if (Qupstream > Qabs) then ! then take out all residual abstraction from upstream inflow
         Qupstream_mod = Qupstream - Qabs
       else ! if residual abstraction is still greater than lateral flow
         Qabs = Qabs - Qupstream ! get residual abstraction after extracting upstream inflow and storage.
         Qupstream_mod = 0._dp ! upstream inflow gets 0 (all is gone to abstracted flow).
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
                     Qupstream_mod,                           &  ! input: total discharge at top of the reach being processed
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

 call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxDW, Qupstream, Qlat, RCHFLX_out(iens,segIndex), verbose, lakeFlag=.false.)

 END SUBROUTINE dfw_rch


 ! *********************************************************************
 ! subroutine: solve diffuisve wave equation
 ! *********************************************************************
 SUBROUTINE diffusive_wave(rch_param,     & ! input: river parameter data structure
                           T0,T1,         & ! input: start and end of the time step
                           Qupstream,    & ! input: discharge from upstream
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
 ! ----------------------------------------------------------------------------------------
 USE globalData, ONLY : nMolecule   ! number of internal nodes for finite difference (including upstream and downstream boundaries)
 USE advection_diffusion, ONLY: solve_ade
 implicit none
 ! Argument variables
 type(RCHPRP), intent(in)        :: rch_param      ! River reach parameter
 real(dp),     intent(in)        :: T0,T1          ! start and end of the time step (seconds)
 real(dp),     intent(in)        :: Qupstream      ! total discharge at top of the reach being processed
 real(dp),     intent(in)        :: Qlat           ! lateral discharge into chaneel [m3/s]
 logical(lgt), intent(in)        :: isHW           ! is this headwater basin?
 type(dwRch),  intent(inout)     :: rstate         ! curent reach states
 type(STRFLX), intent(inout)     :: rflux          ! current Reach fluxes
 logical(lgt), intent(in)        :: verbose        ! reach index to be examined
 integer(i4b), intent(out)       :: ierr           ! error code
 character(*), intent(out)       :: message        ! error message
 ! Local variables
 real(dp)                        :: Qbar           ! 3-point average discharge [m3/s]
 real(dp)                        :: depth          ! flow depth [m]
 real(dp)                        :: ck             ! kinematic wave celerity [m/s]
 real(dp)                        :: dk             ! diffusivity [m2/s]
 real(dp), allocatable           :: Qlocal(:,:)    ! sub-reach & sub-time step discharge at previous and current time step [m3/s]
 real(dp), allocatable           :: Qprev(:)       ! sub-reach discharge at previous time step [m3/s]
 real(dp)                        :: dTsub          ! time inteval for sub time-step [sec]
 real(dp)                        :: pcntReduc      ! flow profile adjustment based on storage [-]
 integer(i4b)                    :: ix,it          ! loop index
 integer(i4b)                    :: ntSub          ! number of sub time-step
 character(len=strLen)           :: fmt1           ! format string
 character(len=strLen)           :: cmessage       ! error message from subroutine

 ierr=0; message='diffusive_wave/'

 ntSub = 1  ! number of sub-time step

 associate(S         => rch_param%R_SLOPE,    & ! channel slope
           n         => rch_param%R_MAN_N,    & ! manning n
           bt        => rch_param%R_WIDTH,    & ! channel bottom width
           bankDepth => rch_param%R_DEPTH,    & ! bankfull depth
           zc        => rch_param%SIDE_SLOPE, & ! channel side slope
           zf        => rch_param%FLDP_SLOPE, & ! floodplain slope
           bankVol   => rch_param%R_STORAGE,  & ! bankful volume
           L         => rch_param%RLENGTH)      ! channel length

 if (.not. isHW .or. hw_drain_point==top_reach) then

   if (L > min_length_route) then

   allocate(Qprev(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! initialize previous time step flow
   Qprev(1:nMolecule%DW_ROUTE) = rstate%molecule%Q     ! flow state at previous time step

   if (verbose) then
     write(iulog,'(A,1X,G12.5)') ' length [m]        =',L
     write(iulog,'(A,1X,G12.5)') ' slope [-]         =',S
     write(iulog,'(A,1X,G12.5)') ' channel width [m] =',bt
     write(iulog,'(A,1X,G12.5)') ' manning coef [-]  =',n
   end if

   ! time-step adjustment so Courant number is less than 1
   dTsub = dt/ntSub

   if (verbose) then
     write(iulog,'(A,1X,I3,A,1X,G12.5)') ' No. sub timestep=',nTsub,' sub time-step [sec]=',dTsub
   end if

   allocate(Qlocal(1:nMolecule%DW_ROUTE, 0:1), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   do it = 1, nTsub
     Qbar = (Qupstream+Qprev(1)+Qprev(nMolecule%DW_ROUTE-1))/3.0 ! 3 point average discharge [m3/s]
     depth = flow_depth(abs(Qbar), bt, zc, S, n, zf=zf, bankDepth=bankDepth) ! compute flow depth as normal depth (a function of flow)
     ck    = celerity(abs(Qbar), depth, bt, zc, S, n, zf=zf, bankDepth=bankDepth)
     dk    = diffusivity(abs(Qbar), depth, bt, zc, S, n, zf=zf, bankDepth=bankDepth)

     call solve_ade(L,                  & ! input: river parameter data structure
                    nMolecule%DW_ROUTE, & ! input: number of sub-segments
                    dTsub,              & ! input: time_step [sec]
                    Qupstream,          & ! input: quantity from upstream [unit of quantity]
                    ck,                 & ! input: velocity [m/s]
                    dk,                 & ! input: diffusivity [m2/s]
                    Qlat,               & ! input: lateral quantity into chaneel [unit of quantity]
                    Qprev,              & ! input: quantity at previous time step [unit of quantity]
                    Qlocal,             & ! inout: quantity soloved at current time step [unit of quantity]
                    verbose,            & ! input: reach index to be examined
                    ierr,message)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end do

   ! For very low flow condition, outflow - inflow may exceed current storage, so limit outflow and adjust flow profile
   if (abs(Qlocal(nMolecule%DW_ROUTE-1,1))>0._dp) then
     pcntReduc = min((rflux%ROUTE(idxDW)%REACH_VOL(1)/dt + Qlocal(1,1) *0.999)/Qlocal(nMolecule%DW_ROUTE-1,1), 1._dp)
     Qlocal(2:nMolecule%DW_ROUTE,1) = Qlocal(2:nMolecule%DW_ROUTE,1)*pcntReduc
   end if

   rflux%ROUTE(idxDW)%REACH_VOL(1) = rflux%ROUTE(idxDW)%REACH_VOL(1) + (Qupstream - Qlocal(nMolecule%DW_ROUTE-1,1))*dt

   ! if reach volume exceeds flood threshold volume, excess water is flooded volume.
   if (rflux%ROUTE(idxDW)%REACH_VOL(1) > bankVol) then
     rflux%ROUTE(idxDW)%FLOOD_VOL(1) = rflux%ROUTE(idxDW)%REACH_VOL(1) - bankVol  ! floodplain volume == overflow volume
   else
     rflux%ROUTE(idxDW)%FLOOD_VOL(1) = 0._dp
   end if
   ! compute surface water height [m]
   rflux%ROUTE(idxDW)%REACH_ELE = water_height(rflux%ROUTE(idxDW)%REACH_VOL(1)/L, bt, zc, zf=zf, bankDepth=bankDepth)

   ! store final outflow in data structure
   rflux%ROUTE(idxDW)%REACH_Q = Qlocal(nMolecule%DW_ROUTE-1,1) + Qlat

   ! update state
   rstate%molecule%Q = Qlocal(:,1)

   else ! length < min_length_route: length is short enough to just pass upstream to downstream
     rflux%ROUTE(idxDW)%REACH_Q = Qupstream + Qlat
     rstate%molecule%Q(1:nMolecule%DW_ROUTE) = 0._dp
     rstate%molecule%Q(nMolecule%DW_ROUTE)   = rflux%ROUTE(idxDW)%REACH_Q

     rflux%ROUTE(idxDW)%REACH_VOL(0) = 0._dp
     rflux%ROUTE(idxDW)%REACH_VOL(1) = 0._dp
     rflux%ROUTE(idxDW)%FLOOD_VOL(1) = 0._dp
     rflux%ROUTE(idxDW)%REACH_ELE    = 0._dp
   end if
 else ! if head-water and pour runnof to the bottom of reach

   rflux%ROUTE(idxDW)%REACH_Q = Qlat

   rflux%ROUTE(idxDW)%REACH_VOL(0) = 0._dp
   rflux%ROUTE(idxDW)%REACH_VOL(1) = 0._dp
   rflux%ROUTE(idxDW)%FLOOD_VOL(1) = 0._dp
   rflux%ROUTE(idxDW)%REACH_ELE    = 0._dp

   rstate%molecule%Q(1:nMolecule%DW_ROUTE) = 0._dp
   rstate%molecule%Q(nMolecule%DW_ROUTE)   = rflux%ROUTE(idxDW)%REACH_Q

   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif

 endif

 end associate

 END SUBROUTINE diffusive_wave

END MODULE dfw_route_module
