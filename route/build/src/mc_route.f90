MODULE mc_route_module

! muskingum-cunge routing

USE nrtype
USE dataTypes,     ONLY: STRFLX          ! fluxes in each reach
USE dataTypes,     ONLY: STRSTA          ! state in each reach
USE dataTypes,     ONLY: RCHTOPO         ! Network topology
USE dataTypes,     ONLY: RCHPRP          ! Reach parameter
USE dataTypes,     ONLY: mcRCH           ! MC specific state data structure
USE public_var,    ONLY: iulog           ! i/o logical unit number
USE public_var,    ONLY: realMissing     ! missing value for real number
USE public_var,    ONLY: integerMissing  ! missing value for integer number
USE public_var,    ONLY: desireId        ! ID or reach where detailed reach state is print in log
USE public_var,    ONLY: dt              ! simulation time step [sec]
USE public_var,    ONLY: qmodOption      ! qmod option (use 1==direct insertion)
USE public_var,    ONLY: hw_drain_point  ! headwater catchment pour point (top_reach==1 or bottom_reach==2)
USE public_var,    ONLY: is_flux_wm      ! logical water management components fluxes should be read
USE public_var,    ONLY: min_length_route! minimum reach length for routing to be performed.
USE globalData,    ONLY: idxMC           ! routing method index for muskingum method
USE water_balance, ONLY: comp_reach_wb   ! compute water balance error
USE base_route,    ONLY: base_route_rch  ! base (abstract) reach routing method class
USE hydraulic,     ONLY: flow_depth      ! function to compute Normal flow depth based on manning equation [m]
USE hydraulic,     ONLY: water_height    ! function to compute water height bsed on flow area [m]
USE hydraulic,     ONLY: celerity        ! function to compute celerity [m/s2]
USE hydraulic,     ONLY: Btop            ! function to compute flow top width [m]
USE data_assimilation, ONLY: direct_insertion ! qmod option (use 1==direct insertion)

implicit none

private
public::mc_route_rch

integer(i4b), parameter :: top_reach=1
integer(i4b), parameter :: bottom_reach=2

type, extends(base_route_rch) :: mc_route_rch
 CONTAINS
   procedure, pass :: route => mc_rch
end type mc_route_rch

CONTAINS

 ! *********************************************************************
 ! subroutine: perform muskingum-cunge routing for one segment
 ! *********************************************************************
 SUBROUTINE mc_rch(this,           & ! mc_route_rch object to bound this procedure
                   iEns, segIndex, & ! input: index of runoff ensemble to be processed
                   T0,T1,          & ! input: start and end of the time step
                   NETOPO_in,      & ! input: reach topology data structure
                   RPARAM_in,      & ! input: reach parameter data structure
                   RCHSTA_out,     & ! inout: reach state data structure
                   RCHFLX_out,     & ! inout: reach flux data structure
                   ierr, message)    ! output: error control

 implicit none
 ! Argument variables
 class(mc_route_rch)                       :: this              ! mc_route_rch object to bound this procedure
 integer(i4b),  intent(in)                 :: iEns              ! runoff ensemble to be routed
 integer(i4b),  intent(in)                 :: segIndex          ! segment where routing is performed
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

 ierr=0; message='mc_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHID == desireId) verbose = .true.

 ! get discharge coming from upstream
 nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
 isHW = .true.
 q_upstream = 0.0_dp

 Qabs = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initial water abstraction (positive) or injection (negative)
 RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initialize actual water abstraction

 ! update volume at previous time step
 RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(0) = RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1)

 if (nUps>0) then ! this hru is not headwater
   isHW = .false.
   do iUps = 1,nUps
     if (.not. NETOPO_in(segIndex)%goodBas(iUps)) cycle ! skip upstream reach which does not any flow due to zero total contributory areas
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxMC)%REACH_Q
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
 endif

 RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_INFLOW = q_upstream ! total inflow from the upstream reaches

 ! Water management - water injection or abstraction (irrigation or industrial/domestic water usage)
 ! For water abstraction, water is extracted from the following priorities:
 ! 1. existing storage(REACH_VOL(0), 2. upstream inflow , 3 lateral flow (BASIN_QR)
 if((RCHFLX_out(iens,segIndex)%REACH_WM_FLUX /= realMissing).and.(is_flux_wm)) then
   if (Qabs > 0) then ! positive == abstraction
     if (RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1)/dt > Qabs) then ! take out all abstraction from strorage
       RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1) - Qabs*dt
     else ! if inital abstraction is greater than volume
       Qabs = Qabs - RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1)/dt ! get residual Qabs after extracting from strorage
       RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1) = 0._dp ! voluem gets 0
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
           RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX - Qabs
         end if
       end if
     end if
   else ! negative == injection
     Qlat = Qlat - Qabs
   endif
 endif

 if(verbose)then
   write(iulog,'(2A)') new_line('a'), '** Check muskingum-cunge routing **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,1X,I12,1X,G12.5)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps), &
             RCHFLX_out(iens, iRch_ups)%ROUTE(idxMC)%REACH_Q
     enddo
   end if
   write(iulog,'(A,1X,G12.5)') ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1)=',RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
 endif

 ! solve muskingum-cunge alogorithm
 call muskingum_cunge(RPARAM_in(segIndex),                     & ! input: parameter at segIndex reach
                      T0,T1,                                   & ! input: start and end of the time step
                      q_upstream_mod,                          & ! input: total discharge at top of the reach being processed
                      Qlat,                                    & ! input: lateral flow [m3/s]
                      isHW,                                    & ! input: is this headwater basin?
                      RCHSTA_out(iens,segIndex)%MC_ROUTE,      & ! inout:
                      RCHFLX_out(iens,segIndex),               & ! inout: updated fluxes at reach
                      verbose,                                 & ! input: reach index to be examined
                      ierr, cmessage)                            ! output: error control
 if(ierr/=0)then
   write(message, '(A,1X,I10,1X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage); return
 endif

 if(verbose)then
   write(iulog,'(A,1X,G12.5)') ' RCHFLX_out(iens,segIndex)%REACH_Q=', RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_Q
 endif

 if (RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1) < 0._dp) then
   write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE VOLUME = ', RCHFLX_out(iens,segIndex)%ROUTE(idxMC)%REACH_VOL(1), &
         'at ', NETOPO_in(segIndex)%REACHID
 end if

 if (qmodOption==1) then
   call direct_insertion(iens, segIndex, & ! input: reach index
                         idxMC,          & ! input: routing method id for diffusive wave routing
                         NETOPO_in,      & ! input: reach topology data structure
                         RCHSTA_out,     & ! inout: reach state data structure
                         RCHFLX_out,     & ! inout: reach fluxes datq structure
                         ierr, cmessage)   ! output: error control
   if(ierr/=0)then
     write(message,'(A,X,I12,X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage); return
   endif
 end if

 if (qmodOption==0) then ! check reach water balance only if data assimilation is off
   call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxMC, q_upstream, Qlat, RCHFLX_out(iens,segIndex), verbose, lakeFlag=.false.)
 end if

 END SUBROUTINE mc_rch

 ! *********************************************************************
 ! subroutine: solve muskingum equation
 ! *********************************************************************
 SUBROUTINE muskingum_cunge(rch_param,     & ! input: river parameter data structure
                            T0,T1,         & ! input: start and end of the time step
                            q_upstream,    & ! input: discharge from upstream
                            Qlat,          & ! input: lateral discharge into chaneel [m3/s]
                            isHW,          & ! input: is this headwater basin?
                            rstate,        & ! inout: reach state at a reach
                            rflux,         & ! inout: reach flux at a reach
                            verbose,       & ! input: reach index to be examined
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

 ! -- EBK 06/26/2023 -- comment out isnan check, doesn't seem to be needed
 ! Use of shr_infnan_isnan will require changes to the standalone build, and
 ! this version is required to work on all compilers.
 !use shr_infnan_mod, only : isnan => shr_infnan_isnan
 implicit none
 ! Argument variables
 type(RCHPRP), intent(in)                 :: rch_param    ! River reach parameter
 real(dp),     intent(in)                 :: T0,T1        ! start and end of the time step (seconds)
 real(dp),     intent(in)                 :: q_upstream   ! total discharge at top of the reach being processed
 real(dp),     intent(in)                 :: Qlat         ! lateral discharge into chaneel [m3/s]
 logical(lgt), intent(in)                 :: isHW         ! is this headwater basin?
 type(mcRCH),  intent(inout)              :: rstate       ! curent reach states
 type(STRFLX), intent(inout)              :: rflux        ! current Reach fluxes
 logical(lgt), intent(in)                 :: verbose      ! reach index to be examined
 integer(i4b), intent(out)                :: ierr         ! error code
 character(*), intent(out)                :: message      ! error message
 ! Local variables
 real(dp)                                 :: theta        ! dT/L
 real(dp)                                 :: X            ! X-factor in descreterized kinematic wave
 real(dp)                                 :: Q(0:1,0:1)   ! discharge at computational molecule
 real(dp)                                 :: Qbar         ! 3-point average discharge [m3/s]
 real(dp)                                 :: depth        ! flow depth [m]
 real(dp)                                 :: topWidth     ! flow top width [m]
 real(dp)                                 :: ck           ! kinematic wave celerity [m/s]
 real(dp)                                 :: Cn           ! Courant number [-]
 real(dp)                                 :: dTsub        ! time inteval for sut time-step [sec]
 real(dp)                                 :: C0,C1,C2     ! muskingum parameters
 real(dp)                                 :: pcntReduc    ! flow profile adjustment based on storage [-]
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

   theta = dt/L    ! [s/m]

   ! compute total flow rate and flow area at upstream end at current time step
   Q(1,0) = q_upstream

   if (verbose) then
     write(iulog,'(A,1X,G12.5)') ' length [m]        =',rch_param%RLENGTH
     write(iulog,'(A,1X,G12.5)') ' slope [-]         =',rch_param%R_SLOPE
     write(iulog,'(A,1X,G12.5)') ' channel width [m] =',rch_param%R_WIDTH
     write(iulog,'(A,1X,G12.5)') ' manning coef [-]  =',rch_param%R_MAN_N
     write(iulog,'(A)')          ' Initial 3 point discharge [m3/s]: '
     write(iulog,'(3(A,1X,G12.5))') ' Qin(t-1) Q(0,0)=',Q(0,0),' Qin(t) Q(1,0)=',Q(1,0),' Qout(t-1) Q(0,1)=',Q(0,1)
   end if

   ! first, using 3-point average in computational molecule, check Cournat number is less than 1, otherwise subcycle within one time step
   Qbar = (Q(0,0)+Q(1,0)+Q(0,1))/3.0  ! average discharge [m3/s]

   if (Qbar>0._dp) then
     depth = flow_depth(abs(Qbar), bt, zc, S, n, zf=zf, bankDepth=bankDepth) ! compute flow depth as normal depth (a function of flow)
     ck   = celerity(abs(Qbar), depth, bt, zc, S, n, zf=zf, bankDepth=bankDepth)
     Cn   = ck*theta                    ! Courant number [-]

     ! time-step adjustment so Courant number is less than 1
     ntSub = 1
     dTsub = dt
     if (Cn>1.0_dp) then
       ntSub = ceiling(dt/L*ck)
       dTsub = dt/ntSub
     end if
     if (verbose) then
       write(iulog,'(A,1X,I3,A,1X,G12.5)') ' No. sub timestep=',nTsub,' sub time-step [sec]=',dTsub
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
       if (Qbar>0._dp) then
         depth    = flow_depth(abs(Qbar), bt, zc, S, n, zf=zf, bankDepth=bankDepth) ! compute flow depth as normal depth (a function of flow)
         topWidth = Btop(depth, bt, zc, zf=zf, bankDepth=bankDepth) ! top width at water level [m] (rectangular channel)
         ck       = celerity(abs(Qbar), depth, bt, zc, S, n, zf=zf, bankDepth=bankDepth)

         X = 0.5*(1.0 - Qbar/(topWidth*S*ck*L))         ! X factor for descreterized kinematic wave equation
         Cn = ck*dTsub/L                                         ! Courant number [-]

         C0 = (-X+Cn*(1-Y))/(1-X+Cn*(1-Y))
         C1 = (X+Cn*Y)/(1-X+Cn*(1-Y))
         C2 = (1-X-Cn*Y)/(1-X+Cn*(1-Y))

         QoutLocal(ix) = C0* QinLocal(ix)+ C1* QinLocal(ix-1)+ C2* QoutLocal(ix-1)
         QoutLocal(ix) = max(0.0_dp, QoutLocal(ix))
       else
         QoutLocal(ix) = 0._dp
       end if

       ! -- EBK 06/26/2023 -- comment out isnan check, doesn't seem to be needed.
       !if (isnan(QoutLocal(ix))) then
       !  ierr=10; message=trim(message)//'QoutLocal is Nan; activate vodose for this segment for diagnosis';return
       !end if

       if (verbose) then
         write(iulog,'(A,I3,1X,A,G12.5,1X,A,G12.5)') '   sub time-step= ',ix,'Courant number= ',Cn, 'Q= ',QoutLocal(ix)
       end if
     end do

     Q(1,1) = sum(QoutLocal(1:nTsub))/real(nTsub,kind=dp)
   else
     Q(1,1) = 0._dp
   end if
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
 if (abs(Q(1,1))>0._dp) then
   pcntReduc = min((rflux%ROUTE(idxMC)%REACH_VOL(1)/dt+Q(1,0))*0.999/Q(1,1), 1._dp)
   Q(1,1) = Q(1,1) *pcntReduc
 end if
 rflux%ROUTE(idxMC)%REACH_VOL(1) = rflux%ROUTE(idxMC)%REACH_VOL(1) + (Q(1,0)-Q(1,1))*dt
 end if

 ! if reach volume exceeds flood threshold volume, excess water is flooded volume.
 if (rflux%ROUTE(idxMC)%REACH_VOL(1) > bankVol) then
   rflux%ROUTE(idxMC)%FLOOD_VOL(1) = rflux%ROUTE(idxMC)%REACH_VOL(1) - bankVol  ! floodplain volume == overflow volume
 else
   rflux%ROUTE(idxMC)%FLOOD_VOL(1) = 0._dp
 end if
 ! compute surface water height [m]
 rflux%ROUTE(idxMC)%REACH_ELE = water_height(rflux%ROUTE(idxMC)%REACH_VOL(1)/L, bt, zc, zf=zf, bankDepth=bankDepth)

 ! add catchment flow
 rflux%ROUTE(idxMC)%REACH_Q = Q(1,1)+ Qlat

 if (verbose) then
   write(iulog,'(A,1X,G12.5)') ' Qout(t)=',Q(1,1)
 endif

 ! save inflow (index 1) and outflow (index 2) at current time step
 rstate%molecule%Q(1) = Q(1,0)
 rstate%molecule%Q(2) = Q(1,1)

 end associate

 END SUBROUTINE muskingum_cunge

END MODULE mc_route_module
