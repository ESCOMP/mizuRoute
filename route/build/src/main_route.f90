MODULE main_route_module

USE nrtype                                      ! variable types, etc.
USE dataTypes,            ONLY: STRSTA          ! state in each reach
USE dataTypes,            ONLY: STRFLX          ! fluxes in each reach
USE dataTypes,            ONLY: RCHTOPO         ! Network topology
USE dataTypes,            ONLY: RCHPRP          ! Reach parameter
USE dataTypes,            ONLY: runoff          ! runoff data type
USE dataTypes,            ONLY: subbasin_omp    ! mainstem+tributary data structures
USE obs_data,             ONLY: gageObs
USE globalData,           ONLY: routeMethods    ! Active routing method IDs
USE public_var,           ONLY: is_lake_sim     ! lake simulation flag
USE public_var,           ONLY: integerMissing  ! missing integer value
USE public_var,           ONLY: tracer          ! logical whether or not tracer is on
USE process_remap_module, ONLY: basin2reach     ! remap HRU variable to reach for volume
USE process_remap_module, ONLY: basin2reach_mass! remap HRU variable to reach for mass concentration
USE basinUH_module,       ONLY: IRF_route_basin ! perform UH convolution for basin routing
USE tracer_module,        ONLY: constituent_rch ! constituent routing per segment

implicit none

private
public::main_route

CONTAINS

 ! ******
 ! public subroutine: main HRU/reach routing routines
 ! ************************
 SUBROUTINE main_route(basinRunoff_in, &  ! input: basin (i.e.,HRU) runoff (m/s)
                       basinEvapo_in,  &  ! input: basin (i.e.,HRU) evaporation (m/s)
                       basinPrecip_in, &  ! input: basin (i.e.,HRU) precipitation (m/s)
                       basinSolute_in, &  ! input: basin (i.e.,HRU) mass flux of constituent(s) (mg/s/m2)
                       reachflux_in,   &  ! input: reach (i.e.,reach) flux (m3/s)
                       reachvol_in,    &  ! input: reach (i.e.,reach) target volume for lakes (m3)
                       ixRchProcessed, &  ! input: indices of reach to be routed
                       river_basin,    &  ! input: OMP basin decomposition
                       NETOPO_in,      &  ! input: reach topology data structure
                       RPARAM_in,      &  ! input: reach parameter data structure
                       RCHFLX_out,     &  ! inout: reach flux data structure
                       RCHSTA_out,     &  ! inout: reach state data structure
                       gage_obs_data,  &  ! inout: gauge data
                       ierr, message,  &  ! output: error control
                       enforce_min_runoff) ! optional input: true->enforce min. runoff, false->allow any runoff (e.g., negative)
   ! Details:
   ! Given HRU (basin) runoff, perform hru routing (optional) to get reach runoff, and then channel routing
   ! This routine is called from each distributed core
   ! Restriction:
   ! 1. Reach order in NETOPO_in, RPARAM_in, RCHFLX_out, RCHSTA_out must be in the same orders
   ! 2. Process a list of reach indices (in terms of NETOPO_in etc.) given by ixRchProcessed
   ! 3. basinRunoff_in is given in the order of NETOPO_in(:)%HRUIX.

   USE globalData, ONLY: TSEC                    ! beginning/ending of simulation time step [sec]
   USE globalData, ONLY: simDatetime             ! current model datetime
   USE globalData, ONLY: rch_routes              !
   USE globalData, ONLY: nTracer                 ! number of active tracers
   USE public_var, ONLY: doesBasinRoute          ! hillslope routing option
   USE public_var, ONLY: is_flux_wm              ! logical whether or not fluxes should be passed
   USE public_var, ONLY: is_vol_wm               ! logical whether or not target volume should be passed
   USE public_var, ONLY: qmodOption              ! options for streamflow modification (DA)

   implicit none
   ! argument variables
   real(dp),           allocatable, intent(in)    :: basinRunoff_in(:)    ! basin (i.e.,HRU) runoff (m/s)
   real(dp),           allocatable, intent(in)    :: basinEvapo_in(:)     ! basin (i.e.,HRU) evaporation (m/s)
   real(dp),           allocatable, intent(in)    :: basinPrecip_in(:)    ! basin (i.e.,HRU) precipitation (m/s)
   real(dp),           allocatable, intent(in)    :: basinSolute_in(:,:)  ! basin (i.e.,HRU) constituent(s) (mg/s/m2)
   real(dp),           allocatable, intent(in)    :: reachflux_in(:)      ! reach (i.e.,reach) flux (m3/s)
   real(dp),           allocatable, intent(in)    :: reachvol_in(:)       ! reach (i.e.,reach) target volume for lakes (m3)
   integer(i4b),       allocatable, intent(in)    :: ixRchProcessed(:)    ! indices of reach to be routed
   type(subbasin_omp), allocatable, intent(in)    :: river_basin(:)       ! OMP basin decomposition
   type(RCHTOPO),      allocatable, intent(in)    :: NETOPO_in(:)         ! River Network topology
   type(RCHPRP),       allocatable, intent(inout) :: RPARAM_in(:)         ! River reach parameter
   type(STRFLX),                    intent(inout) :: RCHFLX_out(:)        ! Reach fluxes (space [reaches]) for decomposed domains
   type(STRSTA),                    intent(inout) :: RCHSTA_out(:)        ! reach state data structure
   type(gageObs),                   intent(inout) :: gage_obs_data        ! gauge observation data
   integer(i4b),                    intent(out)   :: ierr                 ! error code
   character(len=strLen),           intent(out)   :: message              ! error message
   logical(lgt),       optional   , intent(in)    :: enforce_min_runoff   ! enforce minimum threshold runoff or not
   ! local variables
   character(len=strLen)                          :: cmessage             ! error message of downwind routine
   real(dp)                                       :: qobs                 ! observed discharge at a time step and site
   real(dp),           allocatable                :: reachRunoff_local(:) ! reach runoff (m3/s)
   real(dp),           allocatable                :: reachSolute_local(:,:) ! reach constituent(s) (mg/s)
   real(dp),           allocatable                :: reachEvapo_local(:)  ! reach evaporation (m3/s)
   real(dp),           allocatable                :: reachPrecip_local(:) ! reach precipitation (m3/s)
   integer(i4b)                                   :: nSeg                 ! number of reach to be processed
   integer(i4b)                                   :: iSeg                 ! index of reach
   integer(i4b)                                   :: iTracer              ! index of tracers
   integer(i4b)                                   :: ix,jx                ! loop index
   integer(i4b), allocatable                      :: reach_ix(:)
   logical(lgt)                                   :: enfrc_min_ro         ! local variable of enforce_min_runoff
   integer(i4b), parameter                        :: no_mod=0
   integer(i4b), parameter                        :: direct_insert=1

   ierr=0; message = "main_routing/"

  ! optional: enforcing accesptable minimum flow or not (default: yes==true)
  enfrc_min_ro = .true.
  if (present(enforce_min_runoff))then
    enfrc_min_ro = enforce_min_runoff
  end if

  ! number of reaches to be processed
  nSeg = size(ixRchProcessed)

  ! 1. subroutine: map the basin runoff/evapo/precip to the stream network in m3/s
  ! runoff to stream netwrok in (m3/s)
  allocate(reachRunoff_local(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [reachRunoff_local]'; return; endif

  ! Initialize water-management flux (water take, lake volume threshold for release)
  if (is_flux_wm) then
    do iSeg = 1,nSeg
      RCHFLX_out(ixRchProcessed(iSeg))%REACH_WM_FLUX = reachflux_in(iSeg)  ! added or subtracted stremflow for each reach
    end do
  else
    RCHFLX_out(:)%REACH_WM_FLUX = 0._dp
  end if
  if (is_vol_wm .and. is_lake_sim) then
    do iSeg = 1,nSeg
      RCHFLX_out(ixRchProcessed(iSeg))%REACH_WM_VOL = reachvol_in(iSeg)   ! target volume for the lakes
    end do
  else
    RCHFLX_out(:)%REACH_WM_VOL = 0._dp
  end if

  select case(qmodOption)
    case(no_mod) ! do nothing
    case(direct_insert)
      ! read gage observation [m3/s] at current time
      jx = gage_obs_data%time_ix(simDatetime(1))

      if (jx/=integerMissing) then ! there is observation at this time step
        call gage_obs_data%read_obs(ierr, cmessage, index_time=jx)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        ! put qmod at right reach
        reach_ix = gage_obs_data%link_ix()
        do ix=1,size(reach_ix) ! size(reach_ix) is number of gauges
          if (reach_ix(ix)==integerMissing) cycle
          qobs = gage_obs_data%get_obs(tix=1, six=ix)
          if ((qobs/=qobs) .or. (qobs<0)) cycle
          RCHFLX_out(reach_ix(ix))%Qobs = qobs
          RCHFLX_out(reach_ix(ix))%Qelapsed = 0
        end do
      else
        RCHFLX_out(:)%Qelapsed = RCHFLX_out(:)%Qelapsed + 1 ! change only gauge point
      end if
    case default
      ierr=1; message=trim(message)//"Error: qmodOption invalid"; return
  end select

  ! map the basin runoff to the stream network...
  call basin2reach(basinRunoff_in,     & ! input: basin runoff (m/s)
                   NETOPO_in,          & ! input: reach topology
                   RPARAM_in,          & ! input: reach parameter
                   reachRunoff_local,  & ! output: reach runoff (m3/s)
                   ierr, cmessage,     & ! output: error control
                   ixRchProcessed,     & ! optional input: indices of reach to be routed
                   enfrc_min_ro        & ! optional input: flag to enforce minimum basin runoff to runoffMin=1.e-15 [m/s]
                   )
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if(tracer) then
     allocate(reachSolute_local(nSeg,nTracer), source=0._dp, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [reachSolute_local]'; return; endif
     do iTracer=1,nTracer
       call basin2reach_mass(basinSolute_in(:,iTracer),     & ! input: basin total constitunet (mg/s/m2)
                             NETOPO_in,          & ! input: reach topology
                             RPARAM_in,          & ! input: reach parameter
                             reachSolute_local(:,iTracer),  & ! output:total constituent going to reach (mg/s)
                             ierr, cmessage,     & ! output: error control
                             ixRchProcessed)       ! optional input: indices of reach to be routed
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end do
   end if

  if (is_lake_sim) then

    allocate(reachEvapo_local(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [reachEvapo_local]'; return; endif

    ! evaporation to stream netwrok in (m3/s)
     call basin2reach(basinEvapo_in,      & ! input: basin runoff (m/s)
                      NETOPO_in,          & ! input: reach topology
                      RPARAM_in,          & ! input: reach parameter
                      reachEvapo_local,   & ! output: reach Evapo (m3/s)
                      ierr, cmessage,     & ! output: error control
                      ixRchProcessed)       ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     allocate(reachPrecip_local(nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [reachPrecip_local]'; return; endif

     ! precipitation to stream netwrok in (m3/s)
     call basin2reach(basinPrecip_in,     & ! input: basin runoff (m/s)
                      NETOPO_in,          & ! input: reach topology
                      RPARAM_in,          & ! input: reach parameter
                      reachPrecip_local,  & ! output: reach Precip (m3/s)
                      ierr, cmessage,     & ! output: error control
                      ixRchProcessed)       ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   endif

   ! 2. subroutine: basin route
   if (doesBasinRoute == 1) then
     ! instantaneous runoff volume (m3/s) and solute mass (mg/s) to data structure
     do iSeg = 1,nSeg
       RCHFLX_out(ixRchProcessed(iSeg))%BASIN_QI = reachRunoff_local(iSeg)
       if (tracer) then
         if (RCHFLX_out(ixRchProcessed(iSeg))%BASIN_QI>0) then ! this may cause mass inbalance between input and output. so may need to feed flag to land surface model
           do iTracer=1,nTracer
             RCHFLX_out(ixRchProcessed(iSeg))%BASIN_solute_inst(iTracer) = reachSolute_local(iSeg,iTracer)
           end do
         else
           RCHFLX_out(ixRchProcessed(iSeg))%BASIN_solute_inst = 0._dp
         endif
       end if
     enddo
     ! perform Basin routing
     call IRF_route_basin(NETOPO_in,         &  ! input:  reach topology
                          RCHFLX_out,        &  ! inout:  reach flux data structure
                          ierr, cmessage,    &  ! output: error controls
                          ixRchProcessed)       ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else
     ! no basin routing required (handled outside mizuRoute))
     do iSeg = 1,nSeg
       RCHFLX_out(ixRchProcessed(iSeg))%BASIN_QR(0) = RCHFLX_out(iSeg)%BASIN_QR(1)        ! streamflow from previous step
       RCHFLX_out(ixRchProcessed(iSeg))%BASIN_QR(1) = reachRunoff_local(iSeg)             ! streamflow (m3/s)
     end do

     if (tracer) then
       do iSeg = 1,nSeg
         if (RCHFLX_out(ixRchProcessed(iSeg))%BASIN_QR(1)>0) then
           do iTracer=1,nTracer
             RCHFLX_out(ixRchProcessed(iSeg))%BASIN_solute(iTracer) = reachSolute_local(iSeg,iTracer)     ! total constituent going to reach (mg/s)
           end do
         else
           RCHFLX_out(ixRchProcessed(iSeg))%BASIN_solute = 0._dp
         endif
       end do
     end if
   end if

   ! allocating precipitation and evaporation for
   if (is_lake_sim) then
     do iSeg = 1,nSeg
       RCHFLX_out(ixRchProcessed(iSeg))%Basinevapo  = reachEvapo_local(iSeg)  ! Evaporation pass to reach flux (m/s)
       RCHFLX_out(ixRchProcessed(iSeg))%BasinPrecip = reachPrecip_local(iSeg) ! precipitation pass to reach flux (m/s)
     end do
   endif

   ! 3. subroutine: river reach routing
   do ix=1,size(routeMethods)
     call route_network(rch_routes(ix)%rch_route, &  ! input: instantiated routing object
                        routeMethods(ix),         &  ! input: routing method index
                        river_basin,              &  ! input: river basin data type
                        TSEC(1), TSEC(2),         &  ! input: start and end of the time step since simulation start [sec]
                        NETOPO_in,                &  ! input: reach topology data structure
                        RPARAM_in,                &  ! input: reach parameter
                        RCHSTA_out,               &  ! inout: reach state data structure
                        RCHFLX_out,               &  ! inout: reach flux data structure
                        ierr, cmessage,           &  ! output: error controls
                        ixRchProcessed)              ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end do

 END SUBROUTINE main_route

  ! *********************************************************************
  ! private subroutine: route throughout river network
  ! *********************************************************************
  SUBROUTINE route_network(rch_route,            & ! input: instantiated routing object
                           idRoute,              & ! input: routing method id
                           river_basin,          & ! input: river basin information (mainstem, tributary outlet etc.)
                           T0,T1,                & ! input: start and end of the time step since simulation start [sec]
                           NETOPO_in,            & ! input: reach topology data structure
                           RPARAM_in,            & ! input: reach parameter data structure
                           RCHSTA_out,           & ! inout: reach state data structure
                           RCHFLX_out,           & ! inout: reach flux data structure
                           ierr,message,         & ! output: error control
                           ixSubRch)               ! optional input: subset of reach indices to be processed

    USE perf_mod,          ONLY: t_startf,t_stopf    ! timing start/stop
    USE lake_route_module, ONLY: lake_route          ! lake route module
    USE base_route,        ONLY: base_route_rch      !
    USE model_utils,       ONLY: handle_err
    USE public_var,        ONLY: accumRunoff
    USE public_var,        ONLY: kinematicWaveTracking
    USE public_var,        ONLY: impulseResponseFunc
    USE public_var,        ONLY: muskingumCunge
    USE public_var,        ONLY: kinematicWave
    USE public_var,        ONLY: diffusiveWave
    USE globalData,        ONLY: idxSUM, idxKWT, idxIRF, &
                                 idxMC, idxKW, idxDW

    implicit none
    ! Argument variables
    class(base_route_rch), intent(in),    allocatable :: rch_route
    integer(i4b),          intent(in)                 :: idRoute              ! routing method id
    type(subbasin_omp),    intent(in),    allocatable :: river_basin(:)       ! river basin information (mainstem, tributary outlet etc.)
    real(dp),              intent(in)                 :: T0,T1                ! start and end of the time step (seconds)
    type(RCHTOPO),         intent(in),    allocatable :: NETOPO_in(:)         ! River Network topology
    type(RCHPRP),          intent(inout), allocatable :: RPARAM_in(:)         ! River reach parameter
    type(STRSTA),          intent(inout)              :: RCHSTA_out(:)        ! reach state data
    type(STRFLX),          intent(inout)              :: RCHFLX_out(:)        ! Reach fluxes (space [reaches]) for decomposed domains
    integer(i4b),          intent(out)                :: ierr                 ! error code
    character(*),          intent(out)                :: message              ! error message
    integer(i4b),          intent(in), optional       :: ixSubRch(:)          ! subset of reach indices to be processed
    ! local variables
    character(len=strLen)                             :: cmessage             ! error message for downwind routine
    logical(lgt),                      allocatable    :: doRoute(:)           ! logical to indicate which reaches are processed
    integer(i4b)                                      :: idxRoute             ! routing method index
    integer(i4b)                                      :: nOrder               ! number of stream order
    integer(i4b)                                      :: nTrib                ! number of tributary basins
    integer(i4b)                                      :: nSeg                 ! number of reaches in the network
    integer(i4b)                                      :: iSeg, jSeg           ! loop indices - reach
    integer(i4b)                                      :: iTrib                ! loop indices - branch
    integer(i4b)                                      :: ix                   ! loop indices stream order

    ierr=0; message='route_network/'

    select case(idRoute)
      case(accumRunoff);           idxRoute = idxSUM
      case(kinematicWaveTracking); idxRoute = idxKWT
      case(impulseResponseFunc);   idxRoute = idxIRF
      case(muskingumCunge);        idxRoute = idxMC
      case(kinematicWave);         idxRoute = idxKW
      case(diffusiveWave);         idxRoute = idxDW
      case default
        message=trim(message)//'routing method id expect digits 0-5. Check <outOpt> in control file'; ierr=81; return
    end select

    nSeg = size(RCHFLX_out)

    ! number of reach check
    if (size(NETOPO_in)/=nSeg) then
      ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
    endif

    allocate(doRoute(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating space for [doRoute]'; return; endif

    if (present(ixSubRch))then
      doRoute(:)=.false.
      doRoute(ixSubRch) = .true. ! only subset of reaches are on
    else
      doRoute(:)=.true. ! every reach is on
    endif

    nOrder = size(river_basin)

    call t_startf('route_network') ! timing start

    ! routing through river network
    do ix = 1, nOrder
      nTrib=size(river_basin(ix)%branch)
!$OMP PARALLEL DO schedule(dynamic,1)                   & ! chunk size of 1
!$OMP          private(jSeg, iSeg)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(T0,T1)                            & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(doRoute)                          & ! data array shared
!$OMP          shared(NETOPO_in)                        & ! data structure shared
!$OMP          shared(RPARAM_in)                        & ! data structure shared
!$OMP          shared(RCHSTA_out)                       & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(ix,idxRoute)                      & ! indices shared
!$OMP          firstprivate(nTrib)
      do iTrib = 1,nTrib
        do iSeg = 1,river_basin(ix)%branch(iTrib)%nRch
          jSeg = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
          if (.not. doRoute(jSeg)) cycle
          if ((NETOPO_in(jseg)%islake).and.(is_lake_sim).and.idxRoute/=idxSUM) then
            call lake_route(jSeg,          & ! input: reach indices
                            idxRoute,      & ! input: routing method index
                            NETOPO_in,     & ! input: reach topology data structure
                            RPARAM_in,     & ! input: reach parameter data structure
                            RCHFLX_out,    & ! inout: reach flux data structure
                            ierr,cmessage)   ! output: error control
          else
            call rch_route%route(jSeg,           & ! input: array indices
                                 T0,T1,          & ! input: start and end of the time step
                                 NETOPO_in,      & ! input: reach topology data structure
                                 RPARAM_in,      & ! input: reach parameter data structure
                                 RCHSTA_out,     & ! inout: reach state data structure
                                 RCHFLX_out,     & ! inout: reach flux data structure
                                 ierr,cmessage)    ! output: error control
          end if
          if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
          if (tracer .and. idxRoute/=idxSUM) then
            call constituent_rch(jSeg,          & ! input: index of reach to be processed
                                 idxRoute,      & ! input: routing method index
                                 NETOPO_in,     & ! input: reach topology data structure
                                 RPARAM_in,     & ! input: reach parameter data structure
                                 RCHSTA_out,    & ! inout: reach state data structure
                                 RCHFLX_out,    & ! inout: reach flux data structure
                                 ierr, message)   ! output: error control
            if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
          end if
        end do ! reach index
      end do ! tributary
!$OMP END PARALLEL DO
    end do ! basin loop

    call t_stopf('route_network')

  END SUBROUTINE route_network

END MODULE main_route_module
