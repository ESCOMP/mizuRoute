MODULE main_route_module

USE nrtype                                   ! variable types, etc.
! data structures
USE dataTypes, ONLY: STRSTA                  ! state in each reach
USE dataTypes, ONLY: STRFLX                  ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO                 ! Network topology
USE dataTypes, ONLY: RCHPRP                  ! Reach parameter
USE dataTypes, ONLY: runoff                  ! runoff data type
USE dataTypes, ONLY: subbasin_omp            ! mainstem+tributary data structures
! mapping HRU runoff to reach
USE process_remap_module, ONLY: basin2reach
! subroutines: basin routing
USE basinUH_module, ONLY: IRF_route_basin    ! perform UH convolution for basin routing
! subroutines: river routing
USE accum_runoff_module, ONLY: accum_runoff  ! upstream flow accumulation
USE irf_route_module,    ONLY: irf_route     ! unit hydrograph (impulse response function) routing method
USE kwt_route_module,    ONLY: kwt_route     ! lagrangian kinematic wave routing method
USE kw_route_module,     ONLY: kw_route      ! kinematic wave routing method
USE mc_route_module,     ONLY: mc_route      ! muskingum-cunge routing method
USE dfw_route_module,    ONLY: dfw_route     ! diffusive wave routing method

implicit none

private
public::main_route

CONTAINS

 ! ******
 ! public subroutine: main HRU/reach routing routines
 ! ************************
 SUBROUTINE main_route(iens,           &  ! input: ensemble index
                       basinRunoff_in, &  ! input: basin (i.e.,HRU) runoff (m/s)
                       basinEvapo_in,  &  ! input: basin (i.e.,HRU) evaporation (m/s)
                       basinPrecip_in, &  ! input: basin (i.e.,HRU) precipitation (m/s)
                       reachflux_in,   &  ! input: reach (i.e.,reach) flux (m3/s)
                       reachvol_in,    &  ! input: reach (i.e.,reach) target volume for lakes (m3)
                       ixRchProcessed, &  ! input: indices of reach to be routed
                       river_basin,    &  ! input: OMP basin decomposition
                       NETOPO_in,      &  ! input: reach topology data structure
                       RPARAM_in,      &  ! input: reach parameter data structure
                       ixDesire,       &  ! input: index of verbose reach
                       RCHFLX_out,     &  ! inout: reach flux data structure
                       RCHSTA_out,     &  ! inout: reach state data structure
                       ierr, message)     ! output: error control
   ! Details:
   ! Given HRU (basin) runoff, perform hru routing (optional) to get reach runoff, and then channel routing
   ! Restriction:
   ! 1. Reach order in NETOPO_in, RPARAM_in, RCHFLX_out, RCHSTA_out must be in the same orders
   ! 2. Process a list of reach indices (in terms of NETOPO_in etc.) given by ixRchProcessed
   ! 3. basinRunoff_in is given in the order of NETOPO_in(:)%HRUIX.

   USE public_var, ONLY: doesBasinRoute
   USE public_var, ONLY: accumRunoff
   USE public_var, ONLY: impulseResponseFunc
   USE public_var, ONLY: kinematicWaveTracking
   USE public_var, ONLY: kinematicWave
   USE public_var, ONLY: muskingumCunge
   USE public_var, ONLY: diffusiveWave
   USE globalData, ONLY: onRoute                 ! logical to indicate which routing method(s) is on
   USE globalData, ONLY: TSEC                    ! beginning/ending of simulation time step [sec]
   USE public_var, ONLY: is_lake_sim             ! logical whether or not lake should be simulated
   USE public_var, ONLY: is_flux_wm              ! logical whether or not fluxes should be passed
   USE public_var, ONLY: is_vol_wm               ! logical whether or not target volume should be passed

   implicit none
   ! argument variables
   integer(i4b),                    intent(in)    :: iens                 ! ensemble member
   real(dp),           allocatable, intent(in)    :: basinRunoff_in(:)    ! basin (i.e.,HRU) runoff (m/s)
   real(dp),           allocatable, intent(in)    :: basinEvapo_in(:)     ! basin (i.e.,HRU) evaporation (m/s)
   real(dp),           allocatable, intent(in)    :: basinPrecip_in(:)    ! basin (i.e.,HRU) precipitation (m/s)
   real(dp),           allocatable, intent(in)    :: reachflux_in(:)      ! reach (i.e.,reach) flux (m3/s)
   real(dp),           allocatable, intent(in)    :: reachvol_in(:)       ! reach (i.e.,reach) target volume for lakes (m3)
   integer(i4b),       allocatable, intent(in)    :: ixRchProcessed(:)    ! indices of reach to be routed
   type(subbasin_omp), allocatable, intent(in)    :: river_basin(:)       ! OMP basin decomposition
   type(RCHTOPO),      allocatable, intent(in)    :: NETOPO_in(:)         ! River Network topology
   type(RCHPRP),       allocatable, intent(inout) :: RPARAM_in(:)         ! River reach parameter
   integer(i4b),                    intent(in)    :: ixDesire             ! index of the reach for verbose output
   type(STRFLX),       allocatable, intent(inout) :: RCHFLX_out(:,:)      ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   type(STRSTA),       allocatable, intent(inout) :: RCHSTA_out(:,:)      ! reach state data structure
   integer(i4b),                    intent(out)   :: ierr                 ! error code
   character(len=strLen),           intent(out)   :: message              ! error message
   ! local variables
   character(len=strLen)                          :: cmessage             ! error message of downwind routine
   real(dp)                                       :: T0,T1                ! beginning/ending of simulation time step [sec]
   real(dp),           allocatable                :: reachRunoff_local(:) ! reach runoff (m/s)
   real(dp),           allocatable                :: reachEvapo_local(:)  ! reach evaporation (m/s)
   real(dp),           allocatable                :: reachPrecip_local(:) ! reach precipitation (m/s)
   integer(i4b)                                   :: nSeg                 ! number of reach to be processed
   integer(i4b)                                   :: iSeg                 ! index of reach

   ierr=0; message = "main_routing/"

  ! define the start and end of the time step
  T0=TSEC(0); T1=TSEC(1)

  ! number of reaches to be processed
  nSeg = size(ixRchProcessed)

  ! 1. subroutine: map the basin runoff/evapo/precip to the stream network in m3/s
  ! runoff to stream netwrok in (m3/s)
  allocate(reachRunoff_local(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [reachRunoff_local]'; return; endif

  ! passing of the water management fluxes and lake target vol if presence
  if (is_flux_wm) then
    do iSeg = 1,nSeg
      RCHFLX_out(iens,ixRchProcessed(iSeg))%REACH_WM_FLUX  =  reachflux_in(iSeg)  ! added or subtracted stremflow for each reach
    end do
  end if
  if (is_vol_wm.and.is_lake_sim) then
    do iSeg = 1,nSeg
      RCHFLX_out(iens,ixRchProcessed(iSeg))%REACH_WM_VOL   =  reachvol_in(iSeg)   ! target volume for the lakes
    end do
  end if

  ! 1. subroutine: map basin runoff to river network HRUs
  ! map the basin runoff to the stream network...
  call basin2reach(basinRunoff_in,     & ! input: basin runoff (m/s)
                   NETOPO_in,          & ! input: reach topology
                   RPARAM_in,          & ! input: reach parameter
                   reachRunoff_local,  & ! output: reach runoff (m3/s)
                   ierr, cmessage,     & ! output: error control
                   ixRchProcessed)       ! optional input: indices of reach to be routed
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (is_lake_sim) then

    ! allocate evaporation
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

     ! allocate precipitaiton
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
     ! instantaneous runoff volume (m3/s) to data structure
     do iSeg = 1,nSeg
       RCHFLX_out(iens,ixRchProcessed(iSeg))%BASIN_QI = reachRunoff_local(iSeg)
     enddo
     ! perform Basin routing
     call IRF_route_basin(iens,              &  ! input:  ensemble index
                          NETOPO_in,         &  ! input:  reach topology
                          RCHFLX_out,        &  ! inout:  reach flux data structure
                          ierr, cmessage,    &  ! output: error controls
                          ixRchProcessed)       ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else
     ! no basin routing required (handled outside mizuRoute))
     do iSeg = 1,nSeg
       RCHFLX_out(iens,ixRchProcessed(iSeg))%BASIN_QR(0) = RCHFLX_out(iens,iSeg)%BASIN_QR(1)   ! streamflow from previous step
       RCHFLX_out(iens,ixRchProcessed(iSeg))%BASIN_QR(1) = reachRunoff_local(iSeg)             ! streamflow (m3/s)
     end do
   end if

   ! allocating precipitation and evaporation for
   if (is_lake_sim) then
     do iSeg = 1,nSeg
       RCHFLX_out(iens,ixRchProcessed(iSeg))%Basinevapo  = reachEvapo_local(iSeg)  ! Evaporation pass to reach flux (m/s)
       RCHFLX_out(iens,ixRchProcessed(iSeg))%BasinPrecip = reachPrecip_local(iSeg) ! precipitation pass to reach flux (m/s)
     end do
   endif

   ! 3. subroutine: river reach routing
   ! perform upstream flow accumulation
   if (onRoute(accumRunoff)) then
     call accum_runoff(iens,              &  ! input: ensemble index
                       river_basin,       &  ! input: river basin data type
                       ixDesire,          &  ! input: index of verbose reach
                       NETOPO_in,         &  ! input: reach topology data structure
                       RCHFLX_out,        &  ! inout: reach flux data structure
                       ierr, cmessage,    &  ! output: error controls
                       ixRchProcessed)       ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   endif

   if (onRoute(kinematicWaveTracking)) then
     call kwt_route(iens,                 & ! input: ensemble index
                    river_basin,          & ! input: river basin data type
                    T0,T1,                & ! input: start and end of the time step
                    ixDesire,             & ! input: index of verbose reach
                    NETOPO_in,            & ! input: reach topology data structure
                    RPARAM_in,            & ! input: reach parameter data structure
                    RCHSTA_out,           & ! inout: reach state data structure
                    RCHFLX_out,           & ! inout: reach flux data structure
                    ierr,cmessage,        & ! output: error control
                    ixRchProcessed)         ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   if (onRoute(impulseResponseFunc)) then
     call irf_route(iens,                & ! input: ensemble index
                    river_basin,         & ! input: river basin data type
                    ixDesire,            & ! input: index of verbose reach
                    NETOPO_in,           & ! input: reach topology data structure
                    RPARAM_in,           & ! input: reach parameter data structure
                    RCHFLX_out,          & ! inout: reach flux data structure
                    ierr,cmessage,       & ! output: error control
                    ixRchProcessed)        ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   if (onRoute(kinematicWave)) then
     call kw_route(iens,                 & ! input: ensemble index
                   river_basin,          & ! input: river basin data type
                   T0,T1,                & ! input: start and end of the time step
                   ixDesire,             & ! input: index of verbose reach
                   NETOPO_in,            & ! input: reach topology data structure
                   RPARAM_in,            & ! input: reach parameter data structure
                   RCHSTA_out,           & ! inout: reach state data structure
                   RCHFLX_out,           & ! inout: reach flux data structure
                   ierr,cmessage,        & ! output: error control
                   ixRchProcessed)         ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   if (onRoute(muskingumCunge)) then
     call mc_route(iens,                 & ! input: ensemble index
                   river_basin,          & ! input: river basin data type
                   T0,T1,                & ! input: start and end of the time step
                   ixDesire,             & ! input: index of verbose reach
                   NETOPO_in,            & ! input: reach topology data structure
                   RPARAM_in,            & ! input: reach parameter data structure
                   RCHSTA_out,           & ! inout: reach state data structure
                   RCHFLX_out,           & ! inout: reach flux data structure
                   ierr,cmessage,        & ! output: error control
                   ixRchProcessed)         ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   if (onRoute(diffusiveWave)) then
     call dfw_route(iens,                 & ! input: ensemble index
                    river_basin,          & ! input: river basin data type
                    T0,T1,                & ! input: start and end of the time step
                    ixDesire,             & ! input: index of verbose reach
                    NETOPO_in,            & ! input: reach topology data structure
                    RPARAM_in,            & ! input: reach parameter data structure
                    RCHSTA_out,           & ! inout: reach state data structure
                    RCHFLX_out,           & ! inout: reach flux data structure
                    ierr,cmessage,        & ! output: error control
                    ixRchProcessed)         ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

 END SUBROUTINE main_route

END MODULE main_route_module
