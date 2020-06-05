MODULE main_route_module

! variable types
USE nrtype                                   ! variable types, etc.

! data structures
USE dataTypes, ONLY: STRSTA                  ! state in each reach
USE dataTypes, ONLY: STRFLX                  ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO                 ! Network topology
USE dataTypes, ONLY: RCHPRP                  ! Reach parameter
USE dataTypes, ONLY: runoff                  ! runoff data type
USE dataTypes, ONLY: subbasin_omp            ! mainstem+tributary data structures

! mapping HRU runoff to reach
USE remapping, ONLY: basin2reach

! subroutines: basin routing
USE basinUH_module, ONLY: IRF_route_basin    ! perform UH convolution for basin routing

! subroutines: river routing
USE accum_runoff_module, ONLY: accum_runoff  ! upstream flow accumulation
USE kwt_route_module,    ONLY: kwt_route     ! Lagrangian kinematic wave routing method
USE kwe_route_module,    ONLY: kwe_route     ! Eulerian kinematic wave routing method
USE irf_route_module,    ONLY: irf_route     ! unit hydrograph (impulse response function) routing method

implicit none

private

public::main_route

CONTAINS

 ! ******
 ! public subroutine: main HRU/reach routing routines
 ! ************************
 SUBROUTINE main_route(&
                       ! input
                       iens,           &  ! ensemble index
                       basinRunoff_in, &  ! basin (i.e.,HRU) runoff (m/s)
                       ixRchProcessed, &  ! indices of reach to be routed
                       river_basin,    &  ! OMP basin decomposition
                       NETOPO_in,      &  ! reach topology data structure
                       RPARAM_in,      &  ! reach parameter data structure
                       ixDesire,       &  ! index of verbose reach
                       ! inout
                       RCHFLX_out,     &  ! reach flux data structure
                       RCHSTA_out,     &  ! reach state data structure
                       ! output: error handling
                       ierr, message)     ! output: error control

   ! Details:
   ! Given HRU (basin) runoff, perform hru routing (optional) to get reach runoff, and then channel routing
   ! Restriction:
   ! 1. Reach order in NETOPO_in, RPARAM_in, RCHFLX_out, RCHSTA_out must be in the same orders
   ! 2. Process a list of reach indices (in terms of NETOPO_in etc.) given by ixRchProcessed
   ! 3. basinRunoff_in is given in the order of NETOPO_in(:)%HRUIX.

   USE public_var, ONLY: routOpt
   USE public_var, ONLY: doesBasinRoute
   USE public_var, ONLY: doesAccumRunoff
   USE public_var, ONLY: allRoutingMethods
   USE public_var, ONLY: kinematicWave
   USE public_var, ONLY: kinematicWaveEuler
   USE public_var, ONLY: impulseResponseFunc
   USE globalData, ONLY: TSEC                    ! beginning/ending of simulation time step [sec]

   implicit none

   ! input
   integer(i4b),                    intent(in)    :: iens                 ! ensemble member
   real(dp),           allocatable, intent(in)    :: basinRunoff_in(:)    ! basin (i.e.,HRU) runoff (m/s)
   integer(i4b),       allocatable, intent(in)    :: ixRchProcessed(:)    ! indices of reach to be routed
   type(subbasin_omp), allocatable, intent(in)    :: river_basin(:)       ! OMP basin decomposition
   type(RCHTOPO),      allocatable, intent(in)    :: NETOPO_in(:)         ! River Network topology
   type(RCHPRP),       allocatable, intent(in)    :: RPARAM_in(:)         ! River reach parameter
   integer(i4b),                    intent(in)    :: ixDesire             ! index of the reach for verbose output
   ! inout
   type(STRFLX),       allocatable, intent(inout) :: RCHFLX_out(:,:)      ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   type(STRSTA),       allocatable, intent(inout) :: RCHSTA_out(:,:)      ! reach state data structure
   ! output
   integer(i4b),                    intent(out)   :: ierr                 ! error code
   character(len=strLen),           intent(out)   :: message              ! error message
   ! local variables
   character(len=strLen)                          :: cmessage             ! error message of downwind routine
   real(dp)                                       :: T0,T1                ! beginning/ending of simulation time step [sec]
   real(dp),           allocatable                :: reachRunoff_local(:) ! reach runoff (m/s)
   integer(i4b)                                   :: nSeg                 ! number of reach to be processed
   integer(i4b)                                   :: iSeg                 ! index of reach

   ! initialize errors
   ierr=0; message = "main_routing/"

  ! define the start and end of the time step
  T0=TSEC(0); T1=TSEC(1)

  ! number of reaches to be processed
  nSeg = size(ixRchProcessed)

  allocate(reachRunoff_local(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [reachRunoff_local]'; return; endif

  ! 1. subroutine: map basin runoff to river network HRUs
  ! map the basin runoff to the stream network...
  call basin2reach(basinRunoff_in,     & ! input: basin runoff (m/s)
                   NETOPO_in,          & ! input: reach topology
                   RPARAM_in,          & ! input: reach parameter
                   reachRunoff_local,  & ! output: reach runoff (m3/s)
                   ierr, cmessage,     & ! output: error control
                   ixRchProcessed)       ! optional input: indices of reach to be routed
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! 2. subroutine: basin route
  if (doesBasinRoute == 1) then
    ! instantaneous runoff volume (m3/s) to data structure
    do iSeg = 1,nSeg
     RCHFLX_out(iens,ixRchProcessed(iSeg))%BASIN_QI = reachRunoff_local(iSeg)
    enddo
    ! perform Basin routing
    call IRF_route_basin(iens,              &  ! input:  ensemble index
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

  ! 3. subroutine: river reach routing
   ! perform upstream flow accumulation
   if (doesAccumRunoff == 1) then
     call accum_runoff(iens,              &  ! input: ensemble index
                       river_basin,       &  ! input: river basin data type
                       ixDesire,          &  ! input: index of verbose reach
                       NETOPO_in,         &  ! input: reach topology data structure
                       RCHFLX_out,        &  ! inout: reach flux data structure
                       ierr, cmessage,    &  ! output: error controls
                       ixRchProcessed)       ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   endif

   ! perform Lagrangian KWT routing
   if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
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
   endif

   ! perform Eulerian KWT routing
   if (routOpt==kinematicWaveEuler) then
    call kwe_route(iens,                 & ! input: ensemble index
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

   ! perform IRF routing
   if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
    call irf_route(iens,                & ! input: ensemble index
                   river_basin,         & ! input: river basin data type
                   ixDesire,            & ! input: index of verbose reach
                   NETOPO_in,           & ! input: reach topology data structure
                   RCHFLX_out,          & ! inout: reach flux data structure
                   ierr,cmessage,       & ! output: error control
                   ixRchProcessed)        ! optional input: indices of reach to be routed
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   endif

 END SUBROUTINE main_route

END MODULE main_route_module
