module main_route_module

! variable types
USE nrtype                                                  ! variable types, etc.

! data structures
USE dataTypes,           only : KREACH                      ! collection of particles in a given reach
USE dataTypes,           only : STRFLX                      ! fluxes in each reach
USE dataTypes,           only : RCHTOPO                     ! Network topology
USE dataTypes,           only : RCHPRP                      ! Reach parameter
USE dataTypes,           only : runoff                      ! runoff data type
USE dataTypes,           only : basin                       ! river basin data type

! subroutines: general utility
USE nr_utility_module,   only : findIndex                   ! find index within a vector

! mapping HRU runoff to reach
USE remapping,           only : basin2reach
! subroutines: basin routing
USE basinUH_module,      only : IRF_route_basin             ! perform UH convolution for basin routing

! subroutines: river routing
USE accum_runoff_module, only : accum_runoff                ! upstream flow accumulation
!USE kwt_route_module,    only : kwt_route                  ! kinematic wave routing method
!USE irf_route_module,    only : irf_route                  ! unit hydrograph (impulse response function) routing method
USE kwt_route_module,    only : kwt_route => kwt_route_orig ! kinematic wave routing method
USE irf_route_module,    only : irf_route => irf_route_orig ! river network unit hydrograph method

! utilities
USE nr_utility_module, ONLY : arth

implicit none

private

public::serial_route
public::main_route

contains

 ! ******
 ! public subroutine:
 ! ************************
 subroutine serial_route(iens, ierr, message)
  ! shared data
  USE public_var
  USE dataTypes,  only : KREACH           ! derived data type
  USE globalData, only : NETOPO           ! entire river reach netowrk topology structure
  USE globalData, only : RPARAM           ! entire river reach parameter structure
  USE globalData, only : RCHFLX           ! entire reach flux structure
  USE globalData, only : KROUTE           ! entire river reach kwt sate structure
  USE globalData, only : runoff_data      ! runoff data structure
  USE globalData, only : nRch             ! number of reaches in the whoel river network

  implicit none
  ! input variables
  integer(i4b),          intent(in)  :: iens
  ! Output variables
  integer(i4b),          intent(out) :: ierr
  character(len=strLen), intent(out) :: message ! error message
  ! Local variables
  integer(i4b), allocatable          :: ixRchProcessed(:)    ! reach indice list to be processed
  integer(i4b)                       :: mainstem=2           ! basin indicator (mainstem = 2 )
  character(len=strLen)              :: cmessage             ! error message of downwind routine

  ierr=0; message='serial_route/'

  allocate(ixRchProcessed(nRch), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating arrays for []ixRchProcessed'; return; endif
  ixRchProcessed = arth(1,1,nRch)

  ! Perform routing
  call main_route(iens,                    &  ! ensemble index
                  runoff_data%basinRunoff, &  ! basin (i.e.,HRU) runoff (m/s)
                  ixRchProcessed,          &  ! indices of reach to be routed
                  mainstem,                &  ! basinType (1-> tributary, 2->mainstem)
                  NETOPO,                  &  ! reach topology data structure
                  RPARAM,                  &  ! reach parameter data structure
                  ! inout
                  RCHFLX,                  &  ! reach flux data structure
                  KROUTE,                  &  ! reach state data structure
                  ! output: error handling
                  ierr, cmessage)     ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine serial_route


 ! ******
 ! public subroutine: main HRU/reach routing routines
 ! ************************
 subroutine main_route(&
                       ! input
                       iens,           &  ! ensemble index
                       basinRunoff_in, &  ! basin (i.e.,HRU) runoff (m/s)
                       ixRchProcessed, &  ! indices of reach to be routed
                       river_basin,    &  ! OMP basin decomposition
                       basinType,      &  ! basinType (1-> tributary, 2->mainstem)
                       NETOPO_in,      &  ! reach topology data structure
                       RPARAM_in,      &  ! reach parameter data structure
                       ! inout
                       RCHFLX_out,     &  ! reach flux data structure
                       KROUTE_out,     &  ! reach state data structure
                       ! output: error handling
                       ierr, message)     ! output: error control
   ! Details:
   ! Given HRU (basin) runoff, perform hru routing (optional) to get reach runoff, and then channel routing
   ! Restriction:
   ! 1. Reach order in NETOPO_in, RPARAM_in, RCHFLX_out, KROUTE_out must be in the same orders
   ! 2. Process a list of reach indices (in terms of NETOPO_in etc.) given by ixRchProcessed
   ! 3. basinRunoff_in is given in the order of NETOPO_in(:)%HRUIX.

   ! shared data
   USE public_var, only : routOpt
   USE public_var, only : doesBasinRoute
   USE public_var, only : allRoutingMethods
   USE public_var, only : kinematicWave
   USE public_var, only : impulseResponseFunc
   USE globalData, only : TSEC                    ! beginning/ending of simulation time step [sec]
   USE globalData, only : ixPrint                 ! desired reach index to be on-screen print

   implicit none

   ! input
   integer(i4b),               intent(in)    :: iens                 ! ensemble member
   real(dp),      allocatable, intent(in)    :: basinRunoff_in(:)    ! basin (i.e.,HRU) runoff (m/s)
   integer(i4b),  allocatable, intent(in)    :: ixRchProcessed(:)    ! indices of reach to be routed
   type(basin),   allocatable, intent(in)    :: river_basin(:)       ! OMP basin decomposition
   integer(i4b),               intent(in)    :: basinType            ! basinType (1-> tributary, 2->mainstem)
   type(RCHTOPO), allocatable, intent(in)    :: NETOPO_in(:)         ! River Network topology
   type(RCHPRP),  allocatable, intent(in)    :: RPARAM_in(:)         ! River reach parameter
   ! inout
   TYPE(STRFLX),  allocatable, intent(inout) :: RCHFLX_out(:,:)      ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   TYPE(KREACH),  allocatable, intent(inout) :: KROUTE_out(:,:)      ! reach state data structure
   ! output
   integer(i4b),               intent(out)   :: ierr                 ! error code
   character(len=strLen),      intent(out)   :: message              ! error message
   ! local variables
   character(len=strLen)                     :: cmessage             ! error message of downwind routine
   real(dp)                                  :: T0,T1                ! beginning/ending of simulation time step [sec]
   real(dp),      allocatable                :: reachRunoff_local(:) ! reach runoff (m/s)
   integer(i4b)                              :: nSeg                 ! number of reach to be processed
   integer(i4b)                              :: iSeg                 ! index of reach

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
  call basin2reach(&
                  ! input
                  basinRunoff_in,     & ! basin runoff (m/s)
                  NETOPO_in,          & ! reach topology
                  RPARAM_in,          & ! reach parameter
                  ! output
                  reachRunoff_local,  & ! intent(out): reach runoff (m3/s)
                  ierr, cmessage,     & ! intent(out): error control
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
   call accum_runoff(iens,              &  ! input: ensemble index
                     ixPrint,           &  ! input: index of verbose reach
                     NETOPO_in,         &  ! input: reach topology data structure
                     RCHFLX_out,        &  ! inout: reach flux data structure
                     ierr, cmessage,    &  ! output: error controls
                     ixRchProcessed)       ! optional input: indices of reach to be routed
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! perform KWT routing
   if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
    call kwt_route(iens,                 & ! input: ensemble index
                   river_basin,          & ! input: river basin data type
                   T0,T1,                & ! input: start and end of the time step
                   basinType,            & ! input: basinType (0-> tributary, 1->mainstem)
                   ixPrint,              & ! input: index of the desired reach
                   NETOPO_in,            & ! input: reach topology data structure
                   RPARAM_in,            & ! input: reach parameter data structure
                   KROUTE_out,           & ! inout: reach state data structure
                   RCHFLX_out,           & ! inout: reach flux data structure
                   ierr,cmessage,        & ! output: error control
                   ixRchProcessed)         ! optional input: indices of reach to be routed
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   endif

   ! perform IRF routing
   if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
    call irf_route(iens,                & ! input: ensemble index
                   river_basin,         & ! input: river basin data type
                   ixPrint,             & ! input: index of the desired reach
                   NETOPO_in,           & ! input: reach topology data structure
                   RCHFLX_out,          & ! inout: reach flux data structure
                   ierr,cmessage,       & ! output: error control
                   ixRchProcessed)        ! optional input: indices of reach to be routed
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   endif

 end subroutine main_route


end module main_route_module
