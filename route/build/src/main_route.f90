module main_route_module

USE nrtype                                                 ! variable types, etc.
! mapping HRU runoff to reach
USE remapping,           only : basin2reach
! subroutines: basin routing
USE basinUH_module,      only : IRF_route_basin            ! perform UH convolution for basin routing
! subroutines: river routing
USE accum_runoff_module, only : accum_runoff               ! upstream flow accumulation
USE kwt_route_module,    only : kwt_route                  ! lagrangian kinematic wave routing method
USE irf_route_module,    only : irf_route                  ! unit hydrograph (impulse response function) routing method
USE dfw_route_module,    only : dfw_route                  ! diffusive wave routing method
USE kw_route_module,     only : kw_route                   ! kinematic wave routing method
USE mc_route_module,     only : mc_route                   ! muskingum-cunge routing method

implicit none

private

public::main_route

contains

 ! ******
 ! public subroutine: main HRU/reach routing routines
 ! ************************
 subroutine main_route(iens,           &  ! input:  ensemble index
                       ierr, message)     ! output: error control
   ! Details:
   ! Given HRU (basin) runoff, perform hru routing (optional) to get reach runoff, and then channel routing ! Restriction:
   ! 1. Reach order in NETOPO, RPARAM, RCHFLX, RCHSTA must be in the same orders
   ! 2. Process a list of reach indices (in terms of NETOPO etc.) given by ixRchProcessed
   ! 3. basinRunoff_in is given in the order of NETOPO(:)%HRUIX.

   USE public_var, ONLY: doesBasinRoute
   USE public_var, ONLY: accumRunoff
   USE public_var, ONLY: impulseResponseFunc
   USE public_var, ONLY: kinematicWaveTracking
   USE public_var, ONLY: kinematicWave
   USE public_var, ONLY: muskingumCunge
   USE public_var, ONLY: diffusiveWave
   USE globalData, ONLY: onRoute                 ! logical to indicate which routing method(s) is on
   USE globalData, ONLY: TSEC                    ! beginning/ending of simulation time step [sec]
   USE globalData, ONLY: ixPrint                 ! desired reach index to be on-screen print

   USE globalData, ONLY: NETOPO           ! entire river reach netowrk topology structure
   USE globalData, ONLY: RPARAM           ! entire river reach parameter structure
   USE globalData, ONLY: RCHFLX           ! entire reach flux structure
   USE globalData, ONLY: RCHSTA           ! entire river reach restart structure
   USE globalData, ONLY: runoff_data      ! runoff data structure
   USE globalData, ONLY: river_basin      ! OMP basin decomposition
   USE globalData, ONLY: nRch             ! number of reaches in the whoel river network

   implicit none

   ! input
   integer(i4b),               intent(in)    :: iens                 ! ensemble member
   ! output
   integer(i4b),               intent(out)   :: ierr                 ! error code
   character(len=strLen),      intent(out)   :: message              ! error message
   ! local variables
   character(len=strLen)                     :: cmessage             ! error message of downwind routine
   real(dp)                                  :: T0,T1                ! beginning/ending of simulation time step [sec]
   real(dp),      allocatable                :: reachRunoff_local(:) ! reach runoff (m/s)
   integer(i4b)                              :: iSeg                 ! index of reach
   ! timing variables
   integer*8                                 :: cr, startTime, endTime ! rate, start and end time stamps
   real(dp)                                  :: elapsedTime            ! elapsed time for the process

  CALL system_clock(count_rate=cr)

  ! initialize errors
  ierr=0; message = "main_routing/"

  ! define the start and end of the time step
  T0=TSEC(0); T1=TSEC(1)

  allocate(reachRunoff_local(nRch), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [reachRunoff_local]'; return; endif

  ! 1. subroutine: map basin runoff to river network HRUs
  ! map the basin runoff to the stream network...
  call system_clock(startTime)
  call basin2reach(runoff_data%basinRunoff, & ! input:  basin runoff (m/s)
                   NETOPO,                  & ! input:  reach topology
                   RPARAM,                  & ! input:  reach parameter
                   reachRunoff_local,       & ! output: reach runoff (m3/s)
                   ierr, cmessage)            ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  call system_clock(endTime)
  elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
  write(*,"(A,1PG15.7,A)") '      elapsed-time [basin2reach] = ', elapsedTime, ' s'

  ! 2. subroutine: basin route
  if (doesBasinRoute == 1) then
    call system_clock(startTime)
    ! perform Basin routing
    call IRF_route_basin(iens,              &  ! input:  ensemble index
                         reachRunoff_local, &  ! input:  instantaneous Runoff volume [m3/s] flowing into reach
                         RCHFLX,            &  ! inout:  reach flux data structure
                         ierr, cmessage)       ! output: error controls
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [IRF_route_basin] = ', elapsedTime, ' s'
  else
    ! no basin routing required (handled outside mizuRoute))
    do iSeg = 1,nRch
    RCHFLX(iens,iSeg)%BASIN_QR(0) = RCHFLX(iens,iSeg)%BASIN_QR(1)   ! streamflow from previous step
    RCHFLX(iens,iSeg)%BASIN_QR(1) = reachRunoff_local(iSeg)         ! streamflow (m3/s)
    end do
  end if

  ! 3. subroutine: river reach routing
  if (onRoute(accumRunoff)) then
    call system_clock(startTime)
    call accum_runoff(iens,              &  ! input: ensemble index
                      river_basin,       &  ! input: river basin data type
                      ixPrint,           &  ! input: index of verbose reach
                      NETOPO,            &  ! input: reach topology data structure
                      RCHFLX,            &  ! inout: reach flux data structure
                      ierr, cmessage)       ! output: error controls
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [accum_runoff] = ', elapsedTime, ' s'
  endif

  if (onRoute(kinematicWaveTracking)) then
    call system_clock(startTime)
    call kwt_route(iens,                 & ! input: ensemble index
                   river_basin,          & ! input: river basin data type
                   T0,T1,                & ! input: start and end of the time step
                   ixPrint,              & ! input: index of the desired reach
                   NETOPO,               & ! input: reach topology data structure
                   RPARAM,               & ! input: reach parameter data structure
                   RCHSTA,               & ! inout: reach state data structure
                   RCHFLX,               & ! inout: reach flux data structure
                   ierr,cmessage)          ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [kwt_route] = ', elapsedTime, ' s'
  end if

  if (onRoute(impulseResponseFunc)) then
    call system_clock(startTime)
    call irf_route(iens,                & ! input: ensemble index
                   river_basin,         & ! input: river basin data type
                   ixPrint,             & ! input: index of the desired reach
                   NETOPO,              & ! input: reach topology data structure
                   RPARAM,              & ! input: reach parameter data structure
                   RCHFLX,              & ! inout: reach flux data structure
                   ierr,cmessage)         ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [irf_route] = ', elapsedTime, ' s'
  endif

  if (onRoute(kinematicWave)) then
    call system_clock(startTime)
    call kw_route(iens,                & ! input: ensemble index
                  river_basin,         & ! input: river basin data type
                  T0,T1,               & ! input: start and end of the time step
                  ixPrint,             & ! input: index of the desired reach
                  NETOPO,              & ! input: reach topology data structure
                  RPARAM,              & ! input: reach parameter data structure
                  RCHSTA,              & ! inout: reach state data structure
                  RCHFLX,              & ! inout: reach flux data structure
                  ierr,cmessage)         ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [kw_route] = ', elapsedTime, ' s'
  endif

  if (onRoute(muskingumCunge)) then
    call system_clock(startTime)
    call mc_route(iens,                & ! input: ensemble index
                  river_basin,         & ! input: river basin data type
                  T0,T1,               & ! input: start and end of the time step
                  ixPrint,             & ! input: index of the desired reach
                  NETOPO,              & ! input: reach topology data structure
                  RPARAM,              & ! input: reach parameter data structure
                  RCHSTA,              & ! inout: reach state data structure
                  RCHFLX,              & ! inout: reach flux data structure
                  ierr,cmessage)         ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [mc_route] = ', elapsedTime, ' s'
  endif

  if (onRoute(diffusiveWave)) then
    call system_clock(startTime)
    call dfw_route(iens,                & ! input: ensemble index
                   river_basin,         & ! input: river basin data type
                   T0,T1,               & ! input: start and end of the time step
                   ixPrint,             & ! input: index of the desired reach
                   NETOPO,              & ! input: reach topology data structure
                   RPARAM,              & ! input: reach parameter data structure
                   RCHSTA,              & ! inout: reach state data structure
                   RCHFLX,              & ! inout: reach flux data structure
                   ierr,cmessage)         ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [dfw_route] = ', elapsedTime, ' s'
  endif

 end subroutine main_route

end module main_route_module
