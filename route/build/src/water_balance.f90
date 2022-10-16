MODULE water_balance

USE nrtype
! data type
USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO        ! Network topology
USE dataTypes, ONLY: subbasin_omp   ! mainstem+tributary data structures
! global parameters
USE public_var, ONLY: iulog         ! i/o logical unit number
USE public_var, ONLY: dt            ! simulation time step
! subroutines: general
USE perf_mod,   ONLY: t_startf,t_stopf ! timing start/stop

implicit none

private
public:: accum_wb
public:: comp_reach_wb

CONTAINS

  ! ---------------------------------------------------------------------------------------
  ! Public subroutine: compute upstream area water balance for each reach
  ! ---------------------------------------------------------------------------------------
  SUBROUTINE accum_wb(iens,          & ! input: index of runoff ensemble to be processed
                      ixRoute,       & ! input: index of routing method
                      river_basin,   & ! input: river basin information (mainstem, tributary outlet etc.)
                      ixDesire,      & ! input: ReachID to be checked by on-screen printing
                      NETOPO_in,     & ! input: reach topology data structure
                      RCHFLX_out,    & ! inout: reach flux data structure
                      ierr, message, & ! output: error controls
                      ixSubRch)        ! optional input: subset of reach indices to be processed
  ! ----------------------------------------------------------------------------------------
  ! -- Description:
  ! Computed drainage basin water balance for each reach
  ! for each reach, WB = sum(dVol) - (sum(Qlat) - sum(Qtake) - Qout)
  ! sum() is computed over all the upstream catchments
  !
  ! ----------------------------------------------------------------------------------------

  implicit none
  ! argument variables
  integer(i4b),                    intent(in)    :: iens            ! input: runoff ensemble index
  integer(i4b),                    intent(in)    :: ixRoute         ! input: routing method index
  type(subbasin_omp), allocatable, intent(in)    :: river_basin(:)  ! input: river basin information (mainstem, tributary outlet etc.)
  integer(i4b),                    intent(in)    :: ixDesire        ! input: index of the reach for verbose output
  type(RCHTOPO),      allocatable, intent(in)    :: NETOPO_in(:)    ! input: River Network topology
  type(STRFLX),                    intent(inout) :: RCHFLX_out(:,:) ! inout: Reach fluxes (ensembles, space [reaches]) for decomposed domains
  integer(i4b),                    intent(out)   :: ierr            ! output: error code
  character(*),                    intent(out)   :: message         ! output: error message
  integer(i4b),       optional,    intent(in)    :: ixSubRch(:)     ! optional input: subset of reach indices to be processed
  ! local variables
  integer(i4b)                                   :: nSeg            ! number of segments in the network
  integer(i4b)                                   :: nTrib           ! number of tributaries
  integer(i4b)                                   :: nDom            ! number of domains defined by e.g., stream order, tributary/mainstem
  integer(i4b)                                   :: iSeg, jSeg      ! reach segment indices
  integer(i4b)                                   :: iTrib, ix       ! loop indices
  logical(lgt), allocatable                      :: doRoute(:)      ! logical to indicate which reaches are processed
  character(len=strLen)                          :: cmessage        ! error message from subroutines

  ierr=0; message='accum_wb/'

  if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
   ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
  endif

  nSeg = size(RCHFLX_out(iens,:))

  allocate(doRoute(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for [doRoute]'; return; endif

  ! if a subset of reaches is processed
  if (present(ixSubRch))then
    doRoute(:)=.false.
    doRoute(ixSubRch) = .true. ! only subset of reaches are on
  ! if all the reaches are processed
  else
    doRoute(:)=.true. ! every reach is on
  endif

  nDom = size(river_basin)

  call t_startf('route/accum-wb')

  do ix = 1,nDom
    ! 1. Route tributary reaches (parallel)
    ! compute the sum of all upstream runoff at each point in the river network
    nTrib=size(river_basin(ix)%branch)

!$OMP PARALLEL DO schedule(dynamic,1)         &
!$OMP          private(jSeg, iSeg)            & ! private for a given thread
!$OMP          private(ierr, cmessage)        & ! private for a given thread
!$OMP          shared(river_basin)            & ! data structure shared
!$OMP          shared(doRoute)                & ! data array shared
!$OMP          shared(NETOPO_in)              & ! data structure shared
!$OMP          shared(RCHFLX_out)             & ! data structure shared
!$OMP          shared(ix, iens, ixDesire)     & ! indices shared
!$OMP          firstprivate(nTrib)
    do iTrib = 1,nTrib
      do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
        jSeg = river_basin(ix)%branch(iTrib)%segIndex(iSeg)

        if (.not. doRoute(jSeg)) cycle

        call comp_wb(iens, ixRoute, jSeg, ixDesire, NETOPO_in, RCHFLX_out)

      end do
    end do
!$OMP END PARALLEL DO

  end do

  call t_stopf('route/accum-wb')

  CONTAINS

    ! *********************************************************************
    ! subroutine: perform water balance error for a upstream basin
    ! *********************************************************************
    SUBROUTINE comp_wb(iens,       &    ! input: index of runoff ensemble to be processed
                       ixRoute,    &    ! input: index of routing method
                       segIndex,   &    ! input: index of reach to be processed
                       ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                       NETOPO_in,  &    ! input: reach topology data structure
                       RCHFLX_out)      ! inout: reach flux data structure
    implicit none
    ! Argument variables
    integer(i4b), intent(in)                 :: iens           ! runoff ensemble to be routed
    integer(i4b), intent(in)                 :: ixRoute        ! input: routing method index
    integer(i4b), intent(in)                 :: segIndex       ! segment where routing is performed
    integer(i4b), intent(in)                 :: ixDesire       ! index of the reach for verbose output
    type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)   ! River Network topology
    type(STRFLX), intent(inout)              :: RCHFLX_out(:,:)! Reach fluxes (ensembles, space [reaches]) for decomposed domains
    ! Local variables
    real(dp)                                 :: dVol           ! volume change [m3]
    real(dp)                                 :: Qlat           ! lateral flow [m3]
    real(dp)                                 :: precp          ! precipitation into reach [m3]
    real(dp)                                 :: Qtake          ! water take from reach [m3]
    real(dp)                                 :: evapo          ! evaporation out of reach[m3]
    real(dp)                                 :: Qout           ! out flow from reach [m3]
    real(dp)                                 :: WBupstream     ! water balance error for basin including all upstream areas [m3]
    integer(i4b)                             :: nUps           ! number of upstream segment
    integer(i4b)                             :: iUps           ! upstream reach index
    integer(i4b)                             :: iRch_ups       ! index of upstream reach in NETOPO

    ! identify number of upstream segments of the reach being processed
    nUps = size(NETOPO_in(segIndex)%UREACHI)

    ! local water balance
    dVol  = RCHFLX_out(iens,segIndex)%ROUTE(ixRoute)%REACH_VOL(1) &
           -RCHFLX_out(iens,segIndex)%ROUTE(ixRoute)%REACH_VOL(0)         ! volume change for this reach [m3]
    Qlat  = RCHFLX_out(iens,segIndex)%BASIN_QR(1) *dt                     ! later flow for this reach [m3]
    precp = RCHFLX_out(iens,segIndex)%basinprecip *dt                     ! precp for this reach [m3]
    Qtake = -1._dp *RCHFLX_out(iens,segIndex)%REACH_WM_FLUX *dt           ! Water abstraction from this reach [m3]
    evapo = -1._dp *RCHFLX_out(iens,segIndex)%basinevapo *dt              ! evaporation from this reach [m3]
    Qout  = -1._dp *RCHFLX_out(iens,segIndex)%ROUTE(ixRoute)%REACH_Q *dt  ! out flow from this reach [m3]

    RCHFLX_out(iens,segIndex)%ROUTE(ixRoute)%WBupstream = dVol - (Qlat+precp+Qtake+evapo+Qout)

    WBupstream = 0._dp
    if (nUps>0) then
      do iUps = 1,nUps
        iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
        WBupstream = WBupstream + RCHFLX_out(iens,iRch_ups)%ROUTE(ixRoute)%WBupstream &
                                - RCHFLX_out(iens,iRch_ups)%ROUTE(ixRoute)%REACH_Q *dt
      end do
    endif

    RCHFLX_out(iens,segIndex)%ROUTE(ixRoute)%WBupstream = RCHFLX_out(iens,segIndex)%ROUTE(ixRoute)%WBupstream &
                                                        + WBupstream
    ! check
    if(segIndex==ixDesire)then
      write(iulog,'(A,1PG15.7)') '  Total Upstream WBerr [m3] = ', RCHFLX_out(iens,segIndex)%ROUTE(ixRoute)%WBupstream
    endif

    END SUBROUTINE comp_wb

  END SUBROUTINE accum_wb

  ! *********************************************************************
  ! public subroutine: compute reach water balance
  ! *********************************************************************
  SUBROUTINE comp_reach_wb(ixRoute,    &     ! input: index of routing method
                           Qupstream,  &     ! input: inflow from upstream
                           RCHFLX_in,  &     ! inout: reach flux data structure
                           doCheck)

  ! Descriptions
  ! Compute water balance per simulation time step and reach/lake
  ! During a time step dt = t1 - t0 [sec]
  ! dVol     [m3]   = Vol(t1) - Vol(t0)
  ! flux_in  [m3/s] = inflow(dt) + lateral(dt) + Prec(dt)
  ! flux_out [m3/s] = outflow(dt) + abstract(dt) + Evap(dt)
  !
  ! Compare dVol vs (flux_in - flux_out)*dt

  implicit none
  ! Argument variables:
  integer(i4b), intent(in)                 :: ixRoute        ! input: routing method index
  real(dp),     intent(in)                 :: Qupstream      ! input: total inflow from upstream reaches
  type(STRFLX), intent(inout)              :: RCHFLX_in      ! inout: Reach fluxes data structure
  logical(lgt), intent(in)                 :: doCheck        ! input: reach index to be examined
  ! Local variables:
  real(dp)                                 :: dVol           !
  real(dp)                                 :: Qin            !
  real(dp)                                 :: Qlateral       !
  real(dp)                                 :: Qout           !
  real(dp)                                 :: precip         !
  real(dp)                                 :: evapo          !
  real(dp)                                 :: Qtake          !

  ! volume change
  dVol     = RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(1)-RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(0)
  ! in flux
  Qin      = Qupstream *dt
  Qlateral = RCHFLX_in%BASIN_QR(1) *dt
  precip   = RCHFLX_in%basinprecip *dt
  ! out flux
  Qout     = -1._dp *RCHFLX_in%ROUTE(ixRoute)%REACH_Q *dt
  Qtake    = -1._dp *RCHFLX_in%REACH_WM_FLUX *dt
  evapo    = -1._dp *RCHFLX_in%basinevapo *dt

  RCHFLX_in%ROUTE(ixRoute)%WB = dVol - (Qin + Qlateral + precip + Qout + Qtake + evapo)

  if (doCheck) then
    write(iulog,'(A,1PG15.7)') '  WBerr [m3]        = ', RCHFLX_in%ROUTE(ixRoute)%WB
    write(iulog,'(A,1PG15.7)') '  Vol at t0 [m3]    = ', RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(0)
    write(iulog,'(A,1PG15.7)') '  Vol at t1 [m3]    = ', RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(1)
    write(iulog,'(A,1PG15.7)') '  dVol [m3]         = ', dVol
    write(iulog,'(A,1PG15.7)') '  inflow [m3]       = ', Qin
    write(iulog,'(A,1PG15.7)') '  lateral flow [m3] = ', Qlateral
    write(iulog,'(A,1PG15.7)') '  precip [m3]       = ', precip
    write(iulog,'(A,1PG15.7)') '  outflow [m3]      = ', Qout
    write(iulog,'(A,1PG15.7)') '  abstraction [m3]  = ', Qtake
    write(iulog,'(A,1PG15.7)') '  evaporation [m3]  = ', evapo
  endif

  END SUBROUTINE comp_reach_wb

END MODULE water_balance
