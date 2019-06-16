module irf_route_module

!numeric type
USE nrtype
! data type
USE dataTypes,          only : STRFLX         ! fluxes in each reach
USE dataTypes,          only : RCHTOPO        ! Network topology
! global parameters
USE public_var,         only : realMissing    ! missing value for real number
USE public_var,         only : integerMissing ! missing value for integer number
USE globalData,         only : nThreads          ! number of threads used for openMP
! utilities
USE time_utils_module,  only : elapsedSec     ! calculate the elapsed time

! privary
implicit none
private

public::irf_route
public::irf_route_orig

contains

 ! *********************************************************************
 ! subroutine: perform network UH routing
 ! *********************************************************************
 subroutine irf_route(iEns,          &  ! input: index of runoff ensemble to be processed
                      river_basin,   &  ! input: river basin information (mainstem, tributary outlet etc.)
                      ixDesire,      &  ! input: reachID to be checked by on-screen pringing
                      NETOPO_in,     &  ! input: reach topology data structure
                      RCHFLX_out,    &  ! inout: reach flux data structure
                      ierr, message, &  ! output: error control
                      ixSubRch)         ! optional input: subset of reach indices to be processed

 ! global routing data
 USE dataTypes,  only : subbasin_omp   ! mainstem+tributary data structures

 implicit none
 ! Input
 integer(i4b),       intent(in)                  :: iEns                ! runoff ensemble to be routed
 type(subbasin_omp), intent(in),    allocatable  :: river_basin(:)      ! river basin information (mainstem, tributary outlet etc.)
 integer(i4b),       intent(in)                  :: ixDesire            ! index of the reach for verbose output ! Output
 type(RCHTOPO),      intent(in),    allocatable  :: NETOPO_in(:)        ! River Network topology
 ! inout
 TYPE(STRFLX),       intent(inout), allocatable  :: RCHFLX_out(:,:)     ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output variables
 integer(i4b),       intent(out)                 :: ierr                ! error code
 character(*),       intent(out)                 :: message             ! error message
 ! input (optional)
 integer(i4b),       intent(in), optional        :: ixSubRch(:)         ! subset of reach indices to be processed
 ! Local variables to
 character(len=strLen)                           :: cmessage            ! error message from subroutine
 logical(lgt),                      allocatable  :: doRoute(:)          ! logical to indicate which reaches are processed
 integer(i4b)                                    :: nTrib               ! number of tributary basins
 integer(i4b)                                    :: nSeg                ! number of reaches in the network
 integer(i4b)                                    :: iSeg, jSeg          ! loop indices - reach
 integer(i4b)                                    :: iTrib, iStem        ! loop indices - tributary, mainstem
 ! variables needed for timing
 integer(i4b)                                    :: omp_get_thread_num
 integer(i4b), allocatable                       :: ixThread(:)         ! thread id
 integer*8,    allocatable                       :: openMPend(:)        ! time for the start of the parallelization section
 integer*8,    allocatable                       :: timeTribStart(:)    ! time Tributaries start
 real(dp),     allocatable                       :: timeTrib(:)         ! time spent on each Tributary
 integer*8                                       :: cr                  ! rate
 integer*8                                       :: endTrib             ! date/time for the start and end of the initialization
 integer*8                                       :: startMain           ! date/time for the start and end of the initialization
 integer*8                                       :: startTime,endTime   ! date/time for the start and end of the initialization
 real(dp)                                        :: elapsedTime         ! elapsed time for the process

 ! initialize error control
 ierr=0; message='irf_route/'
 call system_clock(count_rate=cr)

 ! number of reach check
 if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
  ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
 endif

 nSeg = size(RCHFLX_out(iens,:))

 allocate(doRoute(nSeg), stat=ierr)

 ! Initialize CHEC_IRF to False.
 RCHFLX_out(iEns,:)%CHECK_IRF=.False.

 if (present(ixSubRch))then
  doRoute(:)=.false.
  doRoute(ixSubRch) = .true. ! only subset of reaches are on
 else
  doRoute(:)=.true. ! every reach is on
 endif

 nTrib=size(river_basin(1)%tributary)

 allocate(ixThread(nTrib), openMPend(nTrib), timeTrib(nTrib), timeTribStart(nTrib), stat=ierr)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to allocate space for Trib timing'; return; endif
 timeTrib(:) = realMissing
 ixThread(:) = integerMissing

 call system_clock(startTime)

 ! 1. Route tributary reaches (parallel)
!$OMP PARALLEL default(none)                            &
!$OMP          private(jSeg, iSeg)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(doRoute)                          & ! data array shared
!$OMP          shared(NETOPO_in)                        & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(iEns, ixDesire)                   & ! indices shared
!$OMP          shared(openMPend, nThreads)              & ! timing variables shared
!$OMP          shared(timeTribStart)                    & ! timing variables shared
!$OMP          shared(timeTrib)                         & ! timing variables shared
!$OMP          shared(ixThread)                         & ! thread id array shared
!$OMP          firstprivate(nTrib)

!$OMP DO schedule(dynamic,1)
  do iTrib = 1,nTrib
!$    ixThread(iTrib) = omp_get_thread_num()
    call system_clock(timeTribStart(iTrib))
    do iSeg=1,river_basin(1)%tributary(iTrib)%nRch
      jSeg = river_basin(1)%tributary(iTrib)%segIndex(iSeg)
      if (.not. doRoute(jSeg)) cycle
      call segment_irf(iEns, jSeg, ixDesire, NETOPO_IN, RCHFLX_out, ierr, cmessage)
!      if(ierr/=0)then; ixmessage(iTrib)=trim(message)//trim(cmessage); exit; endif
    end do
    call system_clock(openMPend(iTrib))
    timeTrib(iTrib) = real(openMPend(iTrib)-timeTribStart(iTrib), kind(dp))
  end do
!$OMP END DO
!$OMP END PARALLEL

  call system_clock(endTrib)
  elapsedTime = real(endTrib-startTime, kind(dp))/real(cr)
  write(*,"(A,1PG15.7,A)") '  total elapsed tributary = ', elapsedTime, ' s'

!  write(*,'(a)') 'iTrib nSeg ixThread nThreads StartTime EndTime'
!  do iTrib=1,nTrib
!    write(*,'(4(i5,1x),2(I20,1x))') iTrib, river_basin(1)%tributary(iTrib)%nRch, ixThread(iTrib), nThreads, timeTribStart(iTrib), openMPend(iTrib)
!  enddo
  deallocate(ixThread, timeTrib, timeTribStart, stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to deallocate space for Trib timing'; return; endif

   ! 2. Route mainstems (serial)
   call system_clock(startMain)
   if (allocated(river_basin(1)%mainstem)) then
     do iStem=1,size(river_basin(1)%mainstem)
       do iSeg=1,river_basin(1)%mainstem(iStem)%nRch
         jSeg = river_basin(1)%mainstem(iStem)%segIndex(iSeg)
         if (.not. doRoute(jSeg)) cycle
         call segment_irf(iEns, jSeg, ixDesire, NETOPO_IN, RCHFLX_out, ierr, message)
         if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
       end do
     end do
   endif

 call system_clock(endTime)
 elapsedTime = real(endTime-startMain, kind(dp))/real(cr)
 write(*,"(A,1PG15.7,A)") '  total elapsed Mainstem = ', elapsedTime, ' s'
 elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
 write(*,"(A,1PG15.7,A)") '  total elapsed entire = ', elapsedTime, ' s'

 end subroutine irf_route


 ! *********************************************************************
 ! subroutine: perform one segment route UH routing
 ! *********************************************************************
 subroutine segment_irf(&
                        ! input
                        iEns,       &    ! input: index of runoff ensemble to be processed
                        segIndex,   &    ! input: index of runoff ensemble to be processed
                        ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                        NETOPO_in,  &    ! input: reach topology data structure
                        ! inout
                        RCHFLX_out, &    ! inout: reach flux data structure
                        ! output
                        ierr, message)   ! output: error control

 implicit none
 ! Input
 INTEGER(I4B), intent(IN)                 :: iEns           ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)                 :: segIndex       ! segment where routing is performed
 INTEGER(I4B), intent(IN)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)   ! River Network topology
 ! inout
 TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! Output
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables to
 type(STRFLX), allocatable                :: uprflux(:)     ! upstream Reach fluxes
 INTEGER(I4B)                             :: nUps           ! number of upstream segment
 INTEGER(I4B)                             :: iUps           ! upstream reach index
 INTEGER(I4B)                             :: iRch_ups       ! index of upstream reach in NETOPO
 INTEGER(I4B)                             :: ntdh           ! number of time steps in IRF
 character(len=strLen)                    :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='segment_irf/'

 ! route streamflow through the river network
  if (.not.allocated(RCHFLX_out(iens,segIndex)%QFUTURE_IRF))then

   ntdh = size(NETOPO_in(segIndex)%UH)

   allocate(RCHFLX_out(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX_out(iens,segIndex)%QFUTURE_IRF'; return; endif

   RCHFLX_out(iens,segIndex)%QFUTURE_IRF(:) = 0._dp

  end if

  ! identify number of upstream segments of the reach being processed
  nUps = size(NETOPO_in(segIndex)%UREACHI)

  allocate(uprflux(nUps), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': uprflux'; return; endif

  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      uprflux(iUps) = RCHFLX_out(iens,iRch_ups)
    end do
  endif

  ! perform river network UH routing
  call conv_upsbas_qr(NETOPO_in(segIndex)%UH,    &    ! input: reach unit hydrograph
                      uprflux,                   &    ! input: upstream reach fluxes
                      RCHFLX_out(iens,segIndex), &    ! inout: updated fluxes at reach
                      ierr, message)                  ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! Check True since now this reach now routed
  RCHFLX_out(iEns,segIndex)%CHECK_IRF=.True.

  ! check
  if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   print*, 'RCHFLX_out(iens,segIndex)%BASIN_QR(1),RCHFLX_out(iens,segIndex)%REACH_Q_IRF = ', &
            RCHFLX_out(iens,segIndex)%BASIN_QR(1),RCHFLX_out(iens,segIndex)%REACH_Q_IRF
  endif

 end subroutine segment_irf


 ! *********************************************************************
 ! subroutine: Compute delayed runoff from the upstream segments
 ! *********************************************************************
 subroutine conv_upsbas_qr(&
                           ! input
                           reach_uh,   &    ! input: reach unit hydrograph
                           rflux_ups,  &    ! input: upstream reach fluxes
                           rflux,      &    ! input: input flux at reach
                           ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Convolute routed basisn flow volume at top of each of the upstream segment at one time step and at each segment
 !
 ! ----------------------------------------------------------------------------------------

 implicit none
 ! Input
 real(dp),     intent(in)               :: reach_uh(:)  ! reach unit hydrograph
 type(STRFLX), intent(in)               :: rflux_ups(:) ! upstream Reach fluxes
 ! inout
 type(STRFLX), intent(inout)            :: rflux        ! current Reach fluxes
 ! Output
 integer(i4b), intent(out)              :: ierr         ! error code
 character(*), intent(out)              :: message      ! error message
 ! Local variables to
 real(dp)                               :: q_upstream   ! total discharge at top of the reach being processed
 INTEGER(I4B)                           :: nTDH         ! number of UH data
 INTEGER(I4B)                           :: iTDH         ! index of UH data (i.e.,future time step)
 INTEGER(I4B)                           :: nUps         ! number of all upstream segment
 INTEGER(I4B)                           :: iUps         ! loop indices for u/s reaches

 ! initialize error control
 ierr=0; message='conv_upsbas_qr/'

 ! identify number of upstream segments of the reach being processed
 nUps = size(rflux_ups)

 q_upstream = 0.0_dp
 if(nUps>0)then
   do iUps = 1,nUps
     ! Find out total q at top of a segment
     q_upstream = q_upstream + rflux_ups(iUps)%REACH_Q_IRF
   end do
 endif

 ! place a fraction of runoff in future time steps
 nTDH = size(reach_uh) ! identify the number of future time steps of UH for a given segment
 do iTDH=1,nTDH
   rflux%QFUTURE_IRF(iTDH) = rflux%QFUTURE_IRF(iTDH) &
                             + reach_uh(iTDH)*q_upstream
 enddo

 ! Add local routed flow
 rflux%REACH_Q_IRF = rflux%QFUTURE_IRF(1) + rflux%BASIN_QR(1)

 ! move array back   use eoshift
 rflux%QFUTURE_IRF=eoshift(rflux%QFUTURE_IRF,shift=1)

 end subroutine conv_upsbas_qr

   ! *********************************************************************
   ! subroutine: perform network UH routing
   ! *********************************************************************
   subroutine irf_route_orig(iEns,         &  ! input: index of runoff ensemble to be processed
                             river_basin,  &  ! input: river basin information (mainstem, tributary outlet etc.)
                             ixDesire,     &  ! input: reachID to be checked by on-screen pringing
                             NETOPO_in,    &  ! input: reach topology data structure
                             RCHFLX_out,   &  ! inout: reach flux data structure
                             ierr, message,&  ! output: error control
                             ixSubRch)        ! optional input: subset of reach indices to be processed
   ! ----------------------------------------------------------------------------------------
   ! Purpose:
   !
   !   Convolute routed basisn flow volume at top of each of the upstream segment at one time step and at each segment
   !
   ! ----------------------------------------------------------------------------------------

   ! global routing data
   USE dataTypes,         only : subbasin_omp  ! mainstem+tributary data strucuture

   implicit none
   ! Input
   integer(I4B), intent(in)                     :: iEns              ! runoff ensemble to be routed
   type(subbasin_omp), intent(in), allocatable  :: river_basin(:)    ! river basin information (mainstem, tributary outlet etc.)
   integer(I4B), intent(in)                     :: ixDesire          ! index of the reach for verbose output
   type(RCHTOPO),intent(in),    allocatable     :: NETOPO_in(:)      ! River Network topology
   ! inout
   TYPE(STRFLX), intent(inout), allocatable     :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   ! Output
   integer(i4b), intent(out)                    :: ierr              ! error code
   character(*), intent(out)                    :: message           ! error message
   ! input (optional)
   integer(i4b), intent(in),   optional         :: ixSubRch(:)       ! subset of reach indices to be processed
   ! Local variables to
   INTEGER(I4B)                                 :: nSeg              ! number of reach segments in the network
   INTEGER(I4B)                                 :: iSeg, jSeg        ! reach segment index
   logical(lgt), allocatable                    :: doRoute(:)        ! logical to indicate which reaches are processed
   character(len=strLen)                        :: cmessage          ! error message from subroutine
   integer*8                                    :: startTime,endTime ! date/time for the start and end of the initialization
   real(dp)                                     :: elapsedTime       ! elapsed time for the process

   ! initialize error control
   ierr=0; message='irf_route_orig/'

   call system_clock(startTime)

   ! check
   if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
    ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
   endif

   ! Initialize CHEC_IRF to False.
   RCHFLX_out(iEns,:)%CHECK_IRF=.False.

   nSeg = size(RCHFLX_out(iens,:))

   allocate(doRoute(nSeg), stat=ierr)

   if (present(ixSubRch))then
    doRoute(:)=.false.
    doRoute(ixSubRch) = .true. ! only subset of reaches are on
   else
    doRoute(:)=.true. ! every reach is on
   endif

   ! route streamflow through the river network
   do iSeg=1,nSeg

    jSeg = NETOPO_in(iSeg)%RHORDER

    if (.not. doRoute(jSeg)) cycle

    call segment_irf(iEns, jSeg, ixDesire, NETOPO_IN, RCHFLX_out, ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end do

   call system_clock(endTime)
   elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!   write(*,"(A,1PG15.7,A)") '  total elapsed entire = ', elapsedTime, ' s'

 end subroutine irf_route_orig

end module irf_route_module

