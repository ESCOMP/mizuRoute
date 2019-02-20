module irf_route_module

!numeric type
USE nrtype
! data type
USE dataTypes,          only : STRFLX         ! fluxes in each reach
! global parameters
USE public_var,         only : realMissing    ! missing value for real number
USE public_var,         only : integerMissing ! missing value for integer number
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
 subroutine irf_route(&
                      ! input
                      iEns,       &    ! input: index of runoff ensemble to be processed
                      river_basin,&    ! input: river basin information (mainstem, tributary outlet etc.)
                      ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                      ! output
                      ierr, message)   ! output: error control

 ! global routing data
 USE globalData, only : RCHFLX ! routing fluxes
 USE dataTypes,  only : basin  ! river basin data type

 implicit none
 ! Input
 integer(i4b), intent(in)               :: iEns                  ! runoff ensemble to be routed
 type(basin),  intent(in), allocatable  :: river_basin(:)        ! river basin information (mainstem, tributary outlet etc.)
 integer(i4b), intent(in)               :: ixDesire              ! index of the reach for verbose output ! Output
 ! output variables
 integer(i4b), intent(out)              :: ierr                  ! error code
 character(*), intent(out)              :: message               ! error message
 ! Local variables to
 integer(i4b)                           :: nOuts                 ! number of outlets
 integer(i4b)                           :: nTrib                 ! number of tributary basins
 integer(i4b)                           :: nStem                 ! number of mainstem in each level
 integer(i4b)                           :: iRch, jRch            ! loop indices - reach
 integer(i4b)                           :: iOut                  ! loop indices - basin outlet
 integer(i4b)                           :: iTrib                 ! loop indices - tributary
 integer(i4b)                           :: iLevel                ! loop indices - mainstem level
 integer(i4b)                           :: iStem                 ! loop inidces - mainstem
 integer(i4b)                           :: maxLevel,minLevel     ! max. and min. mainstem levels
 character(len=strLen)                  :: cmessage              ! error message from subroutine
 ! variables needed for timing
 integer(i4b)                           :: nThreads              ! number of threads
 integer(i4b)                           :: omp_get_num_threads   ! get the number of threads
 integer(i4b)                           :: omp_get_thread_num
 integer(i4b), allocatable              :: ixThread(:)           ! thread id
 integer*8,    allocatable              :: openMPend(:)          ! time for the start of the parallelization section
 integer*8,    allocatable              :: timeTribStart(:)      ! time Tributaries start
 real(dp),     allocatable              :: timeTrib(:)           ! time spent on each Tributary
 integer*8                              :: endTrib               ! date/time for the start and end of the initialization
 integer*8                              :: startTime,endTime     ! date/time for the start and end of the initialization
 integer*8                              :: startMain,endMain     ! date/time for the start and end of the initialization
 real(dp)                               :: elapsedTime           ! elapsed time for the process
 real(dp)                               :: elapsedTrib           ! elapsed time for the process
 real(dp)                               :: elapsedMain           ! elapsed time for the process

 ! initialize error control
 ierr=0; message='irf_route/'

 ! Initialize CHEC_IRF to False.
 RCHFLX(iEns,:)%CHECK_IRF=.False.

 ! Number of Outlets
 nOuts = size(river_basin)

 do iOut=1,nOuts

  maxLevel = 1
  do iLevel=1,size(river_basin(iOut)%level)
    if (.not. allocated(river_basin(iOut)%level(iLevel)%mainstem)) cycle
    if (iLevel > maxLevel) maxLevel = iLevel
  end do
  minLevel = size(river_basin(iOut)%level)
  do iLevel=maxLevel,1,-1
    if (.not. allocated(river_basin(iOut)%level(iLevel)%mainstem)) cycle
    if (iLevel < minLevel) minLevel = iLevel
  end do

  nTrib=size(river_basin(iOut)%tributary)

  allocate(ixThread(nTrib), openMPend(nTrib), timeTrib(nTrib), timeTribStart(nTrib), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to allocate space for Trib timing'; return; endif
  timeTrib(:) = realMissing
  ixThread(:) = integerMissing

  ! 1. Route tributary reaches (parallel)
  call system_clock(startTime)
!$OMP PARALLEL default(none)                            &
!$OMP          private(jRch, iRch)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(iEns, iOut, ixDesire)             & ! indices shared
!$OMP          shared(openMPend, nThreads)              & ! timing variables shared
!$OMP          shared(timeTribStart)                    & ! timing variables shared
!$OMP          shared(timeTrib)                         & ! timing variables shared
!$OMP          shared(ixThread)                         & ! thread id array shared
!$OMP          firstprivate(nTrib)

  nThreads = 1
!$ nThreads = omp_get_num_threads()

!$OMP DO schedule(dynamic,1)
  do iTrib = 1,nTrib
!$    ixThread(iTrib) = omp_get_thread_num()
    call system_clock(timeTribStart(iTrib))
    do iRch=1,river_basin(iOut)%tributary(iTrib)%nRch
      jRch = river_basin(iOut)%tributary(iTrib)%segIndex(iRch)
      call segment_irf(iEns, jRch, ixDesire, ierr, cmessage)
!      if(ierr/=0)then; ixmessage(iTrib)=trim(message)//trim(cmessage); exit; endif
    end do
    call system_clock(openMPend(iTrib))
    timeTrib(iTrib) = real(openMPend(iTrib)-timeTribStart(iTrib), kind(dp))
  end do
!$OMP END DO
!$OMP END PARALLEL

  call system_clock(endTrib)
  elapsedTrib = real(endTrib-startTime, kind(dp))/10e8_dp
  write(*,"(A,1PG15.7,A)") '  total elapsed trib = ', elapsedTrib, ' s'

!  write(*,'(a)') 'iTrib nSeg ixThread nThreads StartTime EndTime'
!  do iTrib=1,nTrib
!    write(*,'(4(i5,1x),2(I20,1x))') iTrib, river_basin(iOut)%tributary(iTrib)%nRch, ixThread(iTrib), nThreads, timeTribStart(iTrib), openMPend(iTrib)
!  enddo
!  deallocate(ixThread, timeTrib, timeTribStart, stat=ierr)
!  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to deallocate space for Trib timing'; return; endif

   ! 2. Route mainstems (serial)
   call system_clock(startMain)
   do iLevel=maxLevel,minLevel,-1
     nStem = size(river_basin(iOut)%level(iLevel)%mainstem)
!$OMP PARALLEL default(none)                            &
!$OMP          private(jRch, iRch)                      & ! private for a given thread
!$OMP          private(iStem)                           & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(iLevel)                           & ! private for a given thread
!$OMP          shared(iEns, iOut, ixDesire)             & ! indices shared
!$OMP          firstprivate(nStem)
!$OMP DO schedule(dynamic,1)
     do iStem=1,nStem
       do iRch=1,river_basin(iOut)%level(iLevel)%mainstem(iStem)%nRch
         jRch = river_basin(iOut)%level(iLevel)%mainstem(iStem)%segIndex(iRch)
         call segment_irf(iEns, jRch, ixDesire, ierr, cmessage)
!         if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
       end do
     end do
!$OMP END DO
!$OMP END PARALLEL
   end do

   call system_clock(endMain)
   elapsedMain = real(endMain-startMain, kind(dp))/10e8_dp
   write(*,"(A,1PG15.7,A)") '  total elapsed mainstem = ', elapsedMain, ' s'

 end do
 call system_clock(endTime)
 elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
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
                        ! output
                        ierr, message)   ! output: error control

 ! global routing data
 USE globalData, only : RCHFLX ! routing fluxes
 USE globalData, only : NETOPO ! Network topology

 implicit none
 ! Input
 INTEGER(I4B), intent(IN)               :: iEns           ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)               :: segIndex       ! segment where routing is performed
 INTEGER(I4B), intent(IN)               :: ixDesire       ! index of the reach for verbose output
 ! Output
 integer(i4b), intent(out)              :: ierr           ! error code
 character(*), intent(out)              :: message        ! error message
 ! Local variables to
 type(STRFLX), allocatable              :: uprflux(:)     ! upstream Reach fluxes
 INTEGER(I4B)                           :: nUps           ! number of upstream segment
 INTEGER(I4B)                           :: iUps           ! upstream reach index
 INTEGER(I4B)                           :: iRch_ups       ! index of upstream reach in NETOPO
 INTEGER(I4B)                           :: ntdh           ! number of time steps in IRF
 character(len=strLen)                  :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='segment_irf/'

 ! route streamflow through the river network
  if (.not.allocated(RCHFLX(iens,segIndex)%QFUTURE_IRF))then

   ntdh = size(NETOPO(segIndex)%UH)

   allocate(RCHFLX(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX(iens,segIndex)%QFUTURE_IRF'; return; endif

   RCHFLX(iens,segIndex)%QFUTURE_IRF(:) = 0._dp

  end if

  ! identify number of upstream segments of the reach being processed
  nUps = size(NETOPO(segIndex)%UREACHI)

  allocate(uprflux(nUps), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': uprflux'; return; endif

  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      uprflux(iUps) = RCHFLX(iens,iRch_ups)
    end do
  endif

  ! perform river network UH routing
  call conv_upsbas_qr(NETOPO(segIndex)%UH,   &    ! input: reach unit hydrograph
                      uprflux,               &    ! input: upstream reach fluxes
                      RCHFLX(iens,segIndex), &    ! inout: updated fluxes at reach
                      ierr, message)              ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! Check True since now this reach now routed
  RCHFLX(iEns,segIndex)%CHECK_IRF=.True.

  ! check
  if(NETOPO(segIndex)%REACHIX == ixDesire)then
   print*, 'RCHFLX(iens,segIndex)%BASIN_QR(1),RCHFLX(iens,segIndex)%REACH_Q_IRF = ', &
            RCHFLX(iens,segIndex)%BASIN_QR(1),RCHFLX(iens,segIndex)%REACH_Q_IRF
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
   subroutine irf_route_orig(&
                             ! input
                             iEns,       &    ! input: index of runoff ensemble to be processed
                             river_basin,&    ! input: river basin information (mainstem, tributary outlet etc.)
                             ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                             ! output
                             ierr, message)   ! output: error control
   ! ----------------------------------------------------------------------------------------
   ! Purpose:
   !
   !   Convolute routed basisn flow volume at top of each of the upstream segment at one time step and at each segment
   !
   ! ----------------------------------------------------------------------------------------

   ! global routing data
   USE globalData, only : RCHFLX ! routing fluxes
   USE globalData, only : NETOPO ! Network topology
   USE dataTypes,  only : basin  ! river basin data type

   implicit none
   ! Input
   INTEGER(I4B), intent(IN)               :: iEns           ! runoff ensemble to be routed
   type(basin),  intent(in), allocatable  :: river_basin(:) ! river basin information (mainstem, tributary outlet etc.)
   INTEGER(I4B), intent(IN)               :: ixDesire       ! index of the reach for verbose output
   ! Output
   integer(i4b), intent(out)              :: ierr           ! error code
   character(*), intent(out)              :: message        ! error message
   ! Local variables to
   INTEGER(I4B)                           :: nRch           ! number of reach segments in the network
   INTEGER(I4B)                           :: iRch           ! reach segment index
   INTEGER(I4B)                           :: jRch           ! reach segment to be routed
   character(len=strLen)                  :: cmessage       ! error message from subroutine
   integer*8                              :: startTime,endTime ! date/time for the start and end of the initialization
   real(dp)                               :: elapsedTime       ! elapsed time for the process

   ! initialize error control
   ierr=0; message='irf_route_orig/'

   ! Initialize CHEC_IRF to False.
   RCHFLX(iEns,:)%CHECK_IRF=.False.

   nRch=size(NETOPO)

   elapsedTime = 0._dp
   call system_clock(startTime)
   ! route streamflow through the river network
   do iRch=1,nRch

    jRch = NETOPO(iRch)%RHORDER

    call segment_irf(iEns, jRch, ixDesire, ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end do
   call system_clock(endTime)
   elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
   write(*,"(A,1PG15.7,A)") '  total elapsed entire = ', elapsedTime, ' s'

 end subroutine irf_route_orig

end module irf_route_module

