module irf_route_module

!numeric type
USE nrtype
! data type
USE dataTypes,          only : STRFLX         ! fluxes in each reach
USE dataTypes,          only : RCHTOPO        ! Network topology
USE dataTypes,          only : RCHPRP         ! Reach parameter
USE dataTypes,          only : irfRCH         ! irf specific state data structure 
! global parameters
USE public_var,         only : realMissing    ! missing value for real number
USE public_var,         only : integerMissing ! missing value for integer number
USE globalData,         only : nThreads       ! number of threads used for openMP
USE globalData,         only : idxIRF          ! index of IRF method 
! subroutines: general
USE model_finalize, ONLY : handle_err

implicit none

private

public::irf_route

contains

 ! *********************************************************************
 ! subroutine: perform network UH routing
 ! *********************************************************************
 subroutine irf_route(iEns,          &  ! input: index of runoff ensemble to be processed
                      river_basin,   &  ! input: river basin information (mainstem, tributary outlet etc.)
                      ixDesire,      &  ! input: reachID to be checked by on-screen pringing
                      NETOPO_in,     &  ! input: reach topology data structure
                      RPARAM_in,     &  ! input: reach parameter data structure
                      RCHFLX_out,    &  ! inout: reach flux data structure
                      ierr, message, &  ! output: error control
                      ixSubRch)         ! optional input: subset of reach indices to be processed

 USE dataTypes,      ONLY : subbasin_omp   ! mainstem+tributary data structures

 implicit none
 ! Input
 integer(i4b),       intent(in)                  :: iEns                ! runoff ensemble to be routed
 type(subbasin_omp), intent(in),    allocatable  :: river_basin(:)      ! river basin information (mainstem, tributary outlet etc.)
 integer(i4b),       intent(in)                  :: ixDesire            ! index of the reach for verbose output ! Output
 type(RCHTOPO),      intent(in),    allocatable  :: NETOPO_in(:)        ! River Network topology
 type(RCHPRP),       intent(in),    allocatable  :: RPARAM_in(:)        ! River reach parameter
 ! inout
 TYPE(STRFLX),       intent(inout), allocatable  :: RCHFLX_out(:,:)     ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output variables
 integer(i4b),       intent(out)                 :: ierr                ! error code
 character(*),       intent(out)                 :: message             ! error message
 ! input (optional)
 integer(i4b),       intent(in), optional        :: ixSubRch(:)         ! subset of reach indices to be processed
 ! Local variables
 character(len=strLen)                           :: cmessage            ! error message from subroutine
 logical(lgt),                      allocatable  :: doRoute(:)          ! logical to indicate which reaches are processed
 integer(i4b)                                    :: nOrder              ! number of stream order
 integer(i4b)                                    :: nTrib               ! number of tributary basins
 integer(i4b)                                    :: nSeg                ! number of reaches in the network
 integer(i4b)                                    :: iSeg, jSeg          ! loop indices - reach
 integer(i4b)                                    :: iTrib               ! loop indices - branch
 integer(i4b)                                    :: ix                  ! loop indices stream order
 ! variables needed for timing
 !integer(i4b)                                    :: omp_get_thread_num
 !integer(i4b), allocatable                       :: ixThread(:)         ! thread id
 !integer*8,    allocatable                       :: openMPend(:)        ! time for the start of the parallelization section
 !integer*8,    allocatable                       :: timeTribStart(:)    ! time Tributaries start
 !real(dp),     allocatable                       :: timeTrib(:)         ! time spent on each Tributary

 ierr=0; message='irf_route/'

 ! number of reach check
 if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
  ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
 endif

 nSeg = size(NETOPO_in)

 allocate(doRoute(nSeg), stat=ierr)

 if (present(ixSubRch))then
  doRoute(:)=.false.
  doRoute(ixSubRch) = .true. ! only subset of reaches are on
 else
  doRoute(:)=.true. ! every reach is on
 endif

 nOrder = size(river_basin)


 do ix = 1,nOrder

   nTrib=size(river_basin(ix)%branch)

!  allocate(ixThread(nTrib), openMPend(nTrib), timeTrib(nTrib), timeTribStart(nTrib), stat=ierr)
!  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to allocate space for Trib timing'; return; endif
!  timeTrib(:) = realMissing
!  ixThread(:) = integerMissing

  ! 1. Route tributary reaches (parallel)
!$OMP PARALLEL DO schedule(dynamic,1)                   &
!$OMP          private(jSeg, iSeg)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(doRoute)                          & ! data array shared
!$OMP          shared(NETOPO_in)                        & ! data structure shared
!$OMP          shared(RPARAM_in)                        & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(ix, iEns, ixDesire)               & ! indices shared
!$OMP          firstprivate(nTrib)
!!$OMP          shared(openMPend, nThreads)              & ! timing variables shared
!!$OMP          shared(timeTribStart)                    & ! timing variables shared
!!$OMP          shared(timeTrib)                         & ! timing variables shared
!!$OMP          shared(ixThread)                         & ! thread id array shared
   trib:do iTrib = 1,nTrib
!!$    ixThread(iTrib) = omp_get_thread_num()
!    call system_clock(timeTribStart(iTrib))
     seg:do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
       jSeg = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
       if (.not. doRoute(jSeg)) cycle
       call segment_irf(iEns, jSeg, ixDesire, NETOPO_IN, RPARAM_in, RCHFLX_out, ierr, cmessage)
       if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
     end do seg
!    call system_clock(openMPend(iTrib))
!    timeTrib(iTrib) = real(openMPend(iTrib)-timeTribStart(iTrib), kind(dp))
   end do trib
!$OMP END PARALLEL DO

!  write(*,'(a)') 'iTrib nSeg ixThread nThreads StartTime EndTime'
!  do iTrib=1,nTrib
!    write(*,'(4(i5,1x),2(I20,1x))') iTrib, river_basin(iOut)%branch(iTrib)%nRch, ixThread(iTrib), nThreads, timeTribStart(iTrib), openMPend(iTrib)
!  enddo
!  deallocate(ixThread, openMPend, timeTrib, timeTribStart, stat=ierr)
!  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to deallocate space for Trib timing'; return; endif

 end do ! basin loop

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
                        RPARAM_in,  &    ! input: reach parameter data structure
                        ! inout
                        RCHFLX_out, &    ! inout: reach flux data structure
                        ! output
                        ierr, message)   ! output: error control

 implicit none
 ! Input
 integer(i4b), intent(in)                 :: iEns           ! runoff ensemble to be routed
 integer(i4b), intent(in)                 :: segIndex       ! segment where routing is performed
 integer(i4b), intent(in)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)   ! River Network topology
 type(RCHPRP), intent(in),    allocatable :: RPARAM_in(:)   ! River reach parameter
 ! inout
 type(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! Output
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables to
 real(dp)                                 :: q_upstream     ! total discharge at top of the reach being processed
 integer(i4b)                             :: nUps           ! number of upstream segment
 integer(i4b)                             :: iUps           ! upstream reach index
 integer(i4b)                             :: iRch_ups       ! index of upstream reach in NETOPO
 integer(i4b)                             :: ntdh           ! number of time steps in IRF
 integer(i4b)                             :: itdh           ! loop index for unit hydrograph
 character(len=strLen)                    :: fmt1           ! format string
 character(len=strLen)                    :: cmessage       ! error message from subroutine

 ierr=0; message='segment_irf/'

 ! initialize future discharge array at first time
  if (.not.allocated(RCHFLX_out(iens,segIndex)%QFUTURE_IRF))then

   ntdh = size(NETOPO_in(segIndex)%UH)

   allocate(RCHFLX_out(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX_out(iens,segIndex)%QFUTURE_IRF'; return; endif

   RCHFLX_out(iens,segIndex)%QFUTURE_IRF(:) = 0._dp

  end if

  ! get discharge coming from upstream
  nUps = size(NETOPO_in(segIndex)%UREACHI)
  q_upstream = 0.0_dp
  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxIRF)%REACH_Q
    end do
  endif

  ! perform UH convolution
  call conv_upsbas_qr(NETOPO_in(segIndex)%UH,    &    ! input: reach unit hydrograph
                      q_upstream,                &    ! input: total discharge at top of the reach being processed
                      RCHFLX_out(iens,segIndex), &    ! inout: updated fluxes at reach
                      RPARAM_in(segIndex)%QTAKE, &    ! input: abstraction(-)/injection(+) [m3/s]
                      RPARAM_in(segIndex)%MINFLOW, &  ! input: minimum environmental flow [m3/s]
                      ierr, message)                  ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! check
  if(segIndex==ixDesire)then
    ntdh = size(NETOPO_in(segIndex)%UH)
    write(fmt1,'(A,I5,A)') '(A, 1X',ntdh,'(1X,F20.7))'
    write(*,'(2a)') new_line('a'),'** Check Impulse Response Function routing **'
    write(*,'(a,x,I10,x,I10)') ' Reach index & ID       =', segIndex, NETOPO_in(segIndex)%REACHID
    write(*,fmt1)              ' Unit-Hydrograph        =', (NETOPO_in(segIndex)%UH(itdh), itdh=1,ntdh)
    write(*,'(a)')             ' * total discharge from upstream(q_upstream) [m3/s], local area discharge [m3/s], and Final discharge [m3/s]:'
    write(*,'(a,x,F15.7)')     ' q_upstream             =', q_upstream
    write(*,'(a,x,F15.7)')     ' RCHFLX_out%BASIN_QR(1) =', RCHFLX_out(iens,segIndex)%BASIN_QR(1)
    write(*,'(a,x,F15.7)')     ' RCHFLX_out%ROUTE%REACH_Q =', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q
  endif

 end subroutine segment_irf


 ! *********************************************************************
 ! subroutine: Compute delayed runoff from the upstream segments
 ! *********************************************************************
 subroutine conv_upsbas_qr(reach_uh,   &    ! input: reach unit hydrograph
                           q_upstream, &    ! input:
                           rflux,      &    ! input: input flux at reach
                           Qtake,      &    ! input: abstraction(-)/injection(+) [m3/s]
                           Qmin,       &    ! input: minimum environmental flow [m3/s]
                           ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Details: Convolute runoff volume of upstream at one reach at one time step
 ! ----------------------------------------------------------------------------------------

 USE public_var, ONLY: dt

 implicit none
 ! Input
 real(dp),     intent(in)               :: reach_uh(:)  ! reach unit hydrograph
 real(dp),     intent(in)               :: q_upstream   ! total discharge at top of the reach being processed
 real(dp),     intent(in)               :: Qtake        ! abstraction(-)/injection(+) [m3/s]
 real(dp),     intent(in)               :: Qmin         ! minimum environmental flow [m3/s]
 ! inout
 type(STRFLX), intent(inout)            :: rflux        ! current Reach fluxes
 ! Output
 integer(i4b), intent(out)              :: ierr         ! error code
 character(*), intent(out)              :: message      ! error message
 ! Local variables to
 real(dp)                               :: QupMod       ! modified total discharge at top of the reach being processed
 real(dp)                               :: Qabs         ! maximum allowable water abstraction rate [m3/s]
 real(dp)                               :: Qmod         ! abstraction rate to be taken from outlet discharge [m3/s]
 integer(i4b)                           :: ntdh         ! number of UH data (i.e., number of future time step
 integer(i4b)                           :: itdh         ! index of UH data

 ierr=0; message='conv_upsbas_qr/'

 ! Q injection, add at top of reach
 QupMod = q_upstream
 if (Qtake>0) then
   QupMod = QupMod+ Qtake
 end if

 ! place a fraction of runoff in future time steps
 ntdh = size(reach_uh) ! number of future time steps of UH for a given segment
 do itdh=1,ntdh
   rflux%QFUTURE_IRF(itdh) = rflux%QFUTURE_IRF(itdh)+ reach_uh(itdh)*QupMod
 enddo

 ! compute volume in reach
 rflux%ROUTE(idxIRF)%REACH_VOL(0) = rflux%ROUTE(idxIRF)%REACH_VOL(1)
 rflux%ROUTE(idxIRF)%REACH_VOL(1) = rflux%ROUTE(idxIRF)%REACH_VOL(0) + (QupMod - rflux%QFUTURE_IRF(1))/dt

 ! Add local routed flow at the bottom of reach
 rflux%ROUTE(idxIRF)%REACH_Q = rflux%QFUTURE_IRF(1) + rflux%BASIN_QR(1)

 ! Q abstraction
 ! Compute actual abstraction (Qabs) m3/s - values should be negative
 ! Compute abstraction (Qmod) m3 taken from outlet discharge (REACH_Q)
 ! Compute REACH_Q subtracted from Qmod abstraction
 ! Compute REACH_VOL subtracted from total abstraction minus abstraction from outlet discharge
 if (Qtake<0) then
   Qabs               = max(-(rflux%ROUTE(idxIRF)%REACH_VOL(1)/dt+rflux%ROUTE(idxIRF)%REACH_Q), Qtake)
   Qmod               = min(rflux%ROUTE(idxIRF)%REACH_VOL(1) + Qabs*dt, 0._dp)
   rflux%ROUTE(idxIRF)%REACH_Q      = max(rflux%ROUTE(idxIRF)%REACH_Q + Qmod/dt, Qmin)
   rflux%ROUTE(idxIRF)%REACH_VOL(1) = rflux%ROUTE(idxIRF)%REACH_VOL(1) + (Qabs*dt - Qmod)
 end if

 ! move array back   use eoshift
 rflux%QFUTURE_IRF=eoshift(rflux%QFUTURE_IRF,shift=1)

 rflux%QFUTURE_IRF(ntdh) = 0._dp

 end subroutine conv_upsbas_qr

end module irf_route_module
