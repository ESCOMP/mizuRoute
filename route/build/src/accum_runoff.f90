MODULE accum_runoff_module

USE nrtype
! data type
USE dataTypes,      ONLY: STRFLX         ! fluxes in each reach
USE dataTypes,      ONLY: RCHTOPO        ! Network topology
USE dataTypes,      ONLY: subbasin_omp   ! mainstem+tributary data structures
! global data
USE public_var,     ONLY: iulog          ! i/o logical unit number
USE globalData,     ONLY: idxSUM         ! index of accumulation method
! subroutines: general
USE model_finalize, ONLY : handle_err

implicit none

private
public::accum_runoff

CONTAINS

 ! ---------------------------------------------------------------------------------------
 ! Public subroutine main driver for basin routing
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE accum_runoff(iEns,          & ! input: index of runoff ensemble to be processed
                         river_basin,   & ! input: river basin information (mainstem, tributary outlet etc.)
                         ixDesire,      & ! input: ReachID to be checked by on-screen printing
                         NETOPO_in,     & ! input: reach topology data structure
                         RCHFLX_out,    & ! inout: reach flux data structure
                         ierr, message, & ! output: error controls
                         ixSubRch)        ! optional input: subset of reach indices to be processed
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 ! Accumulate all the upstream delayed runoff for each reach
 ! mostly for checking if routed runoff at each reach outlet preserve total upstream runoff.
 !
 ! ----------------------------------------------------------------------------------------

 implicit none
 ! input
 integer(i4b),       intent(in)                 :: iens            ! runoff ensemble index
 type(subbasin_omp), intent(in),    allocatable :: river_basin(:)  ! river basin information (mainstem, tributary outlet etc.)
 integer(i4b),       intent(in)                 :: ixDesire        ! index of the reach for verbose output
 type(RCHTOPO),      intent(in),    allocatable :: NETOPO_in(:)    ! River Network topology
 ! inout
 TYPE(STRFLX),       intent(inout), allocatable :: RCHFLX_out(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output
 integer(i4b),       intent(out)                :: ierr            ! error code
 character(*),       intent(out)                :: message         ! error message
 ! input (optional)
 integer(i4b),       intent(in),   optional     :: ixSubRch(:)     ! subset of reach indices to be processed
 ! local variables
 integer(i4b)                                   :: nSeg            ! number of segments in the network
 integer(i4b)                                   :: nTrib           ! number of tributaries
 integer(i4b)                                   :: nDom            ! number of domains defined by e.g., stream order, tributary/mainstem
 integer(i4b)                                   :: iSeg, jSeg      ! reach segment indices
 integer(i4b)                                   :: iTrib, ix       ! loop indices
 logical(lgt), allocatable                      :: doRoute(:)      ! logical to indicate which reaches are processed
 character(len=strLen)                          :: cmessage        ! error message from subroutines

 ierr=0; message='accum_runoff/'

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
!$OMP          shared(ix, iEns, ixDesire)     & ! indices shared
!$OMP          firstprivate(nTrib)
   do iTrib = 1,nTrib
     do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
       jSeg = river_basin(ix)%branch(iTrib)%segIndex(iSeg)

       if (.not. doRoute(jSeg)) cycle

       call accum_qupstream(iens, jSeg, ixDesire, NETOPO_in, RCHFLX_out, ierr, cmessage)
       if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))

     end do
   end do
!$OMP END PARALLEL DO

 end do

 END SUBROUTINE accum_runoff

 ! *********************************************************************
 ! subroutine: perform accumulate immediate upstream flow
 ! *********************************************************************
 SUBROUTINE accum_qupstream(iEns,       &    ! input: index of runoff ensemble to be processed
                            segIndex,   &    ! input: index of reach to be processed
                            ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                            NETOPO_in,  &    ! input: reach topology data structure
                            RCHFLX_out, &    ! inout: reach flux data structure
                            ierr, message)   ! output: error control
 implicit none
 ! Input
 integer(i4b), intent(in)                 :: iEns           ! runoff ensemble to be routed
 integer(i4b), intent(in)                 :: segIndex       ! segment where routing is performed
 integer(i4b), intent(in)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)   ! River Network topology
 ! inout
 TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! Output
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables to
 real(dp)                                 :: q_upstream     ! upstream Reach fluxes
 integer(i4b)                             :: nUps           ! number of upstream segment
 integer(i4b)                             :: iUps           ! upstream reach index
 integer(i4b)                             :: iRch_ups       ! index of upstream reach in NETOPO
 character(len=strLen)                    :: fmt1,fmt2      ! format string

 ierr=0; message='accum_qupstream/'

 ! identify number of upstream segments of the reach being processed
 nUps = size(NETOPO_in(segIndex)%UREACHI)

 RCHFLX_out(iEns,segIndex)%ROUTE(idxSUM)%REACH_Q = RCHFLX_out(iEns,segIndex)%BASIN_QR(1)

 q_upstream = 0._dp
 if (nUps>0) then

   do iUps = 1,nUps
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     q_upstream = q_upstream + RCHFLX_out(iens,iRch_ups)%ROUTE(idxSUM)%REACH_Q
   end do

   RCHFLX_out(iEns,segIndex)%ROUTE(idxSUM)%REACH_Q = RCHFLX_out(iEns,segIndex)%ROUTE(idxSUM)%REACH_Q + q_upstream

 endif

 ! check
 if(segIndex == ixDesire)then
   write(fmt1,'(A,I5,A)') '(A,1X',nUps,'(1X,I10))'
   write(fmt2,'(A,I5,A)') '(A,1X',nUps,'(1X,F20.7))'
   write(iulog,'(2a)') new_line('a'),'** Check upstream discharge accumulation **'
   write(iulog,'(a,x,I10,x,I10)') ' Reach index & ID =', segIndex, NETOPO_in(segIndex)%REACHID
   write(iulog,'(a)')             ' * upstream reach index (NETOPO_in%UREACH) and discharge (uprflux) [m3/s] :'
   write(iulog,fmt1)              ' UREACHK =', (NETOPO_in(segIndex)%UREACHK(iUps), iUps=1,nUps)
   write(iulog,fmt2)              ' prflux  =', (RCHFLX_out(iens,NETOPO_in(segIndex)%UREACHI(iUps))%ROUTE(idxSUM)%REACH_Q, iUps=1,nUps)
   write(iulog,'(a)')             ' * local area discharge (RCHFLX_out%BASIN_QR(1)) and final discharge (RCHFLX_out%ROUTE(idxSUM)%REACH_Q) [m3/s] :'
   write(iulog,'(a,x,F15.7)')     ' RCHFLX_out%BASIN_QR(1) =', RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
   write(iulog,'(a,x,F15.7)')     ' RCHFLX_out%ROUTE(idxSUM)%REACH_Q =', RCHFLX_out(iens,segIndex)%ROUTE(idxSUM)%REACH_Q
 endif

 END SUBROUTINE accum_qupstream

END MODULE accum_runoff_module
