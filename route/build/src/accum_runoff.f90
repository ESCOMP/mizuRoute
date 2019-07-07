MODULE accum_runoff_module

USE nrtype
USE public_var

! data type
USE dataTypes,          only : STRFLX         ! fluxes in each reach
USE dataTypes,          only : RCHTOPO        ! Network topology

implicit none

private

public::accum_runoff

CONTAINS

 ! ---------------------------------------------------------------------------------------
 ! Public subroutine main driver for basin routing
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE accum_runoff(iEns,          & ! input: index of runoff ensemble to be processed
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
 integer(i4b)                                   :: iSeg,jSeg       ! reach segment indices
 logical(lgt), allocatable                      :: doRoute(:)      ! logical to indicate which reaches are processed
 character(len=strLen)                          :: cmessage        ! error message from subroutines
 integer*8                                      :: cr                   ! rate
 integer*8                                      :: startTime,endTime    ! date/time for the start and end of the initialization
 real(dp)                                       :: elapsedTime          ! elapsed time for the process

 ierr=0; message='accum_runoff/'
 call system_clock(count_rate=cr)
 call system_clock(startTime)

 ! check
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

 ! compute the sum of all upstream runoff at each point in the river network
 do iSeg=1,nSeg

   jSeg = NETOPO_in(iSeg)%RHORDER

   if (.not. doRoute(jSeg)) cycle

   call accum_qupstream(iens, jSeg, ixDesire, NETOPO_in, RCHFLX_out, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end do ! looping through stream segments

 call system_clock(endTime)
 elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
 !write(*,"(A,1PG15.7,A)") '  elapsed-time [routing/accum] = ', elapsedTime, ' s'

 END SUBROUTINE accum_runoff

 ! *********************************************************************
 ! subroutine: perform accumulate immediate upstream flow
 ! *********************************************************************
 subroutine accum_qupstream(iEns,       &    ! input: index of runoff ensemble to be processed
                            segIndex,   &    ! input: index of reach to be processed
                            ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                            NETOPO_in,  &    ! input: reach topology data structure
                            RCHFLX_out, &    ! inout: reach flux data structure
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
 real(dp), allocatable                    :: uprflux(:)     ! upstream Reach fluxes
 integer(i4b)                             :: nUps           ! number of upstream segment
 integer(i4b)                             :: iUps           ! upstream reach index
 integer(i4b)                             :: iRch_ups       ! index of upstream reach in NETOPO
 character(len=strLen)                    :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='accum_qupstream/'

 ! identify number of upstream segments of the reach being processed
 nUps = size(NETOPO_in(segIndex)%UREACHI)

 RCHFLX_out(iEns,segIndex)%UPSTREAM_QI = RCHFLX_out(iEns,segIndex)%BASIN_QR(1)

 if (nUps>0) then

   allocate(uprflux(nUps), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': uprflux'; return; endif

   do iUps = 1,nUps
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     uprflux(iUps) = RCHFLX_out(iens,iRch_ups)%UPSTREAM_QI
   end do

   RCHFLX_out(iEns,segIndex)%UPSTREAM_QI = RCHFLX_out(iEns,segIndex)%UPSTREAM_QI + sum(uprflux)

 endif

 ! check
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
  associate( ixUpstream => NETOPO_in(segIndex)%RCHLIST)
  print*, 'ixUpstream = ', NETOPO_in(ixUpstream(1:size(ixUpstream)))%REACHIX
  print*, 'idUpstream = ', NETOPO_in(ixUpstream(1:size(ixUpstream)))%REACHID
  end associate
  print*, 'RCHFLX_out%UPSTREAM_QI = ', RCHFLX_out(iens,segIndex)%UPSTREAM_QI
 endif

 end subroutine accum_qupstream

END MODULE accum_runoff_module
