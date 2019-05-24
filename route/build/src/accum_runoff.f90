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
 SUBROUTINE accum_runoff(&
                         iEns,          & ! input: index of runoff ensemble to be processed
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

 USE nr_utility_module, ONLY : arth

 implicit none
 ! input
 integer(i4b), intent(in)                 :: iens           ! runoff ensemble index
 integer(i4b), intent(in)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)   ! River Network topology
 ! inout
 TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! input (optional)
 integer(i4b), intent(in),   optional     :: ixSubRch(:)     ! subset of reach indices to be processed
 ! local variables
 integer(i4b)                             :: nSeg           ! number of segments in the network
 integer(i4b)                             :: iSeg, jSeg     ! reach segment indices
 integer(i4b), allocatable                :: ixRch(:)       ! a list of reach indices to be processed
 integer(i4b)                             :: nUpstream      ! number of all the upstream reach segments for a segment
 integer(i4b), allocatable                :: iUpstream(:)   ! indices of all the upstream reach segments for a segment
 real(dp),     allocatable                :: qUpstream(:)   ! runoff at all the upstream reach segments for a segment
 character(len=strLen)                    :: cmessage       ! error message from subroutines

 ierr=0; message='accum_runoff/'

 ! check
 if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
  ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
 endif

 ! if a subset of reaches is processed
 if (present(ixSubRch))then
  nSeg=size(ixSubRch)
  allocate(ixRch(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for [ixRch]'; return; endif
  ixRch = ixSubRch
 ! if all the reaches are processed
 else
  nSeg = size(RCHFLX_out(iens,:))
  allocate(ixRch(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for [ixRch]'; return; endif
  ixRch = arth(1,1,nSeg)
 endif

 ! compute the sum of all upstream runoff at each point in the river network
 do iSeg=1,nSeg

   jSeg = ixRch(iSeg)

   ! identify how many reaches are upstream
   nUpstream = size(NETOPO_in(jSeg)%RCHLIST)

   ! allocate space for upstream vectors
   allocate(iUpstream(nUpstream), qUpstream(nUpstream), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': iUpstream or qUpstream'; return; endif

   ! get indices for all reaches upstream
   iUpstream(1:nUpstream) = NETOPO_in(jSeg)%RCHLIST(1:nUpstream)
   ! get streamflow for all reaches upstream
   qUpstream(1:nUpstream) = RCHFLX_out(iEns,iUpstream(1:nUpstream))%BASIN_QR(1)
   ! get mean streamflow
   RCHFLX_out(iEns,jSeg)%UPSTREAM_QI = sum(qUpstream)

   if(NETOPO_in(jSeg)%REACHIX == ixDesire)then
    print*, 'ixUpstream = ', NETOPO_in(iUpstream(1:nUpstream))%REACHIX
    print*, 'idUpstream = ', NETOPO_in(iUpstream(1:nUpstream))%REACHID
    print*, 'qUpstream = ', qUpstream
   endif

   ! deallocate space for upstream vectors
   deallocate(iUpstream, qUpstream, stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': iUpstream or qUpstream'; return; endif

 end do ! looping through stream segments

 END SUBROUTINE accum_runoff


END MODULE accum_runoff_module
