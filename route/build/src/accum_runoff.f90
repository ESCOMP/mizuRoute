MODULE accum_runoff_module

USE nrtype
USE public_var

! data type
USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO        ! Network topology
! subroutines: general
USE perf_mod,  ONLY: t_startf,t_stopf ! timing start/stop

implicit none

private

public::accum_runoff

CONTAINS

 ! ---------------------------------------------------------------------------------------
 ! Public subroutine main driver for basin routing
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE accum_runoff(iEns,          & ! input: index of runoff ensemble to be processed
                         river_basin,   & ! input: river basin information (mainstem, tributary outlet etc.)
                         ixDesire,      & ! input: index of verbose reach
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

 USE dataTypes,   ONLY: subbasin_omp   ! mainstem+tributary data structures
 USE model_utils, ONLY: handle_err

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

 nDom = size(river_basin)

 call t_startf('route/accum-runoff')

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

 end do ! looping through stream segments

 call t_stopf('route/accum-runoff')

 END SUBROUTINE accum_runoff

 ! *********************************************************************
 ! subroutine: perform accumulate immediate upstream flow
 ! *********************************************************************
 subroutine accum_qupstream(iEns,       &    ! input: index of runoff ensemble to be processed
                            segIndex,   &    ! input: index of reach to be processed
                            ixDesire,   &    ! input: index of verbose reach
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
  write(iulog,*) 'CHECK ACCUM_RUNOFF'
  write(iulog,*) ' UREACHK, uprflux = ', (NETOPO_in(segIndex)%UREACHK(iUps), uprflux(iUps), iUps=1,nUps)
  write(iulog,*) ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1) = ', RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
  write(iulog,*) ' RCHFLX_out%UPSTREAM_QI = ', RCHFLX_out(iens,segIndex)%UPSTREAM_QI
 endif

 end subroutine accum_qupstream

END MODULE accum_runoff_module
