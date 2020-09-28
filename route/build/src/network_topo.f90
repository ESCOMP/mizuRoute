MODULE network_topo

! data types
USE nrtype,    only : i4b,dp,lgt
USE nrtype,    only : strLen         ! length of a string
USE dataTypes, only : var_ilength    ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength    ! double precision type: var(:)%dat

! metadata on data structures
USE globalData, only : meta_struct  ! structure information
USE globalData, only : meta_HRU     ! HRU properties
USE globalData, only : meta_HRU2SEG ! HRU-to-segment mapping
USE globalData, only : meta_SEG     ! stream segment properties
USE globalData, only : meta_NTOPO   ! network topology

! named variables
USE public_var, only : true,false   ! named integers for true/false
USE public_var, only : verySmall    ! a very small value
USE public_var, only : integerMissing

! named variables
USE var_lookup,only:ixStruct, nStructures  ! index of data structures
USE var_lookup,only:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup,only:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup,only:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup,only:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology

! external utilities
USE nr_utility_module, ONLY: findIndex     ! Num. Recipies utilities
USE nr_utility_module, ONLY: indexx        ! Num. Recipies utilities
USE nr_utility_module, ONLY: arth          ! Num. Recipies utilities

! privacy
implicit none

private
public :: hru2segment    ! compute correspondence between HRUs and segments
public :: up2downSegment ! mapping between upstream and downstream segments
public :: reachOrder     ! define the processing order
public :: streamOrdering ! get stream order for each reach
public :: reach_list     ! get a list of reaches above each reach
public :: reach_mask     ! get a mask that defines all segments above a given segment
public :: reach_mask_orig     ! get a mask that defines all segments above a given segment

contains

 ! *********************************************************************
 ! new subroutine: compute correspondence between HRUs and segments
 ! *********************************************************************
 subroutine hru2segment(&
                        ! input
                        nHRU,       &   ! input: number of HRUs
                        nSeg,       &   ! input: number of stream segments
                        ! input-output: data structures
                        structHRU,     & ! ancillary data for HRUs
                        structSeg,     & ! ancillary data for stream segments
                        structHRU2seg, & ! ancillary data for mapping hru2basin
                        structNTOPO,   & ! ancillary data for network toopology
                        ! output
                        total_hru,  &   ! output: total number of HRUs that drain into any segments
                        ierr, message)  ! output: error control
 implicit none
 ! input variables
 integer(i4b), intent(in)                      :: nHRU              ! number of HRUs
 integer(i4b), intent(in)                      :: nSeg              ! number of stream segments
 ! input-output: data structures
 type(var_dlength), intent(inout), allocatable :: structHRU(:)      ! HRU properties
 type(var_dlength), intent(inout), allocatable :: structSeg(:)      ! stream segment properties
 type(var_ilength), intent(inout), allocatable :: structHRU2seg(:)  ! HRU-to-segment mapping
 type(var_ilength), intent(inout), allocatable :: structNTOPO(:)    ! network topology
 ! output variables
 integer(i4b), intent(out)                     :: total_hru         ! total number of HRUs that drain into any segments
 integer(i4b), intent(out)                     :: ierr              ! error code
 character(*), intent(out)                     :: message           ! error message
 ! local variables
 logical(lgt),parameter          :: checkMap=.true.   ! flag to check the mapping
 character(len=strLen)           :: cmessage          ! error message of downwind routine
 integer(i4b)                    :: hruIndex          ! index of HRU (from another data structure)
 integer(i4b)                    :: iHRU              ! index of HRU
 integer(i4b)                    :: iSeg              ! index of stream segment
 integer(i4b)                    :: segId(nSeg)       ! unique identifier of the HRU
 integer(i4b)                    :: hruSegId(nHRU)    ! unique identifier of the segment where HRU drains
 integer(i4b)                    :: segHRUix(nHRU)    ! index of segment where HRU drains
 integer(i4b)                    :: nHRU2seg(nSeg)    ! number of HRUs that drain into a given segment
 real(dp)                        :: totarea           ! total area of all HRUs feeding into a given stream segment (m2)
 !integer*8                       :: cr, startTime, endTime

 ! initialize error control
 ierr=0; message='hru2segment/'
 !call system_clock(count_rate=cr)

 !print*, 'PAUSE: start of '//trim(message); read(*,*)

 ! initialize timing
 !call system_clock(startTime)

 ! ---------- get the index of the stream segment that a given HRU drains into ------------------------------

 ! get input vectors
 do iSeg=1,nSeg
  segId(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1)
 end do
 do iHRU=1,nHRU
  hruSegId(iHRU) = structHRU2seg(iHRU)%var(ixHRU2seg%hruSegId)%dat(1)
 end do

 call downReachIndex(&
                     ! input
                     nHRU,          & ! number of upstream elements
                     nSeg,          & ! number of stream segments
                     segId,         & ! unique identifier of the stream segments
                     hruSegId,      & ! unique identifier of the segment where water drains
                     ! output
                     segHRUix,      & ! index of downstream stream segment
                     nHRU2seg,      & ! number of elements that drain into each segment
                     ierr,cmessage)   ! error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! populate data structure
 do iHRU=1,nHRU
  structHRU2seg(iHRU)%var(ixHRU2seg%hruSegIndex)%dat(1) = segHRUix(iHRU)
 end do

 ! get the total number of HRUs that drain into any segments
 total_hru = sum(nHRU2seg)

 ! ---------- allocate space for the mapping structures -----------------------------------------------------

 ! loop through stream segments
 do iSeg=1,nSeg
  ! allocate space (number of elements that drain into each segment)
  allocate(structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat( nHRU2seg(iSeg) ), &
           structNTOPO(iSeg)%var(ixNTOPO%hruContribId)%dat( nHRU2seg(iSeg) ), &
           structSEG(  iSeg)%var(ixSEG%hruArea       )%dat( nHRU2seg(iSeg) ), &
           structSEG(  iSeg)%var(ixSEG%weight        )%dat( nHRU2seg(iSeg) ), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for hru2seg structure component'; return; endif
  ! initialize the number of HRUs
  structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1) = 0
 end do

 !call system_clock(endTime)
 !print*, 'timing: allocate space = ', real(endTime-startTime,kind(dp))/real(cr)

 ! ---------- populate structure components for HRU-2-Segment mapping ---------------------------------------

 ! loop through HRUs
 do iHRU=1,nHRU

  ! identify the index of the stream segment that the HRU drains into
  iSeg = segHRUix(iHRU)

  ! if there is no stream segment associated with current hru
  if (iSeg == integerMissing) cycle

  ! associate variables in data structure
  associate(nContrib       => structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1),      & ! contributing HRUs
            hruContribIx   => structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat, & ! index of contributing HRU
            hruContribId   => structNTOPO(iSeg)%var(ixNTOPO%hruContribId)%dat  ) ! unique ids of contributing HRU

  ! increment the HRU counter
  nContrib = nContrib + 1

  ! populate structure components
  hruContribIx(nContrib)   = iHRU
  hruContribId(nContrib)   = structHRU2seg(iHRU)%var(ixHRU2SEG%HRUid)%dat(1)

  ! end associations
  end associate

  ! save the HRU index
  structHRU2seg(iHRU)%var(ixHRU2SEG%HRUindex)%dat(1) = iHRU

 end do ! looping through HRUs

 ! check
 if(checkMap)then
  do iSeg=1,nSeg
   if(nHRU2seg(iSeg)/=structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1))then
    message=trim(message)//'problems identifying the HRUs draining into stream segment'
    ierr=20; return
   endif
  end do
 endif

 ! get timing
 !call system_clock(endTime)
 !print*, 'timing: populate structure components = ', real(endTime-startTime,kind(dp))/real(cr)

 ! ---------- compute additional variables ------------------------------------------------------------------

 ! loop through segments
 do iSeg=1,nSeg

  ! skip segments with no HRU drainage
  if(structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1)==0) cycle

  ! copy the HRU areas to the contrinuting HRU structure
  do iHRU=1,structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1)
   hruIndex = structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat(iHRU)
   structSEG(iSeg)%var(ixSEG%hruArea)%dat(iHRU) = structHRU(hruIndex)%var(ixHRU%area)%dat(1)
  end do

  ! compute the weights
  totarea = sum(structSEG(iSeg)%var(ixSEG%hruArea)%dat)
  structSEG(iSeg)%var(ixSEG%weight)%dat(:) = structSEG(iSeg)%var(ixSEG%hruArea)%dat(:) / totarea

 end do  ! (looping thru stream segments)

 ! get timing
 !call system_clock(endTime)
 !print*, 'timing: compute HRU weights = ', real(endTime-startTime,kind(dp))/real(cr)
 !print*, 'PAUSE: end of '//trim(message); read(*,*)

 end subroutine hru2segment


 ! *********************************************************************
 ! new subroutine: mapping between upstream and downstream segments
 ! *********************************************************************
 subroutine up2downSegment(&
                           ! input
                           nRch,         & ! input: number of stream segments
                           ! input-output: data structures
                           structNTOPO,  & ! ancillary data for network toopology
                           ! output
                           total_upseg,  & ! output: sum of immediate upstream segments
                           ierr, message)  ! output (error control)
 implicit none
 ! input variables
 integer(i4b)      , intent(in)                 :: nRch             ! number of reaches
 ! input-output: data structures
 type(var_ilength) , intent(inout), allocatable :: structNTOPO(:)   ! network topology
 ! output variables
 integer(i4b)      , intent(out)                :: total_upseg      ! sum of immediate upstream segments
 integer(i4b)      , intent(out)                :: ierr             ! error code
 character(*)      , intent(out)                :: message          ! error message
 ! local variables
 logical(lgt),parameter          :: checkMap=.true.     ! flag to check the mapping
 logical(lgt)                    :: checkDownID         ! flag to invalid downstream id
 character(len=strLen)           :: cmessage            ! error message of downwind routine
 integer(i4b)                    :: iRch                ! reach index
 integer(i4b)                    :: ixDownRch           ! index of the downstream reach
 integer(i4b)                    :: segId(nRch)         ! unique identifier of the stream segments
 integer(i4b)                    :: downSegId(nRch)     ! unique identifier of the downstream segment
 integer(i4b)                    :: downIndex(nRch)     ! index of downstream stream segment
 integer(i4b)                    :: nUpstream(nRch)     ! number of elements that drain into each segment
 integer(i4b)                    :: mUpstream(nRch)     ! number of elements that drain into each segment
 ! initialize error control
 ierr=0; message='up2downSegment/'

 ! populate reach index
 do iRch = 1,nRch
  structNTOPO(iRch)%var(ixNTOPO%segIndex)%dat(1) = iRch
 enddo

 ! ---------- define the index of the downstream reach ID ----------------------------------------------------

 checkDownID = .false.

 ! get the segid and downstream segment
 do iRch=1,nRch
  segId(iRch)     = structNTOPO(iRch)%var(ixNTOPO%segId)%dat(1)
  downSegId(iRch) = structNTOPO(iRch)%var(ixNTOPO%downSegId)%dat(1)

  ! check topology
  if (segId(iRch) == downSegId(iRch)) then
    checkDownID=.true.
    write(*,'(a,i0,a,i0)') 'Reach-ID= ', segId(iRch), ' Downstream-reach-ID= ', downSegId(iRch)
  endif

 end do

 if (checkDownID) then
  ierr=10; write(message,'(a,i0,a)') trim(message)//'reach ID and downstream ID are identical for above reaches!!!'
  return
 endif

 ! define the index of the downstream reach ID
 call downReachIndex(&
                     ! input
                     nRch,          & ! number of upstream elements
                     nRch,          & ! number of stream segments
                     segId,         & ! unique identifier of the stream segments
                     downSegId,     & ! unique identifier of the downstream segment
                     ! output
                     downIndex,     & ! index of downstream stream segment
                     nUpstream,     & ! number of elements that drain into each segment
                     ierr,cmessage)   ! error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the total number of HRUs that drain into any segments
 total_upseg = sum(nUpstream)

 ! ---------- allocate space for the number of upstream reaches ---------------------------------------------

 ! loop through the reaches
 do iRch=1,nRch
  allocate(structNTOPO(iRch)%var(ixNTOPO%upSegIds    )%dat( nUpstream(iRch) ), &
           structNTOPO(iRch)%var(ixNTOPO%upSegIndices)%dat( nUpstream(iRch) ), &
           structNTOPO(iRch)%var(ixNTOPO%goodBasin   )%dat( nUpstream(iRch) ), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for upstream reaches'; return; endif
 end do

 ! ---------- populate data structures for the upstream reaches ----------------------------------------------

 ! initialize the number of upstream elements in each reach
 mUpstream(:)=0

 ! loop through the reaches
 do iRch=1,nRch

  ! identify the index of the downstream segment
  ixDownRch = downIndex(iRch)
  if(ixDownRch == integerMissing) cycle

  ! increment the number of upstream elements in the downstream segment
  mUpstream(ixDownRch) = mUpstream(ixDownRch)+1
  if(mUpstream(ixDownRch)>nUpstream(ixDownRch))then
   message=trim(message)//'upstream index exceeds dimension'
   ierr=20; return
  endif

  ! populate the structure components
  structNTOPO(ixDownRch)%var(ixNTOPO%upSegIndices)%dat( mUpstream(ixDownRch) ) = structNTOPO(iRch)%var(ixNTOPO%segIndex)%dat(1)
  structNTOPO(ixDownRch)%var(ixNTOPO%upSegIds    )%dat( mUpstream(ixDownRch) ) = structNTOPO(iRch)%var(ixNTOPO%segId   )%dat(1)

 end do  ! looping through reaches

 ! set missing values to -1
 ! NOTE: check if the -1 is special and if the replacement with integer missing is necessary
 where(downIndex==integerMissing) downIndex=-1

 ! populate data structures
 do iRch=1,nRch
  structNTOPO(iRch)%var(ixNTOPO%downSegIndex)%dat(1) = downIndex(iRch)
 end do

 end subroutine up2downSegment

 ! *********************************************************************
 ! new subroutine: define index of downstream reach
 ! *********************************************************************
 subroutine downReachIndex(&
                           ! input
                           nUp,          & ! number of upstream elements
                           nSeg,         & ! number of stream segments
                           segId,        & ! unique identifier of the stream segments
                           downId,       & ! unique identifier of the segment where water drains
                           ! output
                           downSegIndex, & ! index of downstream stream segment
                           nElement2Seg, & ! number of elements that drain into each segment
                           ierr,message)
 ! external modules
 !USE nr_utility_module, ONLY: indexx  ! Num. Recipies utilities
 implicit none
 ! input variables
 integer(i4b), intent(in)        :: nUp             ! number of upstream elements
 integer(i4b), intent(in)        :: nSeg            ! number of stream segments
 integer(i4b), intent(in)        :: segId(:)        ! unique identifier of the stream segments
 integer(i4b), intent(in)        :: downId(:)       ! unique identifier of the segment where water drains
 ! output variables
 integer(i4b), intent(out)       :: downSegIndex(:) ! index of downstream stream segment
 integer(i4b), intent(out)       :: nElement2Seg(:) ! number of elements that drain into each segment
 integer(i4b), intent(out)       :: ierr            ! error code
 character(*), intent(out)       :: message         ! error message
 ! local variables
 integer(i4b)                    :: iUp                 ! index of upstream element
 integer(i4b)                    :: iSeg,jSeg           ! index of stream segment
 integer(i4b)                    :: rankSegId           ! ranked Id of the stream segment
 integer(i4b)                    :: rankDownId          ! ranked Id of the downstream stream segment
 integer(i4b)                    :: rankSeg(nSeg)       ! rank index of each segment in the nRch vector
 integer(i4b)                    :: rankDownSeg(nUp)    ! rank index of each downstream stream segment
 logical(lgt),parameter          :: checkLink=.false.   ! flag to check the links
 logical(lgt),parameter          :: checkMap=.true.     ! flag to check the mapping
 integer(i4b),parameter          :: nProgress=100000    ! print every nProgress step

 ! initialize error control
 ierr=0; message='downReachIndex/'

 ! initialize output
 nElement2Seg(:) = 0
 downSegIndex(:) = integerMissing

 ! rank the ids of the stream segments
 call indexx(segId, rankSeg)

 ! rank the ids of the downstream stream segments
 call indexx(downId, rankDownSeg)

 iSeg=1  ! second counter
 ! loop through the upstream elements
 do iUp=1,nUp

  ! print progress
  if(mod(iUp,nProgress)==0) print*, 'Getting downstream link for reach: ', iUp, nUp

  ! get Ids for the up2seg mapping vector
  rankDownId = downId( rankDownSeg(iUp) )

  ! check if part of the network
  if (rankDownId<=0) cycle ! upstream element does not drain into any stream segment (closed basin or coastal HRU)

  ! keep going until found the index
  do jSeg=iSeg,nSeg ! normally a short loop

   ! get Id of the stream segment
   rankSegId = segId( rankSeg(jSeg) )

   ! check
   if(jSeg>iSeg+100 .and. checkLink)then
    print*, 'iUp, iSeg, jSeg, rankDownId, rankSegId, rankDownSeg(iUp), rankSeg(jSeg) = ', &
             iUp, iSeg, jSeg, rankDownId, rankSegId, rankDownSeg(iUp), rankSeg(jSeg)
    if(jSeg>iSeg+200)then
     print*, trim(message)//'PAUSE : '; read(*,*)
    endif
   endif

   ! check if the basin is missing
   if(rankSegId>rankDownId) exit ! exit the segment loop

   ! define the index where we have a match
   if(rankDownId==rankSegId)then

    ! identify the index of the segment that the HRU drains into
    downSegIndex( rankDownSeg(iUp) ) = rankSeg(jSeg)
    nElement2Seg( rankSeg(jSeg)    ) = nElement2Seg( rankSeg(jSeg) ) + 1

    ! check if we should increment the stream segment
    ! NOTE: we can have multiple upstream elements draining into the same segment
    !        --> in this case, we want to keep the segment the same
    if(iUp<nUp .and. jSeg<nSeg)then
     if(downId( rankDownSeg(iUp+1) ) >= segId( rankSeg(jSeg+1) ) ) iSeg=jSeg+1
    endif

    ! identified the segment so exit the segment loop and evaluate the next upstream element
    exit

   endif  ! match between the upstream drainage segment and the stream segment
  end do  ! skipping segments that have no input

 end do  ! looping through upstream elements

 ! check
 if(checkMap)then
  do iUp=1,nUp
   if(downId(iUp)<=0 .or. downSegIndex(iUp)==integerMissing) cycle
   if(downId(iUp) /= segId( downSegIndex(iUp) ) )then
    message=trim(message)//'problems identifying the index of the stream segment that a given HRU drains into'
    ierr=20; return
   endif
  end do
 endif

 end subroutine downReachIndex

 ! *********************************************************************
 ! subroutine: define processing order for the individual
 !                 stream segments in the river network
 ! *********************************************************************
 subroutine REACHORDER(NRCH,         &   ! input:        number of reaches
                       structNTOPO,  &   ! input:output: network topology
                       ierr, message)    ! output:       error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Defines the processing order for the individual stream segments in the river network
 !
 ! ----------------------------------------------------------------------------------------
 IMPLICIT NONE
 ! input variables
 INTEGER(I4B), INTENT(IN)               :: NRCH            ! number of stream segments
 type(var_ilength) , intent(inout)      :: structNTOPO(:)  ! network topology structure
 ! output variables
 integer(i4b), intent(out)              :: ierr            ! error code
 character(*), intent(out)              :: message         ! error message
 ! local variables
 INTEGER(I4B)                           :: IRCH,JRCH,KRCH  ! loop through reaches
 INTEGER(I4B)                           :: IUPS            ! loop through upstream reaches
 INTEGER(I4B)                           :: ICOUNT          ! counter for the gutters
 INTEGER(I4B)                           :: NASSIGN         ! # reaches currently assigned
 LOGICAL(LGT),DIMENSION(:),ALLOCATABLE  :: RCHFLAG         ! TRUE if reach is processed
 INTEGER(I4B)                           :: NUPS            ! number of upstream reaches
 INTEGER(I4B)                           :: UINDEX          ! upstream reach index
 integer(i4b)                           :: jCount          ! counter
 ! initialize error control
 ierr=0; message='reachorder/'
 ! ----------------------------------------------------------------------------------------
 NASSIGN = 0
 ALLOCATE(RCHFLAG(NRCH),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for RCHFLAG'; return; endif
 RCHFLAG(1:NRCH) = .FALSE.
 ! ----------------------------------------------------------------------------------------
 ICOUNT=0
 jCount=0
 DO  ! do until all reaches are assigned
  NASSIGN = 0
  DO IRCH=1,NRCH
   ! print progress
   if(jCount>0)then
    if(mod(jCount,100000)==0) print*, 'Getting order for reach: ', count(RCHFLAG), nRch
   endif
   ! check if the reach is assigned yet
   IF(RCHFLAG(IRCH)) THEN
    NASSIGN = NASSIGN + 1
    CYCLE
   ENDIF
   ! climb upstream as far as possible
   JRCH = IRCH    ! the first reach under investigation
   DO  ! do until get to a "most upstream" reach that is not assigned
    nUps = size(structNTOPO(jRch)%var(ixNTOPO%upSegIds)%dat)
    !print*, 'iRch, nUps = ', iRch, nUps
    IF (NUPS.GE.1) THEN     ! (if NUPS = 0, then it is a first-order basin)
     KRCH = JRCH   ! the reach under investigation
     ! loop through upstream reaches
     DO IUPS=1,NUPS
      jCount = jCount+1
      uIndex = structNTOPO(jRch)%var(ixNTOPO%upSegIndices)%dat(iUps)  ! index of the upstream reach
      !print*, 'jRch, uIndex = ', jRch, uIndex, structNTOPO(jRch)%var(ixNTOPO%upSegIds)%dat
      ! check if the reach is NOT assigned
      IF (.NOT.RCHFLAG(UINDEX)) THEN
       JRCH = UINDEX
       EXIT    ! exit IUPS loop
      END IF  ! if the reach is assigned
     END DO  ! (looping through upstream reaches)
     ! check if all upstream reaches are already assigned (only if KRCH=JRCH)
     IF (JRCH.EQ.KRCH) THEN
      ! assign JRCH
      ICOUNT=ICOUNT+1
      RCHFLAG(JRCH) = .TRUE.
      structNTOPO(iCount)%var(ixNTOPO%rchOrder)%dat(1) = jRch
      EXIT
     ENDIF
     CYCLE   ! if jrch changes, keep looping (move upstream)
    ELSE    ! if the reach is a first-order basin
     ! assign JRCH
     ICOUNT=ICOUNT+1
     RCHFLAG(JRCH) = .TRUE.
     structNTOPO(iCount)%var(ixNTOPO%rchOrder)%dat(1) = jRch
     EXIT
    ENDIF
   END DO   !  climbing upstream (do-forever)
  END DO   ! looping through reaches
  IF (NASSIGN.EQ.NRCH) EXIT
 END DO  ! do forever (do until all reaches are assigned)
 DEALLOCATE(RCHFLAG,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for RCHFLAG'; return; endif
 ! test
 !print*, 'reachOrder = ', (structNTOPO(jRch)%var(ixNTOPO%rchOrder)%dat(1), jRch=1,nRch)
 ! ----------------------------------------------------------------------------------------
 end subroutine reachOrder

 ! *********************************************************************
 ! public subroutine: compute Shreve stream order for each segment
 ! *********************************************************************
 subroutine streamOrdering(nSeg,         &   ! input:        number of reaches
                           structNTOPO,  &   ! input:output: network topology
                           ierr, message)    ! output:       error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Defines the shreve stream order for the individual stream segments in the river network
 !   Starting 1 at most upstream segments, sum all orders of immediate upstream reaches.
 !
 !   Variables that needs to be computed
 !   ixNTOPO%upSegIndices, ixNTOPO%rchOrder
 !
 ! ----------------------------------------------------------------------------------------
 implicit none
 ! input variables
 integer(i4b), intent(in)          :: nSeg           ! number of stream segments
 type(var_ilength) , intent(inout) :: structNTOPO(:) ! network topology structure
 ! output variables
 integer(i4b), intent(out)         :: ierr           ! error code
 character(*), intent(out)         :: message        ! error message
 ! local variables
 integer(i4b)                      :: iSeg,jSeg,kSeg ! loop through reaches
 integer(i4b)                      :: iUps           ! loop through immeidate upstream reach
 integer(i4b)                      :: strmOrder      ! maximum stream order
 integer(i4b)                      :: nImmediate     ! number of immediate upstream reaches
 logical(lgt)                      :: increment      ! increment order or not

 ierr=0; message='streamOrdering/'

 do iSeg = 1,nSeg

   jSeg = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)

   nImmediate = size(structNTOPO(jSeg)%var(ixNTOPO%upSegIndices)%dat)

   ! Initialize streamorder as 1 for most upstream reach, otherwise 0
   if (nImmediate<1) then

     structNTOPO(jSeg)%var(ixNTOPO%streamOrder)%dat(1) = 1

   elseif (nImmediate==1) then

     kSeg = structNTOPO(jSeg)%var(ixNTOPO%upSegIndices)%dat(nImmediate)
     structNTOPO(jSeg)%var(ixNTOPO%streamOrder)%dat(1) = structNTOPO(kSeg)%var(ixNTOPO%streamOrder)%dat(1)

   else
     ! take order from one of upstream reach for initialization
     kSeg = structNTOPO(jSeg)%var(ixNTOPO%upSegIndices)%dat(1)
     strmOrder = structNTOPO(kSeg)%var(ixNTOPO%streamOrder)%dat(1)
     increment = .true.
     do iUps = 2,nImmediate
       kSeg = structNTOPO(jSeg)%var(ixNTOPO%upSegIndices)%dat(iUps)
       ! if the other upstream reaches are different, assign higher one and no need to increment order
       if (structNTOPO(kSeg)%var(ixNTOPO%streamOrder)%dat(1) /= strmOrder) then
         increment = .False.
         if (structNTOPO(kSeg)%var(ixNTOPO%streamOrder)%dat(1) > strmOrder) then
           strmOrder = structNTOPO(kSeg)%var(ixNTOPO%streamOrder)%dat(1)
         endif
       endif
     enddo

     structNTOPO(jSeg)%var(ixNTOPO%streamOrder)%dat(1) = strmOrder
     if (increment) then
       structNTOPO(jSeg)%var(ixNTOPO%streamOrder)%dat(1) = structNTOPO(jSeg)%var(ixNTOPO%streamOrder)%dat(1) + 1
     endif

   endif

  enddo

 end subroutine streamOrdering


 ! *********************************************************************
 ! new subroutine: identify all reaches above the current reach
 ! *********************************************************************
 SUBROUTINE REACH_LIST(&
                       ! input
                       NRCH,        & ! Number of reaches
                       doReachList, & ! flag to compute the list of upstream reaches
                       structNTOPO, & ! Network topology
                       ! output
                       structSEG,   & ! Reach properties
                       NTOTAL,      & ! Total number of upstream reaches for all reaches
                       ierr,message)  ! Error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Generates a list of all reaches upstream of each reach (used to compute total runoff
 !     at each point in the river network)
 !
 ! ----------------------------------------------------------------------------------------
 USE nr_utility_module, ONLY : arth                          ! Num. Recipies utilities
 IMPLICIT NONE
 ! input variables
 INTEGER(I4B)      , INTENT(IN)                 :: NRCH            ! number of stream segments
 logical(lgt)      , intent(in)                 :: doReachList     ! flag to compute the list of upstream reaches
 type(var_ilength) , intent(inout), allocatable :: structNTOPO(:)  ! network topology structure
 ! output variables
 type(var_dlength) , intent(inout)              :: structSEG(:)    ! reach properties structure
 integer(i4b)      , intent(out)                :: NTOTAL          ! total number of upstream reaches for all reaches
 integer(i4b)      , intent(out)                :: ierr            ! error code
 character(*)      , intent(out)                :: message         ! error message
 ! local variables
 INTEGER(I4B)                                   :: IRCH,JRCH,KRCH  ! loop through reaches
 logical(lgt)                                   :: processedReach(nRch)  ! flag to define if reaches are processed
 integer(i4b)                                   :: nImmediate      ! number of immediate upstream reaches
 integer(i4b)                                   :: nUpstream       ! total number of upstream reaches
 integer(i4b)                                   :: iUps            ! index of upstream reaches
 integer(i4b)                                   :: iPos            ! position in vector
 integer(i4b),parameter                         :: nProgress=100000! print every nProgress step
 ! ----------------------------------------------------------------------------------------
 message='REACH_LIST/'

 ! check if the list of upstream reaches is desired
 if(.not.doReachList)then
  ! allocate zero length vector in each reach
  do iRch=1,nRch
   allocate(structNTOPO(iRch)%var(ixNTOPO%allUpSegIndices)%dat(0), stat=ierr)
   if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for RCHLIST'; return; endif
  end do
  ! early retrun
  NTOTAL=0
  return
 endif

 ! initialize
 nTotal = 0                   ! total number of upstream reaches
 processedReach(:) = .false.  ! check that we processed the reach already

 ! Loop through reaches
 do kRch=1,nRch

  ! ---------- identify reach in the ordered vector ----------------------------------------

  ! print progress
  if(mod(kRch,nProgress)==0) print*, 'Getting list of all upstream reaches: kRch, nRch = ', kRch, nRch

  ! NOTE: Reaches are ordered
  !        -->  kRch cannpt be processed until all upstream reaches are processed
  iRch = structNTOPO(kRch)%var(ixNTOPO%rchOrder)%dat(1)  ! reach index
  processedReach(iRch)=.true.

  ! ---------- allocate space for the list of all upstream reaches  ------------------------

  ! get the number of upstream reaches
  nUpstream  = 1  ! count the reach itself
  nImmediate = size(structNTOPO(iRch)%var(ixNTOPO%upSegIndices)%dat)

  if(nImmediate > 0)then
   do iUps=1,nImmediate ! get the upsteam segments

    ! get upstream reach
    jRch = structNTOPO(iRch)%var(ixNTOPO%upSegIndices)%dat(iUps)
    if(.not.processedReach(jRch))then
     write(message,'(a,i0,a)') trim(message)//'expect reach index ', jRch, ' to be processed already'
     ierr=20; return
    endif

    ! update size
    nUpstream = nUpstream + size(structNTOPO(jRch)%var(ixNTOPO%allUpSegIndices)%dat)

   end do  ! looping through immediate upstream reaches
  endif   ! if upstream segments exist

  ! allocate space
  allocate(structNTOPO(iRch)%var(ixNTOPO%allUpSegIndices)%dat(nUpstream), stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for RCHLIST'; return; endif

  ! update total upstream
  nTotal = nTotal + nUpstream

  ! ---------- add the list of upstream reaches to the current reach -----------------------

  ! initialize the upstream area
  structSEG(iRch)%var(ixSEG%upsArea)%dat(1) = 0._dp

  ! add the current segment
  iPos=1
  structNTOPO(iRch)%var(ixNTOPO%allUpSegIndices)%dat(iPos) = iRch

  ! add the upstream segments, if they exist
  if(nImmediate > 0)then
   do iUps=1,nImmediate ! get the upstream segments

    ! get the list of reaches
    jRch      = structNTOPO(iRch)%var(ixNTOPO%upSegIndices)%dat(iUps)
    nUpstream = size(structNTOPO(jRch)%var(ixNTOPO%allUpSegIndices)%dat)
    structNTOPO(iRch)%var(ixNTOPO%allUpSegIndices)%dat(iPos+1:iPos+nUpstream) = &
    structNTOPO(jRch)%var(ixNTOPO%allUpSegIndices)%dat(     1:     nUpstream)
    iPos = iPos + nUpstream

    ! get the upstream area (above the top of the reach)
    structSEG(iRch)%var(ixSEG%upsArea)%dat(1) = structSEG(iRch)%var(ixSEG%upsArea)%dat(1) + &
                                                structSEG(jRch)%var(ixSEG%totalArea)%dat(1)


   end do  ! looping through immediate upstream segments
  endif   ! if immediate upstream segments exist

  ! ---------- compute the upstream area ---------------------------------------------------

  ! compute the local basin area
  if(structNTOPO(iRch)%var(ixNTOPO%nHRU)%dat(1) > 0)then
   structSEG(iRch)%var(ixSEG%basArea)%dat(1) = sum(structSEG(iRch)%var(ixSEG%hruArea)%dat)
  else
   structSEG(iRch)%var(ixSEG%basArea)%dat(1) = 0._dp
  endif

  ! compute the total area
  structSEG(iRch)%var(ixSEG%totalArea)%dat(1) = structSEG(iRch)%var(ixSEG%basArea)%dat(1) + &
                                                structSEG(iRch)%var(ixSEG%upsArea)%dat(1)

  ! ---------- define good basins ---------------------------------------------------------

  ! "goodBasin" == Drainage area at the bottom of the immediate reaches is greater than threshold, "verySmall" value.

  if(nImmediate > 0)then
   do iUps=1,nImmediate ! get the upstream segments
    structNTOPO(iRch)%var(ixNTOPO%goodBasin)%dat(iUps) = merge(true, false, structSEG(iRch)%var(ixSEG%totalArea)%dat(1) > verySmall)
   end do
  endif

 end do  ! looping through reaches


 END SUBROUTINE REACH_LIST

 ! *********************************************************************
 ! new subroutine: identify all reaches above a given reach
 ! *********************************************************************
 SUBROUTINE reach_mask(&
                       ! input
                       outletId,      &  ! input: outlet reach id
                       structNTOPO,   &  ! input: network topology structures
                       structSeg,     &  ! input: Reach property structures
                       nHRU,          &  ! input: number of HRUs
                       nRch,          &  ! input: number of reaches
                       ! output: updated dimensions
                       tot_hru,       &  ! input+output: total number of all the upstream hrus for all stream segments
                       tot_upseg,     &  ! input+output: sum of immediate upstream segments
                       tot_upstream,  &  ! input+output: total number of upstream reaches for all reaches
                       tot_uh,        &  ! input+output: total number of unit hydrograph dimensions
                       ! output: dimension masks
                       ixHRU_desired, &  ! output: indices of desired hrus
                       ixSeg_desired, &  ! output: indices of desired reaches
                       ! output: error control
                       ierr, message )   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Generates a list of all reaches upstream of a given reach
 !
 ! ----------------------------------------------------------------------------------------
 USE nrtype
 USE nr_utility_module, ONLY : arth                                 ! Num. Recipies utilities
 IMPLICIT NONE
 ! input variables
 integer(i4b)      , intent(in)                :: outletId          ! id of the outlet reach
 type(var_ilength) , intent(inout)             :: structNTOPO(:)    ! network topology structure
 type(var_dlength) , intent(in)                :: structSeg(:)      ! stream segment properties
 integer(i4b)      , intent(in)                :: nHRU              ! number of HRUs
 integer(i4b)      , intent(in)                :: nRch              ! number of reaches
 ! input+output: updated dimensions
 integer(i4b)      , intent(inout)             :: tot_hru           ! total number of all the upstream hrus for all stream segments
 integer(i4b)      , intent(inout)             :: tot_upseg         ! sum of immediate upstream segments
 integer(i4b)      , intent(inout)             :: tot_upstream      ! total number of upstream reaches for all reaches
 integer(i4b)      , intent(inout)             :: tot_uh            ! total number of unit hydrograph dimensions
 ! output: dimension masks
 integer(i4b)      , intent(out) , allocatable :: ixHRU_desired(:)  ! indices of desired hrus
 integer(i4b)      , intent(out) , allocatable :: ixSeg_desired(:)  ! indices of desired reaches
 ! output: error control
 integer(i4b)      , intent(out)               :: ierr              ! error code
 character(*)      , intent(out)               :: message           ! error message
 ! ----------------------------------------------------------------------------------------
 ! general local variables
 integer(i4b)                                  :: nHRU_desire       ! num desired HRUs
 integer(i4b)                                  :: nRch_desire       ! num desired reaches
 INTEGER(I4B)                                  :: iRch              ! loop through reaches
 integer(i4b)                                  :: ixHRU_map(nHRU)   ! mapping indices for HRUs
 integer(i4b)                                  :: idSeg_vec(nRch)   ! unique ids for stream segments
 integer(i4b)                    , allocatable :: ixRank(:)         ! indices of a sorted vector
 integer(i4b)                                  :: ixDesire          ! index of desired reach
 integer(i4b)                                  :: jxDesire          ! index of desired reach
 ! ----------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='reach_mask/'

 ! check if we actually want the mask
 if(outletId<0)then

  ! ---------- case 1: no mask desired ---------------------------------------------------------------------------------------------------

  ! allocate space
  allocate(ixHRU_desired(nHRU), ixSeg_desired(nRch), stat=ierr)
  if(ierr/=0) message=trim(message)//'unable to allocate space for the vectors of desired reaches'

  ! include every index
  ixHRU_desired = arth(1,1,nHRU)
  ixSeg_desired = arth(1,1,nRch)
  return

 ! mask is desired
 else

  ! ---------- case 2: mask desired ------------------------------------------------------------------------------------------------------

  ! initialize the HRU map
  ixHRU_map(:) = integerMissing

  ! find the index of the reach
  do iRch = 1,nRch
   idSeg_vec(iRch) = structNTOPO(iRch)%var(ixNTOPO%segId)%dat(1)
  end do
  ixDesire = findIndex(idSeg_vec,outletId,integerMissing)
  if(ixDesire==integerMissing)then
   message=trim(message)//'unable to find index of desired reach id'
   ierr=20; return
  endif

  ! set the downstream reach to missing
  structNTOPO(ixDesire)%var(ixNTOPO%downSegId)%dat(1)    = integerMissing
  structNTOPO(ixDesire)%var(ixNTOPO%downSegIndex)%dat(1) = integerMissing

  ! allocate space for the number of stream segments
  nRch_desire = size(structNTOPO(ixDesire)%var(ixNTOPO%allUpSegIndices)%dat)
  allocate(ixSeg_desired(nRch_desire), ixRank(nRch_desire), stat=ierr)
  if(ierr/=0) message=trim(message)//'unable to allocate space for the vectors of desired reaches'

  ! sort the list of upstream reaches
  call indexx(structNTOPO(ixDesire)%var(ixNTOPO%allUpSegIndices)%dat, ixRank)

  ! loop through desired stream segments
  do iRch=1,nRch_desire

   ! get the desired stream segments
   ixSeg_desired(iRch) = structNTOPO(ixDesire)%var(ixNTOPO%allUpSegIndices)%dat( ixRank(iRch) )

   ! get the desired HRUs
   jxDesire = ixSeg_desired(iRch)
   if(structNTOPO(jxDesire)%var(ixNTOPO%nHRU)%dat(1)>0)then
    ixHRU_map( structNTOPO(jxDesire)%var(ixNTOPO%hruContribIx)%dat(:) ) = &
               structNTOPO(jxDesire)%var(ixNTOPO%hruContribIx)%dat(:)
   endif

  end do  ! looping through desired stream segments

  ! allocate space for the number of HRUs
  nHRU_desire = count(ixHRU_map/=integerMissing)
  allocate(ixHRU_desired(nHRU_desire), stat=ierr)
  if(ierr/=0) message=trim(message)//'unable to allocate space for the vectors of desired HRUs'

  ! pack the desired indices into the HRU vector
  ixHRU_desired(:) = pack(ixHRU_map, ixHRU_map/=integerMissing)

  tot_hru       = 0  ! total number of all the upstream hrus for all stream segments
  tot_upseg     = 0  ! sum of immediate upstream segments
  tot_upstream  = 0  ! total number of upstream reaches for all reaches
  tot_uh        = 0  ! total number of unit hydrograph dimensions
  ! get the updated dimensions
  do iRch=1,nRch_desire
   tot_hru      = tot_hru      + structNTOPO( ixSeg_desired(iRch) )%var(ixNTOPO%nHRU)%dat(1)
   tot_upseg    = tot_upseg    + size(structNTOPO( ixSeg_desired(iRch) )%var(ixNTOPO%upSegIds)%dat)
   tot_upstream = tot_upstream + size(structNTOPO( ixSeg_desired(iRch) )%var(ixNTOPO%allUpSegIndices)%dat)
   tot_uh       = tot_uh       + size(structSeg(   ixSeg_desired(iRch) )%var(ixSEG%timeDelayHist)%dat)
  end do

 endif  ! if the mask is desired

 end subroutine reach_mask

 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================
 ! ====================================================================================================

 ! ROUTINES BELOW HERE ARE NOT USED
 ! THEY ARE INCLUDED IN THE MODULE SINCE THE CONCEPTS COULD BE USEFUL

 ! *********************************************************************
 ! new subroutine: identify all reaches above a given reach
 ! *********************************************************************
 SUBROUTINE REACH_MASK_ORIG(&
                            ! input
                            desireId,      &  ! input: reach index
                            structNTOPO,   &  ! input: network topology structures
                            nHRU,          &  ! input: number of HRUs
                            nRch,          &  ! input: number of reaches
                            ! output: updated dimensions
                            tot_hru,       &  ! input+output: total number of all the upstream hrus for all stream segments
                            tot_upseg,     &  ! input+output: sum of immediate upstream segments
                            tot_upstream,  &  ! input+output: total number of upstream reaches for all reaches
                            ! output: dimension masks
                            ixHRU_desired, &  ! output: indices of desired hrus
                            ixSeg_desired, &  ! output: indices of desired reaches
                            ! output: error control
                            ierr, message )   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Generates a list of all reaches upstream of a given reach
 !
 ! ----------------------------------------------------------------------------------------
 USE nrtype
 USE nr_utility_module, ONLY : arth                                 ! Num. Recipies utilities
 IMPLICIT NONE
 ! input variables
 integer(i4b)      , intent(in)                :: desireId          ! id of the desired reach
 type(var_ilength) , intent(in)                :: structNTOPO(:)    ! network topology structure
 integer(i4b)      , intent(in)                :: nHRU              ! number of HRUs
 integer(i4b)      , intent(in)                :: nRch              ! number of reaches
 ! input+output: updated dimensions
 integer(i4b)      , intent(inout)             :: tot_hru           ! total number of all the upstream hrus for all stream segments
 integer(i4b)      , intent(inout)             :: tot_upseg         ! sum of immediate upstream segments
 integer(i4b)      , intent(inout)             :: tot_upstream      ! total number of upstream reaches for all reaches
 ! output: dimension masks
 integer(i4b)      , intent(out) , allocatable :: ixHRU_desired(:)  ! indices of desired hrus
 integer(i4b)      , intent(out) , allocatable :: ixSeg_desired(:)  ! indices of desired reaches
 ! output: error control
 integer(i4b)      , intent(out)               :: ierr              ! error code
 character(*)      , intent(out)               :: message           ! error message
 ! ----------------------------------------------------------------------------------------
 ! general local variables
 integer(i4b)      , parameter                 :: startReach=-9999  ! starting reach to test
 logical(lgt)                                  :: checkList         ! flag to check the list
 integer(i4b)                                  :: nHRU_desire       ! num desired HRUs
 integer(i4b)                                  :: nRch_desire       ! num desired reaches
 INTEGER(I4B)                                  :: IRCH,JRCH,KRCH    ! loop through reaches
 logical(lgt)                                  :: isTested(nRch)    ! flags to define that the reach has been tested
 logical(lgt)                                  :: isDesired(nRch)   ! flags to define that the reach is desired
 logical(lgt)                                  :: isWithinBasin     ! flag to denote that the reaches are within the basin
 integer(i4b)                                  :: ixHRU_map(nHRU)   ! mapping indices for HRUs
 integer(i4b)                                  :: ixSeg_map(nRch)   ! mapping indices for stream segments
 character(len=strLen)                         :: cmessage          ! error message of downwind routine
 ! structure in a linked list
 TYPE NODE_STRUCTURE
  INTEGER(I4B)                                 :: IX                ! index of downstream reach
  TYPE(NODE_STRUCTURE), POINTER                :: NEXT              ! next node in linked list
 END TYPE NODE_STRUCTURE
 ! a structure object in a linked list
 TYPE(NODE_STRUCTURE), POINTER                 :: NODE              ! node in linked list
 TYPE(NODE_STRUCTURE), POINTER                 :: HPOINT            ! Head pointer in linked list
 ! vectors of downstream reaches
 integer(i4b)                                  :: nDown             ! number d/s reaches
 integer(i4b)                                  :: ixDesire          ! index of desired reach
 integer(i4b),allocatable                      :: ixDownstream(:)   ! indices of downstream reaches
 ! ----------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='REACH_MASK/'

 ! check if we actually want the mask
 if(desireId<0)then
  ! allocate space
  allocate(ixHRU_desired(nHRU), ixSeg_desired(nRch), stat=ierr)
  if(ierr/=0) message=trim(message)//'unable to allocate space for the vectors of desired reaches'
  ! include every index
  ixHRU_desired = arth(1,1,nHRU)
  ixSeg_desired = arth(1,1,nRch)
  return
 endif

 ! initialize vectors
 ixHRU_map(:) = integerMissing  ! set to the HRU index if the HRU is desired
 ixSeg_map(:) = integerMissing  ! set to the reach index if the reach is desired
 isDesired(:) = .false.         ! .true. if reach is desired
 isTested(:)  = .false.         ! .true. if we have processed a reach already

 ! loop through all reaches
 DO KRCH=1,nRch

  ! print progress
  !if(mod(kRch,100000)==0) print*, 'Getting reach subset: kRch, nRch = ', kRch, nRch

  ! check the list
  checkList = (structNTOPO(kRch)%var(ixNTOPO%segId)%dat(1)==startReach)
  if(checkList) print*, 'Checking the list for reach ', startReach

  ! skip if we have already tested krch
  if(isTested(kRch)) cycle

  ! ---------- get a vector of reaches downstream of a given reach -----------------------------------------------------------------------

  ! initialise the reach array
  nDown = 0        ! initialize the number of upstream reaches
  NULLIFY(HPOINT)  ! set pointer to a linked list to NULL

  ! initialize desired index
  ixDesire      = integerMissing
  isWithinBasin = .false.

  ! ensure take streamflow from surrounding basin (a reach is upstream of itself!)
  ! but only if reach is one of the selected
  CALL ADD2LIST(KRCH,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! climb down the network
  IRCH = KRCH
  DO  ! (do-forever)

   ! index of downstream reach
   JRCH = structNTOPO(iRch)%var(ixNTOPO%downSegIndex)%dat(1)
   if(checkList) print*, 'Downstream id = ', structNTOPO(iRch)%var(ixNTOPO%downSegId)%dat(1)

   ! check if reached the outlet
   if(jRch<=0) exit         ! negative = missing, which means that the reach is the outlet

   ! check if the downstream reach is already tested
   ! NOTE: this means that all reaches downstream of this reach have also been tested
   if(isTested(jRch))then
    isWithinBasin = isDesired(jRch)   ! reach is within basin if isDesired(jRch) is true
    exit ! can exit here and accept the full linked list
   endif

   ! add downstream reach to the list
   CALL ADD2LIST(JRCH,ierr,cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! move the index downstream
   IRCH = JRCH

  END DO  ! do forever (go to the outlet)

  ! accept the full linked list
  if(isWithinBasin) ixDesire = nDown

  ! check
  if(checkList)then
   print*, trim(message)//'pause: '; read(*,*)
  endif

  ! ---------- extract the vector from the linked list -----------------------------------------------------------------------------------

  ! allocate space
  allocate(ixDownstream(nDown),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for downstream reaches'; return; endif

  ! set node to first node in the list of reaches
  node => HPOINT

  jRch = nDown
  ! keep going while node still points to something
  do while (associated(node))

   ! extract the reach index from the linked list
   ixDownstream(jRch) = node%ix

   ! if desired reach is contained in downstream vector, then set the index in downstream vector
   if(structNTOPO(node%ix)%var(ixNTOPO%segId)%dat(1)==desireId) ixDesire=jRch

   ! point to the next node
   node=>node%next

   ! remove the previous node
   deallocate(HPOINT, STAT=IERR)
   if(ierr/=0)then; ierr=20; message=trim(message)//'problem deallocating space for linked list'; return; endif

   ! reset head pointer
   HPOINT=>node
   jRch = jRch-1

  end do  ! while

  ! set the flags to define that these reaches have been tested
  isTested(ixDownstream) = .true.

  ! set the flags to denote the reach is desired
  if(ixDesire/=integerMissing)then
   isDesired(ixDownstream(1:ixDesire)) = .true.
   ixSeg_map(ixDownstream(1:ixDesire)) = ixDownstream(1:ixDesire)
  endif  ! if the reach is desired

  ! deallocate space
  deallocate(ixDownstream,stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for downstream reaches'; return; endif

 END DO  ! assess each reach

 ! ---------- check for errors ----------------------------------------------------------------------------------------------------------

 ! get the counts
 nRch_desire = count(isDesired)

 ! check that we processed all reaches
 if(count(isTested)/=nRch)then
  message=trim(message)//'did not process all reaches'
  ierr=20; return
 endif

 ! check that the desired reach exists
 if(nRch_desire==0)then
  message=trim(message)//'desired reach does not exist: NOTE: reach Ids must be >0'
  ierr=20; return
 endif

 ! ---------- get the subset of indices -------------------------------------------------------------------------------------------------

 ! allocate space
 allocate(ixSeg_desired(nRch_desire), stat=ierr)
 if(ierr/=0) message=trim(message)//'unable to allocate space for the vectors of desired reaches'

 ! get the desired stream segments
 ixSeg_desired(:) = pack(ixSeg_map, isDesired)

 ! get the desired HRUs
 do iRch=1,nRch_desire
  ixDesire = ixSeg_desired(iRch)
  if(structNTOPO(ixDesire)%var(ixNTOPO%nHRU)%dat(1)>0)then
   ixHRU_map( structNTOPO(ixDesire)%var(ixNTOPO%hruContribIx)%dat(:) ) = &
              structNTOPO(ixDesire)%var(ixNTOPO%hruContribIx)%dat(:)
  endif
 end do  ! looping through desired stream segments

 ! allocate space for the number of HRUs
 nHRU_desire = count(ixHRU_map/=integerMissing)
 allocate(ixHRU_desired(nHRU_desire), stat=ierr)
 if(ierr/=0) message=trim(message)//'unable to allocate space for the vectors of desired HRUs'

 ! pack the desired indices into the HRU vector
 ixHRU_desired(:) = pack(ixHRU_map, ixHRU_map/=integerMissing)

 ! --------- update the dimension lengths ---------------------------------------------------------------------------------------------


 tot_hru       = 0_i4b  ! total number of all the upstream hrus for all stream segments
 tot_upseg     = 0_i4b  ! sum of immediate upstream segments
 tot_upstream  = 0_i4b  ! total number of upstream reaches for all reaches
 ! get the updated dimensions
 do iRch=1,nRch_desire
  tot_hru      = tot_hru      + structNTOPO( ixSeg_desired(iRch) )%var(ixNTOPO%nHRU)%dat(1)
  tot_upseg    = tot_upseg    + size(structNTOPO( ixSeg_desired(iRch) )%var(ixNTOPO%upSegIds)%dat)
  tot_upstream = tot_upstream + size(structNTOPO( ixSeg_desired(iRch) )%var(ixNTOPO%allUpSegIndices)%dat)
 end do

 ! check
 !  print*, 'ixHRU_desired = ', ixHRU_desired
 !  print*, 'ixSeg_desired = ', ixSeg_desired
 !  do iRch=1,nRch_desire
 !   print*, isDesired(ixSeg_desired(iRch)), ixSeg_desired(iRch), ixSeg_map(ixSeg_desired(iRch)), &
 !    structNTOPO( ixSeg_desired(iRch) )%var(ixNTOPO%nHRU)%dat(:), structNTOPO( ixSeg_desired(iRch) )%var(ixNTOPO%hruContribIx)%dat(:)
 !  end do
 !  print*, trim(message)//'PAUSE: '; read(*,*)

 ! ----------------------------------------------------------------------------------------
 ! ----------------------------------------------------------------------------------------
 CONTAINS

  ! for a given upstream reach, add a downstream reach index to the list
  subroutine add2list(ixDown,ierr,message)
  integer(i4b) ,intent(in)         :: ixDown   ! downstream reach index
  integer(i4b), intent(out)        :: ierr     ! error code
  character(*), intent(out)        :: message  ! error message
  message='ADD2LIST/'

  ! increment number of downstream reaches
  nDown = nDown+1

  ! allocate node to be added to the list
  allocate(node, stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for linked list'; return; endif

  ! add downstream reach index
  node%ix = ixDown

  ! insert downstream reach at the beginning of existing list
  node%next => HPOINT
  HPOINT    => node

  ! nullify pointer
  NULLIFY(node)

  END SUBROUTINE ADD2LIST

 ! ----------------------------------------------------------------------------------------
 END SUBROUTINE REACH_MASK_ORIG

end module network_topo
