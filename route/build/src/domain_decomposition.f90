MODULE domain_decomposition_module

USE nrtype
USE public_var
USE pfafstetter_module
USE dataTypes,         ONLY: var_clength   ! character type:   var(:)%dat
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
USE dataTypes,         ONLY: basin         !
USE dataTypes,         ONLY: reach         !
USE var_lookup,        ONLY: ixPFAF        ! index of variables for the pfafstetter code
USE var_lookup,        ONLY: ixNTOPO       ! index of variables for the netowork topolgy
USE nr_utility_module, ONLY: indexx        ! Num. Recipies utilities
USE nr_utility_module, ONLY: indexTrue     ! Num. Recipies utilities
USE nr_utility_module, ONLY: arth          ! Num. Recipies utilities

implicit none

private

public :: classify_river_basin
public :: classify_river_basin_mpi

contains

 subroutine classify_river_basin_mpi(nSeg, structPFAF, structNTOPO, maxSegs, ierr, message)
   USE globalData, only : ixDomain     ! count of domains
   ! Domain decomposition wrapper
   implicit none
   ! Input variables
   integer(i4b),                   intent(in)  :: nSeg                   ! number of stream segments
   type(var_clength), allocatable, intent(in)  :: structPFAF(:)          ! pfafstetter code
   type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)         ! network topology
   integer(i4b),                   intent(in)  :: maxSegs                ! threshold number of tributary reaches
   ! Output variables
   integer(i4b),                   intent(out) :: ierr
   character(len=strLen),          intent(out) :: message                ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage               ! error message from subroutine
   character(len=32)                           :: pfafs(nSeg)            ! pfaf_codes for all the segment
   character(len=32)                           :: pfafOutlet             ! pfaf_codes for a outlet segment
   character(len=32)                           :: pfafCommon             ! common pfaf_codes over the entire basin
   character(len=32), allocatable              :: pfafOutlets(:)         ! pfaf_codes for outlet segments
   integer(i4b)                                :: level                  ! number of digits of common pfaf codes given pfaf code at outlet reach
   integer(i4b)                                :: downIndex(nSeg)        ! downstream segment index for all the segments
   integer(i4b)                                :: segIndex(nSeg)         ! reach index for all the segments
   integer(i4b)                                :: nOuts                  ! number of outlets
   integer(i4b)                                :: iSeg, iOut             ! loop indices

   ierr=0; message='classify_river_basin_mpi/'

   ixDomain = 0

   forall(iSeg=1:nSeg) pfafs(iSeg)     = structPFAF(iSeg)%var(ixPFAF%code)%dat(1)
   forall(iSeg=1:nSeg) downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)
   forall(iSeg=1:nSeg) segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)

   ! Number of outlets
   nOuts=count(downIndex<0)

   allocate(pfafOutlets(nOuts), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [pfasOutlets]'; return; endif

   ! Outlet information - pfaf code
   pfafOutlets = pack(pfafs, downIndex<0)

   do iOut = 1,nOuts
     pfafOutlet = adjustl(pfafOutlets(iOut))

     ! Identify pfaf level for a river basin given outlet pfaf
     call get_common_pfaf(pfafs, pfafOutlet, level, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! perform domain decomposition for a river basin
     pfafCommon = trim(pfafOutlet(1:level))   ! to be removed
     call decomposition(pfafs, segIndex, downIndex, pfafCommon, maxSegs, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end do

 end subroutine classify_river_basin_mpi


 recursive subroutine decomposition(pfafs, segIndex, downIndex, pfafCode, maxSegs, ierr, message)
   implicit none
   ! Input variables
   character(len=32),        intent(in)  :: pfafs(:)        ! pfaf_code list
   integer(i4b),             intent(in)  :: downIndex(:)    ! downstream segment index for all the segments
   integer(i4b),             intent(in)  :: segIndex(:)     ! reach index for all the segments
   character(len=32),        intent(in)  :: pfafCode        ! common pfaf codes within a river basin
   integer(i4b)                          :: maxSegs         ! maximum number of reaches in basin
   ! Output variables
   integer(i4b),             intent(out) :: ierr
   character(len=strLen),    intent(out) :: message         ! error message
   ! Local variables
   character(len=strLen)                 :: cmessage        ! error message from subroutine
   integer(i4b)                          :: iPfaf,iSeg      ! loop indices
   integer(i4b)                          :: nMatch          ! count for pfafs matching with pfafCode
   integer(i4b)                          :: nSeg            ! number of stream segments
   logical(lgt),     allocatable         :: isMatch(:)      ! logical to indicate xxxxx
   character(len=32)                     :: pfafCodeTmp     ! copy of pfafCode
   character(len=32),allocatable         :: subPfafs(:)     ! subset of pfaf_codes
   integer(i4b),     allocatable         :: subSegIndex(:)  ! subset of segment indices
   integer(i4b),     allocatable         :: subDownIndex(:) ! subset of segment indices
   character(len=32)                     :: subPfaf         ! pfaf_code appended by next level
   character(len=32)                     :: pfaf            ! a pfaf_code
   character(len=1)                      :: cPfaf           ! character digit

   ierr=0; message='decomposition/'

   pfafCodeTmp = trim(pfafCode)
   nSeg = size(pfafs)

   do iPfaf = 1,9
     write(cPfaf, '(I1)') iPfaf
     pfafCodeTmp = trim(pfafCodeTmp)//cPfaf

     if (.not. allocated(isMatch)) then
       allocate(isMatch(nSeg), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating [isMatch]'; return; endif
     endif
     isMatch(:) = .false.

     do iSeg = 1,size(pfafs)
       pfaf = adjustl(pfafs(iSeg))
       subPfaf = pfaf(1:len(trim(pfafCodeTmp)))
       if (subPfaf == pfafCodeTmp) then
         isMatch(iSeg) = .true.
       end if
     end do
     nMatch = count(isMatch)

     print*,'pfafCode, number of matches = ', trim(pfafCodeTmp), nMatch

     if (nMatch==0) then
       pfafCodeTmp = pfafCodeTmp(1:len(trim(pfafCodeTmp))-1)  ! decrement
       cycle
     end if

     if (nMatch < maxSegs) then
       print*,'Aggregate: nMatch less than nThreh = ', maxSegs
       allocate(subPfafs(nMatch), subDownIndex(nMatch), subSegIndex(nMatch), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating [subPfafs, subDownIndex, subSegIndex]'; return; endif
       subPfafs     = pack(pfafs, isMatch)
       subSegIndex  = pack(segIndex, isMatch)
       subDownIndex = pack(downIndex, isMatch)
       call aggregate(pfafCodeTmp, subPfafs, subSegIndex, subDownIndex, ierr, cmessage) ! populate reach classification data structure
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
       deallocate(subPfafs, subDownIndex, subSegIndex, stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem deallocating [subPfafs, subDownIndex, subSegIndex]'; return; endif
     else
       print*,'Disaggregate: nMatch more than maxSegs = ', maxSegs
       call decomposition(pfafs, segIndex, downIndex, pfafCodeTmp, maxSegs, ierr, cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end if

     pfafCodeTmp = pfafCodeTmp(1:len(trim(pfafCodeTmp))-1) ! decrement

   end do

 end subroutine decomposition


 subroutine aggregate(pfafCode, pfafs, segIndex, downIndex, ierr, message)
   USE globalData, only : domains      ! domain data structure
   USE globalData, only : ixDomain     ! count of domains
   implicit none
   ! Input variables
   character(len=32),          intent(in)  :: pfafCode                  ! Common pfaf codes within a river basin
   character(len=32),          intent(in)  :: pfafs(:)                  ! pfaf code for all the reaches
   integer(i4b),               intent(in)  :: downIndex(:)              ! downstream reach index for all the reaches
   integer(i4b),               intent(in)  :: segIndex(:)               ! reach index for all the reaches
   ! Output variables
   integer(i4b),               intent(out) :: ierr                      ! error code
   character(len=strLen),      intent(out) :: message                   ! error message
   ! Local variables
   character(len=strLen)                   :: cmessage                  ! error message from subroutine
   character(len=32)                       :: mainCode                  ! mainstem code
   character(len=32)                       :: pfaf                      ! a pfaf code for a reach
   integer(i4b)                            :: iSeg, jSeg, iTrib         ! loop indices
   integer(i4b)                            :: nSeg, nTrib, nMainstem    ! number of reaches, tributaries, and mainstems, respectively
   integer(i4b)                            :: cc                        ! counter
   integer(i4b)                            :: pfafLen
   integer(i4b)                            :: mainCodeLen
   integer(i4b)                            :: pfaf_old, pfaf_new
   integer(i4b),allocatable                :: ixSubset(:)               ! subset indices based on logical array from global index array
   integer(i4b),allocatable                :: downIndexNew(:)           !
   character(len=32)                       :: subPfaf
   character(len=32),allocatable           :: trib_outlet_pfafs(:)
   logical(lgt)                            :: isInterbasin, isHeadwater
   logical(lgt),allocatable                :: isClosed(:)
   logical(lgt),allocatable                :: isTribOutlet(:)
   logical(lgt),allocatable                :: isMainstem(:)
   logical(lgt),allocatable                :: isMainstem2d(:,:)

   ierr=0; message='aggregate/'

   ! Initialization
   ixDomain = ixDomain + 1
   nSeg = size(pfafs)
   allocate(domains(ixDomain)%pfaf(nSeg), domains(ixDomain)%nRch(nSeg), domains(ixDomain)%segIndex(nSeg), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [domains%pfaf, domains%nSeg]'; return; endif
   domains(ixDomain)%pfaf(:)  = '-999'
   domains(ixDomain)%nRch(:)  = -999
   domains(ixDomain)%segIndex = segIndex

   ! get pfaf old and pfaf new
   pfafLen = len(trim(pfafCode))
   read(pfafCode(pfafLen-1:pfafLen-1),'(I1)') pfaf_old
   read(pfafCode(pfafLen:pfafLen),'(I1)') pfaf_new

   ! check if an interbasin
   isInterbasin = (mod(pfaf_new, 2)==1)
   isHeadwater  = (mod(pfaf_old, 2)==0 .and. pfaf_new==9)

   print*, 'basin, isInterbasin, isHeadwater = ',pfafCode, isInterbasin, isHeadWater

   ! get closed basin segments
   allocate(isClosed(nSeg), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [isClosed]'; return; endif
   call find_closed_basin(pfafCode, pfafs, isClosed, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (isInterbasin .and. (.not. isHeadwater)) then  ! if a river reach is in inter-basin and not headwater
     ! populate mainstem segments
     call find_mainstems(pfafCode, pfafs, isMainstem, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     nMainstem = count(isMainstem)
     do iSeg = 1, nSeg
       if (isMainstem(iSeg)) then
         domains(ixDomain)%pfaf(iSeg) = '-'//trim(pfafCode)
         domains(ixDomain)%nRch(iSeg) = nMainstem
       end if
     end do

     ! Populate tributary segments
     ! Get logical array indicating tributary outlet segments (size = nSeg)
     isMainstem2d =  spread(isMainstem,2,1)
     allocate(downIndexNew(nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [downIndexNew]'; return; endif
     downIndexNew=-999
     do iSeg = 1, nSeg
       do jSeg = 1, nSeg
         if (downIndex(iSeg) == segIndex(jSeg)) then
           downIndexNew(iSeg) = jSeg
           exit
         end if
       end do
     end do
     call lgc_tributary_outlet(isMainstem2d, downIndexNew, isTribOutlet, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     ! Get indices of tributary outlet segments (size = number of tributaries)
     nTrib = count(isTribOutlet)

     allocate(ixSubset(nTrib), trib_outlet_pfafs(nTrib), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [ixSubset, trib_outlet_pfafs]'; return; endif
     ixSubset = pack(arth(1,1,nSeg),isTribOutlet)
     trib_outlet_pfafs = pfafs(ixSubset)
     deallocate(ixSubset, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem deallocating [ixSubset]'; return; endif

     do iTrib = 1,nTrib             ! loop through each tributary
       mainCode = mainstem_code(trim(trib_outlet_pfafs(iTrib)))
       mainCodeLen = len(trim(mainCode))
       cc = 0
       do iSeg = 1, nSeg
         pfaf = adjustl(pfafs(iSeg))
         subPfaf = trim(pfaf(1:mainCodeLen))
         if (subPfaf == mainCode) then
           cc = cc + 1
           domains(ixDomain)%pfaf(iSeg) = mainCode
         end if
       end do
       !domains(ixDomain)%nRch(iSegList) = cc
     end do
   else    ! a river reach is in tributary basin or headwater
     domains(ixDomain)%pfaf(:) = trim(pfafCode)
     domains(ixDomain)%nRch(:) = nSeg
   end if

   allocate(ixSubset(count(isClosed)), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [ixSubset]'; return; endif
   ixSubset = pack(arth(1,1,nSeg),isClosed)
   domains(ixDomain)%pfaf(ixSubset) = '-777'
   domains(ixDomain)%nRch(ixSubset) = -777
   deallocate(ixSubset, stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem deallocating [ixSubset]'; return; endif

 end subroutine aggregate


 subroutine classify_river_basin(nSeg, structPFAF, structNTOPO, river_basin, maxSegs, ierr, message)
  ! Identify tributary basin and mainstems using river network data and pfafstetter code
  ! Output: return populated basin dataType
  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: nSeg                   ! number of stream segments
  type(var_clength), allocatable, intent(in)  :: structPFAF(:)          ! pfafstetter code
  type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)         ! network topology
  integer(i4b)                  , intent(in)  :: maxSegs                ! maximum reach number in each tributary
  ! Output variables
  type(basin),       allocatable, intent(out) :: river_basin(:)         ! river basin data
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                ! error message
  ! Local variables
  type(reach), allocatable                    :: tmpMainstem(:)         ! temporary reach structure for mainstem
  character(len=strLen)                       :: cmessage               ! error message from subroutine
  character(len=32)                           :: pfafs(nSeg)            ! pfaf_codes for all the segment
  character(len=32), allocatable              :: pfaf_out(:)            ! pfaf_codes for outlet segments
  character(len=32)                           :: upPfaf                 ! pfaf_code for one upstream segment
  character(len=32)                           :: outMainstemCode
  logical(lgt),      allocatable              :: mainstems(:,:)         ! logical to indicate segment is mainstem at each level
  logical(lgt),      allocatable              :: updated_mainstems(:,:) ! logical to indicate segment is mainstem at each level
  logical(lgt),      allocatable              :: lgc_trib_outlet(:)     ! logical to indicate segment is outlet of tributary
  logical(lgt)                                :: done                   ! logical
  integer(i4b)                                :: maxLevel=20
  integer(i4b)                                :: downIndex(nSeg)        ! downstream segment index for all the segments
  integer(i4b)                                :: segIndex(nSeg)         ! reach index for all the segments
  integer(i4b)                                :: segOrder(nSeg)         ! reach order for all the segments
  integer(i4b)                                :: rankSegOrder(nSeg)     ! ranked reach order for all the segments
  integer(i4b), allocatable                   :: segIndex_out(:)        ! index for outlet segment
  integer(i4b)                                :: level                  ! manstem level
  integer(i4b)                                :: nMains,nUpSegMain      ! number of mainstems in a level, and number of segments in a mainstem
  integer(i4b)                                :: nOuts                  ! number of outlets
  integer(i4b)                                :: nUpSegs                ! number of upstream segments for a specified segment
  integer(i4b)                                :: nTrib                  ! Number of Tributary basins
  integer(i4b), allocatable                   :: nSegTrib(:)            ! Number of segments in each tributary basin
  integer(i4b), allocatable                   :: nSegMain(:)            ! number of mainstem only segments excluding tributary segments
  integer(i4b), allocatable                   :: rankTrib(:)            ! index for ranked tributary based on num. of upstream reaches
  integer(i4b), allocatable                   :: rankMain(:)            ! index for ranked mainstem based on num. of mainstem reaches
  integer(i4b), allocatable                   :: segOrderTrib(:)        ! index for tributary upstream reach order
  integer(i4b), allocatable                   :: segMain(:)             ! index for segment in one mainstem
  integer(i4b), allocatable                   :: segOrderMain(:)        ! index for ordered segment in one mainstem
  integer(i4b), allocatable                   :: trPos(:)               ! tributary outlet indices
  integer(i4b), allocatable                   :: msPos(:)               ! mainstem indices
  integer(i4b)                                :: iSeg,jSeg,iOut         ! loop indices
  integer(i4b)                                :: iTrib,jTrib            ! loop indices
  integer(i4b)                                :: iMain,jMain            ! loop indices

  ierr=0; message='classify_river_basin/'

  forall(iSeg=1:nSeg) pfafs(iSeg)     = structPFAF(iSeg)%var(ixPFAF%code)%dat(1)
  forall(iSeg=1:nSeg) downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)
  forall(iSeg=1:nSeg) segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
  forall(iSeg=1:nSeg) segOrder(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)

  call indexx(segOrder,rankSegOrder)

  ! Number of outlets
  nOuts=count(downIndex<0)

  allocate(river_basin(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating river_basin'; return; endif
  allocate(pfaf_out(nOuts), segIndex_out(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating [pfaf_out, segIndex_out]'; return; endif

  ! Outlet information - pfaf code and segment index
  pfaf_out = pack(pfafs, downIndex<0)
  segIndex_out = pack(segIndex, downIndex<0)

  do iOut = 1,nOuts

    print*, 'working on outlet:'//trim(adjustl(pfaf_out(iOut)))

    allocate(river_basin(iOut)%level(maxLevel), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%mainstem'; return; endif

    river_basin(iOut)%outIndex = segIndex_out(iOut)

!    call system_clock(startTime)
    ! Identify pfaf level given a river network
    call get_common_pfaf(pfafs, pfaf_out(iOut), level, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    call system_clock(endTime)
!    elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!    write(*,"(A,1PG15.7,A)") '       elapsed-time [get_common_pfaf] = ', elapsedTime, ' s'

!    call system_clock(startTime)
    ! Identify mainstem segments at all levels (up to maxLevel)
    call lgc_mainstems(pfafs, pfaf_out(iOut), maxLevel, mainstems, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    call system_clock(endTime)
!    elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!    write(*,"(A,1PG15.7,A)") '      elapsed-time [lgc_mainstems] = ', elapsedTime, ' s'

    ! Initial assignment of mainstem segments
    allocate(updated_mainstems(size(mainstems,1),size(mainstems,2)), stat=ierr)
    updated_mainstems = .false.
    updated_mainstems(:,level) = mainstems(:,level)

    ! Identify the lowest level mainstem segment index
    call indexTrue(mainstems(:,level), msPos)

    allocate(river_basin(iOut)%level(level)%mainstem(1), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%mainstem(level)%mainstem'; return; endif
    allocate(river_basin(iOut)%level(level)%mainstem(1)%segIndex(size(msPos)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(:)%mainstem(:)%segIndex'; return; endif
    allocate(segOrderMain(size(msPos)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating segOrderMain'; return; endif

    ! Compute reach order for mainstem segment
    call indexx(rankSegOrder(msPos), segOrderMain)

    river_basin(iOut)%level(level)%mainstem(1)%segIndex(1:size(msPos))  = msPos(segOrderMain)
    river_basin(iOut)%level(level)%mainstem(1)%nRch = size(msPos)

!    call system_clock(startTime)
    ! Identify tributary outlets into a mainstem at the lowest level
    !  i.e. the segment that is not on any mainstems AND flows into any mainstem segments
    call lgc_tributary_outlet(mainstems(:,:level), downIndex, lgc_trib_outlet, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    call system_clock(endTime)
!    elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!    write(*,"(A,1PG15.7,A)") '      elapsed-time [lgc_tributary_outlet] = ', elapsedTime, ' s'

    deallocate(mainstems, segIndex_out, segOrderMain, stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem deallocating mainstems'; return; endif

    do
      level = level + 1
      print*, 'Mainstem Level = ', level

      ! number of tributary basins
      nTrib = count(lgc_trib_outlet)

      allocate(nSegTrib(nTrib),rankTrib(nTrib), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating nSegTrib'; return; endif

      ! Extract array elements with only tributary outlet (keep indices in master array
      call indexTrue(lgc_trib_outlet, trPos)

      ! number of tributary segments
      do iTrib=1,nTrib
        nSegTrib(iTrib) = size(structNTOPO(trPos(iTrib))%var(ixNTOPO%allUpSegIndices)%dat)
      end do

      ! rank the tributary outlets based on number of upstream segments
      call indexx(nSegTrib, rankTrib)

      ! Identify mainstems of large tributaries (# of segments > maxSeg) and update updated_mainstems.
      done=.true.
      nMains = 0
      do iTrib=nTrib,1,-1
        if (nSegTrib(rankTrib(iTrib)) > maxSegs) then
           write(*,'(A,A,I5)') 'Exceed maximum number of segments: ', pfafs(trPos(rankTrib(iTrib))), nSegTrib(rankTrib(iTrib))
           nMains=nMains+1
           done=.false.
        else
           exit
        endif
      end do

      if (done) then ! if no mainstem/tributary updated, update tributary reach info, and then exist loop

!        print*, 'Main/Trib update finished, Now Update tributary reach ....'

        allocate(river_basin(iOut)%tributary(nTrib), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating river_basin%tributary'; return; endif

        do iTrib=nTrib,1,-1

!          call system_clock(startTime)

          jTrib = nTrib - iTrib + 1

!          if (mod(jTrib,1000)==0) print*, 'jTrib = ',jTrib

          allocate(river_basin(iOut)%tributary(jTrib)%segIndex(nSegTrib(rankTrib(iTrib))), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(:)%tributary%segIndex'; return; endif
          allocate(segOrderTrib(nSegTrib(rankTrib(iTrib))), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating segOrderTrib'; return; endif

          ! compute reach order for tributary segments
          associate( segIndexTrib => structNTOPO(trPos(rankTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat)
          call indexx(rankSegOrder(segIndexTrib), segOrderTrib)
          river_basin(iOut)%tributary(jTrib)%segIndex(:) = segIndexTrib(segOrderTrib)
          end associate

          ! compute number of segments in each tributary
          river_basin(iOut)%tributary(jTrib)%nRch = nSegTrib(rankTrib(iTrib))

          deallocate(segOrderTrib, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem deallocating segOrderTrib'; return; endif

!          if (mod(jTrib,1000)==0) then
!          call system_clock(endTime)
!          elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!          write(*,"(A,1PG15.7,A)") ' total elapsed-time [Update Trib. reach] = ', elapsedTime, ' s'
!          endif

        end do

        exit

      else ! if mainstem/tributary updated, store mainstem reach info at current level, update lgc_trib_outlet and go onto next level

        allocate(river_basin(iOut)%level(level)%mainstem(nMains), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%level(level)%mainstem'; return; endif
        allocate(tmpMainstem(nMains), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating tmpMainstem'; return; endif
        allocate(nSegMain(nMains), rankMain(nMains), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating [nSegMain,rankMain]'; return; endif

        do iTrib=nTrib,nTrib-nMains+1,-1

!          if (mod(nTrib-iTrib+1,10)==0) print*, 'iTrib = ', nTrib-iTrib+1
!          call system_clock(startTime)

          ! Outlet mainstem code
          outMainstemCode = mainstem_code(pfafs(trPos(rankTrib(iTrib))))

          associate(upSegIndex => structNTOPO(trPos(rankTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat)

          nUpSegs = size(upSegIndex)

          allocate(segMain(nUpSegs), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating [segMain]'; return; endif

          nUpSegMain = 0
          do jSeg=1,nUpSegs
            upPfaf = pfafs(upSegIndex(jSeg))
            if (trim(mainstem_code(upPfaf)) /= trim(outMainstemCode)) cycle
            updated_mainstems(upSegIndex(jSeg),level) = .true.
            nUpSegMain = nUpSegMain+1
            segMain(nUpSegMain) = upSegIndex(jSeg)
          end do

          allocate(tmpMainstem(nTrib-iTrib+1)%segIndex(nUpSegMain), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating tmpMainstem(nUpSegMain)%segIndex'; return; endif
          allocate(segOrderMain(nUpSegMain), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating segOrderMain'; return; endif

          call indexx(rankSegOrder(segMain(1:nUpSegMain)), segOrderMain)
          tmpMainstem(nTrib-iTrib+1)%segIndex(:) = segMain(segOrderMain)
          nSegMain(nTrib-iTrib+1) = nUpSegMain

          end associate

          deallocate(segMain, segOrderMain, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem deallocating [segMain, segOrderMain]'; return; endif

!          if (mod(nTrib-iTrib+1,10)==0) then
!          call system_clock(endTime)
!          elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!          write(*,"(A,1PG15.7,A)") ' total elapsed-time [update trib each loop] = ', elapsedTime, ' s'
!          endif

        end do ! tributary loop

        ! rank the mainstem based on number of upstream segments (excluding tributaries)
        call indexx(nSegMain, rankMain)
        ! populate mainstem segment component in river basin structure
        do iMain=1,nMains
          jMain=rankMain(nMains-iMain+1)
          allocate(river_basin(iOut)%level(level)%mainstem(iMain)%segIndex(nSegMain(jMain)), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%level(level)%mainstem(iMain)%segIndex'; return; endif
          river_basin(iOut)%level(level)%mainstem(iMain)%segIndex(:) = tmpMainstem(jMain)%segIndex(:)
          river_basin(iOut)%level(level)%mainstem(iMain)%nRch = nSegMain(jMain)
        enddo

        ! update lgc_trib_outlet based on added mainstem
        call lgc_tributary_outlet(updated_mainstems(:,:level), downIndex, lgc_trib_outlet, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        deallocate(tmpMainstem, nSegMain, rankMain, stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem deallocating [tmpMainstem,nSegMain,rankMain]'; return; endif
        deallocate(nSegTrib, rankTrib, stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem deallocating [nSegTrib, rankTrib]'; return; endif

      endif
    end do ! end of tributary update loop

  end do ! outlet loop

 end subroutine classify_river_basin


end module domain_decomposition_module
