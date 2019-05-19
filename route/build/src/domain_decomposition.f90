MODULE domain_decomposition

! External modules (general modules)
! numeric types
USE nrtype
! derived data types
USE dataTypes,         ONLY: var_clength   ! character type:   var(:)%dat
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
! variable indices
USE var_lookup,        ONLY: ixPFAF        ! index of variables for the pfafstetter code
USE var_lookup,        ONLY: ixNTOPO       ! index of variables for the netowork topolgy
! General utilities
USE nr_utility_module, ONLY: indexx        ! sorting array
USE nr_utility_module, ONLY: indexTrue     ! index at only true in array
USE nr_utility_module, ONLY: arth          ! generate sequential array
! updated and saved data
USE public_var

implicit none

private

public :: classify_river_basin_omp

contains

 ! ***************************************************************
 ! public subroutine: Main routine for OMP domain decomposition
 ! ***************************************************************
 subroutine classify_river_basin_omp(nSeg,          & ! input: number of reaches in the entire river network
                                     structPFAF,    & ! input: pfafstetter code data structure
                                     structNTOPO,   & ! input: river network data structure
                                     river_basin,   & ! output: tributary/mainstem domain decomposition
                                     ierr, message)
  ! Details:
  ! Identify tributary reaches and mainstems using river network topology and pfafstetter code
  ! for each basin. Basins are identified based on outlet reach in topology (downSegIndex < 0)

  ! The following data struct components need to be populated
  !   structPFAF(:)%var(ixPFAF%code)%dat(1)
  !   structNTOPO(:)%var(ixNTOPO%segIndex)%dat(1)
  !   structNTOPO(:)%var(ixNTOPO%downSegIndex)%dat(1)
  !   structNTOPO(:)%var(ixNTOPO%rchOrder)%dat(1)
  !   structNTOPO(:)%var(ixNTOPO%allUpSegIndices)%dat(:)
  !
  ! Populate the domain data structure
  !   river_basin(ixOut)%outIndex                                   : oulet reach index for "ixOut" basin
  !   river_basin(ixOut)%level(ixLvl)%mainstem(ixMain)%segIndex(:)  : reach indice for "ixLvl" level mainstem "ixMain" of "ixOut" basin
  !                                                   %nRch         : number of reaches for "ixLvl" level mainstem "ixMain" of "ixOut" basin
  !   river_basin(ixOut)%tributary(ixTrib)%segIndex(:)              : reach indice for "ixTrib" tributary of "ixOut" basin
  !                                       %nRch                     : number of reaches for "ixTrib" tributary of "ixOut" basin

  ! External modules
  ! derive data types
  USE dataTypes,          only: basin                ! basin data structure
  USE dataTypes,          only: reach                !
  USE public_var,         only: maxSegs
  USE public_var,         only: maxLevel
  USE public_var,         only: pfafMissing
  ! pfafstetter routines
  USE pfafstetter_module, only: get_common_pfaf
  USE pfafstetter_module, only: lgc_mainstems
  USE pfafstetter_module, only: lgc_tributary_outlet
  USE pfafstetter_module, only: mainstem_code
  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: nSeg                   ! number of stream segments
  type(var_clength), allocatable, intent(in)  :: structPFAF(:)          ! pfafstetter code
  type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)         ! network topology
  ! Output variables
  type(basin),       allocatable, intent(out) :: river_basin(:)         ! river basin data
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                ! error message
  ! Local variables
  type(reach), allocatable                    :: tmpMainstem(:)         ! temporary reach structure for mainstem
  character(len=strLen)                       :: cmessage               ! error message from subroutine
  character(len=32)                           :: pfafs(nSeg)            ! pfaf_codes for all the segment
  character(len=32), allocatable              :: pfafOutlets(:)         ! pfaf_codes for outlet segments
  character(len=32)                           :: upPfaf                 ! pfaf_code for one upstream segment
  character(len=32)                           :: outMainstemCode
  logical(lgt),      allocatable              :: mainstems(:,:)         ! logical to indicate segment is mainstem at each level
  logical(lgt),      allocatable              :: updated_mainstems(:,:) ! logical to indicate segment is mainstem at each level
  logical(lgt),      allocatable              :: lgc_trib_outlet(:)     ! logical to indicate segment is outlet of tributary
  logical(lgt),      allocatable              :: isValid(:)             ! ogical to indicate reach with vlid pfaf code (non "0")
  logical(lgt)                                :: done                   ! logical
  integer(i4b)                                :: downIndex(nSeg)        ! downstream segment index for all the segments
  integer(i4b)                                :: segIndex(nSeg)         ! reach index for all the segments
  integer(i4b)                                :: segOrder(nSeg)         ! reach order for all the segments
  integer(i4b)                                :: rankSegOrder(nSeg)     ! ranked reach order for all the segments
  integer(i4b), allocatable                   :: ixUpSeg(:)             ! list of upstream reach indices at a given outlet
  integer(i4b), allocatable                   :: ixSubset(:)            ! subset indices based on logical array from global index array
  integer(i4b), allocatable                   :: ixOutlets(:)           ! index for outlet segment
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

!  integer(i4b)                                :: i1,i2,jLevel,level1,dangle
!  integer(i4b)                                :: nUps
!  integer(i4b), allocatable                   :: idxList(:)            ! mainstem indices
!  integer(i4b)                                :: nSeg_tmp
!  integer(i4b)                                :: MainOutIndex
!  character(len=32)                           :: levelArray(nSeg)        !
!  integer(i4b)                                :: tribNseg(nSeg)         !
!  integer(i4b)                                :: tribIndex(nSeg)         !

  integer*8                                   :: cr, startTime, endTime ! rate, start and end time stamps
  real(dp)                                    :: elapsedTime            ! elapsed time for the process

  CALL system_clock(count_rate=cr)

  ierr=0; message='classify_river_basin_omp/'

  forall(iSeg=1:nSeg) pfafs(iSeg)     = structPFAF(iSeg)%var(ixPFAF%code)%dat(1)
  forall(iSeg=1:nSeg) downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)
  forall(iSeg=1:nSeg) segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
  forall(iSeg=1:nSeg) segOrder(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)

  ! sorting reach processing order
  call indexx(segOrder,rankSegOrder)

  ! Number of outlets
  nOuts=count(downIndex<0)

  allocate(river_basin(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating river_basin'; return; endif
  allocate(pfafOutlets(nOuts), ixOutlets(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating [pfafOutlets, ixOutlets]'; return; endif

  ! Outlet information - pfaf code and segment index
  pfafOutlets = pack(pfafs, downIndex<0)
  ixOutlets   = pack(segIndex, downIndex<0)

  ! Process basin by basin
  do iOut = 1,nOuts

    print*, 'working on outlet:'//trim(adjustl(pfafOutlets(iOut)))

    allocate(river_basin(iOut)%level(maxLevel), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%mainstem'; return; endif

    ! basin outlet reach index
    river_basin(iOut)%outIndex = ixOutlets(iOut)

    ! get all the upstream reach indices for this outlet
    ! check if upstream segment include reaches with pCode==0. if so remove them.
    associate (ixUpSeg_tmp => structNTOPO(ixOutlets(iOut))%var(ixNTOPO%allUpSegIndices)%dat)
    if (allocated(isValid)) deallocate(isValid)
    allocate(isValid(size(ixUpSeg_tmp)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating [isValid]'; return; endif
    isValid =.false.
    do iSeg = 1,size(ixUpSeg_tmp)
      if (trim(adjustl(pfafs(ixUpSeg_tmp(iSeg))))/=pfafMissing) isValid(iSeg)=.true.
    end do
    if (all(.not.isValid)) cycle
    if (allocated(ixUpSeg)) deallocate(ixUpSeg)
    allocate(ixUpSeg(count(isValid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating [ixUpSeg]'; return; endif
    call indexTrue(isValid, ixSubset)
    ixUpSeg = ixUpSeg_tmp(ixSubset)
    deallocate(ixSubset, stat=ierr)

    ! Special case: small basin (number of reaches < maxSegs
    ! Put all the small basin reaches under tributary data structures
    if (size(ixUpSeg)<maxSegs) then

      allocate(river_basin(iOut)%tributary(1), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%tributary'; return; endif
      allocate(river_basin(iOut)%tributary(1)%segIndex(size(ixUpSeg)), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%tributary(1)%segIndex'; return; endif
      allocate(segOrderTrib(size(ixUpSeg)), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating segOrderTrib'; return; endif

      ! Compute reach order for only small basin
      call indexx(rankSegOrder(ixUpSeg), segOrderTrib)

      river_basin(iOut)%tributary(1)%segIndex(:) = ixUpSeg(segOrderTrib)
      river_basin(iOut)%tributary(1)%nRch        = size(ixUpSeg)

      deallocate(segOrderTrib, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem deallocating segOrderTrib'; return; endif

      cycle

    end if
    end associate

    call system_clock(startTime)
    ! Identify pfaf level given a river network
    call get_common_pfaf(pfafs(ixUpSeg), pfafOutlets(iOut), level, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '       elapsed-time [get_common_pfaf] = ', elapsedTime, ' s'

!    call system_clock(startTime)
    ! Identify mainstem segments at all levels (up to maxLevel)
    call lgc_mainstems(pfafs(ixUpSeg), pfafOutlets(iOut), maxLevel, mainstems, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    call system_clock(endTime)
!    elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!    write(*,"(A,1PG15.7,A)") '      elapsed-time [lgc_mainstems] = ', elapsedTime, ' s'

    ! Initial assignment of mainstem segments
    allocate(updated_mainstems(nSeg,size(mainstems,2)), stat=ierr)
    updated_mainstems = .false.
    updated_mainstems(ixUpSeg,level) = mainstems(:,level)

    ! Identify the lowest level mainstem segment index
    call indexTrue(updated_mainstems(:,level), msPos)

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
    call lgc_tributary_outlet(updated_mainstems(:,:level), downIndex, lgc_trib_outlet, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    call system_clock(endTime)
!    elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!    write(*,"(A,1PG15.7,A)") '      elapsed-time [lgc_tributary_outlet] = ', elapsedTime, ' s'

    deallocate(mainstems, segOrderMain, stat=ierr)
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
          if(ierr/=0)then; message=trim(message)//'problem deallocating [segOrderTrib]'; return; endif

!          if (mod(jTrib,1000)==0) then
!          call system_clock(endTime)
!          elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!          write(*,"(A,1PG15.7,A)") ' total elapsed-time [Update Trib. reach] = ', elapsedTime, ' s'
!          endif

        end do

        deallocate(nSegTrib, rankTrib, stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem deallocating [nSegTrib, rankTrib]'; return; endif

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

 end subroutine classify_river_basin_omp

! ------------------------------------------------------------------
! START PRINT OUT
! ------------------------------------------------------------------
!     1. Print out number of upstream reaches for each tributaries
!     do iTrib=1,nTrib
!        nUpSegs=size(river_basin(iOut)%tributary(iTrib)%segOrder)
!        print*,'iTrib, nUpSegs = ', iTrib, nUpSegs
!     end do

!     2. Print out for each mainstem level
!     write(*,'(A)') 'level mainstem nUps'
!     do jLevel =1,maxLevel
!       nMains = 0
!       if (allocated(river_basin(iOut)%level(jLevel)%mainstem)) then
!         nMains = size(river_basin(iOut)%level(jLevel)%mainstem)
!         allocate(nSegMain(nMains))
!         do iMain=1,nMains
!           nSeg_tmp = size(river_basin(iOut)%level(jLevel)%mainstem(iMain)%segIndex)
!           associate(idx => river_basin(iOut)%level(jLevel)%mainstem(iMain)%segIndex(nSeg_tmp))
!           nSegMain(iMain) = size(structNTOPO(idx)%var(ixNTOPO%allUpSegIndices)%dat)
!           end associate
!           write(*,'(I2,x,I3,x,I8)') jLevel, iMain, nSegMain(iMain)
!         enddo
!         deallocate(nSegMain)
!       endif
!     enddo

!     3. Print mainstem code if in one of the mainstems for each reach
!     3.1. assign mainstem code to the reach that belong to that mainstem
!     levelArray = '-999'
!     do jLevel = 1,maxLevel
!       if (allocated(river_basin(iOut)%level(jLevel)%mainstem)) then
!         do iMain=1,size(river_basin(iOut)%level(jLevel)%mainstem)
!           nSeg_tmp = size(river_basin(iOut)%level(jLevel)%mainstem(iMain)%segIndex)
!           allocate(idxList(nSeg_tmp), stat=ierr)
!           idxList(1:nSeg_tmp) = river_basin(iOut)%level(jLevel)%mainstem(iMain)%segIndex(1:nSeg_tmp)
!           do iSeg = 1,size(idxList)
!             levelArray(idxList(iSeg)) = trim(mainstem_code(pfafs(idxList(iSeg))))
!           enddo
!           deallocate(idxList, stat=ierr)
!         enddo
!       endif
!     end do
!     3.2. compute number of tributary reaches and assign the tributary indix to the reaches within that tribuary
!     tribNseg = -999
!     tribIndex = -999
!     do iTrib=1,nTrib
!       nUpSegs=size(river_basin(iOut)%tributary(iTrib)%segIndex)
!       tribNseg(river_basin(iOut)%tributary(iTrib)%segIndex) = nUpSegs
!       tribIndex(river_basin(iOut)%tributary(iTrib)%segIndex) = iTrib
!     end do
!     do iSeg=1,nSeg
!       write(*,'(I10,x,A,x,I5,x,I5)') structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1), levelArray(iSeg), tribIndex(iSeg), tribNseg(iSeg)
!     end do

!     4. Print mainstem code of the basin that reach belongs to
!     levelArray = '-999'
!     do jLevel =maxLevel,1,-1
!       if (allocated(river_basin(iOut)%level(jLevel)%mainstem)) then
!         do iMain=1,size(river_basin(iOut)%level(jLevel)%mainstem)
!           nSeg_tmp = size(river_basin(iOut)%level(jLevel)%mainstem(iMain)%segIndex)
!           MainOutIndex = river_basin(iOut)%level(jLevel)%mainstem(iMain)%segIndex(nSeg_tmp)
!           nUps = size(structNTOPO(MainOutIndex)%var(ixNTOPO%allUpSegIndices)%dat(:))
!           allocate(idxList(nUps), stat=ierr)
!           idxList(1:nUps) = structNTOPO(MainOutIndex)%var(ixNTOPO%allUpSegIndices)%dat(1:nUps)
!           do iSeg = 1,size(idxList)
!             if (trim(levelArray(idxList(iSeg))) == '-999') then
!               levelArray(idxList(iSeg)) = trim(mainstem_code(pfafs(MainOutIndex)))
!             endif
!           enddo
!           deallocate(idxList, stat=ierr)
!         enddo
!       endif
!     end do
!     do iSeg=1,nSeg
!       write(*,'(I10,x,A)') structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1), levelArray(iSeg)
!     end do

!     5. Print mainstem level or tributary dangling logic if not in mainstem for each reach
!     do iSeg=1,nSeg
!       level1=-999
!       dangle=0
!       do jLevel =1,maxLevel
!         if (allocated(river_basin(iOut)%level(jLevel)%mainstem)) then
!           do iMain=1,size(river_basin(iOut)%level(jLevel)%mainstem)
!             i1 = findIndex(river_basin(iOut)%level(jLevel)%mainstem(iMain)%segIndex, iSeg, -999)
!             if (i1 /= -999) then
!               level1 = jLevel
!               exit
!             end if
!           enddo
!         endif
!         if (level1/=-999) exit
!       end do
!       ! to get dangle
!       if (level1==-999) exit
!         do iTrib=1,nTrib
!           nUpSegs=size(river_basin(iOut)%tributary(iTrib)%segIndex)
!           if (iSeg == river_basin(iOut)%tributary(iTrib)%segIndex(nUpSegs)) then
!             dangle=1
!             exit
!           end if
!         end do
!       endif
!       write(*,'(I10,x,I4,x,I1)') structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1), level1, dangle
!     end do
!     stop
! ------------------------------------------------------------------
! END PRINT OUT
! ------------------------------------------------------------------

end module domain_decomposition
