MODULE pfafstetter_module
USE nrtype
USE public_var
USE dataTypes,  only : var_clength         ! character type:   var(:)%dat
USE dataTypes,  only : var_ilength         ! integer type:     var(:)%dat
USE dataTypes,  only : basin               ! integer type:     var(:)%dat
USE var_lookup, only : ixPFAF              ! index of variables for the pfafstetter code
USE var_lookup, only : ixNTOPO             ! index of variables for the pfafstetter code
USE nr_utility_module, ONLY: indexx        ! Num. Recipies utilities
USE nr_utility_module, ONLY: indexTrue     ! Num. Recipies utilities
USE nr_utility_module, ONLY: findIndex     ! Num. Recipies utilities

implicit none

private

public :: classify_river_basin
public :: subset_order

interface char2int
  module procedure :: char2int_1d
  module procedure :: char2int_2d
end interface

contains

 subroutine classify_river_basin(nSeg, structPFAF, structNTOPO, river_basin, ierr, message)
  ! Process river network data using pfafstetter code and network topology
  ! Output: return populated basin dataType
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
  character(len=strLen)                       :: cmessage               ! error message from subroutine
  character(len=32)                           :: pfafs(nSeg)            ! pfaf_codes for all the segment
  character(len=32), allocatable              :: pfaf_out(:)            ! pfaf_codes for outlet segments
  character(len=32)                           :: upPfaf                 ! pfaf_code for one upstream segment
  character(len=32)                           :: outMainstemCode
  logical(lgt),      allocatable              :: mainstems(:,:)         ! logical to indicate segment is mainstem at each level
  logical(lgt),      allocatable              :: updated_mainstems(:,:) ! logical to indicate segment is mainstem at each level
  logical(lgt),      allocatable              :: lgc_trib_outlet(:)     ! logical to indicate segment is outlet of tributary
  logical(lgt)                                :: done                   ! logical
  integer(i4b)                                :: maxSegs=100
  integer(i4b)                                :: minSegs=20             ! if number of segments in a basin, flag it as minor (do not break up into mainstems and tributaries)
  integer(i4b)                                :: maxLevel=10
  integer(i4b)                                :: downIndex(nSeg)        ! downstream segment index for all the segments
  integer(i4b)                                :: segIndex(nSeg)         ! reach index for all the segments
  integer(i4b)                                :: order(nSeg)            ! reach order for all the segments
  integer(i4b), allocatable                   :: segIndex_out(:)        ! index for outlet segment
  integer(i4b), allocatable                   :: tributary_outlet(:)    ! index for tributary outlet reach
  integer(i4b)                                :: level                  ! manstem level
  integer(i4b)                                :: nOuts                  ! number of outlets
  integer(i4b)                                :: nUpSegs                ! number of upstream segments for a specified segment
  integer(i4b)                                :: tot_trib               ! Number of Tributary basins
  integer(i4b), allocatable                   :: upSegIndex(:)          ! upstream segment indices
  integer(i4b), allocatable                   :: nTrib(:)               ! Number of segments in each tributary basin
  integer(i4b), allocatable                   :: rank_nTrib(:)
  integer(i4b), allocatable                   :: trPos(:)
  integer(i4b), allocatable                   :: msPos(:)
  integer(i4b)                                :: iSeg,jSeg,iOut,iTrib   ! loop indices
  integer(i4b)                                :: i1,i2,jLevel,level1,dangle

  ierr=0; message='classify_river_basin/'

  forall(iSeg=1:nSeg) pfafs(iSeg)     = structPFAF(iSeg)%var(ixPFAF%code)%dat(1)
  forall(iSeg=1:nSeg) downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)
  forall(iSeg=1:nSeg) segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
  forall(iSeg=1:nSeg) order(iSeg)     = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)

  ! Number of outlets
  nOuts=count(downIndex<0)

  allocate(river_basin(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating river_basin'; return; endif

  allocate(pfaf_out(nOuts), segIndex_out(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating pfaf_out or segIndex_out'; return; endif

  ! Outlet information - pfaf code and segment index
  pfaf_out = pack(pfafs, downIndex<0)
  segIndex_out = pack(segIndex, downIndex<0)

 ! river_basin(1:nOuts)%isMajor  = .true.

  do iOut = 1,nOuts

    print*, 'working on outlet:'//trim(adjustl(pfaf_out(iOut)))

    allocate(river_basin(iOut)%mainstem(maxLevel), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%mainstem'; return; endif

    river_basin(iOut)%outIndex = segIndex_out(iOut)

 !   if (size(structNTOPO(segIndex_out(iOut))%var(ixNTOPO%allUpSegIndices)%dat) < minSegs) then
 !     river_basin(iOut)%isMajor  = .false.
 !     cycle
 !   endif

    ! Identify pfaf level given a river network
    call get_common_pfaf(pfafs, pfaf_out(iOut), level, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Identify mainstem segments at all levels (up to maxLevel)
    call lgc_mainstems(pfafs, pfaf_out(iOut), maxLevel, mainstems, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Initial assignment of mainstem segments
    allocate(updated_mainstems(size(mainstems,1),size(mainstems,2)), stat=ierr)
    updated_mainstems = .false.
    updated_mainstems(:,level) = mainstems(:,level)

    ! Identify the lowest level mainstem segment index
    call indexTrue(mainstems(:,level), msPos)
    allocate(river_basin(iOut)%mainstem(level)%segIndex(size(msPos)), &
             river_basin(iOut)%mainstem(level)%segOrder(size(msPos)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(:)%mainstem(:)%segIndex or segOrder'; return; endif
    river_basin(iOut)%mainstem(level)%segIndex(1:size(msPos))  = msPos(1:size(msPos))
    ! Compute reach order for mainstem segment
    call subset_order(order, &
                      river_basin(iOut)%mainstem(level)%segIndex, &
                      river_basin(iOut)%mainstem(level)%segOrder, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    river_basin(iOut)%mainstem(level)%nRch = size(river_basin(iOut)%mainstem(level)%segIndex)

    ! Identify tributary outlets into a mainstem at the lowest level
    !  i.e. the segment that is not on any mainstems AND flows into any mainstem segments
    call lgc_tributary_outlet(mainstems(:,:level), downIndex, lgc_trib_outlet, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    deallocate(mainstems, segIndex_out, stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem deallocating mainstems'; return; endif

    do
      level = level + 1

      ! number of tributary basins
      tot_trib = count(lgc_trib_outlet)

      allocate(nTrib(tot_trib),rank_nTrib(tot_trib), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating nTrib'; return; endif

      ! Extract array elements with only tributary outlet (keep indices in master array
      call indexTrue(lgc_trib_outlet, trPos)

      ! number of tributary segments
      do iTrib=1,tot_trib
        nTrib(iTrib) = size(structNTOPO(trPos(iTrib))%var(ixNTOPO%allUpSegIndices)%dat)
      end do

      ! rank the tributary outlets based on number of upstream segments
      call indexx(nTrib, rank_nTrib)

      ! Identify mainstems of large tributaries (# of segments > maxSeg) and update updated_mainstems.
      done=.true.
      do iTrib=tot_trib,1,-1
        ! nUpSegs = nTrib(rank_nTrib(iTrib))
        if (nTrib(rank_nTrib(iTrib)) > maxSegs) then
          write(*,'(A, A,I)') 'Exceed maximum number of segments: ', pfafs(trPos(rank_nTrib(iTrib))), nTrib(rank_nTrib(iTrib))

          ! Outlet mainstem code
          outMainstemCode = mainstem_code(pfafs(trPos(rank_nTrib(iTrib))))

          nUpSegs = size(structNTOPO(trPos(rank_nTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat)
          allocate(upSegIndex(nUpSegs), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating upSegIndex'; return; endif

          upSegIndex(:) = structNTOPO(trPos(rank_nTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat
          do jSeg=1,nUpSegs
            upPfaf = pfafs(upSegIndex(jSeg))
            if (trim(mainstem_code(upPfaf)) /= trim(outMainstemCode)) cycle
            updated_mainstems(upSegIndex(jSeg),level) = .true.
          end do

          deallocate(upSegIndex, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem deallocating nSegIndex'; return; endif

          done=.false.

        else ! if all the tributaries have segments less than maxSegs, exit
          exit
        endif
      end do ! tributary loop

      if (done) then ! if no mainstem/tributary updated, update tributary reach info, and then exist loop
        ! get tributary outlet segment index
        allocate(tributary_outlet(tot_trib), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating ributary'; return; endif
        tributary_outlet = trPos

        allocate(river_basin(iOut)%tributary(tot_trib), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%tributary'; return; endif

        do iTrib = 1,tot_trib
          allocate(river_basin(iOut)%tributary(iTrib)%segIndex(nTrib(iTrib)), &
                   river_basin(iOut)%tributary(iTrib)%segOrder(nTrib(iTrib)), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(:)%tributary(:)%segIndex or %segOrder'; return; endif
          ! compute reach index for tributary segments
          river_basin(iOut)%tributary(iTrib)%segIndex(:) = structNTOPO(tributary_outlet(iTrib))%var(ixNTOPO%allUpSegIndices)%dat
          ! compute reach order for tributary segments
          call subset_order(order, &
                            river_basin(iOut)%tributary(iTrib)%segIndex, &
                            river_basin(iOut)%tributary(iTrib)%segOrder, ierr, cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
          river_basin(iOut)%tributary(iTrib)%nRch = size(river_basin(iOut)%tributary(iTrib)%segIndex)
        end do
        exit
      else ! if mainstem/tributary updated, store mainstem reach info at current level, lgc_trib_outlet and go onto next level
        ! Compute reach index for mainstem segments
        call indexTrue(updated_mainstems(:,level), msPos)
        allocate(river_basin(iOut)%mainstem(level)%segIndex(size(msPos)), &
                 river_basin(iOut)%mainstem(level)%segOrder(size(msPos)), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(:)%mainstem(:)%segIndex or segOrder'; return; endif
        river_basin(iOut)%mainstem(level)%segIndex(1:size(msPos)) = msPos(1:size(msPos))
        ! Compute reach order for mainstem segments
        call subset_order(order, &
                          river_basin(iOut)%mainstem(level)%segIndex, &
                          river_basin(iOut)%mainstem(level)%segOrder, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        river_basin(iOut)%mainstem(level)%nRch = size(river_basin(iOut)%mainstem(level)%segIndex)

        ! update lgc_trib_outlet based on added mainstem
        call lgc_tributary_outlet(updated_mainstems(:,:level), downIndex, lgc_trib_outlet, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        deallocate(nTrib, rank_nTrib, stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem deallocating nTrib or rank_nTrib'; return; endif
      endif
    end do

!    do iSeg=1,nSeg
!      level1=-999
!      dangle=0
!      do jLevel =1,maxLevel
!        if (allocated(river_basin(iOut)%mainstem(jLevel)%segIndex)) then
!          i1 = findIndex(river_basin(iOut)%mainstem(jlevel)%segIndex, iSeg, -999)
!          if (i1 /= -999) then
!            level1 = jLevel
!            exit
!          end if
!        endif
!      end do
!
!      do iTrib=1,tot_trib
!        nUpSegs=size(river_basin(iOut)%tributary(iTrib)%segOrder)
!        associate(outIndex => river_basin(iOut)%tributary(iTrib)%segOrder(nUpSegs))
!        if (iSeg == river_basin(iOut)%tributary(iTrib)%segIndex(outIndex)) then
!          dangle=1
!          exit
!        end if
!        end associate
!      end do
!
!      write(*,'(A,I,I)') pfafs(iSeg), level1, dangle
!    end do

  end do ! outlet loop

 end subroutine classify_river_basin


 subroutine subset_order(order, subIndices, new_order, ierr, message)
  ! Input variables
  integer(i4b),              intent(in)       :: order(:)
  integer(i4b),              intent(in)       :: subIndices(:)
  ! Output variables
  integer(i4b), allocatable, intent(out)      :: new_order(:)
  integer(i4b),              intent(out)      :: ierr
  character(len=strLen),     intent(out)      :: message          ! error message
  ! local variables
  integer(i4b)                                :: iSeg
  integer(i4b), allocatable                   :: rank_order(:)
  integer(i4b), allocatable                   :: tmp_order(:)

  ierr=0; message='subset_order/'

  allocate(rank_order(size(order)), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating rank_orderx'; return; endif

  allocate(tmp_order(size(subIndices)), new_order(size(subIndices)), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating tmp_order or new_order'; return; endif
  ! rank order
  call indexx(order,rank_order)

  ! Extract ranked order for a given index list
  forall(iSeg=1:size(subIndices)) tmp_order(iSeg) = rank_order(subIndices(iSeg))

  ! Sort extracted ranked order in ascending order ==> new_order
  call indexx(tmp_order, new_order)

 end subroutine subset_order


 subroutine lgc_tributary_outlet(mainstem, ixDown, isDangle, ierr, message)
  ! Return 1D logical array thats indicate tributary outlet
  implicit none
  ! Input variables
  logical(lgt),          intent(in)      :: mainstem(:,:)
  integer(i4b),          intent(in)      :: ixDown(:)
  ! Output variables
  logical(lgt), allocatable, intent(out) :: isDangle(:)
  integer(i4b), intent(out)              :: ierr
  character(len=strLen),intent(out)      :: message          ! error message
  ! Local variables
  logical(lgt)                           :: isMainstem(size(ixDown))
  integer(i4b)                           :: iSeg

  ierr=0; message='lgc_tributary_outlet/'

  isMainstem = (count(mainstem,2)/=0)

  allocate(isDangle(size(ixDown)), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating isDangle'; return; endif

  isDangle = .false.
  do iSeg=1,size(isMainstem)
    if (ixDown(iSeg) > 0) then
      isDangle(iSeg) = (.not.isMainstem(iSeg) .and. isMainstem(ixDown(iSeg)))
    endif
  enddo

 end subroutine lgc_tributary_outlet


 subroutine lgc_mainstems(pfafs, pfaf_out, nLevel, mainstem, ierr, message)
  ! Return 2D logical array (nseg x nlevel) to indicate which level of mainstem is assigned to each segment
  implicit none
  ! Input variables
  character(*), intent(in)               :: pfafs(:)         ! input: pfafstetter code list
  character(*), intent(in)               :: pfaf_out         ! input: pfafstetter code at outlet
  integer(i4b), intent(in)               :: nLevel           ! input:  from comLevel
  ! output variables
  logical(lgt), allocatable, intent(out) :: mainstem(:,:)
  integer(i4b), intent(out)              :: ierr
  character(len=strLen),intent(out)      :: message          ! error message
  ! Local variables
  integer(i4b), allocatable              :: pfaf_int(:,:)
  integer(i4b)                           :: comLevel
  integer(i4b)                           :: iLevel
  integer(i4b)                           :: iSeg
  integer(i4b)                           :: idx
  character(len=strLen)                  :: cmessage       ! error message from subroutine
  character(len=2)                       :: x1             ! string converted from integer

  ierr=0; message='lgc_mainstems/'

  call get_common_pfaf(pfafs, pfaf_out, comLevel, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (nLevel .le. comLevel) then
    write(x1,'(I2)') comLevel ! converting integer to string using a 'internal file'
    ierr=10; message=trim(message)//'nLevel has to be greater than comLevel:'//trim(x1); return
  endif
  if (nLevel .ge. len(pfafs)) then;
    write(x1,'(I2)') len(pfafs) ! converting integer to string using a 'internal file'
    ierr=10; message=trim(message)//'nLevel has to be smaller than max. pfafs size:'//trim(x1); return
  endif

  allocate(mainstem(size(pfafs),nLevel), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating mainstem'; return; endif

  call char2int(pfafs, pfaf_int, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  mainstem = .false.

  do iSeg = 1, size(pfafs)

    if ( all(pfaf_int(iSeg,:)==-1) ) cycle

    do iLevel = size(pfaf_int(iSeg,:)),comLevel+1,-1
      if (pfaf_int(iSeg,iLevel) == -1) cycle
      if (mod(pfaf_int(iSeg,iLevel),2)==0) exit
    enddo

    idx = size(pfaf_int(iSeg,:iLevel))

    if ( idx == nLevel .and. pfaf_int(iSeg,idx+1)/=-1) then
      mainstem(iSeg, idx) = .true.
    else if ( idx .lt. nLevel ) then
      mainstem(iSeg, idx) = .true.
    endif

  enddo

 end subroutine lgc_mainstems


 character(len=32) function mainstem_code(pfaf)
  ! mainstem: mainstem code (trimd  n last consecutive odd digits)
  implicit none
  ! input variables
  character(*), intent(in)       :: pfaf          ! a pfafstetter code
  ! locat variables
  integer(i4b)                   :: pfaf_int      ! integer of one fafstetter code digit
  integer(i4b)                   :: iDigit        ! index of pfafstetter code digit
  integer(i4b)                   :: nDigits       ! number of pfafstetter code digit

  mainstem_code = '-999'   ! indicating segment is not manstem

  if (trim(adjustl(pfaf)) == '-999') then
   return
  end if

  nDigits = len(pfaf)

  do iDigit = nDigits, 1, -1

    if (pfaf(iDigit:iDigit)=="") cycle !if character digit is empty, go next

    read(pfaf(iDigit:iDigit),'(i)') pfaf_int

    if (mod(pfaf_int, 2) == 0) exit

    if (iDigit /= 1) then
      mainstem_code = trim(adjustl(pfaf(1:iDigit-1)))
    else
      mainstem_code = trim(adjustl(pfaf(1:iDigit)))
    end if

  end do

 end function mainstem_code

 logical(lgt) function isUpstream(pfaf_a, pfaf_b)

  ! check if "pfaf_a" is an upstream segment of "pfaf_b"
  ! return true, if pfaf_a is upstream of pfaf_b

  implicit none

  ! Input variables
  character(*), intent(in)       :: pfaf_a         ! a pfafstetter code for seg a
  character(*), intent(in)       :: pfaf_b         ! a pfafstetter code for seg b
  ! Local variables
  integer(i4b)                   :: nDigits       ! number of pfafstetter code digit
  integer(i4b)                   :: pfaf_a_int    ! integer of one fafstetter code digit
  integer(i4b)                   :: pfaf_b_int    ! integer of one fafstetter code digit
  integer(i4b)                   :: iDigit        ! index of pfafstetter code digit
  integer(i4b)                   :: nth           ! number of leading digits that two pfaf codes match

  nDigits = min(len(pfaf_a), len(pfaf_b))

  ! Find first nth digits that match
  nth = 0
  do iDigit =1,nDigits
    if ( pfaf_a(iDigit:iDigit) == pfaf_b(iDigit:iDigit) ) then
      nth = nth + 1
    else
      exit
    endif
  end do

  isUpstream = .true.

  ! if none of digit matches, A and B are not the same river system
  if (nth == 0) then
    isUpstream = .false.
  else
    if (len(pfaf_b) /= nth .and. len(pfaf_a) /= nth) then
      read(pfaf_a(nth+1:nth+1),'(i)') pfaf_a_int
      read(pfaf_b(nth+1:nth+1),'(i)') pfaf_b_int
      if (pfaf_b_int .gt. pfaf_a_int) then
        isUpstream = .false.
      else
        do iDigit = nth+1, len(pfaf_b)
          read(pfaf_b(iDigit:iDigit),'(i)') pfaf_b_int
          if ( mod(pfaf_b_int, 2) == 0) then
            isUpstream = .false.
            exit
          endif
        end do
      endif
    elseif (len(pfaf_b) /= nth .and. len(pfaf_a) == nth) then
      isUpstream = .false.
    end if
  end if

 end function isUPstream


 subroutine get_common_pfaf(pfafs, pfaf_out, comLevel, ierr, message)
  ! Return level of pfaf code from left side where integers are common to the ones of a given outlet segment
  implicit none
  ! Input variables
  character(*), intent(in)               :: pfafs(:)
  character(*), intent(in)               :: pfaf_out
  ! output variables
  integer(i4b), intent(out)              :: comLevel
  integer(i4b), intent(out)              :: ierr
  character(len=strLen),intent(out)      :: message          ! error message
  ! Local variables
  integer(i4b)                           :: iLevel
  integer(i4b), allocatable              :: pfaf_int(:,:)
  integer(i4b), allocatable              :: pfaf_out_int(:)
  character(len=strLen)                  :: cmessage         ! error message from subroutine

  ierr=0; message='get_common_pfaf/'

  call char2int(pfafs, pfaf_int, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call char2int(pfaf_out, pfaf_out_int, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  do iLevel = 1,len(pfafs)
   if ( any(pfaf_int(:,iLevel) /= pfaf_out_int(iLevel) .and. pfaf_int(:,iLevel) /= -1) ) exit
  end do

  comLevel = iLevel - 1

 end subroutine get_common_pfaf


 subroutine char2int_1d(char_array, int_array, ierr, message)
  ! Convert character array to one digit integer array
  ! if character array is '-9999' or '0', int_array(:) = (/-1,...-1/)

  implicit none
  ! Input variables
  character(*), intent(in)                :: char_array
  ! output variables
  integer(i4b), allocatable, intent(out)  :: int_array(:)
  integer(i4b), intent(out)               :: ierr
  character(len=strLen),intent(out)       :: message          ! error message
  ! local variables
  character(len=strLen)                   :: string
  integer(i4b)                            :: str_len
  integer(i4b)                            :: iChr

  ierr=0; message='char2int_1d/'

  allocate(int_array(len(char_array)), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating int_array'; return; endif

  int_array = -1

  str_len = len(trim(adjustl(char_array)))
  string = adjustl(char_array)

  if (trim(string) == '-9999' .or. trim(string) == '0') then
    int_array(1) = -1
  else
    do iChr =1,str_len
      read(string(iChr:iChr),'(i)') int_array(iChr)
    end do
  endif

  end subroutine char2int_1d

 subroutine char2int_2d(char_array, int_array, ierr, message)
  ! convert character array to one digit integer array
  ! if character array is '-9999' or '0', int_array(i,:) = (/-1,...-1/)

  implicit none
  ! Input variables
  character(*), intent(in)                :: char_array(:)
  ! output variables
  integer(i4b), allocatable, intent(out)  :: int_array(:,:)
  integer(i4b), intent(out)               :: ierr
  character(len=strLen),intent(out)       :: message          ! error message
  ! local variables
  character(len=strLen)                   :: string
  integer(i4b)                            :: str_len
  integer(i4b)                            :: iSize
  integer(i4b)                            :: iChr

  ierr=0; message='char2int_2d/'

  allocate(int_array(size(char_array),len(char_array)), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating int_array'; return; endif

  int_array = -1

  do iSize = 1, size(char_array)
    str_len = len(trim(adjustl(char_array(iSize))))
    string = adjustl(char_array(iSize))
    if (trim(string) == '-9999' .or. trim(string) == '0') then
      int_array(iSize, 1) = -1
    else
      do iChr =1,str_len
        read(string(iChr:iChr),'(i)') int_array(iSize, iChr)
      end do
    end if
  end do

  end subroutine char2int_2d

 subroutine group_array(narray, group, ierr, message)
  ! Group array elements based on number
  implicit none
  ! Input variables
  integer(i4b),          intent(in)  :: narray(:)      ! number of stream segments
  ! Output variables
  integer(i4b),          intent(out) :: group(:)           ! process group
  integer(i4b),          intent(out) :: ierr
  character(len=strLen), intent(out) :: message            ! error message
  ! Local variables
  integer(i4b),allocatable           :: rank_array(:)      ! number of stream segments
  integer(i4b)                       :: maxNgroup=36
  integer(i4b)                       :: arraySize
  integer(i4b)                       :: maxElem
  integer(i4b)                       :: cc
  integer(i4b)                       :: iGroup
  integer(i4b)                       :: ii                  ! loop index

  ierr=0; message='group_array/'

  arraySize = size(narray)
  if (arraySize >= maxNgroup) arraySize=maxNgroup

  maxElem = sum(narray)/arraySize

  allocate(rank_array(size(narray)),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating rank_array'; return; endif

  call indexx(narray,rank_array)

  iGroup = 1
  cc = 0
  do ii=size(narray),1,-1

    group(rank_array(ii)) = iGroup
    cc = cc + narray(rank_array(ii))
    if (cc .gt. maxElem)then
       iGroup=iGroup+1
       cc=0
     endif

  end do

  if (iGroup .gt. arraySize) then
   print*,'yes'
  end if

 end subroutine group_array

end module pfafstetter_module

!      subroutine group_river_network(nSeg, structPFAF, structNTOPO, ierr, message)
!       !
!       implicit none
!       ! Input variables
!       integer(i4b),                   intent(in)  :: nSeg                   ! number of stream segments
!       type(var_clength), allocatable, intent(in)  :: structPFAF(:)          ! pfafstetter code
!       type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)         ! network topology
!       ! output variables
!       integer(i4b),                   intent(out) :: ierr
!       character(len=strLen),          intent(out) :: message                ! error message
!       ! Local variables
!       character(len=strLen)                       :: cmessage               ! error message from subroutine
!       type(basin), allocatable                    :: river_basin(:)         ! river basin data
!       !integer(i4b)                                :: order(nSeg)            ! process order
!       integer(i4b)                                :: tot_trib               ! Number of Tributary basins
!       integer(i4b)                                :: nOuts
!       integer(i4b)                                :: iSeg,jSeg,iOut         ! loop indices
!       integer(i4b)                                :: segIndex_trib          ! index of tributary outlet segment
!       !integer(i4b), allocatable                   :: new_order(:)
!       integer(i4b), allocatable                   :: tribGroup(:)
!       integer(i4b), allocatable                   :: nTrib(:)
!       !integer(i4b), allocatable                   :: msGroup(:,:)
!
!       ierr=0; message='group_river_network/'
!
!       call classify_river_basin(nSeg, structPFAF, structNTOPO, river_basin, ierr, message)
!       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!
!       nOuts = size(river_basin)
!
!       do iOut=1,nOuts
!         if (river_basin(iOut)%isMajor) then
!           ! Total number of tributary
!           tot_trib=size(river_basin(iOut)%tributary_outlet)
!
!           allocate(nTrib(tot_trib),tribGroup(tot_trib), stat=ierr)
!           if(ierr/=0)then; message=trim(message)//'problem allocating nTrib'; return; endif
!
!           do iSeg=1,tot_trib
!             segIndex_trib=river_basin(iOut)%tributary_outlet(iSeg)
!             nTrib(iSeg) = size(structNTOPO(segIndex_trib)%var(ixNTOPO%allUpSegIndices)%dat)
!           end do
!           ! Group tributary
!           call group_array(nTrib, tribGroup, ierr, message)
!
!           deallocate(nTrib, stat=ierr)
!           if(ierr/=0)then; message=trim(message)//'problem deallocating nTrib'; return; endif
!
!      !     print*, tribGroup
!
!         end if
!
!       end do
!
!       !update reachorder for each tributaries
!     !  forall(iSeg=1:nSeg) order(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)
!
!     !  do iSeg=1,tot_trib
!     !    segIndex_trib=river_basin(iOut)%tributary_outlet(iSeg)
!     !    nTrib(iSeg) = size(structNTOPO(segIndex_trib)%var(ixNTOPO%allUpSegIndices)%dat)
!     !
!     !    allocate(upSegIndex(nTrib(iSeg)), new_order(nTrib(iSeg)), stat=ierr)
!     !    if(ierr/=0)then; message=trim(message)//'problem allocating upSegIndex'; return; endif
!     !
!     !    ! Extract upstream segment indices
!     !    upSegIndex(:) = structNTOPO(segIndex_trib)%var(ixNTOPO%allUpSegIndices)%dat
!     !
!     !    call subset_order(order, upSegIndex, new_order, ierr, cmessage)
!     !    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!     !
!     !    print*, '------------------------- '
!     !    print*, 'segment: ', trim(pfafs(segIndex_trib))
!     !    print*, 'No. upstream seg= ', nTrib(iSeg)
!     !    do jSeg=1,size(upSegIndex)
!     !      print*, jSeg, trim(pfafs(upSegIndex(jSeg))), new_order(jSeg)
!     !    enddo
!     !
!     !    deallocate(upSegIndex, new_order, stat=ierr)
!     !
!     !  end do
!      stop
!
!      end subroutine


! JUNK LINES
! print out mainstem, and tributary outlet
  !  do iSeg=1,nSeg
  !    level1=-999
  !    dangle=0
  !    do jLevel =1,maxLevel
  !      if (updated_mainstems(iSeg,jLevel)) then
  !        level1 = jLevel
  !        exit
  !      end if
  !    end do
  !   if (lgc_trib_outlet(iSeg)) dangle=1
  !   write(*,'(A,I,I)') pfafs(iSeg), level1, dangle
  !  end do
