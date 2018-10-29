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

interface char2int
  module procedure :: char2int_1d
  module procedure :: char2int_2d
end interface

contains

 subroutine classify_river_basin(nSeg, structPFAF, structNTOPO, river_basin, ierr, message)
  ! Identify tributary basin and mainstems using river network data and pfafstetter code
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
  integer(i4b)                                :: maxLevel=10
  integer(i4b)                                :: downIndex(nSeg)        ! downstream segment index for all the segments
  integer(i4b)                                :: segIndex(nSeg)         ! reach index for all the segments
  integer(i4b)                                :: segOrder(nSeg)         ! reach order for all the segments
  integer(i4b), allocatable                   :: segIndex_out(:)        ! index for outlet segment
  integer(i4b)                                :: level                  ! manstem level
  integer(i4b)                                :: nOuts                  ! number of outlets
  integer(i4b)                                :: nUpSegs                ! number of upstream segments for a specified segment
  integer(i4b)                                :: nTrib                  ! Number of Tributary basins
  integer(i4b), allocatable                   :: nSegTrib(:)            ! Number of segments in each tributary basin
  integer(i4b), allocatable                   :: rankSegTrib(:)         ! index for ranked tributary based on num. of upstream reaches
  integer(i4b), allocatable                   :: segOrderTrib(:)        ! index for tributary upstream reach order
  integer(i4b), allocatable                   :: trPos(:)               ! tributary outlet indices
  integer(i4b), allocatable                   :: msPos(:)               ! mainstem indices
  integer(i4b)                                :: iSeg,jSeg,iOut         ! loop indices
  integer(i4b)                                :: iTrib,jTrib            ! loop indices
  integer(i4b)                                :: i1,i2,jLevel,level1,dangle

  ierr=0; message='classify_river_basin/'

  forall(iSeg=1:nSeg) pfafs(iSeg)     = structPFAF(iSeg)%var(ixPFAF%code)%dat(1)
  forall(iSeg=1:nSeg) downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)
  forall(iSeg=1:nSeg) segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
  forall(iSeg=1:nSeg) segOrder(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)

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

    allocate(river_basin(iOut)%mainstem(maxLevel), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%mainstem'; return; endif

    river_basin(iOut)%outIndex = segIndex_out(iOut)

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
    call subset_order(segOrder, &
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
      nTrib = count(lgc_trib_outlet)

      allocate(nSegTrib(nTrib),rankSegTrib(nTrib), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating nSegTrib'; return; endif

      ! Extract array elements with only tributary outlet (keep indices in master array
      call indexTrue(lgc_trib_outlet, trPos)

      ! number of tributary segments
      do iTrib=1,nTrib
        nSegTrib(iTrib) = size(structNTOPO(trPos(iTrib))%var(ixNTOPO%allUpSegIndices)%dat)
      end do

      ! rank the tributary outlets based on number of upstream segments
      call indexx(nSegTrib, rankSegTrib)

      ! Identify mainstems of large tributaries (# of segments > maxSeg) and update updated_mainstems.
      done=.true.
      do iTrib=nTrib,1,-1

        if (nSegTrib(rankSegTrib(iTrib)) > maxSegs) then
          write(*,'(A,A,I5)') 'Exceed maximum number of segments: ', pfafs(trPos(rankSegTrib(iTrib))), nSegTrib(rankSegTrib(iTrib))

          ! Outlet mainstem code
          outMainstemCode = mainstem_code(pfafs(trPos(rankSegTrib(iTrib))))

          nUpSegs = size(structNTOPO(trPos(rankSegTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat)

          associate(upSegIndex => structNTOPO(trPos(rankSegTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat)
          do jSeg=1,nUpSegs
            upPfaf = pfafs(upSegIndex(jSeg))
            if (trim(mainstem_code(upPfaf)) /= trim(outMainstemCode)) cycle
            updated_mainstems(upSegIndex(jSeg),level) = .true.
          end do
          end associate

          done=.false.

        else ! if all the tributaries have segments less than maxSegs, exit

          exit

        endif

      end do ! tributary loop

      if (done) then ! if no mainstem/tributary updated, update tributary reach info, and then exist loop

        allocate(river_basin(iOut)%tributary(nTrib), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating river_basin%tributary'; return; endif

        do iTrib=nTrib,1,-1

          jTrib = nTrib - iTrib + 1

          allocate(river_basin(iOut)%tributary(jTrib)%segIndex(nSegTrib(rankSegTrib(iTrib))), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(:)%tributary%segIndex'; return; endif
          allocate(segOrderTrib(nSegTrib(rankSegTrib(iTrib))), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating segOrderTrib'; return; endif

          ! compute reach index for tributary segments
          associate( segIndexTrib => structNTOPO(trPos(rankSegTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat)

          ! compute reach order for tributary segments
          call subset_order(segOrder, segIndexTrib, segOrderTrib, ierr, cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

!          do iSeg = 1,nSegTrib(rankSegTrib(iTrib))
!            river_basin(iOut)%tributary(jTrib)%segIndex(iSeg) = segIndexTrib(segOrderTrib(iSeg))
!          end do
          river_basin(iOut)%tributary(jTrib)%segIndex(:) = segIndexTrib(segOrderTrib)
          end associate

          river_basin(iOut)%tributary(jTrib)%nRch = nSegTrib(rankSegTrib(iTrib))

          deallocate(segOrderTrib, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem deallocating segOrderTrib'; return; endif

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
        call subset_order(segOrder, &
                          river_basin(iOut)%mainstem(level)%segIndex, &
                          river_basin(iOut)%mainstem(level)%segOrder, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        river_basin(iOut)%mainstem(level)%nRch = size(river_basin(iOut)%mainstem(level)%segIndex)

        ! update lgc_trib_outlet based on added mainstem
        call lgc_tributary_outlet(updated_mainstems(:,:level), downIndex, lgc_trib_outlet, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        deallocate(nSegTrib, rankSegTrib, stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem deallocating [nSegTrib, rankSegTrib]'; return; endif

      endif
    end do ! end of tributary update loop

!     do iTrib=1,nTrib
!        associate(jTrib => river_basin(iOut)%tributary_order(iTrib))
!        nUpSegs=size(river_basin(iOut)%tributary(jTrib)%segOrder)
!        print*,'iTrib, jTrib, nUpSegs = ', iTrib, jTrib, nUpSegs
!        end associate
!     end do

!     do iSeg=1,nSeg
!       level1=-999
!       dangle=0
!       do jLevel =1,maxLevel
!         if (allocated(river_basin(iOut)%mainstem(jLevel)%segIndex)) then
!           i1 = findIndex(river_basin(iOut)%mainstem(jlevel)%segIndex, iSeg, -999)
!           if (i1 /= -999) then
!             level1 = jLevel
!             exit
!           end if
!         endif
!       end do
!
!       do iTrib=1,nTrib
!         nUpSegs=size(river_basin(iOut)%tributary(iTrib)%segIndex)
!         if (iSeg == river_basin(iOut)%tributary(iTrib)%segIndex(nUpSegs)) then
!           dangle=1
!           exit
!         end if
!       end do
!       write(*,'(A,I4,x,I1)') pfafs(iSeg), level1, dangle
!     end do
!     stop

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

    read(pfaf(iDigit:iDigit),'(I1)') pfaf_int

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
      read(pfaf_a(nth+1:nth+1),'(I1)') pfaf_a_int
      read(pfaf_b(nth+1:nth+1),'(I1)') pfaf_b_int
      if (pfaf_b_int .gt. pfaf_a_int) then
        isUpstream = .false.
      else
        do iDigit = nth+1, len(pfaf_b)
          read(pfaf_b(iDigit:iDigit),'(I1)') pfaf_b_int
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
      read(string(iChr:iChr),'(I1)') int_array(iChr)
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
        read(string(iChr:iChr),'(I1)') int_array(iSize, iChr)
      end do
    end if
  end do

  end subroutine char2int_2d

end module pfafstetter_module


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
