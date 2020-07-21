MODULE pfafstetter_module

USE nrtype
USE public_var

implicit none

private

public :: lgc_tributary_outlet
public :: find_mainstems
public :: lgc_mainstems
public :: mainstem_code
public :: isUPstream
public :: find_closed_basin
public :: get_common_pfaf

interface char2int
  module procedure :: char2int_1d
  module procedure :: char2int_2d
end interface

contains


 ! ----------------------------------------------------------------------------------------
 ! Tributary related routines/functions
 ! ----------------------------------------------------------------------------------------
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

 ! ----------------------------------------------------------------------------------------
 ! Mainstem related routines/functions
 ! ----------------------------------------------------------------------------------------
 subroutine find_mainstems(pfafCode, pfafs, isMainstem, ierr, message)
  ! Find mainstem reaches whoes mainstem code matches with "pfafCode"
  implicit none
  ! Input variables
  character(*), intent(in)               :: pfafs(:)         ! input: pfafstetter code list
  character(*), intent(in)               :: pfafCode         ! input: pfafstetter code at outlet
  ! output variables
  logical(lgt), allocatable, intent(out) :: isMainstem(:)
  integer(i4b), intent(out)              :: ierr
  character(len=strLen),intent(out)      :: message          ! error message
  ! Local variables
  integer(i4b)                           :: iSeg, iPfaf      ! loop indices
  integer(i4b)                           :: nSeg
  integer(i4b)                           :: pfafLen
  integer(i4b)                           :: pfaf_int
  integer(i4b)                           :: nPfaf            ! number of digits in a pfaf code
  integer(i4b)                           :: cc               ! counter
  character(len=32)                      :: pfaf             ! a pfafstetter code

  ierr=0; message='find_mainstems/'

  ! initialize output
  nSeg = size(pfafs)
  allocate(isMainstem(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating isMainstem'; return; endif
  isMainstem = .false.
  pfafLen=len(trim(pfafCode))

  if (nSeg == 1) then ! simple case of one stream segment
    isMainstem(1) = .true.
  else            ! standard case: 1)find odd sequence after pfafCode length and 2) digit at length of pfafCode equal to last digit of pfafCode
    do iSeg = 1,nSeg
      pfaf = adjustl(pfafs(iSeg))
      nPfaf = len(trim(pfaf))
      cc = 0
      do iPfaf = pfafLen,nPfaf
        read(pfaf(iPfaf:iPfaf),'(I1)') pfaf_int
        if (mod(pfaf_int, 2) == 0) exit
        cc = cc + 1
      enddo
      if (cc /= nPfaf-pfafLen+1) cycle
      if (pfaf(pfafLen:pfafLen) /= pfafCode(pfafLen:pfafLen)) cycle
      isMainstem(iSeg) = .true.
    enddo
  endif

 end subroutine find_mainstems


 subroutine lgc_mainstems(pfafs, pfaf_out, nLevel, mainstem, ierr, message)
  ! Return 2D logical array (nseg x nlevel) to indicate which level of mainstem is assigned to each segment
  implicit none
  ! Input variables
  character(*), intent(in)               :: pfafs(:)         ! input: segment pfafstetter code list
  character(*), intent(in)               :: pfaf_out         ! input: pfafstetter code at outlet segment
  integer(i4b), intent(in)               :: nLevel           ! input: from comLevel
  ! output variables
  logical(lgt), allocatable, intent(out) :: mainstem(:,:)    ! output:
  integer(i4b), intent(out)              :: ierr
  character(len=strLen),intent(out)      :: message          ! error message
  ! Local variables
  integer(i4b), allocatable              :: pfaf_int(:,:)    !
  integer(i4b)                           :: comLevel
  integer(i4b)                           :: iLevel
  integer(i4b)                           :: iSeg
  integer(i4b)                           :: idx
  character(len=strLen)                  :: cmessage         ! error message from subroutine
  character(len=2)                       :: x1               ! string converted from integer

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
  character(*), intent(in)       :: pfaf          ! a segment pfafstetter code
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

!    if (iDigit /= 1) then
!      mainstem_code = trim(adjustl(pfaf(1:iDigit-1)))
!    else
!      mainstem_code = trim(adjustl(pfaf(1:iDigit)))
!    end if
!
  end do

  mainstem_code = trim(adjustl(pfaf(1:iDigit)))

 end function mainstem_code


 ! ----------------------------------------------------------------------------------------
 ! Other pfafstetter decoding routines/functions
 ! ----------------------------------------------------------------------------------------
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


 subroutine find_closed_basin(pfafBasin, pfafs, isClosed, ierr, message)
  ! Given a list of pfafs of river network, identify closed basin reach within the basin defined by pfafBasin
  implicit none
  ! input variables
  character(len=32),    intent(in)   :: pfafBasin       ! a river basin pfaf code
  character(len=32),    intent(in)   :: pfafs(:)        ! segment pfaf codes arrays
  ! output variables
  logical(lgt),         intent(out)  :: isClosed(:)     ! logical to indicate segment is in closed basin
  integer(i4b),         intent(out)  :: ierr
  character(len=strLen),intent(out)  :: message         ! error message
  ! local variables
  logical(lgt)                       :: isInBasin       ! logical to indicate the reach is within the basin
  integer(i4b)                       :: iSeg,iPfaf      ! loop indices
  integer(i4b)                       :: nSeg            ! number of reaches
  integer(i4b)                       :: uniqLen         ! pfafBasin digit length
  character(len=32)                  :: pfaf            ! a pfafstetter code

  ierr=0; message='find_closed_basin/'

  ! initialize output
  nSeg = size(pfafs)
  if (size(isClosed)/=nSeg)then;ierr=10;message=trim(message)//'isClosed array sized is inconsistent with pfaf array size';return;endif
  isClosed(:) = .false.
  ! get the length of the basin ID (=pfafBasin)
  uniqLen = len(trim(pfafBasin))

  ! check if each reach is in closed basin
  do iSeg = 1,nSeg

    pfaf = pfafs(iSeg) ! a pfaf code

    ! check a "pfaf" basin is inside "pfafBasin" basin
    isInBasin = .true.
    do iPfaf = uniqLen,len(trim(pfaf))
      if (pfaf(iPfaf:iPfaf) /= pfafBasin(iPfaf:iPfaf)) then
        isInBasin = .false.
        exit
      endif
    end do

    ! if a "pfaf" basin is inside, check a basin is closed basin (include 0 digit after common pfaf code with "pfafBasin")
    if (isInBasin) then
      do iPfaf = uniqLen,len(trim(pfaf))
        if (pfaf(iPfaf:iPfaf) == '0') then
          isClosed(iSeg) = .true.
          exit
        end if
      end do
    end if

  end do ! seg loop

 end subroutine find_closed_basin


 ! ----------------------------------------------------------------------------------------
 ! character-to-integer routines/functions
 ! ----------------------------------------------------------------------------------------
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
