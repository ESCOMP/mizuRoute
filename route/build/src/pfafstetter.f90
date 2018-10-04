MODULE pfafstetter_module

USE nrtype
USE public_var
USE dataTypes,  only : var_clength         ! character type:   var(:)%dat
USE dataTypes,  only : var_ilength         ! integer type:     var(:)%dat
USE var_lookup, only : ixPFAF,nVarsPFAF    ! index of variables for the pfafstetter code
USE var_lookup, only : ixNTOPO,nVarsNTOPO  ! index of variables for the pfafstetter code
USE nr_utility_module, ONLY: indexx        ! Num. Recipies utilities
USE nr_utility_module, ONLY: indexTrue     ! Num. Recipies utilities

implicit none

private

public :: process_pfaf

interface char2int
  module procedure :: char2int_1d
  module procedure :: char2int_2d
end interface

contains

 subroutine process_pfaf(nSeg, structPFAF, structNTOPO, isDangle, ierr, message)
  ! To obtain outlets of independent tributaries
  ! Output:  isDangle(nSeg),  True: outlet of tributaries
  implicit none
  ! Input variables
  integer(i4b),      intent(in)               :: nSeg             ! number of stream segments
  type(var_clength), intent(in), allocatable  :: structPFAF(:)    ! pfafstetter code
  type(var_ilength), intent(in), allocatable  :: structNTOPO(:)   ! network topology
  ! Output variables
  logical(lgt), allocatable,intent(out)       :: isDangle(:)
  integer(i4b),             intent(out)       :: ierr
  character(len=strLen),    intent(out)       :: message          ! error message
  ! Local variables
  character(len=strLen)                       :: cmessage           ! error message from subroutine
  character(len=32)                           :: pfafs(nSeg)
  character(len=32), allocatable              :: pfaf_out(:)
  character(len=32)                           :: upPfaf
  character(len=32)                           :: outMainstemCode
  logical(lgt), allocatable                   :: mainstems(:,:)     ! store mainstem level
  logical(lgt), allocatable                   :: mainstem(:,:)      ! store mainstem
  integer(i4b)                                :: maxSegs=300
  integer(i4b)                                :: maxLevel=10
  integer(i4b)                                :: downIndex(nseg)
  integer(i4b)                                :: comLevel
  integer(i4b), allocatable                   :: upSegIndex(:)
  integer(i4b)                                :: nOuts
  integer(i4b)                                :: nUpSegs
  integer(i4b)                                :: tot_trib           ! Number of Tributary basins
  integer(i4b), allocatable                   :: nTrib(:)           ! Number of segments in each tributary basin
  integer(i4b), allocatable                   :: rank_nTrib(:)
  integer(i4b), allocatable                   :: pos(:)
  integer(i4b)                                :: iSeg,jSeg,iOut   ! loop indices
  !integer(i4b)                                :: iLevel, jLevel,Level,dangle

  ierr=0; message='process_pfaf/'

  forall(iSeg=1:nSeg) pfafs(iSeg)     = structPFAF(iSeg)%var(ixPFAF%code)%dat(1)
  forall(iSeg=1:nSeg) downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)

  ! Number of outlets
  nOuts=count(downIndex<0)
  allocate(pfaf_out(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating pfaf_out'; return; endif

  ! Outlet pfaf codes
  pfaf_out = pack(pfafs, downIndex<0)

  do iOut = 1,nOuts

    print*, 'working on outlet:'//trim(adjustl(pfaf_out(iOut)))

    ! Identify pfaf level given a river network
    call get_common_pfaf(pfafs, pfaf_out(iOut), comLevel, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Identify mainstem segments at all levels (up to maxLevel)
    call lgc_mainstems(pfafs, pfaf_out(iOut), maxLevel, mainstems, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Initial assignment of mainstem segments
    allocate(mainstem(size(mainstems,1),1), stat=ierr)
    mainstem(:,1) = mainstems(:,comLevel+1)

    ! Identify tributary outlets into a mainstem at the lowest level
    !  i.e. the segment that is not on any mainstems AND flows into any mainstem segments
    call lgc_tributary_outlet(mainstems(:,:comLevel+1), downIndex, isDangle, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    !do

    ! number of tributary basins
    tot_trib = count(isDangle)

    allocate(nTrib(tot_trib),rank_nTrib(tot_trib), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating nTrib'; return; endif

    ! Extract array elements with only tributary outlet (keep indices in master array
    call indexTrue(isDangle,pos)

    ! number of tributary segments
    do iSeg=1,tot_trib
      nTrib(iSeg) = size(structNTOPO(pos(iSeg))%var(ixNTOPO%allUpSegIndices)%dat)
    end do

    ! rank the tributary outlets based on number of upstream segments
    call indexx(nTrib, rank_nTrib)

    ! Identify mainstems of large tributaries (# of segments > maxSeg)
    do iSeg=tot_trib,1,-1

      if (nTrib(rank_nTrib(iSeg)) > maxSegs) then
        write(*,'(A, A,I)') 'Exceed maximum number of segments: ', pfafs(pos(rank_nTrib(iSeg))), nTrib(rank_nTrib(iSeg))

        ! Outlet mainstem code
        outMainstemCode = mainstem_code(pfafs(pos(rank_nTrib(iSeg))))

        ! this tributary outlet is now mainstem
        mainstem(pos(rank_nTrib(iSeg)),1) = .true.

        nUpSegs = size(structNTOPO(pos(rank_nTrib(iSeg)))%var(ixNTOPO%allUpSegIndices)%dat)
        allocate(upSegIndex(nUpSegs), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating upSegIndex'; return; endif

        upSegIndex(:) = structNTOPO(pos(rank_nTrib(iSeg)))%var(ixNTOPO%allUpSegIndices)%dat

        do jSeg=1,nUpSegs
          upPfaf = pfafs(upSegIndex(jSeg))
          if (trim(mainstem_code(upPfaf)) /= trim(outMainstemCode)) cycle
          mainstem(jSeg,1) = .true.
        end do

        ! update isDangle based on added mainstem
        call lgc_tributary_outlet(mainstem, downIndex, isDangle, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      else

        exit

      endif

      deallocate(upSegIndex, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem deallocating nSegs'; return; endif

    end do

    !end do

    stop

    do iSeg=1,tot_trib
      write(*,'(A,L,I,I)') pfafs(pos(iSeg)),isDangle(pos(iSeg)), nTrib(iSeg), nTrib(rank_nTrib(iSeg))
    end do

    print*, trim(message)//'pause: '; read(*,*)

      !if (tot_tributaries .lt. threshSegs) exit
      !if (tot_tributaries .lt. threshSegs) then
      !  print*,iLevel, tot_tributaries
      !  do iSeg=1,nSeg
      !    level=-999
      !    dangle=0
      !    do jLevel =1,iLevel
      !      if (mainstem(iSeg,jLevel)) then
      !        level = jLevel
      !        cycle
      !      end if
      !    end do
      !   if (isDangle(iSeg)) dangle=1
      !   write(*,'(A,I,I)') pfafs(iSeg), level, dangle
      !  end do
      !  exit
      !endif

  end do

 end subroutine process_pfaf


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

  if (pfaf == '-999') then
   return
  end if

  mainstem_code = pfaf

  nDigits = len(pfaf)

  do iDigit = nDigits, 1, -1

    read(pfaf(iDigit:iDigit),'(i)') pfaf_int

    if (mod(pfaf_int, 2) == 0) then
      exit
    endif

    if (iDigit /= 1) then
      mainstem_code = pfaf(1:iDigit-1)
    else
      mainstem_code = pfaf(1:iDigit)
    end if

  end do

 end function mainstem_code


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

end module pfafstetter_module
