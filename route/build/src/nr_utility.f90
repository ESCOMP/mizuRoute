MODULE nr_utility_module
! contains functions that should really be part of the fortran standard, but are not

USE nrtype

implicit none

interface arth
  MODULE PROCEDURE arth_r, arth_d, arth_i4b, arth_i8b
end interface

interface indexx
  module procedure indexx_i4b
  module procedure indexx_i8b
end interface

interface swap
  module procedure swap_i4b
  module procedure swap_i8b
end interface

interface unique
  module procedure unique_i4b
  module procedure unique_i8b
end interface

INTERFACE char2int
  module procedure :: char2int_1d
  module procedure :: char2int_2d
END INTERFACE

private
public::arth
public::indexx
public::findIndex
public::match_index
public::indexTrue
public::unique
public::get_digits
public::char2int

CONTAINS

 ! *************************************************************************************************
 ! * the arth function, used to build a vector of regularly spaced numbers
 ! *************************************************************************************************
 FUNCTION arth_r(first,increment,n)
 implicit none
 REAL(SP), INTENT(IN) :: first,increment
 INTEGER(I4B), INTENT(IN) :: n
 REAL(SP), DIMENSION(n) :: arth_r
 INTEGER(I4B) :: k
 arth_r(1)=first
 if(n>1)then
  do k=2,n
   arth_r(k) = arth_r(k-1) + increment
  end do
 end if
 END FUNCTION arth_r
 ! ------------------------------------------------------------------------------------------------
 FUNCTION arth_d(first,increment,n)
 implicit none
 REAL(DP), INTENT(IN) :: first,increment
 INTEGER(I4B), INTENT(IN) :: n
 REAL(DP), DIMENSION(n) :: arth_d
 INTEGER(I4B) :: k
 arth_d(1)=first
 if(n>1)then
  do k=2,n
   arth_d(k) = arth_d(k-1) + increment
  end do
 end if
 END FUNCTION arth_d
 ! ------------------------------------------------------------------------------------------------
 FUNCTION arth_i4b(first,increment,n)
 implicit none
 INTEGER(I4B), INTENT(IN) :: first,increment,n
 INTEGER(I4B), DIMENSION(n) :: arth_i4b
 INTEGER(I4B) :: k
 arth_i4b(1)=first
 if(n>1)then
  do k=2,n
   arth_i4b(k) = arth_i4b(k-1) + increment
  end do
 end if
 END FUNCTION arth_i4b
 ! ------------------------------------------------------------------------------------------------
 FUNCTION arth_i8b(first,increment,n)
 implicit none
 INTEGER(I8B), INTENT(IN) :: first,increment,n
 INTEGER(I8B), DIMENSION(n) :: arth_i8b
 INTEGER(I8B) :: k
 arth_i8b(1)=first
 if(n>1)then
  do k=2,n
   arth_i8b(k) = arth_i8b(k-1) + increment
  end do
 end if
 END FUNCTION arth_i8b

 ! *************************************************************************************************
 ! * sort function, used to sort numbers in ascending order
 ! *************************************************************************************************
 SUBROUTINE indexx_i4b(arr,index)
 IMPLICIT NONE
 INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
 INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
 INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
 INTEGER(I4B) :: a
 INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
 INTEGER(I4B), DIMENSION(NSTACK) :: istack
 n=size(arr)
 index=arth(1,1,n)
 jstack=0
 l=1
 r=n
 do
     if (r-l < NN) then
         do j=l+1,r
             indext=index(j)
             a=arr(indext)
             do i=j-1,1,-1
                 if (arr(index(i)) <= a) exit
                 index(i+1)=index(i)
             end do
             index(i+1)=indext
         end do
         if (jstack == 0) RETURN
         r=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
     else
         k=(l+r)/2
         call swap(index(k),index(l+1))
         call icomp_xchg(index(l),index(r))
         call icomp_xchg(index(l+1),index(r))
         call icomp_xchg(index(l),index(l+1))
         i=l+1
         j=r
         indext=index(l+1)
         a=arr(indext)
         do
             do
                 i=i+1
                 if (arr(index(i)) >= a) exit
             end do
             do
                 j=j-1
                 if (arr(index(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(index(i),index(j))
         end do
         index(l+1)=index(j)
         index(j)=indext
         jstack=jstack+2
         if (r-i+1 >= j-l) then
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
         else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
         end if
     end if
 end do
 CONTAINS
 ! internal subroutine
 SUBROUTINE icomp_xchg(i,j)
 INTEGER(I4B), INTENT(INOUT) :: i,j
 INTEGER(I4B) :: swp
 if (arr(j) < arr(i)) then
     swp=i
     i=j
     j=swp
 end if
 END SUBROUTINE icomp_xchg
 END SUBROUTINE indexx_i4b
 ! ------------------------------------------------------------------------------------------------
 SUBROUTINE indexx_i8b(arr,index)
 IMPLICIT NONE
 INTEGER(I8B), DIMENSION(:), INTENT(IN) :: arr
 INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
 INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
 INTEGER(I8B) :: a
 INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
 INTEGER(I4B), DIMENSION(NSTACK) :: istack
 n=size(arr)
 index=arth(1,1,n)
 jstack=0
 l=1
 r=n
 do
     if (r-l < NN) then
         do j=l+1,r
             indext=index(j)
             a=arr(indext)
             do i=j-1,1,-1
                 if (arr(index(i)) <= a) exit
                 index(i+1)=index(i)
             end do
             index(i+1)=indext
         end do
         if (jstack == 0) RETURN
         r=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
     else
         k=(l+r)/2
         call swap(index(k),index(l+1))
         call icomp_xchg(index(l),index(r))
         call icomp_xchg(index(l+1),index(r))
         call icomp_xchg(index(l),index(l+1))
         i=l+1
         j=r
         indext=index(l+1)
         a=arr(indext)
         do
             do
                 i=i+1
                 if (arr(index(i)) >= a) exit
             end do
             do
                 j=j-1
                 if (arr(index(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(index(i),index(j))
         end do
         index(l+1)=index(j)
         index(j)=indext
         jstack=jstack+2
         if (r-i+1 >= j-l) then
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
         else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
         end if
     end if
 end do
 CONTAINS
 ! internal subroutine
 SUBROUTINE icomp_xchg(i,j)
 INTEGER(I4B), INTENT(INOUT) :: i,j
 INTEGER(I4B) :: swp
 if (arr(j) < arr(i)) then
     swp=i
     i=j
     j=swp
 end if
 END SUBROUTINE icomp_xchg
 END SUBROUTINE indexx_i8b

 ! ************************************************************************************************
 ! private subroutine
 ! ************************************************************************************************
 SUBROUTINE swap_i4b(a,b)
 INTEGER(I4B), INTENT(INOUT) :: a,b
 INTEGER(I4B) :: dum
 dum=a
 a=b
 b=dum
 END SUBROUTINE swap_i4b
 ! ------------------------------------------------------------------------------------------------
 SUBROUTINE swap_i8b(a,b)
 INTEGER(I8B), INTENT(INOUT) :: a,b
 INTEGER(I8B) :: dum
 dum=a
 a=b
 b=dum
 END SUBROUTINE swap_i8b

 ! ************************************************************************************************
 ! findIndex: find the first index within a vector
 ! ************************************************************************************************
 function findIndex(vector,desiredValue,missingValue)
   ! NOTE: if the index does not exist, returns zero
   !        workaround for (not-yet-implemented) f2008 intrinsic findloc
   implicit none
   ! argument variables
   integer(i4b),intent(in)            :: vector(:)    ! vector to search
   integer(i4b),intent(in)            :: desiredValue ! desired value in the vector
   integer(i4b),intent(in),optional   :: missingValue ! desired missing value if desiredValue is not found
   integer(i4b)                       :: findIndex    ! first index of the desired value in the vector
   ! local variables
   integer(i4b),dimension(1)          :: vecIndex     ! first index of the desired value in the vector (vec of length=1)

   ! check if the value exisits
   if(any(vector==desiredValue))then
     ! get the index: merge provides a vector with 1s where mask is true and 0s otherwise, so maxloc(merge) is the first index of value=1
     vecIndex=maxloc( merge(1, 0, vector==desiredValue) )
   else
     if(present(missingValue))then
       vecIndex=missingValue
     else
       vecIndex=0
     endif
   endif
   findIndex=vecIndex(1)
 end function findIndex

 ! *************************************************************************************************
 ! Return indices of True in TF array
 ! *************************************************************************************************
 subroutine indexTrue(TF,pos)
   implicit none
   ! argument variables
   logical(lgt),intent(in)                :: TF(:)           ! Logical vector (True or False)
   integer(i4b), allocatable, intent(out) :: pos(:)          ! position of "true" conditions
   ! Local variables
   integer(i4b)                           :: npos            ! number of "true" conditions
   integer(i4b)                           :: idx(size(TF))   ! vector of all positions

   idx = arth(1,1,size(TF))           ! Enumerate all positions
   npos  = count(TF)                  ! Count the elements of TF that are .True.
   allocate(pos(npos))
   pos = pack(idx,TF)                 ! With Pack function, verify position of true conditions
 end subroutine indexTrue

 ! *************************************************************************************************
 ! Find unique elements and indices given array
 ! *************************************************************************************************
 SUBROUTINE unique_i4b(array, unq, idx)
   implicit none
   ! argument variables
   integer(i4b),            intent(in)  :: array(:)             ! integer array including duplicated elements
   integer(i4b),allocatable,intent(out) :: unq(:)               ! integer array including unique elements
   integer(i4b),allocatable,intent(out) :: idx(:)               ! integer array including unique element index
   ! local variables
   integer(i4b)                         :: ranked(size(array))  !
   integer(i4b)                         :: unq_tmp(size(array)) !
   logical(lgt)                         :: flg_tmp(size(array)) !
   integer(i4b)                         :: ix                   ! loop index, counter
   integer(i4b)                         :: last_unique          ! last unique element

   flg_tmp = .false.
   call indexx(array, ranked)

   unq_tmp(ranked(1)) = array(ranked(1))
   flg_tmp(ranked(1)) = .true.
   last_unique = array(ranked(1))
   do ix = 2,size(ranked)
     if (last_unique==array(ranked(ix))) cycle
     flg_tmp(ranked(ix)) = .true.
     unq_tmp(ranked(ix)) = array(ranked(ix))
     last_unique = array(ranked(ix))
   end do

   allocate(unq(count(flg_tmp)),idx(count(flg_tmp)))
   idx = pack(arth(1,1,size(array)), flg_tmp)
   unq = unq_tmp(idx)
 END SUBROUTINE unique_i4b
 ! ------------------------------------------------------------------------------------------------
 SUBROUTINE unique_i8b(array, unq, idx)
  implicit none
  ! argument variables
  integer(i8b),            intent(in)  :: array(:)             ! integer array including duplicated elements
  integer(i8b),allocatable,intent(out) :: unq(:)               ! integer array including unique elements
  integer(i4b),allocatable,intent(out) :: idx(:)               ! integer array including unique element index
  ! local variables
  integer(i4b)                         :: ranked(size(array))  !
  integer(i8b)                         :: unq_tmp(size(array)) !
  logical(lgt)                         :: flg_tmp(size(array)) !
  integer(i4b)                         :: ix                   ! loop index, counter
  integer(i8b)                         :: last_unique          ! last unique element

  flg_tmp = .false.
  call indexx(array, ranked)

  unq_tmp(ranked(1)) = array(ranked(1))
  flg_tmp(ranked(1)) = .true.
  last_unique = array(ranked(1))
  do ix = 2,size(ranked)
    if (last_unique==array(ranked(ix))) cycle
    flg_tmp(ranked(ix)) = .true.
    unq_tmp(ranked(ix)) = array(ranked(ix))
    last_unique = array(ranked(ix))
  end do

  allocate(unq(count(flg_tmp)),idx(count(flg_tmp)))

  idx = pack(arth(1,1,size(array)), flg_tmp)
  unq = unq_tmp(idx)

 END SUBROUTINE unique_i8b

 ! *************************************************************************************************
 ! * convert integer to digit array
 ! *************************************************************************************************
 FUNCTION get_digits(num) result(digs)
   implicit none
   ! argument variables
   integer(i4b), intent(in)  :: num
   ! local variables
   integer(i4b), allocatable :: digs(:)
   integer(i4b)              :: num_digits, ix, rem
   if (num==0) then
     allocate(digs(1))
     digs=0
   else
     num_digits = floor(log10(real(num))+1)
     allocate(digs(num_digits))
     rem = num
     do ix = 1, num_digits
        digs(num_digits-ix+1) = rem - (rem/10)*10  ! Take advantage of integer division
        rem = rem/10
     end do
   end if
 END FUNCTION get_digits

 ! *************************************************************************************************
 ! character-to-integer routines/functions
 ! *************************************************************************************************
 SUBROUTINE char2int_1d(char_array, int_array, invalid_value)
   ! Convert character array to one digit integer array
   ! if character array is '-9999' or '0', int_array(:) = [invalid_value,...,invalid_value]
   implicit none
   ! argument variables
   character(*),              intent(in)   :: char_array
   integer(i4b), allocatable, intent(out)  :: int_array(:)
   integer(i4b), optional,    intent(in)   :: invalid_value
   ! local variables
   character(len=strLen)                   :: string
   integer(i4b)                            :: str_len
   integer(i4b)                            :: iChr
   integer(i4b)                            :: invalValue

   if (present(invalid_value)) then
     invalValue = invalid_value
   else
     invalValue = -1
   endif

   allocate(int_array(len(char_array)))
   int_array = invalValue

   string = adjustl(char_array)
   str_len = len(trim(string))

   if (trim(string)=='-9999' .or. trim(string)=='0') then
     int_array(1) = invalValue
   else
     do iChr =1,str_len
       read(string(iChr:iChr),'(I1)') int_array(iChr)
     end do
   endif
 END SUBROUTINE char2int_1d

 SUBROUTINE char2int_2d(char_array, int_array, invalid_value)
   ! convert character array to one digit integer array
   ! if character array is '-9999' or '0', int_array(:) = [invalid_value,...,invalid_value]
   implicit none
   ! Argument variables
   character(*),              intent(in)   :: char_array(:)
   integer(i4b), allocatable, intent(out)  :: int_array(:,:)
   integer(i4b), optional,    intent(in)   :: invalid_value
   ! local variables
   character(len=strLen)                   :: string
   integer(i4b)                            :: str_len
   integer(i4b)                            :: iSize
   integer(i4b)                            :: iChr
   integer(i4b)                            :: invalValue

  if (present(invalid_value)) then
    invalValue = invalid_value
  else
    invalValue = -1
  endif

   allocate(int_array(size(char_array),len(char_array)))
   int_array = invalValue

   do iSize = 1, size(char_array)
     str_len = len(trim(adjustl(char_array(iSize))))
     string = adjustl(char_array(iSize))
     if (trim(string) == '-9999' .or. trim(string) == '0') then
       int_array(iSize, 1) = invalValue
     else
       do iChr =1,str_len
         read(string(iChr:iChr),'(I1)') int_array(iSize, iChr)
       end do
     end if
   end do
 END SUBROUTINE char2int_2d

  ! ************************************************************************************************
  ! match_index: find array1 indix for each array2 element if array2 includes matching element in array1
  ! ************************************************************************************************
  FUNCTION match_index(array1, array2, missingValue) RESULT(index1)
    implicit none
    ! Argument variables:
    integer(i4b), allocatable, intent(in)  :: array1(:)
    integer(i4b), allocatable, intent(in)  :: array2(:)
    integer(i4b), optional,    intent(in)  :: missingValue ! desired missing value if desiredValue is not found
    ! Local variables:
    integer(i4b), allocatable              :: index1(:)
    integer(i4b), allocatable              :: rnkArray1(:)
    integer(i4b), allocatable              :: rnkArray2(:)
    integer(i4b)                           :: ix, jx, begIx

    allocate(index1(size(array2)), rnkArray1(size(array1)), rnkArray2(size(array2)) )

    if(present(missingValue))then
      index1=missingValue
    else
      index1 = -9999
    endif

    call indexx(array1, rnkArray1)
    call indexx(array2, rnkArray2)

    begIx=1
    do ix=1,size(rnkArray2)
      do jx=begIx,size(rnkArray1)
        if (array2(rnkArray2(ix))==array1(rnkArray1(jx))) then
          index1(rnkArray2(ix)) = rnkArray1(jx)
          begIx=jx
          exit
        else if (array2(rnkArray2(ix))<array1(rnkArray1(jx))) then
          begIx=jx
          exit
        end if
      end do
    end do
  END FUNCTION

END MODULE nr_utility_module
