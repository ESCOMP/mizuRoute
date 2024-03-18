MODULE nr_utils

USE nrtype
USE public_var

implicit none

INTERFACE arth
  module procedure arth_r
  module procedure arth_d
  module procedure arth_i4b
END INTERFACE

INTERFACE indexx
  module procedure indexx_i4b
  module procedure indexx_i8b
END INTERFACE

INTERFACE unique
  module procedure unique_i4b
  module procedure unique_i8b
END INTERFACE

INTERFACE sizeo
  module procedure sizeo_i4b
  module procedure sizeo_dp
  module procedure sizeo_sp
END INTERFACE

INTERFACE char2int
  module procedure :: char2int_1d
  module procedure :: char2int_2d
END INTERFACE

INTERFACE match_index
  module procedure :: match_index_i4b, match_index_i8b
END INTERFACE

private
public::arth
public::indexx
public::findIndex
public::match_index
public::indexTrue
public::unique
public::sizeo
public::get_digits
public::char2int

CONTAINS

 ! *************************************************************************************************
 ! * the arth function, used to build a vector of regularly spaced numbers
 ! *************************************************************************************************
 pure FUNCTION arth_r(first,increment,n)
 implicit none
 real(sp), intent(in) :: first,increment
 integer(i4b), intent(in) :: n
 real(sp), dimension(n) :: arth_r
 integer(i4b) :: k
 arth_r(1)=first
 if(n>1)then
  do k=2,n
   arth_r(k) = arth_r(k-1) + increment
  end do
 end if
 END FUNCTION arth_r
 ! ------------------------------------------------------------------------------------------------
 pure FUNCTION arth_d(first,increment,n)
 implicit none
 real(dp), intent(in) :: first,increment
 integer(i4b), intent(in) :: n
 real(dp), dimension(n) :: arth_d
 integer(i4b) :: k
 arth_d(1)=first
 if(n>1)then
  do k=2,n
   arth_d(k) = arth_d(k-1) + increment
  end do
 end if
 END FUNCTION arth_d
 ! ------------------------------------------------------------------------------------------------
 pure FUNCTION arth_i4b(first,increment,n)
 implicit none
 integer(i4b), intent(in) :: first,increment,n
 integer(i4b), dimension(n) :: arth_i4b
 integer(i4b) :: k
 arth_i4b(1)=first
 if(n>1)then
  do k=2,n
   arth_i4b(k) = arth_i4b(k-1) + increment
  end do
 end if
 END FUNCTION arth_i4b

 ! *************************************************************************************************
 ! * sort function, used to sort numbers in ascending order
 ! *************************************************************************************************
 pure SUBROUTINE indexx_i4b(arr,index)
 implicit none
 integer(i4b), dimension(:), intent(IN) :: arr
 integer(i4b), dimension(:), intent(OUT) :: index
 integer(i4b), parameter :: NN=15, NSTACK=50
 integer(i4b) :: a
 integer(i4b) :: n,k,i,j,indext,jstack,l,r
 integer(i4b), dimension(NSTACK) :: istack
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
         if ( k /= l+1 ) call swap(index(k),index(l+1))
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
             if ( i /= j ) call swap(index(i),index(j))
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
   pure SUBROUTINE icomp_xchg(i,j)
   integer(i4b), intent(inout) :: i,j
   integer(i4b) :: swp
   if (arr(j) < arr(i)) then
       swp=i
       i=j
       j=swp
   end if
   END SUBROUTINE icomp_xchg

 END SUBROUTINE indexx_i4b

 pure SUBROUTINE indexx_i8b(arr,index)
 implicit none
 integer(i8b), dimension(:), intent(in) :: arr
 integer(i4b), dimension(:), intent(out) :: index
 integer(i4b), PARAMETER :: NN=15, NSTACK=50
 integer(i8b) :: a
 integer(i4b) :: n,k,i,j,indext,jstack,l,r
 integer(i4b), dimension(NSTACK) :: istack
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
         if ( k /= l+1 ) call swap(index(k),index(l+1))
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
             if ( i /= j ) call swap(index(i),index(j))
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
   pure SUBROUTINE icomp_xchg(i,j)
   integer(i4b), intent(inout) :: i,j
   integer(i4b) :: swp
   if (arr(j) < arr(i)) then
       swp=i
       i=j
       j=swp
   end if
   END SUBROUTINE icomp_xchg
 END SUBROUTINE indexx_i8b

 ! private subroutine
 pure SUBROUTINE swap(a,b)
 integer(i4b), intent(inout) :: a,b
 integer(i4b) :: dum
 dum=a
 a=b
 b=dum
 END SUBROUTINE swap

 ! ************************************************************************************************
 ! findIndex: find the first index within a vector
 ! ************************************************************************************************
 pure function findIndex(vector,desiredValue,missingValue)
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
 pure subroutine indexTrue(TF,pos)
   implicit none
   ! argument variables
   logical(lgt),intent(in)                :: TF(:)           ! Logical vector (True or False)
   integer(i4b), allocatable, intent(out) :: pos(:)          ! position of "true" conditions
   ! Local variable
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
 pure SUBROUTINE unique_i4b(array, unq, idx)
   implicit none
   ! argument variables
   integer(i4b),            intent(in)  :: array(:)             ! integer array including duplicated elements
   integer(i4b),allocatable,intent(out) :: unq(:)               ! integer array including unique elements
   integer(i4b),allocatable,intent(out) :: idx(:)               ! integer array including unique element index
   ! local
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

 ! ---------------------------------
 pure SUBROUTINE unique_i8b(array, unq, idx)
   implicit none
   ! argument variables
   integer(i8b),            intent(in)  :: array(:)             ! integer array including duplicated elements
   integer(i8b),allocatable,intent(out) :: unq(:)               ! integer array including unique elements
   integer(i4b),allocatable,intent(out) :: idx(:)               ! integer array including unique element index
   ! local
   integer(i4b)                         :: ranked(size(array))  !
   integer(i4b)                         :: unq_tmp(size(array)) !
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
 ! * size of array, if not allocated, return zero
 ! *************************************************************************************************
 pure FUNCTION sizeo_i4b(var) RESULT(asize)
   implicit none
   integer(i4b), allocatable, intent(in) :: var(:)
   integer(i4b)                          :: asize
   if (.not.(allocated(var))) then
     asize = 0
   else
     asize = size(var)
   end if
 END FUNCTION sizeo_i4b
 ! ------------------------------------------------------------------------------------------------
 pure FUNCTION sizeo_dp(var) RESULT(asize)
   implicit none
   real(dp), allocatable, intent(in) :: var(:)
   integer(i4b)                      :: asize
   if (.not.(allocated(var))) then
     asize = 0
   else
     asize = size(var)
   end if
 END FUNCTION sizeo_dp
 ! ------------------------------------------------------------------------------------------------
 pure FUNCTION sizeo_sp(var) RESULT(asize)
   implicit none
   real(sp), allocatable, intent(in) :: var(:)
   integer(i4b)                      :: asize
   if (.not.(allocated(var))) then
     asize = 0
   else
     asize = size(var)
   end if
 END FUNCTION sizeo_sp

 ! *************************************************************************************************
 ! * convert integer to digit array
 ! *************************************************************************************************
 pure FUNCTION get_digits(num) result(digs)
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
  FUNCTION match_index_i4b(array1, array2, ierr, message) RESULT(index1)
    implicit none
    ! Argument variables:
    integer(i4b),              intent(in)  :: array1(:)
    integer(i4b),              intent(in)  :: array2(:)
    integer(i4b),              intent(out) :: ierr         ! error code
    character(*),              intent(out) :: message      ! error message
    ! Local variables:
    integer(i4b)                           :: index1(size(array2))
    integer(i4b)                           :: rnkArray1(size(array1))
    integer(i4b)                           :: rnkArray2(size(array2))
    integer(i4b)                           :: ix, jx, begIx

    ierr=0; message='match_index/'

    index1 = integerMissing

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

    ! check
    do ix=1,size(array2)
      if(index1(ix) == integerMissing) cycle
      if(array2(ix) /= array1( index1(ix) ) )then
        write(iulog,'(a,2(1x,I12,1x,I15))') 'ERROR Mapping: ix, ID(ix), index(ix), masterID(index(ix))=', ix, array2(ix), index1(ix), array1(index1(ix))
        message=trim(message)//'unable to find the match'
        ierr=20; return
      endif
    end do

  END FUNCTION match_index_i4b

  FUNCTION match_index_i8b(array1, array2, ierr, message) RESULT(index1)
    implicit none
    ! Argument variables:
    integer(i8b),              intent(in)  :: array1(:)
    integer(i8b),              intent(in)  :: array2(:)
    integer(i4b),              intent(out) :: ierr         ! error code
    character(*),              intent(out) :: message      ! error message
    ! Local variables:
    integer(i4b)                           :: index1(size(array2))
    integer(i4b)                           :: rnkArray1(size(array1))
    integer(i4b)                           :: rnkArray2(size(array2))
    integer(i4b)                           :: ix, jx, begIx

    ierr=0; message='match_index/'

    index1 = integerMissing

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

    ! check
    do ix=1,size(array2)
      if(index1(ix) == integerMissing) cycle
      if(array2(ix) /= array1( index1(ix) ) )then
        write(iulog,'(a,2(1x,I12,1x,I15))') 'ERROR Mapping: ix, ID(ix), index(ix), masterID(index(ix))=', ix, array2(ix), index1(ix), array1(index1(ix))
        message=trim(message)//'unable to find the match'
        ierr=20; return
      endif
    end do

  END FUNCTION match_index_i8b

END MODULE nr_utils
