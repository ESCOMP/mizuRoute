module nr_utility_module
USE nrtype
! contains functions that should really be part of the fortran standard, but are not
implicit none
INTERFACE arth
 MODULE PROCEDURE arth_r, arth_d, arth_i
END INTERFACE
! (everything private unless otherwise specifed)
private
public::arth
public::indexx
public::findIndex
public::indexTrue
contains

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
 FUNCTION arth_i(first,increment,n)
 implicit none
 INTEGER(I4B), INTENT(IN) :: first,increment,n
 INTEGER(I4B), DIMENSION(n) :: arth_i
 INTEGER(I4B) :: k
 arth_i(1)=first
 if(n>1)then
  do k=2,n
   arth_i(k) = arth_i(k-1) + increment
  end do
 end if
 END FUNCTION arth_i

 ! *************************************************************************************************
 ! * sort function, used to sort numbers in ascending order
 ! *************************************************************************************************
 SUBROUTINE indexx(arr,index)
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
 END SUBROUTINE indexx

 ! private subroutine
 SUBROUTINE swap(a,b)
 INTEGER(I4B), INTENT(INOUT) :: a,b
 INTEGER(I4B) :: dum
 dum=a
 a=b
 b=dum
 END SUBROUTINE swap

 ! ************************************************************************************************
 ! * findIndex: find the first index within a vector
 ! ************************************************************************************************
 function findIndex(vector,desiredValue,missingValue)
 ! finds the first index within a vector
 !  -- if the index does not exist, returns zero
 ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc
 implicit none
 ! dummy variables
 integer(i4b),intent(in)            :: vector(:)    ! vector to search
 integer(i4b),intent(in)            :: desiredValue ! desired value in the vector
 integer(i4b),intent(in),optional   :: missingValue ! desired missing value if desiredValue is not found
 integer(i4b)                       :: findIndex    ! first index of the desired value in the vector
 ! local variables
 integer(i4b),dimension(1)          :: vecIndex     ! first index of the desired value in the vector (vec of length=1)

 ! check if the value exisits
 if(any(vector==desiredValue))then

  ! get the index: merge provides a vector with 1s where mask is true and 0s otherwise, so maxloc(merge) is the first index of value=1
  ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc
  vecIndex=maxloc( merge(1, 0, vector==desiredValue) )

 ! value does not exist
 else
  if(present(missingValue))then
   vecIndex=missingValue
  else
   vecIndex=0
  endif
 endif

 ! return function value (extract into a scalar)
 findIndex=vecIndex(1)

 end function findIndex

 subroutine indexTrue(TF,pos)
  ! Return indices of True in TF array
  implicit none
  ! Inlet variables
  logical(lgt),intent(in)                :: TF(:)           ! Logical vector (True or False)
  ! Outlet variables
  integer(i4b), allocatable, intent(out) :: pos(:)          ! position of "true" conditions
  ! Local variable
  integer(i4b)                           :: npos            ! number of "true" conditions
  integer(i4b)                           :: idx(size(TF))   ! vector of all positions

  idx = arth(1,1,size(TF))           ! Enumerate all positions
  npos  = count(TF)                  ! Count the elements of TF that are .True.
  allocate(pos(npos))
  pos = pack(idx,TF)                 ! With Pack function, verify position of true conditions

 end subroutine indexTrue

end module nr_utility_module
