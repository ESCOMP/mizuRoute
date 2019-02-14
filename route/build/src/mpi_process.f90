MODULE mpi_routine_module

USE nrtype
USE public_var
USE mpi
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
USE var_lookup,        ONLY: ixNTOPO       ! index of variables for the netowork topolgy
USE globalData,        ONLY: domains       ! domain data structure - for each domain, pfaf codes and list of segment indices
USE globalData,        ONLY: nDomain       ! count of decomposed domains (tributaries + mainstems)
USE nr_utility_module, ONLY: arth          ! Num. Recipies utilities
USE nr_utility_module, ONLY: indexx        ! Num. Recipies utilities

implicit none

private

public :: commun_reach_param

contains

  ! create mapping array to map global index to local index
  !  ixSub(nSeg) : index that belong to a subbasin in global index
  !  iySub(nSeg) : local indices for a subbasin
  !  ixMap(nSeg) : map indix from global to local
  !
  !  examples
  !  global_ix, iySub, ixSub,  ixMap, idNode
  !    1            1      2       1      0
  !    2            2      7       1      0
  !    3            1      5       2      0
  !    4            2      8       3      0
  !    5            3      9       1      0
  !    6            1      1       1      1
  !    7            2      3       2      1
  !    8            3      4       2      1
  !    9            1      6       3      2
  !   10            2     10       2      2
  !
  !  if global index of interest is 4(=ixGlob), ixMap(ixGlob)=3 meaning element with 4 is located in 3rd in local array.


 subroutine commun_reach_param(nSeg, pid, nNodes, structNTOPO, ierr, message)
  ! 1. convert structNTOPO(:)%segIndex into segIndex(:) array
  !    segIndex(:) is ordered in the node assignment from 0 to nNodes-1
  ! 2.
  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: nSeg                       ! number of total segments
  integer(i4b),                   intent(in)  :: pid                        ! process id (MPI)
  integer(i4b),                   intent(in)  :: nNodes                     ! number of processes (MPI)
  type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)             ! network topology
  ! Output variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                    ! error message
  ! Local variables
  character(len=strLen)                       :: cmessage                   ! error message from subroutine
  integer(i4b)                                :: segId(nSeg)                ! reach id for all the segments
  integer(i4b)                                :: ixSub(nSeg)                ! global index in the order of domains
  integer(i4b)                                :: iySub(nSeg)                ! local reach index
  integer(i4b)                                :: ixMap(nSeg)                ! map global index to local
  integer(i4b)                                :: downIndex(nSeg)            ! downstream reach index for all the segments
  integer(i4b)                                :: ixNode(nSeg)               !
  character(len=32)                           :: pfaf(nSeg)                 !
  integer(i4b)                                :: idNode(nDomain)            ! node id array for each domain
  integer(i4b)                                :: rnkIdNode(nDomain)         ! ranked node id array for each domain
  logical(lgt)                                :: isTrib(nSeg)               ! logical to indicate tributary seg or not
  integer(i4b), allocatable                   :: segId_local(:)             ! reach id for decomposed network
  integer(i4b), allocatable                   :: segIndex_local(:)          ! reach index for decomposed network
  integer(i4b), allocatable                   :: downIndex_local(:)         ! downstream segment index for decomposed network
  integer(i4b)                                :: elem_per_proc(0:nNodes-1)  ! number of segments assigned to each proc (i.e., node)
  integer(i4b)                                :: nSubSeg
  integer(i4b)                                :: num_elems_to_received
  integer(i4b)                                :: start_ix                   ! start indices for each node assignment
  integer(i4b)                                :: end_ix                     ! end indices for each node assignment
  integer(i4b)                                :: iSeg                       ! segment loop indices
  integer(i4b)                                :: myid                       ! process id indices
  integer(i4b)                                :: ix,ixx                     ! loop indices
  integer(i4b)                                :: ix1,ix2                    ! starting index and ending index, respectively
  integer(i4b)                                :: idx                        ! node indix (1, ... , nNodes)
  integer(i4b), parameter                     :: root=0                     ! root node id
  integer(i4b), parameter                     :: send_data_tag=2001
  integer(i4b), parameter                     :: return_data_tag=2002
  integer(i4b)                                :: status(MPI_STATUS_SIZE)

  ierr=0; message='commun_reach_param/'

  if (pid == root) then ! this is a root process
    ! Create segIndex array from domains derived type. The array is sorted from node 0 through nNodes-1
    ! SegIndex Array needs to be contiguous when a chunk is sent to computing node (use sort function...)
    ! start with mainstem domain assigned to root node

    forall(ix=1:nDomain) idNode(ix) = domains(ix)%idNode
    call indexx(idNode,rnkIdNode)
    ix2 = 0
    do ix = 1, nDomain
     ixx = rnkIdNode(ix)
     nSubSeg = size(domains(ixx)%segIndex)
     ix1 = ix2+1
     ix2 = ix1+nSubSeg-1
     ixSub(ix1:ix2)  = domains(ixx)%segIndex(1:nSubSeg)
     iySub(ix1:ix2)  = arth(1,1,nSubSeg)
     isTrib(ix1:ix2) = domains(ixx)%isTrib

     ixNode(ix1:ix2) = domains(ixx)%idNode
     pfaf(ix1:ix2)   = adjustl(trim(domains(ixx)%pfaf))

    end do

    ! create index mapping from global to local
    do iSeg = 1,nSeg
      ixMap(ixSub(iSeg)) = iySub(iSeg)
    enddo

    ! covert component of derived data type to the arrays
    downIndex = -1
    do iSeg = 1,nSeg
     segId(iSeg)      = structNTOPO(ixSub(iSeg))%var(ixNTOPO%segId)%dat(1)
     if (structNTOPO(ixSub(iSeg))%var(ixNTOPO%downSegIndex)%dat(1) > 0)then
       downIndex(iSeg)  = ixMap(structNTOPO(ixSub(iSeg))%var(ixNTOPO%downSegIndex)%dat(1))
     endif
    end do

    ! Compute the number of elements in each node
    elem_per_proc = 0
    do ix = 1,nDomain
     idx = domains(ix)%idNode
     elem_per_proc(idx) = elem_per_proc(idx) + size(domains(ix)%segIndex)
    end do

    ! send a portion of array to each process
    do myid = 1, nNodes-1
     ! Get starting and ending indices
     end_ix   = sum(elem_per_proc(0:myid))
     start_ix = end_ix - elem_per_proc(myid) + 1

     ! Send number of elements (segments)
     call MPI_SEND(elem_per_proc(myid), 1, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)

     ! Send regular 1D array
     call MPI_SEND(ixSub(start_ix), elem_per_proc(myid), MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
     call MPI_SEND(segId(start_ix), elem_per_proc(myid), MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
     call MPI_SEND(downIndex(start_ix), elem_per_proc(myid), MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)

     ! Send lagged 1D array
!     start_ix = myid*elem_per_proc + 1
!     end_ix   = start_ix(myid) + elem_per_proc - 1
!     call MPI_SEND(num_elems_to_send, 1, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)

!     call MPI_SEND(uh(start_ix(myid)), num_elems_to_send, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
!     call MPI_SEND(xxxx(start_ix(myid)), num_elems_to_send, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)

    end do

!    do ix = 1,nSeg
!      print*, segId(ix), ixSub(ix), iySub(ix), ixNode(ix), pfaf(ix)
!    enddo

  else
   ! recieve number of elements to be recieved
   call MPI_RECV(num_elems_to_received, 1, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   ! allocate array local to a node
   allocate(segId_local(num_elems_to_received), segIndex_local(num_elems_to_received), stat=ierr)

   ! recieve a local array and hold it in local arrays
   call MPI_RECV(segId_local, num_elems_to_received, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(segIndex_local, num_elems_to_received, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

!   call MPI_RECV(uh_local, num_elems_to_received, MPI_REAL8, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

!   ! convert arrays into data structures
!   call array_to_dtype()

  endif

 end subroutine commun_reach_param


end module mpi_routine_module

!    print*, 'node=0 : mainstem domain'
!    ix2 = 0
!    do ixx = 1,nDomain
!     if (domains(ixx)%idNode==0 .and. domains(ixx)%pfaf(1:1)=='-') then
!      ix1 = ix2+1
!      ix2 = ix1+size(domains(ixx)%segIndex)-1
!      segIndex(ix1:ix2) = domains(ixx)%segIndex
!      print*,domains(ixx)%pfaf, size(domains(ixx)%segIndex)
!     endif
!    end do
!    ! second, small tributary domain assigned to root node
!    print*, 'node=0 : small tributary domain'
!    do ixx = 1,nDomain
!     if (domains(ixx)%idNode==0 .and. domains(ixx)%pfaf(1:1)/='-') then
!      ix1 = ix2+1
!      ix2 = ix1+size(domains(ixx)%segIndex)-1
!      segIndex(ix1:ix2) = domains(ixx)%segIndex
!      print*,domains(ixx)%pfaf, size(domains(ixx)%segIndex)
!     endif
!    end do
!    ! finally large tributary domain assigned to computing node
!    do ix=1,nNodes-1
!     print*, 'node= ', ix
!     do ixx = 1,nDomain
!      if (domains(ixx)%idNode == ix) then
!        ix1 = ix2+1
!        ix2 = ix1+size(domains(ixx)%segIndex)-1
!        segIndex(ix1:ix2) = domains(ixx)%segIndex
!        print*,domains(ixx)%pfaf, size(domains(ixx)%segIndex)
!       endif
!     end do
!    end do

