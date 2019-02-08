MODULE mpi_routine_module

USE nrtype
USE public_var
USE mpi
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
USE var_lookup,        ONLY: ixNTOPO       ! index of variables for the netowork topolgy
USE globalData,        ONLY: domains       ! domain data structure - for each domain, pfaf codes and list of segment indices
USE globalData,        ONLY: nDomain       ! count of decomposed domains (tributaries + mainstems)
USE globalData,        ONLY: node_id       ! node id for each domain

implicit none

private

public :: commun_reach_param

contains

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
  integer(i4b)                                :: segIndex(nSeg)             ! reach index for all the segments
  integer(i4b)                                :: downIndex(nSeg)            ! downstream reach index for all the segments
  integer(i4b), allocatable                   :: segId_local(:)             ! reach id for decomposed network
  integer(i4b), allocatable                   :: segIndex_local(:)          ! reach index for decomposed network
  integer(i4b), allocatable                   :: downIndex_local(:)         ! downstream segment index for decomposed network
  integer(i4b)                                :: elem_per_proc(0:nNodes-1)  ! number of segments assigned to each proc (i.e., node)
  integer(i4b)                                :: num_elems_to_received
  integer(i4b)                                :: start_ix(nNodes)           ! start indices for each node assignment
  integer(i4b)                                :: end_idx(nNodes)            ! end indices for each node assignment
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
    ! SegIndex Array needs to be contiguous when a chunk is sent to computing node
    ix2 = 0
    do ix=0,nNodes-1
     do ixx = 1,nDomain
      if (node_id(ixx) == ix) then
        ix1 = ix2+1
        ix2 = ix1+size(domains(ixx)%segIndex)-1
        segIndex(ix1:ix2) = domains(ixx)%segIndex
       endif
     end do
    end do
    ! covert derived data type to the rest of the arrays
    do iSeg = 1,nSeg
     segId(iSeg)     = structNTOPO(segIndex(iSeg))%var(ixNTOPO%segId)%dat(1)
     downIndex(iSeg) = structNTOPO(segIndex(iSeg))%var(ixNTOPO%downSegIndex)%dat(1)
    end do
    ! Compute the number of elements in each node
    elem_per_proc = 0
    do ix = 1,nDomain
     idx = node_id(ix)
     elem_per_proc(idx) = elem_per_proc(idx) + size(domains(ix)%segIndex)
    end do

    ! prepare for data structure for small tributaries and mainstems (routed at root node: 0)
    ! split segments assigned to root node into mainstems and tributaries
    ix1 = 1
    ix2 = elem_per_proc(1)
    segIndex(ix1:ix2)

    ! send a portion of array to each process
    do myid = 1, nNodes-1
     ! Get starting and ending indices
     end_idx(myid)   = sum(elem_per_proc(0:myid))
     start_ix(myid) = end_idx(myid) - elem_per_proc(myid) + 1

     ! Send number of elements (segments)
     call MPI_SEND(elem_per_proc(myid), 1, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)

     ! Send regular 1D array
     call MPI_SEND(segIndex(start_ix(myid)), elem_per_proc(myid), MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
     call MPI_SEND(segId(start_ix(myid)), elem_per_proc(myid), MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
     call MPI_SEND(downIndex(start_ix(myid)), elem_per_proc(myid), MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)

     ! Send lagged 1D array
!     start_ix(myid) = myid*elem_per_proc + 1
!     end_idx(myid)   = start_ix(myid) + elem_per_proc - 1
!     call MPI_SEND(num_elems_to_send, 1, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)

!     call MPI_SEND(uh(start_ix(myid)), num_elems_to_send, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
!     call MPI_SEND(xxxx(start_ix(myid)), num_elems_to_send, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)

    end do

  else
   ! recieve number of elements to be recieved
   call MPI_RECV(num_elems_to_received, 1, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   ! allocate array local to a node
   allocate(segId_local(num_elems_to_received), segIndex_local(num_elems_to_received), stat=ierr)

   ! recieve a local array and hold it in local arrays
   call MPI_RECV(segId_local, num_elems_to_received, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(segIndex_local, num_elems_to_received, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   print*, pid, size(segId_local)

!   call MPI_RECV(uh_local, num_elems_to_received, MPI_REAL8, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

!   ! convert arrays into data structures
!   call array_to_dtype()

  endif

 end subroutine commun_reach_param


end module mpi_routine_module
