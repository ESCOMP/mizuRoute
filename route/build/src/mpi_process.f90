MODULE mpi_routine

USE nrtype
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
USE dataTypes,         ONLY: var_dlength   ! double precision type: var(:)%dat
USE dataTypes,         ONLY: var_clength   ! character type:        var(:)%dat

! public data
USE public_var

! named variables
USE var_lookup,        ONLY:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup,        ONLY:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup,        ONLY:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup,        ONLY:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology
USE var_lookup,        ONLY:ixPFAF,   nVarsPFAF    ! index of variables for the pfafstetter code

! utility
USE nr_utility_module, ONLY: arth          ! Num. Recipies utilities
USE nr_utility_module, ONLY: indexx        ! Num. Recipies utilities

! saved/updated global data
USE globalData,        ONLY: domains       ! domain data structure - for each domain, pfaf codes and list of segment indices
USE globalData,        ONLY: nDomain       ! count of decomposed domains (tributaries + mainstems)

USE alloc_data,        ONLY: alloc_struct_scalar

USE popMetadat_module, ONLY: popMetadat    ! populate metadata

USE process_ntopo,     ONLY: augment_ntopo ! compute all the additional network topology (only compute option = on)

USE mpi

implicit none

private

public :: comm_ntopo_data

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


 ! *********************************************************************
 ! public subroutine: send decomposed reach/hru information and populate data structures
 ! *********************************************************************
 subroutine comm_ntopo_data(pid,          & ! input: proc id
                            nNodes,       & ! input: number of procs
                            nSeg,         & ! input: number of stream segments in whole domain
                            nHRU,         & ! input: number of HRUs in whole domain
                            structHRU,    & ! input: ancillary data for HRUs
                            structSEG,    & ! input: ancillary data for stream segments
                            structHRU2seg,& ! input: ancillary data for mapping hru2basin
                            structNTOPO,  & ! input: ancillary data for network toopology
                            structPFAF,   & ! ancillary data for pfafstetter code
                            ierr,message)   ! output: error control

  USE globalData, only: fshape, tscale, velo, diff, mann_n, wscale         ! parameter

  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: pid                      ! process id (MPI)
  integer(i4b),                   intent(in)  :: nNodes                   ! number of processes (MPI)
  integer(i4b),                   intent(in)  :: nSeg                     ! number of total segments
  integer(i4b),                   intent(in)  :: nHRU                     ! number of total hru
  type(var_dlength) , intent(in), allocatable :: structHRU(:)             ! HRU properties
  type(var_dlength) , intent(in), allocatable :: structSeg(:)             ! stream segment properties
  type(var_ilength) , intent(in), allocatable :: structHRU2SEG(:)         ! HRU to SEG mapping
  type(var_ilength) , intent(in), allocatable :: structNTOPO(:)           ! network topology
  type(var_clength) , intent(in), allocatable :: structPFAF(:)            ! pfafstetter code
  ! Output variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                    ! error message
  ! Local variables
  ! data structure for decomposed river network per reach/hru
  type(var_dlength), allocatable              :: structHRU_local(:)         ! ancillary data for HRUs
  type(var_dlength), allocatable              :: structSeg_local(:)         ! ancillary data for stream segments
  type(var_ilength), allocatable              :: structNTOPO_local(:)       ! network topology
  type(var_ilength), allocatable              :: structHRU2seg_local(:)     ! ancillary data for mapping hru2basin
  type(var_clength), allocatable              :: structPFAF_local(:)        ! ancillary data for pfafstetter code
  ! number of data in data strucutes for decomposed river network per reach/hru
  integer(i4b)                                :: tot_upstream_local         ! total number of all of the upstream stream segments for all stream segments
  integer(i4b)                                :: tot_upseg_local            ! total number of immediate upstream segments for all  stream segments
  integer(i4b)                                :: tot_hru_local              ! total number of all the upstream hrus for all stream segments
  integer(i4b)                                :: tot_uh_local               ! total number of unit hydrograph from all the stream segments
  ! indices for decomposed river network per reach/hru
  integer(i4b),      allocatable              :: ixHRU_desired_local(:)     ! indices of desired hrus
  integer(i4b),      allocatable              :: ixSeg_desired_local(:)     ! indices of desired reaches
  ! flat array for decomposed river network per reach/hru
  integer(i4b),      allocatable              :: segId_local(:)             ! reach id for decomposed network
  integer(i4b),      allocatable              :: downSegId_local(:)         ! downstream reach id for decomposed network
  integer(i4b),      allocatable              :: hruId_local(:)             ! hru id array in decomposed network
  integer(i4b),      allocatable              :: hruSegId_local(:)          ! downstream reach id array in decomposed network
  real(dp),          allocatable              :: slope_local(:)             ! reach slope array in decomposed network
  real(dp),          allocatable              :: length_local(:)            ! reach length array in decomposed network
  real(dp),          allocatable              :: area_local(:)              ! hru area in decomposed network
  ! flat array for the entire river network per reach/hru
  integer(i4b)                                :: hruId(nHRU)                ! hru id for all the HRUs
  integer(i4b)                                :: hruSegId(nSeg)             ! hru-to-seg mapping for each hru
  integer(i4b)                                :: segId(nSeg)                ! reach id for all the segments
  integer(i4b)                                :: downIndex(nSeg)            ! downstream reach index for each reach
  integer(i4b)                                :: downSegId(nSeg)            ! downstream reach ID for each reach
  real(dp)                                    :: slope(nSeg)                ! reach slope array for each reach
  real(dp)                                    :: length(nSeg)               ! reach length array for each reach
  real(dp)                                    :: area(nHRU)                 ! hru area for each hru
  integer(i4b)                                :: ixNode(nSeg)               ! node assignment for each reach
  character(len=32)                           :: pfaf(nSeg)                 ! reach pfafcode for each reach
  integer(i4b)                                :: ixSubHRU(nHRU)             ! global HRU index in the order of domains
  integer(i4b)                                :: iySubHRU(nHRU)             ! local HRU index
  integer(i4b)                                :: ixSub(nSeg)                ! global reach index in the order of domains
  integer(i4b)                                :: iySub(nSeg)                ! local reach index
!  integer(i4b)                                :: ixMap(nSeg)                ! map global index to local
  logical(lgt)                                :: isTrib(nSeg)               ! logical to indicate tributary seg or not
  ! flat array for decomposed river network per domain
  integer(i4b)                                :: idNode(nDomain)            ! node id array for each domain
  integer(i4b)                                :: rnkIdNode(nDomain)         ! ranked node id array for each domain
  ! mpi related variables
  integer(i4b)                                :: iSeg,iHru                  ! reach and hru loop indices
  integer(i4b)                                :: ix,ixx                     ! loop indices
  integer(i4b)                                :: seg_per_proc(0:nNodes-1)   ! number of reaches assigned to each proc (i.e., node)
  integer(i4b)                                :: hru_per_proc(0:nNodes-1)   ! number of hrus assigned to each proc (i.e., node)
  integer(i4b)                                :: num_seg_received
  integer(i4b)                                :: num_hru_received
  integer(i4b)                                :: myid                       ! process id indices
  integer(i4b)                                :: ixSeg1,ixSeg2              ! starting index and ending index, respectively, for reach array
  integer(i4b)                                :: ixHru1,ixHru2              ! starting index and ending index, respectively, for HRU array
  integer(i4b)                                :: idx                        ! node indix (1, ... , nNodes)
  integer(i4b), parameter                     :: root=0                     ! root node id
  integer(i4b), parameter                     :: send_data_tag=2001
  integer(i4b), parameter                     :: return_data_tag=2002
  integer(i4b)                                :: status(MPI_STATUS_SIZE)
  character(len=strLen)                       :: cmessage                   ! error message from subroutine

  ierr=0; message='commun_reach_param/'

  if (pid == root) then ! this is a root process

    ! Create segIndex array from domains derived type. The array is sorted from node 0 through nNodes-1
    ! SegIndex Array needs to be contiguous when a chunk is sent to computing node (use sort function...)
    ! start with mainstem domain assigned to root node

    forall(ix=1:nDomain) idNode(ix) = domains(ix)%idNode
    call indexx(idNode,rnkIdNode)

    ixSeg2=0; ixHru2=0
    do ix = 1, nDomain

     ixx = rnkIdNode(ix)
     associate (nSubSeg => size(domains(ixx)%segIndex), nSubHru => size(domains(ixx)%hruIndex) )

     ixSeg1 = ixSeg2+1
     ixSeg2 = ixSeg1+nSubSeg-1
     ixSub(ixSeg1:ixSeg2)  = domains(ixx)%segIndex(1:nSubSeg)   ! global seg index per node
     iySub(ixSeg1:ixSeg2)  = arth(1,1,nSubSeg)                  ! local hru indix per node

     if (nSubHru>0) then
       ixHru1 = ixHru2+1
       ixHru2 = ixHru1+nSubHru-1
       ixSubHRU(ixHru1:ixHru2)  = domains(ixx)%hruIndex(1:nSubHru) ! global hru index per node
       iySubHRU(ixHru1:ixHru2)  = arth(1,1,nSubHru)                ! local hru indix per node
     end if

     isTrib(ixSeg1:ixSeg2) = domains(ixx)%isTrib                 ! if domain is tributary, T otherwise, F
     ixNode(ixSeg1:ixSeg2) = domains(ixx)%idNode                 ! node id
     pfaf(ixSeg1:ixSeg2)   = adjustl(trim(domains(ixx)%pfaf))    ! basin pfaf code

     end associate
    end do

!    ! create index mapping from global to local
!    do iSeg = 1,nSeg
!      ixMap(ixSub(iSeg)) = iySub(iSeg)
!    enddo

    ! covert component of derived data type to the arrays
    ! reach array
!    downIndex = -1
    do iSeg = 1,nSeg
     segId(iSeg)     = structNTOPO(ixSub(iSeg))%var(ixNTOPO%segId)%dat(1)
     downSegId(iSeg) = structNTOPO(ixSub(iSeg))%var(ixNTOPO%downSegId)%dat(1)
     slope(iSeg)     = structSEG(ixSub(iSeg))%var(ixSEG%slope)%dat(1)
     length(iSeg)    = structSEG(ixSub(iSeg))%var(ixSEG%length)%dat(1)
!     if (structNTOPO(ixSub(iSeg))%var(ixNTOPO%downSegIndex)%dat(1) > 0)then
!       downIndex(iSeg)  = ixMap(structNTOPO(ixSub(iSeg))%var(ixNTOPO%downSegIndex)%dat(1))
!     endif
    end do

    ! hru array
    do iHru = 1,nHRU
      hruId(iHru)    = structHRU2SEG(ixSubHRU(iHru))%var(ixHRU2SEG%HRUid)%dat(1)
      hruSegId(iHru) = structHRU2SEG(ixSubHRU(iHru))%var(ixHRU2SEG%hruSegId)%dat(1)
      area(iHru)     = structHRU(ixSubHRU(iHru))%var(ixHRU%area)%dat(1)
    enddo

    ! Compute the number of elements in each node
    seg_per_proc = 0
    hru_per_proc = 0
    do ix = 1,nDomain
     idx = domains(ix)%idNode
     seg_per_proc(idx) = seg_per_proc(idx) + size(domains(ix)%segIndex)
     hru_per_proc(idx) = hru_per_proc(idx) + size(domains(ix)%hruIndex)
    end do

    ! send a portion of array to each process
    do myid = 1, nNodes-1
      ! Get starting and ending indices
      ixSeg2   = sum(seg_per_proc(0:myid))
      ixSeg1 = ixSeg2 - seg_per_proc(myid) + 1

      ixHru2   = sum(hru_per_proc(0:myid))
      ixHru1 = ixHru2 - hru_per_proc(myid) + 1

      ! Send number of elements (segments)
      call MPI_SEND(seg_per_proc(myid), 1, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(hru_per_proc(myid), 1, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      ! Send regular 1D array
      ! reach
      call MPI_SEND(segId(ixSeg1),     seg_per_proc(myid), MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(downSegId(ixSeg1), seg_per_proc(myid), MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(slope(ixSeg1),     seg_per_proc(myid), MPI_REAL8, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(length(ixSeg1),    seg_per_proc(myid), MPI_REAL8, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      ! hru
      call MPI_SEND(hruId(ixHru1),    hru_per_proc(myid), MPI_INT,   myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(hruSegId(ixHru1), hru_per_proc(myid), MPI_INT,   myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(area(ixHru1),     hru_per_proc(myid), MPI_REAL8, myid, send_data_tag, MPI_COMM_WORLD, ierr)
    end do

!    print*, 'ix, hruId, ixSubHRU, iySubHRU, hruSegId'
!    do ix = 1,nHru
!      print*, ix, hruId(ix), ixSubHRU(ix), iySubHRU(ix), hruSegId(ix)
!    end do
!    print*, 'ix, segId, ixSubHRU, iySubHRU, ixNode, pfaf'
!    do ix = 1,nSeg
!      print*, ix, segId(ix), ixSub(ix), iySub(ix), ixNode(ix), pfaf(ix)
!    enddo

  else
   ! recieve number of elements to be recieved
   call MPI_RECV(num_seg_received, 1, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(num_hru_received, 1, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

   ! allocate array local to a node
   allocate(segId_local    (num_seg_received), &
            downSegId_local(num_seg_received), &
            slope_local    (num_seg_received), &
            length_local   (num_seg_received), &
            hruId_local    (num_hru_received), &
            hruSegId_local (num_hru_received), &
            area_local     (num_hru_received), &
            stat=ierr)

   call popMetadat(ierr,cmessage)

   call alloc_struct_scalar(&
                          num_hru_received,      & ! output: number of HRUs
                          num_seg_received,      & ! output: number of stream segments
                          structHRU_local,       & ! inout: ancillary data for HRUs
                          structSeg_local,       & ! inout: ancillary data for stream segments
                          structHRU2seg_local,   & ! inout: ancillary data for mapping hru2basin
                          structNTOPO_local,     & ! inout: ancillary data for network toopology
                          structPFAF_local,      & ! inout: ancillary data for pfafstetter code
                          ierr,cmessage)           ! output: error control

   ! recieve a local array and hold it in local arrays
   ! reach
   call MPI_RECV(segId_local,     num_seg_received, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(downSegId_local, num_seg_received, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(length_local,    num_seg_received, MPI_REAL8, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(slope_local,     num_seg_received, MPI_REAL8, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   ! hru
   call MPI_RECV(hruId_local,     num_hru_received, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(hruSegId_local,  num_hru_received, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(area_local,      num_hru_received, MPI_REAL8, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

   ! Populate data structure
   ! reach
   forall(ix=1:num_seg_received) structNTOPO_local(ix)%var(ixNTOPO%segId)%dat(1)     = segId_local(ix)
   forall(ix=1:num_seg_received) structNTOPO_local(ix)%var(ixNTOPO%downSegId)%dat(1) = downSegId_local(ix)
   forall(ix=1:num_seg_received) structSEG_local  (ix)%var(ixSEG%length)%dat(1)      = length_local(ix)
   forall(ix=1:num_seg_received) structSEG_local  (ix)%var(ixSEG%slope)%dat(1)       = slope_local(ix)
   ! hru
   forall(ix=1:num_hru_received) structHRU2SEG_local(ix)%var(ixHRU2SEG%HRUid)%dat(1)    = hruId_local(ix)
   forall(ix=1:num_hru_received) structHRU2SEG_local(ix)%var(ixHRU2SEG%hruSegId)%dat(1) = hruSegId_local(ix)
   forall(ix=1:num_hru_received) structHRU_local    (ix)%var(ixHRU%area)%dat(1)         = area_local(ix)

   dt = 86400._dp
   fshape=2.5_dp
   tscale=86400._dp

   call augment_ntopo(&
                   ! input: model control
                   num_hru_received,             & ! number of HRUs
                   num_seg_received,             & ! number of stream segments
                   ! inout: populate data structures
                   structHRU_local,              & ! ancillary data for HRUs
                   structSeg_local,              & ! ancillary data for stream segments
                   structHRU2seg_local,          & ! ancillary data for mapping hru2basin
                   structNTOPO_local,            & ! ancillary data for network toopology
                   structPFAF_local,             & ! ancillary data for pfafstetter code
                   ! output
                   tot_hru_local,                & ! total number of all the upstream hrus for all stream segments
                   tot_upseg_local,              & ! total number of all the immediate upstream segments for all stream segments
                   tot_upstream_local,           & ! total number of all the upstream segments for all stream segments
                   tot_uh_local,                 & ! total number of unit hydrograph for all stream segments
                   ixHRU_desired_local,          & ! indices of desired hrus
                   ixSeg_desired_local,          & ! indices of desired reaches
                   ! output: error control
                   ierr, message)

   print*, pid, size(structNTOPO_local)


  endif

 end subroutine comm_ntopo_data


end module mpi_routine

!    print*, 'node=0 : mainstem domain'
!    ixSeg2 = 0
!    do ixx = 1,nDomain
!     if (domains(ixx)%idNode==0 .and. domains(ixx)%pfaf(1:1)=='-') then
!      ixSeg1 = ixSeg2+1
!      ixSeg2 = ixSeg1+size(domains(ixx)%segIndex)-1
!      segIndex(ixSeg1:ixSeg2) = domains(ixx)%segIndex
!      print*,domains(ixx)%pfaf, size(domains(ixx)%segIndex)
!     endif
!    end do
!    ! second, small tributary domain assigned to root node
!    print*, 'node=0 : small tributary domain'
!    do ixx = 1,nDomain
!     if (domains(ixx)%idNode==0 .and. domains(ixx)%pfaf(1:1)/='-') then
!      ixSeg1 = ixSeg2+1
!      ixSeg2 = ixSeg1+size(domains(ixx)%segIndex)-1
!      segIndex(ixSeg1:ixSeg2) = domains(ixx)%segIndex
!      print*,domains(ixx)%pfaf, size(domains(ixx)%segIndex)
!     endif
!    end do
!    ! finally large tributary domain assigned to computing node
!    do ix=1,nNodes-1
!     print*, 'node= ', ix
!     do ixx = 1,nDomain
!      if (domains(ixx)%idNode == ix) then
!        ixSeg1 = ixSeg2+1
!        ixSeg2 = ixSeg1+size(domains(ixx)%segIndex)-1
!        segIndex(ixSeg1:ixSeg2) = domains(ixx)%segIndex
!        print*,domains(ixx)%pfaf, size(domains(ixx)%segIndex)
!       endif
!     end do
!    end do

