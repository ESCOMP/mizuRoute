MODULE mpi_routine

USE mpi

USE nrtype
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
USE dataTypes,         ONLY: var_dlength   ! double precision type: var(:)%dat
USE dataTypes,         ONLY: var_clength   ! character type:        var(:)%dat

! named variables
USE var_lookup,        ONLY:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup,        ONLY:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup,        ONLY:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup,        ONLY:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology
USE var_lookup,        ONLY:ixPFAF,   nVarsPFAF    ! index of variables for the pfafstetter code

! utility
USE nr_utility_module, ONLY: arth                  !
USE nr_utility_module, ONLY: indexx                ! sorted index array
USE nr_utility_module, only: findIndex             ! find index within a vector

implicit none

private

public :: comm_ntopo_data
public :: comm_runoff_data

contains

 ! *********************************************************************
 ! public subroutine: send reach/hru information to tasks and populate data structures
 ! *********************************************************************
 subroutine comm_ntopo_data(pid,          & ! input: proc id
                            nNodes,       & ! input: number of procs
                            nRch_in,      & ! input: number of stream segments in whole domain
                            nHRU_in,      & ! input: number of HRUs that are connected to reaches
                            structHRU,    & ! input: data structure for HRUs
                            structSEG,    & ! input: data structure for stream segments
                            structHRU2seg,& ! input: data structure for mapping hru2basin
                            structNTOPO,  & ! input: data structure for network toopology
                            ixSubHRU,     & ! output: sorted hru index array based on proc assignment
                            ixSubSEG,     & ! output: sorted seg index array based on proc assignment
                            hru_per_proc, & ! number of hrus assigned to each proc (i.e., node)
                            seg_per_proc, & ! number of reaches assigned to each proc (i.e., node)
                            ierr,message)   ! output: error control

  USE popMetadat_module, ONLY: popMetadat           ! populate metadata
  USE alloc_data,        ONLY: alloc_struct
  USE process_ntopo,     ONLY: augment_ntopo        ! compute all the additional network topology (only compute option = on)
  USE process_ntopo,     ONLY: put_data_struct      !
  USE globalData,        ONLY: fshape, tscale       ! routing parameters
  USE globalData,        ONLY: velo, diff           ! routing parameters
  USE globalData,        ONLY: mann_n, wscale       ! routing parameters
  USE globalData,        ONLY: ixDesire             ! desired reach index
  USE globalData,        ONLY: domains              ! domain data structure - for each domain, pfaf codes and list of segment indices
  USE globalData,        ONLY: nDomain              ! count of decomposed domains (tributaries + mainstems)
  USE globalData,        ONLY: RCHFLX, RCHFLX_local
  USE globalData,        ONLY: KROUTE, KROUTE_local
  USE globalData,        ONLY: nEns
  USE public_var

  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: pid                      ! process id (MPI)
  integer(i4b),                   intent(in)  :: nNodes                   ! number of processes (MPI)
  integer(i4b),                   intent(in)  :: nRch_in                  ! number of total segments
  integer(i4b),                   intent(in)  :: nHRU_in                  ! number of total hru
  type(var_dlength), allocatable, intent(in)  :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(in)  :: structSEG(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(in)  :: structHRU2SEG(:)         ! HRU to SEG mapping
  type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)           ! network topology
  ! Output
  integer(i4b),      allocatable, intent(out) :: ixSubHRU(:)              ! global HRU index in the order of domains
  integer(i4b),      allocatable, intent(out) :: ixSubSEG(:)              ! global reach index in the order of domains
  integer(i4b),      allocatable, intent(out) :: hru_per_proc(:)          ! number of hrus assigned to each proc (i.e., node)
  integer(i4b),      allocatable, intent(out) :: seg_per_proc(:)          ! number of reaches assigned to each proc (i.e., node)
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                   ! error message
  ! Local variables
  ! data structure for decomposed river network per reach/hru
  type(var_dlength), allocatable              :: structHRU_local(:)        ! ancillary data for HRUs
  type(var_dlength), allocatable              :: structSEG_local(:)        ! ancillary data for stream segments
  type(var_ilength), allocatable              :: structNTOPO_local(:)      ! network topology
  type(var_ilength), allocatable              :: structHRU2seg_local(:)    ! ancillary data for mapping hru2basin
  type(var_clength), allocatable              :: structPFAF_local(:)       ! ancillary data for pfafstetter code
  ! number of data in data strucutes for decomposed river network per reach/hru
  integer(i4b)                                :: tot_upstream_local        ! total number of all of the upstream stream segments for all stream segments
  integer(i4b)                                :: tot_upseg_local           ! total number of immediate upstream segments for all  stream segments
  integer(i4b)                                :: tot_hru_local             ! total number of all the upstream hrus for all stream segments
  integer(i4b)                                :: tot_uh_local              ! total number of unit hydrograph from all the stream segments
  ! indices for decomposed river network per reach/hru
  integer(i4b),      allocatable              :: ixHRU_desired_local(:)    ! indices of desired hrus
  integer(i4b),      allocatable              :: ixSeg_desired_local(:)    ! indices of desired reaches
  ! flat array for decomposed river network per reach/hru
  integer(i4b),      allocatable              :: segId_local(:)            ! reach id for decomposed network
  integer(i4b),      allocatable              :: downSegId_local(:)        ! downstream reach id for decomposed network
  integer(i4b),      allocatable              :: hruId_local(:)            ! hru id array in decomposed network
  integer(i4b),      allocatable              :: hruSegId_local(:)         ! downstream reach id array in decomposed network
  real(dp),          allocatable              :: slope_local(:)            ! reach slope array in decomposed network
  real(dp),          allocatable              :: length_local(:)           ! reach length array in decomposed network
  real(dp),          allocatable              :: area_local(:)             ! hru area in decomposed network
  ! flat array for the entire river network per reach/hru
  integer(i4b)                                :: hruId(nHRU_in)            ! hru id for all the HRUs
  integer(i4b)                                :: hruSegId(nRch_in)         ! hru-to-seg mapping for each hru
  integer(i4b)                                :: segId(nRch_in)            ! reach id for all the segments
  integer(i4b)                                :: downSegId(nRch_in)        ! downstream reach ID for each reach
  real(dp)                                    :: slope(nRch_in)            ! reach slope array for each reach
  real(dp)                                    :: length(nRch_in)           ! reach length array for each reach
  real(dp)                                    :: area(nHRU_in)             ! hru area for each hru
  integer(i4b)                                :: ixNode(nRch_in)           ! node assignment for each reach
  character(len=32)                           :: pfaf(nRch_in)             ! reach pfafcode for each reach
  integer(i4b)                                :: iySubHRU(nHRU_in)         ! local HRU index
  integer(i4b)                                :: iySubSEG(nRch_in)         ! local reach index
  logical(lgt)                                :: isTrib(nRch_in)           ! logical to indicate tributary seg or not
  ! flat array for decomposed river network per domain (sub-basin)
  integer(i4b)                                :: idNode(nDomain)           ! node id array for each domain
  integer(i4b)                                :: rnkIdNode(nDomain)        ! ranked node id array for each domain
  ! mpi related variables
  integer(i4b)                                :: iSeg,iHru                 ! reach and hru loop indices
  integer(i4b)                                :: ix,ixx                    ! loop indices
  integer(i4b)                                :: num_seg_received
  integer(i4b)                                :: num_hru_received
  integer(i4b)                                :: nHRU_root, nSeg_root      ! number of hrus and reaches in a root proc
  integer(i4b)                                :: myid                      ! process id indices
  integer(i4b)                                :: ixSeg1,ixSeg2             ! starting index and ending index, respectively, for reach array
  integer(i4b)                                :: ixHru1,ixHru2             ! starting index and ending index, respectively, for HRU array
  integer(i4b)                                :: idx                       ! node indix (1, ... , nNodes)
  integer(i4b), parameter                     :: root=0                    ! root node id
  integer(i4b), parameter                     :: send_data_tag=2001
  integer(i4b), parameter                     :: return_data_tag=2002
  integer(i4b)                                :: status(MPI_STATUS_SIZE)
  character(len=strLen)                       :: cmessage                  ! error message from subroutine

  ierr=0; message='commun_ntopo_data/'

  if (pid == root) then ! this is a root process

    allocate(RCHFLX(nEns,nRch_in), KROUTE(nEns,nRch_in), stat=ierr)
    allocate(ixSubHRU(nHRU_in),ixSubSEG(nRch_in), stat=ierr)

    ! Create segIndex array from domains derived type. The array is sorted from node 0 through nNodes-1
    ! SegIndex Array needs to be contiguous when a chunk is sent to computing node (use sort function...)
    ! start with mainstem domain assigned to root node

    forall(ix=1:nDomain) idNode(ix) = domains(ix)%idNode
    call indexx(idNode,rnkIdNode)

    ixSeg2=0; ixHru2=0
    do ix = 1, nDomain

     ixx = rnkIdNode(ix)
     associate (nSubSeg => size(domains(ixx)%segIndex), nSubHru => size(domains(ixx)%hruIndex) )

     ! reach index array in order of node assignment
     ixSeg1 = ixSeg2+1
     ixSeg2 = ixSeg1+nSubSeg-1
     ixSubSEG(ixSeg1:ixSeg2)  = domains(ixx)%segIndex(1:nSubSeg)   ! global seg index per node
     iySubSEG(ixSeg1:ixSeg2)  = arth(1,1,nSubSeg)                  ! local hru indix per node

     ! hru index array in order of node assignment
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

    ! covert component of derived data type to the arrays
    ! reach array
    do iSeg = 1,nRch_in
     segId(iSeg)     = structNTOPO(ixSubSEG(iSeg))%var(ixNTOPO%segId)%dat(1)
     downSegId(iSeg) = structNTOPO(ixSubSEG(iSeg))%var(ixNTOPO%downSegId)%dat(1)
     slope(iSeg)     = structSEG(ixSubSEG(iSeg))%var(ixSEG%slope)%dat(1)
     length(iSeg)    = structSEG(ixSubSEG(iSeg))%var(ixSEG%length)%dat(1)
    end do

    ! hru array
    do iHru = 1,nHRU_in
      hruId(iHru)    = structHRU2SEG(ixSubHRU(iHru))%var(ixHRU2SEG%HRUid)%dat(1)
      hruSegId(iHru) = structHRU2SEG(ixSubHRU(iHru))%var(ixHRU2SEG%hruSegId)%dat(1)
      area(iHru)     = structHRU(ixSubHRU(iHru))%var(ixHRU%area)%dat(1)
    enddo

    ! Compute the number of elements in each node
    allocate(seg_per_proc(0:nNodes), hru_per_proc(0:nNodes), stat=ierr)
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
      call MPI_SEND(length(ixSeg1),    seg_per_proc(myid), MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(slope(ixSeg1),     seg_per_proc(myid), MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      ! hru
      call MPI_SEND(hruId(ixHru1),    hru_per_proc(myid), MPI_INT,   myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(hruSegId(ixHru1), hru_per_proc(myid), MPI_INT,   myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(area(ixHru1),     hru_per_proc(myid), MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
     ! other public/global data
      call MPI_SEND(desireId,1, MPI_INT,   myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(routOpt, 1, MPI_INT,   myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(dt,      1, MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(fshape,  1, MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(tscale,  1, MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(velo,    1, MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(diff,    1, MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(mann_n,  1, MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(wscale,  1, MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)

    end do

    ! reach assigned to root proc
    nSeg_root=seg_per_proc(root)   ! number of seg in root proc
    nHRU_root=hru_per_proc(root)   ! number of hru in root proc

    ! find index of desired reach
    ixDesire = findIndex(segId(1:nSeg_root), desireId, integerMissing)

    call alloc_struct(nHRU_root,             & ! output: number of HRUs
                      nSeg_root,             & ! output: number of stream segments
                      structHRU_local,       & ! inout: ancillary data for HRUs
                      structSEG_local,       & ! inout: ancillary data for stream segments
                      structHRU2seg_local,   & ! inout: ancillary data for mapping hru2basin
                      structNTOPO_local,     & ! inout: ancillary data for network toopology
                      structPFAF_local,      & ! inout: ancillary data for pfafstetter code
                      ierr,cmessage)           ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    allocate(RCHFLX_local(nEns,nSeg_root), KROUTE_local(nEns,nSeg_root), stat=ierr)

    ! Populate data structure
    ! reach
    do ix = 1, nSeg_root
      structNTOPO_local(ix)%var(ixNTOPO%segId)%dat(1)     = segId(ix)
      structNTOPO_local(ix)%var(ixNTOPO%downSegId)%dat(1) = downSegId(ix)
      structSEG_local  (ix)%var(ixSEG%length)%dat(1)      = length(ix)
      structSEG_local  (ix)%var(ixSEG%slope)%dat(1)       = slope(ix)
    end do
    ! hru
    do ix = 1, nHRU_root
      structHRU2SEG_local(ix)%var(ixHRU2SEG%HRUid)%dat(1)    = hruId(ix)
      structHRU2SEG_local(ix)%var(ixHRU2SEG%hruSegId)%dat(1) = hruSegId(ix)
      structHRU_local    (ix)%var(ixHRU%area)%dat(1)         = area(ix)
    end do

    call augment_ntopo(&
                       ! input: model control
                       nHRU_root,                   & ! number of stream segments
                       nSeg_root,                   & ! number of HRUs
                       ! inout: populate data structures
                       structHRU_local,              & ! ancillary data for HRUs
                       structSEG_local,              & ! ancillary data for stream segments
                       structHRU2seg_local,          & ! ancillary data for mapping hru2basin
                       structNTOPO_local,            & ! ancillary data for network toopology
                       ! output
                       tot_hru_local,                & ! total number of all the upstream hrus for all stream segments
                       tot_upseg_local,              & ! total number of all the immediate upstream segments for all stream segments
                       tot_upstream_local,           & ! total number of all the upstream segments for all stream segments
                       tot_uh_local,                 & ! total number of unit hydrograph for all stream segments
                       ixHRU_desired_local,          & ! indices of desired hrus
                       ixSeg_desired_local,          & ! indices of desired reaches
                       ! output: error control
                       ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    print*, 'pid, nHRU_root, nSeg_root, tot_hru_local, tot_upseg_local, tot_upstream_local, tot_uh_local =', pid, nHRU_root, nSeg_root, tot_hru_local, tot_upseg_local, tot_upstream_local, tot_uh_local
    call put_data_struct(nSeg_root, structSEG_local, structNTOPO_local, ierr, cmessage)

!    print*, 'ix, hruId, ixSubHRU, iySubHRU, hruSegId'
!    do ix = 1,nHRU_in
!      print*, ix, hruId(ix), ixSubHRU(ix), iySubHRU(ix), hruSegId(ix)
!    end do
!    print*, 'ix, segId, ixSubSEG, iySubSEG, ixNode, pfaf'
!    do ix = 1,nRch_in
!      print*, ix, segId(ix), ixSubSEG(ix), iySubSEG(ix), ixNode(ix), pfaf(ix)
!    enddo

  else
   call popMetadat(ierr,cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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

   allocate(RCHFLX_local(nEns,num_seg_received), KROUTE_local(nEns,num_seg_received), stat=ierr)

   call alloc_struct(num_hru_received,      & ! output: number of HRUs
                     num_seg_received,      & ! output: number of stream segments
                     structHRU_local,       & ! inout: ancillary data for HRUs
                     structSEG_local,       & ! inout: ancillary data for stream segments
                     structHRU2seg_local,   & ! inout: ancillary data for mapping hru2basin
                     structNTOPO_local,     & ! inout: ancillary data for network toopology
                     structPFAF_local,      & ! inout: ancillary data for pfafstetter code
                     ierr,cmessage)           ! output: error control

   ! recieve a local array and hold it in local arrays
   ! reach
   call MPI_RECV(segId_local,     num_seg_received, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(downSegId_local, num_seg_received, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(length_local,    num_seg_received, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(slope_local,     num_seg_received, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   ! hru
   call MPI_RECV(hruId_local,     num_hru_received, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(hruSegId_local,  num_hru_received, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(area_local,      num_hru_received, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   ! other public/global data
   call MPI_RECV(desireId,1, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(routOpt, 1, MPI_INT,   root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(dt,      1, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(fshape,  1, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(tscale,  1, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(velo,    1, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(diff,    1, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(mann_n,  1, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
   call MPI_RECV(wscale,  1, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

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

   ! find index of desired reach
   ixDesire = findIndex(segId_local, desireId, integerMissing)

   call augment_ntopo(&
                   ! input: model control
                   num_hru_received,             & ! number of HRUs
                   num_seg_received,             & ! number of stream segments
                   ! inout: populate data structures
                   structHRU_local,              & ! ancillary data for HRUs
                   structSEG_local,              & ! ancillary data for stream segments
                   structHRU2seg_local,          & ! ancillary data for mapping hru2basin
                   structNTOPO_local,            & ! ancillary data for network toopology
                   ! output
                   tot_hru_local,                & ! total number of all the upstream hrus for all stream segments
                   tot_upseg_local,              & ! total number of all the immediate upstream segments for all stream segments
                   tot_upstream_local,           & ! total number of all the upstream segments for all stream segments
                   tot_uh_local,                 & ! total number of unit hydrograph for all stream segments
                   ixHRU_desired_local,          & ! indices of desired hrus
                   ixSeg_desired_local,          & ! indices of desired reaches
                   ! output: error control
                   ierr, cmessage)

   print*, 'pid, num_hru_received, num_seg_received, tot_hru_local, tot_upseg_local, tot_upstream_local, tot_uh_local =', pid, num_hru_received, num_seg_received, tot_hru_local, tot_upseg_local, tot_upstream_local, tot_uh_local
   call put_data_struct(num_seg_received, structSEG_local, structNTOPO_local, ierr, cmessage)

  endif

 end subroutine comm_ntopo_data


 ! *********************************************************************
 ! public subroutine: send decomposed hru runoff to tasks and populate data structures
 ! *********************************************************************
 subroutine comm_runoff_data(pid,           & ! input: proc id
                             nNodes,        & ! input: number of procs
                             ixSubHRU,      & ! input: global HRU index in the order of domains
                             hru_per_proc,  & ! input: number of hrus assigned to each proc
                             seg_per_proc,  & ! input: number of hrus assigned to each proc
                             ierr,message)    ! output: error control
  USE public_var
  USE globalData, only : RCHFLX_local           ! reach flux structure
  USE globalData, only : river_basin            ! OMP domain decomposition
  USE globalData, only : ixDesire               ! desired reach index
  USE globalData, only : TSEC                   ! beginning/ending of simulation time step [sec]
  USE globalData, only : runoff_data            !
  USE globalData, only : nHRU                   ! number of HRUs in the whoel river network
  USE globalData, only : nEns                   ! number of ensembles
  ! subroutine: mapping basin runoff to reach runoff
  USE remapping,  only : basin2reach            ! remap runoff from routing basins to routing reaches
  ! subroutines: basin routing
  USE basinUH_module, only : IRF_route_basin    ! perform UH convolution for basin routing
  ! subroutines: river routing
  USE accum_runoff_module, only : accum_runoff  ! upstream flow accumulation
  USE kwt_route_module,    only : kwt_route     ! kinematic wave routing method
  USE irf_route_module,    only : irf_route     ! unit hydrograph (impulse response function) routing method
  !USE kwt_route_module,    only : kwt_route => kwt_route_orig   ! kinematic wave routing method
  !USE irf_route_module,    only : irf_route => irf_route_orig    ! river network unit hydrograph method

  implicit none

  ! input variables
  integer(i4b),             intent(in)  :: pid                      ! process id (MPI)
  integer(i4b),             intent(in)  :: nNodes                   ! number of processes (MPI)
  integer(i4b),allocatable, intent(in)  :: ixSubHRU(:)              ! global HRU index in the order of domains
  integer(i4b),allocatable, intent(in)  :: hru_per_proc(:)          ! number of hrus assigned to each proc (i.e., node)
  integer(i4b),allocatable, intent(in)  :: seg_per_proc(:)          ! number of hrus assigned to each proc (i.e., node)
  ! Output variables
  integer(i4b),             intent(out) :: ierr
  character(len=strLen),    intent(out) :: message                  ! error message
  ! local variables
  real(dp)                              :: basinRunoff_sorted(nHRU) ! sorted basin runoff (m/s) for whole domain
  real(dp),  allocatable                :: basinRunoff_local(:)     ! basin runoff (m/s) for decomposed domains
  real(dp),  allocatable                :: reachRunoff_local(:)     ! reach runoff (m/s) for decomposed domains
  real(dp)                              :: T0, T1                   ! beginning/ending of simulation time step [sec]
  integer(i4b)                          :: num_hru_received
  integer(i4b)                          :: num_seg_received
  integer(i4b)                          :: iens
  integer(i4b)                          :: ixHru1,ixHru2            ! starting index and ending index, respectively, for HRU array
  integer(i4b)                          :: myid                     ! process id indices
  integer(i4b)                          :: iHru                     ! process id indices
  integer(i4b), parameter               :: root=0                   ! root node id
  integer(i4b), parameter               :: send_data_tag=2001
  integer(i4b), parameter               :: return_data_tag=2002
  integer(i4b)                          :: status(MPI_STATUS_SIZE)
  character(len=strLen)                 :: cmessage                 ! error message from subroutine

  ierr=0; message='comm_runoff_data/'

  if (pid == root) then ! this is a root process

    T0=TSEC(0); T1=TSEC(1)

    do iHru = 1,nHRU
      basinRunoff_sorted(iHru) = runoff_data%basinRunoff(ixSubHRU(iHru))
    enddo

    do myid = 1, nNodes-1
      ! Get starting and ending indices
      ixHru2   = sum(hru_per_proc(0:myid))
      ixHru1 = ixHru2 - hru_per_proc(myid) + 1

      ! Send number of elements (hrus and segs)
      call MPI_SEND(hru_per_proc(myid), 1, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      call MPI_SEND(seg_per_proc(myid), 1, MPI_INT, myid, send_data_tag, MPI_COMM_WORLD, ierr)
      ! Send regular 1D array
      call MPI_SEND(basinRunoff_sorted(ixHru1), hru_per_proc(myid), MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
    end do

    allocate(basinRunoff_local(hru_per_proc(root)), stat=ierr)
    allocate(reachRunoff_local(seg_per_proc(root)), stat=ierr)
    basinRunoff_local(1:hru_per_proc(root)) = basinRunoff_sorted(1:hru_per_proc(root))

    ! 1. subroutine: map basin runoff to river network HRUs
    ! map the basin runoff to the stream network...
    call basin2reach(&
                    ! input
                    basinRunoff_local,       & ! intent(in):  basin runoff (m/s)
                    ! output
                    reachRunoff_local,       & ! intent(out): reach runoff (m3/s)
                    ierr, cmessage)            ! intent(out): error control

    ! 2. subroutine: basin route
    if (doesBasinRoute == 1) then
      ! instantaneous runoff volume (m3/s) to data structure
      do iens = 1,nEns
        RCHFLX_local(iens,:)%BASIN_QI = reachRunoff_local(:)
        ! perform Basin routing
        call IRF_route_basin(iens, num_seg_received, ierr, cmessage)
      enddo
    else
      ! no basin routing required (handled outside mizuRoute))
      RCHFLX_local(iens,:)%BASIN_QR(0) = RCHFLX_local(iens,:)%BASIN_QR(1)   ! streamflow from previous step
      RCHFLX_local(iens,:)%BASIN_QR(1) = reachRunoff_local(:)                     ! streamflow (m3/s)
    end if

    ! 3. subroutine: river reach routing
    ! perform upstream flow accumulation
    do iens = 1,nEns
      call accum_runoff(iens,              &    ! input: ensemble index
                        num_seg_received,  &    ! input: number of reaches in the river network
                        ixDesire,          &    ! input: index of verbose reach
                        ierr, cmessage)         ! output: error controls

      ! perform KWT routing
      if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
       call kwt_route(iens,                 & ! input: ensemble index
                      river_basin,          & ! input: river basin data type
                      T0,T1,                & ! input: start and end of the time step
                      ixDesire,             & ! input: index of the desired reach
                      ierr,cmessage)          ! output: error control
      endif

      ! perform IRF routing
      if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
       call irf_route(iens,                 & ! input: ensemble index
                      river_basin,          & ! input: river basin data type
                      ixDesire,             & ! input: index of the desired reach
                      ierr,cmessage)          ! output: error control
      endif
    end do

    ! update model time step bound
    TSEC(0) = TSEC(0) + dt
    TSEC(1) = TSEC(0) + dt

  else

    T0=TSEC(0); T1=TSEC(1)

    call MPI_RECV(num_hru_received, 1, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)
    call MPI_RECV(num_seg_received, 1, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

    allocate(basinRunoff_local(num_hru_received), stat=ierr)
    allocate(reachRunoff_local(num_seg_received), stat=ierr)

    call MPI_RECV(basinRunoff_local, num_hru_received, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

    ! 1. subroutine: map basin runoff to river network HRUs
    ! map the basin runoff to the stream network...
    call basin2reach(&
                    ! input
                    basinRunoff_local,       & ! intent(in):  basin runoff (m/s)
                    ! output
                    reachRunoff_local,       & ! intent(out): reach runoff (m3/s)
                    ierr, cmessage)            ! intent(out): error control

    ! 2. subroutine: basin routing
    if (doesBasinRoute == 1) then
      ! instantaneous runoff volume (m3/s) to data structure
      do iens = 1,nEns
        RCHFLX_local(iens,:)%BASIN_QI = reachRunoff_local(:)
        ! perform Basin routing
        call IRF_route_basin(iens, num_seg_received, ierr, cmessage)
      enddo
    else
      ! no basin routing required (handled outside mizuRoute))
      RCHFLX_local(iens,:)%BASIN_QR(0) = RCHFLX_local(iens,:)%BASIN_QR(1)   ! streamflow from previous step
      RCHFLX_local(iens,:)%BASIN_QR(1) = reachRunoff_local(:)               ! streamflow (m3/s)
    end if

    ! 3. subroutine: river reach routing
    ! perform upstream flow accumulation
    do iens = 1,nEns
      call accum_runoff(iens,              &    ! input: ensemble index
                        num_seg_received,  &    ! input: number of reaches in the river network
                        ixDesire,          &    ! input: index of verbose reach
                        ierr, cmessage)         ! output: error controls

      ! perform KWT routing
      if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
       call kwt_route(iens,                 & ! input: ensemble index
                      river_basin,          & ! input: river basin data type
                      T0,T1,                & ! input: start and end of the time step
                      ixDesire,             & ! input: index of the desired reach
                      ierr,cmessage)          ! output: error control
      endif

      ! perform IRF routing
      if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
       call irf_route(iens,                 & ! input: ensemble index
                      river_basin,          & ! input: river basin data type
                      ixDesire,             & ! input: index of the desired reach
                      ierr,cmessage)          ! output: error control
      endif
    end do

    ! update model time step bound
    TSEC(0) = TSEC(0) + dt
    TSEC(1) = TSEC(0) + dt

  endif ! end of tasks

 end subroutine comm_runoff_data

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




  ! create mapping array to map global index to local index
  !  ixSubSEG(nRch_in) : index that belong to a subbasin in global index
  !  iySubSEG(nRch_in) : local indices for a subbasin
  !  ixMap(nRch_in) : map indix from global to local
  !
  !  examples
  !  global_ix, iySub, ixSubSEG,  ixMap, idNode
  !    1            1         2       1       0
  !    2            2         7       1       0
  !    3            1         5       2       0
  !    4            2         8       3       0
  !    5            3         9       1       0
  !    6            1         1       1       1
  !    7            2         3       2       1
  !    8            3         4       2       1
  !    9            1         6       3       2
  !   10            2        10       2       2
  !
  !  if global index of interest is 4(=ixGlob), ixMap(ixGlob)=3 meaning element with 4 is located in 3rd in local array.
