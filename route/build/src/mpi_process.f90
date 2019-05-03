MODULE mpi_routine

USE mpi

! general public variable
USE public_var,       ONLY: integerMissing

! numeric definition
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

! general utility
USE nr_utility_module, ONLY: indexx               ! sorted index array
USE nr_utility_module, ONLY: arth                 !
USE nr_utility_module, ONLY: findIndex            ! find index within a vector

! MPI utility
USE mpi_mod,           ONLY: shr_mpi_bcast
USE mpi_mod,           ONLY: shr_mpi_gatherV
USE mpi_mod,           ONLY: shr_mpi_scatterV
USE mpi_mod,           ONLY: shr_mpi_allgather

implicit none

integer(i4b),parameter  :: scatter=1  ! communication identifier
integer(i4b),parameter  :: gather=2   ! communication identifier

private

public :: comm_ntopo_data
public :: mpi_route
public :: pass_global_data

contains

 ! *********************************************************************
 ! public subroutine: send reach/hru information to tasks and populate data structures
 ! *********************************************************************
 subroutine comm_ntopo_data(pid,                & ! input: proc id
                            nNodes,             & ! input: number of procs
                            nRch_in,            & ! input: number of stream segments in whole domain
                            nHRU_in,            & ! input: number of HRUs that are connected to reaches
                            structHRU,          & ! input: data structure for HRUs
                            structSEG,          & ! input: data structure for stream segments
                            structHRU2seg,      & ! input: data structure for mapping hru2basin
                            structNTOPO,        & ! input: data structure for network toopology
                            ierr,message)         ! output: error control

  USE public_var
  USE globalData,        ONLY: ixPrint              ! desired reach index
  USE globalData,        ONLY: domains              ! domain data structure - for each domain, pfaf codes and list of segment indices
  USE globalData,        ONLY: nDomain              ! count of decomposed domains (tributaries + mainstems)
  USE globalData,        ONLY: RCHFLX_trib          ! Reach flux data structures (per proc, tributary only)
  USE globalData,        ONLY: KROUTE_trib          ! Reach k-wave data structures (per proc, tributary only)
  USE globalData,        ONLY: NETOPO_trib
  USE globalData,        ONLY: RPARAM_trib
  USE globalData,        ONLY: nEns
  USE globalData,        ONLY: hru_per_proc         ! number of hrus assigned to each proc (size = num of procs+1)
  USE globalData,        ONLY: rch_per_proc         ! number of reaches assigned to each proc (size = num of procs+1)
  USE globalData,        ONLY: ixHRU_order          ! global HRU index in the order of proc assignment
  USE globalData,        ONLY: ixRch_order          ! global reach index in the order of proc assignment
  USE globalData,        ONLY: tribOutlet_per_proc  ! number of tributary outlets per proc (size = nNodes)
  USE globalData,        ONLY: global_ix_comm       ! global reach index at tributary reach outlets to mainstem (size = sum of tributary outlets within entire network)
  USE globalData,        ONLY: local_ix_comm        ! local reach index at tributary reach outlets to mainstem (size = sum of tributary outlets within entire network)
  USE alloc_data,        ONLY: alloc_struct
  USE process_ntopo,     ONLY: augment_ntopo        ! compute all the additional network topology (only compute option = on)
  USE process_ntopo,     ONLY: put_data_struct      !

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
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                   ! error message
  ! Local variables
  ! data structure for decomposed river network data
  type(var_dlength), allocatable              :: structHRU_local(:)        ! ancillary data for HRUs
  type(var_dlength), allocatable              :: structSEG_local(:)        ! ancillary data for stream segments
  type(var_ilength), allocatable              :: structNTOPO_local(:)      ! network topology
  type(var_ilength), allocatable              :: structHRU2seg_local(:)    ! ancillary data for mapping hru2basin
  type(var_clength), allocatable              :: structPFAF_local(:)       ! ancillary data for pfafstetter code
  ! flat array for decomposed river network data
  integer(i4b),      allocatable              :: segId_local(:)            ! reach id for decomposed network
  integer(i4b),      allocatable              :: downSegId_local(:)        ! downstream reach id for decomposed network
  integer(i4b),      allocatable              :: hruId_local(:)            ! hru id array in decomposed network
  integer(i4b),      allocatable              :: hruSegId_local(:)         ! downstream reach id array in decomposed network
  real(dp),          allocatable              :: slope_local(:)            ! reach slope array in decomposed network
  real(dp),          allocatable              :: length_local(:)           ! reach length array in decomposed network
  real(dp),          allocatable              :: area_local(:)             ! hru area in decomposed network
  logical(lgt),      allocatable              :: tribOutlet_local(:)       ! logical to indicate tributary outlet to mainstems
  ! flat array for the entire river network
  integer(i4b)                                :: hruId(nHRU_in)            ! hru id for all the HRUs
  integer(i4b)                                :: hruSegId(nRch_in)         ! hru-to-seg mapping for each hru
  integer(i4b)                                :: segId(nRch_in)            ! reach id for all the segments
  integer(i4b)                                :: downSegId(nRch_in)        ! downstream reach ID for each reach
  real(dp)                                    :: slope(nRch_in)            ! reach slope array for each reach
  real(dp)                                    :: length(nRch_in)           ! reach length array for each reach
  real(dp)                                    :: area(nHRU_in)             ! hru area for each hru
  integer(i4b)                                :: ixNode(nRch_in)           ! node assignment for each reach
  character(len=32)                           :: pfaf(nRch_in)             ! reach pfafcode for each reach
  integer(i4b)                                :: ixLocalSubHRU(nHRU_in)    ! local HRU index
  integer(i4b)                                :: ixLocalSubSEG(nRch_in)    ! local reach index
  logical(lgt),      allocatable              :: tribOutlet(:)             ! logical to indicate tributary outlet to mainstems over entire network
  integer(i4b)                                :: nRch_mainstem             ! number of reaches on the main stem
  integer(i4b)                                :: nHRU_mainstem             ! number of hrus on the main stem
  integer(i4b)                                :: nRch_trib                 ! number of tributary outlets for each proc (scalar)
  integer(i4b),     allocatable               :: nRch_trib_array(:)        ! number of tributary outlets for each proc (array)
  ! flat array for decomposed river network per domain (sub-basin)
  integer(i4b)                                :: idNode(nDomain)           ! node id array for each domain
  integer(i4b)                                :: rnkIdNode(nDomain)        ! ranked node id array for each domain
  integer(i4b)                                :: jHRU,jSeg                 ! ranked indices
  ! miscellaneous
  integer(i4b), allocatable                   :: seq_array(:)
  integer(i4b)                                :: iSeg,iHru                 ! reach and hru loop indices
  integer(i4b)                                :: ix1, ix2
  integer(i4b)                                :: ix,ixx                    ! loop indices
  integer(i4b)                                :: ixSeg1,ixSeg2             ! starting index and ending index, respectively, for reach array
  integer(i4b)                                :: ixHru1,ixHru2             ! starting index and ending index, respectively, for HRU array
  integer(i4b)                                :: idx                       ! node indix (1, ... , nNodes)
  character(len=strLen)                       :: cmessage                  ! error message from subroutine

  ierr=0; message='comm_ntopo_data/'

  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! Part 1: define routing vectors ordered by domain/node
  !  - define the global indices ordered by domain/node
  !  - define the number of reaches/hrus on each processor
  !  - copy the data from the data structures to the ordered routing vectors
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************

  if (pid == root) then ! this is a root process

    ! allocate local and global indices
    allocate(rch_per_proc(-1:nNodes-1), hru_per_proc(-1:nNodes-1), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [rch_per_proc, hru_per_proc]'; return; endif

    allocate(ixHRU_order(nHRU_in),ixRch_order(nRch_in), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ixHRU_order, ixRch_order]'; return; endif

    ! Create segIndex array from domains derived type. The array is sorted from node 0 through nNodes-1
    ! SegIndex Array needs to be contiguous when a chunk is sent to computing node (use sort function...)
    ! start with mainstem domain assigned to root node

    ! domain is a contiguous collection of reaches/HRUs -- multiple domains may be on a single processor

    forall(ix=1:nDomain) idNode(ix) = domains(ix)%idNode ! extracts the processing node from the "domain" data structire
    call indexx(idNode,rnkIdNode) ! rank the processor nodes

    ! loop through the domains
    ixSeg2=0; ixHru2=0 ! last indices of domain chunks
    domain:do ix = 1, nDomain

     ! get the number of stream segments and HRUs in each domain
     ixx = rnkIdNode(ix)
     associate (nSubSeg => size(domains(ixx)%segIndex), nSubHru => size(domains(ixx)%hruIndex) )

     ! define reach index array in order of node assignment
     ixSeg1 = ixSeg2+1             ! start index in the mapping vector
     ixSeg2 = ixSeg1+nSubSeg-1     ! end index in the mapping vector
     ixRch_order(ixSeg1:ixSeg2) = domains(ixx)%segIndex(1:nSubSeg)   ! global seg index per node
     ixLocalSubSEG(ixSeg1:ixSeg2)  = arth(1,1,nSubSeg)                  ! local hru indix per node

     ! define hru index array in order of node assignment
     if (nSubHru>0) then
       ixHru1 = ixHru2+1
       ixHru2 = ixHru1+nSubHru-1
       ixHRU_order(ixHru1:ixHru2)  = domains(ixx)%hruIndex(1:nSubHru) ! global hru index per node
       ixLocalSubHRU(ixHru1:ixHru2)  = arth(1,1,nSubHru)                 ! local hru indix per node
     end if

     ! extra information (debugging)
     ixNode(ixSeg1:ixSeg2)       = domains(ixx)%idNode                 ! node id
     pfaf(ixSeg1:ixSeg2)         = adjustl(trim(domains(ixx)%pfaf))    ! basin pfaf code

     end associate
    end do domain

    ! Count the number of reaches and hrus in each node
    ! index of seg_per_proc and hru_per_proc: -1 -> mainstem, 0 -> small tributaries, 1 through nNodes-1 -> large tributaries
    rch_per_proc = 0
    hru_per_proc = 0
    do ix = 1,nDomain
     idx = domains(ix)%idNode
     rch_per_proc(idx) = rch_per_proc(idx) + size(domains(ix)%segIndex)
     hru_per_proc(idx) = hru_per_proc(idx) + size(domains(ix)%hruIndex)
    end do

    ! define routing vectors ordered by domain/node

    ! reach array
    do iSeg = 1,nRch_in
     jSeg = ixRch_order(iSeg) ! global index, ordered by domain/node
     segId(iSeg)     = structNTOPO(jSeg)%var(ixNTOPO%segId)%dat(1)
     downSegId(iSeg) = structNTOPO(jSeg)%var(ixNTOPO%downSegId)%dat(1)
     slope(iSeg)     = structSEG(  jSeg)%var(ixSEG%slope)%dat(1)
     length(iSeg)    = structSEG(  jSeg)%var(ixSEG%length)%dat(1)
    end do

    ! hru array
    do iHru = 1,nHRU_in
      jHRU = ixHRU_order(iHru)  ! global index, ordered by domain/node
      hruId(iHru)    = structHRU2SEG(jHRU)%var(ixHRU2SEG%HRUid)%dat(1)
      hruSegId(iHru) = structHRU2SEG(jHRU)%var(ixHRU2SEG%hruSegId)%dat(1)
      area(iHru)     = structHRU(    jHRU)%var(ixHRU%area)%dat(1)
    enddo

!    print*, 'ix, segId, ixRch_order, ixLocalSubSEG, ixNode, pfaf'
!    do ix = 1,nRch_in
!      print*, segId(ix), ixRch_order(ix), ixLocalSubSEG(ix), ixNode(ix), pfaf(ix)
!    enddo

  endif  ! if pid==root

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! Part 2: Send the information to individual processors
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************

  ! sends the number of reaches/hrus per proc to all processors
  call shr_mpi_bcast(rch_per_proc, ierr, cmessage)
  call shr_mpi_bcast(hru_per_proc, ierr, cmessage)

  call shr_mpi_bcast(ixRch_order, ierr, cmessage)
  call shr_mpi_bcast(ixHRU_order, ierr, cmessage)

  ! define the number of reaches/hrus on the main stem
  nRch_mainstem = rch_per_proc(-1)
  nHRU_mainstem = hru_per_proc(-1)

  call shr_mpi_scatterV(segId    (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), segId_local,     ierr, cmessage)
  call shr_mpi_scatterV(downSegId(nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), downSegId_local, ierr, cmessage)
  call shr_mpi_scatterV(slope    (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), slope_local,     ierr, cmessage)
  call shr_mpi_scatterV(length   (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), length_local,    ierr, cmessage)

  call shr_mpi_scatterV(hruId    (nHRU_mainstem+1:nHRU_in), hru_per_proc(0:nNodes-1), hruId_local,    ierr, cmessage)
  call shr_mpi_scatterV(hruSegId (nHRU_mainstem+1:nHRU_in), hru_per_proc(0:nNodes-1), hruSegId_local, ierr, cmessage)
  call shr_mpi_scatterV(area     (nHRU_mainstem+1:nHRU_in), hru_per_proc(0:nNodes-1), area_local,     ierr, cmessage)

  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! Part 3: populate local data structures and compute additional ancillary information
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************

  ! allocate space for tributary data structures
  allocate(RCHFLX_trib(nEns,rch_per_proc(pid)), KROUTE_trib(nEns,rch_per_proc(pid)), stat=ierr)

  ! allocate space for local data structures
  call alloc_struct(hru_per_proc(pid),     & ! input: number of HRUs
                    rch_per_proc(pid),     & ! input: number of stream segments
                    structHRU_local,       & ! inout: ancillary data for HRUs
                    structSEG_local,       & ! inout: ancillary data for stream segments
                    structHRU2seg_local,   & ! inout: ancillary data for mapping hru2basin
                    structNTOPO_local,     & ! inout: ancillary data for network toopology
                    structPFAF_local,      & ! inout: ancillary data for pfafstetter code
                    ierr,cmessage)           ! output: error control

  ! Populate local data structures

  ! reach
  do ix = 1,rch_per_proc(pid)
   structNTOPO_local(ix)%var(ixNTOPO%segId)%dat(1)     = segId_local(ix)
   structNTOPO_local(ix)%var(ixNTOPO%downSegId)%dat(1) = downSegId_local(ix)
   structSEG_local  (ix)%var(ixSEG%length)%dat(1)      = length_local(ix)
   structSEG_local  (ix)%var(ixSEG%slope)%dat(1)       = slope_local(ix)
  end do

  ! hru
  do ix=1,hru_per_proc(pid)
   structHRU2SEG_local(ix)%var(ixHRU2SEG%HRUid)%dat(1)    = hruId_local(ix)
   structHRU2SEG_local(ix)%var(ixHRU2SEG%hruSegId)%dat(1) = hruSegId_local(ix)
   structHRU_local    (ix)%var(ixHRU%area)%dat(1)         = area_local(ix)
  end do

  ! find index of desired reach
  if (desireId/=integerMissing) ixPrint = findIndex(segId_local, desireId, integerMissing)

  ! compute additional ancillary infomration
  call augment_ntopo(&
                     ! input: model control
                     hru_per_proc(pid),            & ! number of HRUs
                     rch_per_proc(pid),            & ! number of stream segments
                     ! inout: populate data structures
                     structHRU_local,              & ! ancillary data for HRUs
                     structSEG_local,              & ! ancillary data for stream segments
                     structHRU2seg_local,          & ! ancillary data for mapping hru2basin
                     structNTOPO_local,            & ! ancillary data for network toopology
                     ! output: error control
                     ierr, cmessage)

  ! copy data to routing structres RPARAM_trib and NETOPO_trib
  call put_data_struct(rch_per_proc(pid), structSEG_local, structNTOPO_local, & ! input
                       RPARAM_trib, NETOPO_trib,                              & ! output:
                       ierr, cmessage)

  ! -----------------------------------------------------------------------------
  ! Find "dangling reach", or tributary outlet reaches that link to mainstems
  ! -----------------------------------------------------------------------------
  allocate(tribOutlet_local(rch_per_proc(pid)), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [tribOutlet_local]'; return; endif
  tribOutlet_local = .false.
  do ix=1,rch_per_proc(pid)
    if (structNTOPO_local(ix)%var(ixNTOPO%downSegIndex)%dat(1) == -1) then
      tribOutlet_local(ix) = .true.
    endif
  enddo

  ! ather array for number of tributary outlet reaches per proc
  call shr_mpi_allgather(tribOutlet_local, rch_per_proc(root:nNodes-1), tribOutlet, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! gather array for number of tributary outlet reaches from each proc
  nRch_trib = count(tribOutlet_local) ! number of tributary outlets for each proc
  call shr_mpi_allgather(nRch_trib, 1, nRch_trib_array, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! need to use 0:nNodes-1 bound instead of 1:nNodes
  allocate(tribOutlet_per_proc(0:nNodes-1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [tribOutlet_per_proc]'; return; endif
  tribOutlet_per_proc(0:nNodes-1) = nRch_trib_array

  allocate(global_ix_comm(count(tribOutlet)),seq_array(rch_per_proc(pid)), local_ix_comm(nRch_trib), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [global_ix_comm, loca_ix_comm]'; return; endif

  ! mask non tributary outlet reaches
  ! mainstem part (1:nRch_mainstem) has to be removed from ixRch_order
  ix1 = nRch_mainstem+1; ix2 = sum(rch_per_proc)
  global_ix_comm = pack(ixRch_order(ix1:ix2), tribOutlet)   ! size = number of tributary outlet reaches from all the procs

  seq_array = arth(1,1,rch_per_proc(pid))
  local_ix_comm = pack(seq_array, tribOutlet_local) ! size = number of tributary outlet reaches depending on proc

!deleteme
!   if (pid==0) then
!   print*, 'ix, local_ix_comm, NETOPO_trib%REACHIX, reachId, downstream ID, downstream local ix'
!   do ix =1, nRch_trib
!    print*,  ix, local_ix_comm(ix), NETOPO_trib(local_ix_comm(ix))%REACHIX, NETOPO_trib(local_ix_comm(ix))%REACHID, NETOPO_trib(local_ix_comm(ix))%DREACHK, NETOPO_trib(local_ix_comm(ix))%DREACHI
!   enddo
!   do ix =1, size(global_ix_comm)
!    print*,  ix, global_ix_comm(ix), NETOPO(global_ix_comm(ix))%REACHIX, NETOPO(global_ix_comm(ix))%REACHID, NETOPO(global_ix_comm(ix))%DREACHK, NETOPO(global_ix_comm(ix))%DREACHI
!   enddo
!   do ix =1, size(NETOPO_trib)
!   print*, NETOPO_trib(ix)%REACHID, NETOPO_trib(ix)%DREACHK, NETOPO_trib(ix)%DREACHI, tribOutlet(ix)
!   enddo
!   endif
!deleteme

 end subroutine comm_ntopo_data


 ! *********************************************************************
 ! public subroutine: send decomposed hru runoff to tasks and populate data structures
 ! *********************************************************************
 subroutine mpi_route(pid,           & ! input: proc id
                      nNodes,        & ! input: number of procs
                      iens,          & ! input: ensemble index
                      ierr,message)    ! output: error control
  ! shared data
  USE public_var
  USE dataTypes,  only : KREACH           ! derived data type
  USE globalData, only : NETOPO_trib      ! tributary and mainstem reach netowrk topology structure
  USE globalData, only : NETOPO           ! entire river reach netowrk topology structure
  USE globalData, only : RPARAM_trib      ! tributary and mainstem reach parameter structure
  USE globalData, only : RPARAM           ! entire river reach parameter structure
  USE globalData, only : RCHFLX_trib      ! tributary reach flux structure
  USE globalData, only : RCHFLX           ! entire reach flux structure
  USE globalData, only : KROUTE_trib      ! tributary reach kwt data structure
  USE globalData, only : KROUTE           ! entire river reach kwt sate structure
  USE globalData, only : river_basin      ! OMP domain decomposition
  USE globalData, only : runoff_data      ! runoff data structure
  USE globalData, only : nHRU             ! number of HRUs in the whoel river network
  USE globalData, only : nRch             ! number of reaches in the whoel river network
  USE globalData, only : ixHRU_order      ! global HRU index in the order of proc assignment
  USE globalData, only : ixRch_order      ! global reach index in the order of proc assignment
  USE globalData, only : hru_per_proc     ! number of hrus assigned to each proc (i.e., node)
  USE globalData, only : rch_per_proc     ! number of reaches assigned to each proc (i.e., node)
  USE globalData, only : tribOutlet_per_proc ! number of tributary outlets per proc (size = nNodes)
  USE globalData, only : global_ix_comm   ! global reach index at tributary reach outlets to mainstem (size = sum of tributary outlets in all the procs)
  USE globalData, only : local_ix_comm    ! local reach index at tributary reach outlets to mainstem for each proc (size = sum of tributary outlets in proc)
  ! routing driver
  USE main_route_module, only: main_route ! routing driver

  implicit none

  ! input variables
  integer(i4b),             intent(in)  :: pid                      ! process id (MPI)
  integer(i4b),             intent(in)  :: nNodes                   ! number of processes (MPI)
  integer(i4b),             intent(in)  :: iens                     ! ensemble index
  ! Output variables
  integer(i4b),             intent(out) :: ierr
  character(len=strLen),    intent(out) :: message                  ! error message
  ! local variables
  real(dp)                              :: basinRunoff_sorted(nHRU) ! sorted basin runoff (m/s) for whole domain
  real(dp),     allocatable             :: basinRunoff_local(:)     ! basin runoff (m/s) for tributaries
  integer(i4b), allocatable             :: ixRchProcessed(:)        ! reach indice list to be processed
  integer(i4b)                          :: iHru                     ! loop indices
  integer(i4b)                          :: nSegTrib                 ! number of reaches from one tributary
  integer(i4b)                          :: nSegMain                 ! number of reaches from mainstems
  integer(i4b)                          :: tributary=1              !
  integer(i4b)                          :: mainstem=2               !
  character(len=strLen)                 :: cmessage                 ! error message from subroutine
  ! timing
  integer*8     :: cr, startTime, endTime
  real(dp)      :: rate, elapsedTime

  ierr=0; message='mpi_route/'

  ! Initialize the system_clock
  call system_clock(count_rate=cr)
  rate = real(cr)

  ! Reaches/HRU assigned to root node include BOTH small tributaries and mainstem
  ! First, route "small tributaries" while routing over other bigger tributaries (at slave nodes).

 ! sort the basin runoff in terms of nodes/domains
 if (pid == root) then ! this is a root process
    do iHru = 1,nHRU
      basinRunoff_sorted(iHru) = runoff_data%basinRunoff(ixHRU_order(iHru))
    enddo
  end if

call system_clock(startTime)
  ! Distribute the basin runoff to each process
  call shr_mpi_scatterV(basinRunoff_sorted(hru_per_proc(-1)+1:nHRU), hru_per_proc(0:nNodes-1), basinRunoff_local, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/rate
write(*,"(A,I2,A,1PG15.7,A)") 'pid=',pid,',   elapsed-time [routing/scatter-ro] = ', elapsedTime, ' s'

  ! --------------------------------
  ! Perform tributary routing (for all procs)
  ! --------------------------------
call system_clock(startTime)
  !Idenfity number of tributary reaches for each procs
  nSegTrib = rch_per_proc(pid)
  allocate(ixRchProcessed(nSegTrib), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [ixRchProcessed]'; return; endif

  ! Define processing reach indices in terms of tributary data sets
  ixRchProcessed = arth(1,1,nSegTrib)

  ! Perform routing
  call main_route(iens,              &  ! ensemble index
                  basinRunoff_local, &  ! basin (i.e.,HRU) runoff (m/s)
                  ixRchProcessed,    &  ! indices of reach to be routed
                  river_basin,       &  ! OMP basin decomposition
                  tributary,         &  ! basinType (1-> tributary, 2->mainstem)
                  NETOPO_trib,       &  ! reach topology data structure
                  RPARAM_trib,       &  ! reach parameter data structure
                  ! inout
                  RCHFLX_trib,       &  ! reach flux data structure
                  KROUTE_trib,       &  ! reach state data structure
                  ! output: error handling
                  ierr, message)     ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! make sure that routing at all the procs finished
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/rate
write(*,"(A,I2,A,1PG15.7,A)") 'pid=',pid,',   elapsed-time [routing/tributary-route] = ', elapsedTime, ' s'

  ! --------------------------------
  ! Collect all the tributary flows
  ! --------------------------------
call system_clock(startTime)
  ! flux communication
  call mpi_comm_flux(pid, nNodes,                              &
                     iens,                                     &
                     rch_per_proc(root:nNodes-1),              &
                     ixRch_order(rch_per_proc(root-1)+1:nRch), &
                     arth(1,1,rch_per_proc(pid)),              &
                     gather,                                   & ! input: communication type
                     ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! KWT state communication
  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then

    call mpi_comm_kwt_state(pid, nNodes,         &
                            iens,                &
                            rch_per_proc(root:nNodes-1),              &
                            !tribOutlet_per_proc, & ! input: number of reaches communicate per node (dimension size == number of proc)
                            ixRch_order(rch_per_proc(root-1)+1:nRch), &
                            !global_ix_comm,      & ! input: global reach indices to communicate (dimension size == sum of nRearch)
                            arth(1,1,rch_per_proc(pid)),              &
                            !local_ix_comm,       & ! input: local reach indices per proc (dimension size depends on procs )
                            gather,               & ! communication type
                            ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif

  ! make sure that routing at all the procs finished
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/rate
write(*,"(A,I2,A,1PG15.7,A)") 'pid=',pid,',   elapsed-time [routing/gater-state-flux] = ', elapsedTime, ' s'

  ! --------------------------------
  ! perform mainstem routing
  ! --------------------------------
  if (pid==root) then
call system_clock(startTime)
    ! number of HRUs and reaches from Mainstems
    nSegMain = rch_per_proc(-1)

    if (allocated(ixRchProcessed)) then
      deallocate(ixRchProcessed, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem deallocating array for [ixRchProcessed]'; return; endif
    end if
    allocate(ixRchProcessed(nSegMain), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ixRchProcessed]'; return; endif

    ! Define processing reach indices
    ixRchProcessed = ixRch_order(1:nSegMain)

    call main_route(iens,                    &  ! input: ensemble index
                    runoff_data%basinRunoff, &  ! input: basin (i.e.,HRU) runoff (m/s)
                    ixRchProcessed,          &  ! input: indices of reach to be routed
                    river_basin,             &  ! input: OMP basin decomposition
                    mainstem,                &  ! input: basinType (1-> tributary, 2->mainstem)
                    NETOPO,                  &  ! input: reach topology data structure
                    RPARAM,                  &  ! input: reach parameter data structure
                    RCHFLX,                  &  ! inout: reach flux data structure
                    KROUTE,                  &  ! inout: reach state data structure
                    ierr, message)              ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/rate
write(*,"(A,I2,A,1PG15.7,A)") 'pid=',pid,',   elapsed-time [routing/main_route] = ', elapsedTime, ' s'
  endif ! end of root proc

  ! make sure that routing at all the procs finished
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! --------------------------------
  ! Distribute global states to processors to update states upstream reaches
  ! --------------------------------
  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
call system_clock(startTime)
   call mpi_comm_kwt_state(pid, nNodes,                        & ! input:
                           iens,                               & ! input:
                           tribOutlet_per_proc,                & ! input: number of reaches communicate per node (dimension size == number of proc)
                           global_ix_comm,                     & ! input: global reach indices to communicate (dimension size == sum of nRearch)
                           local_ix_comm,                      & ! input: local reach indices per proc (dimension size depends on procs )
                           scatter,                            & ! input: 1 = scatter
                           ierr, message)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/rate
write(*,"(A,I2,A,1PG15.7,A)") 'pid=',pid,',   elapsed-time [routing/scatter-kwt-state] = ', elapsedTime, ' s'

 endif ! end of kwt option

 end subroutine mpi_route


 ! *********************************************************************
 ! subroutine: fluxes communication
 ! *********************************************************************
 subroutine mpi_comm_flux(pid,          &
                          nNodes,       &
                          iens,         &
                          nReach,       &
                          rchIdxGlobal, &
                          rchIdxLocal,  &
                          commType,     &
                          ierr, message)

  USE public_var,       ONLY: root
  USE globalData,       ONLY: RCHFLX
  USE globalData,       ONLY: RCHFLX_trib

  ! input variables
  integer(i4b),             intent(in)  :: pid                   ! process id (MPI)
  integer(i4b),             intent(in)  :: nNodes                ! number of processes (MPI)
  integer(i4b),             intent(in)  :: iens                  ! ensemble index
  integer(i4b),             intent(in)  :: nReach(0:nNodes-1)    ! number of reaches communicate per node (dimension size == number of proc)
  integer(i4b),             intent(in)  :: rchIdxGlobal(:)       ! reach indices (w.r.t. global) to be transfer (dimension size == sum of nRearch)
  integer(i4b),             intent(in)  :: rchIdxLocal(:)        ! reach indices (w.r.t. local) (dimension size depends on procs )
  integer(i4b),             intent(in)  :: commType              ! communication type 1->scatter, 2->gather otherwise error
  ! output variables
  integer(i4b),             intent(out) :: ierr                  ! error code
  character(len=strLen),    intent(out) :: message               ! error message
  ! local variables
  character(len=strLen)                 :: cmessage              ! error message from a subroutine
  real(dp),     allocatable             :: flux(:,:)             !
  real(dp),     allocatable             :: flux_local(:,:)       !
  real(dp),     allocatable             :: vec_out(:)            ! output vector from mpi gather/scatter routine
  integer(i4b)                          :: nSeg                  ! number of reaches
  integer(i4b)                          :: iSeg, jSeg            ! reach looping index
  integer(i4b)                          :: ix                    ! general looping index
  integer(i4b),parameter                :: nFluxes=6             ! number of flux variables

  ierr=0; message='mpi_comm_flux/'

  ! Number of total reaches to be communicated
  nSeg = sum(nReach)

  if (commType == scatter) then

    if (pid==root) then

      allocate(flux(nSeg, nFluxes), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux]'; return; endif

      do iSeg =1,nSeg ! Loop through tributary reaches
       jSeg = rchIdxGlobal(iSeg)
       flux(iSeg,1) = RCHFLX(iens,jSeg)%BASIN_QR(0)  ! HRU routed flow (previous time step)
       flux(iSeg,2) = RCHFLX(iens,jSeg)%BASIN_QR(1)  ! HRU routed flow (current time step)
       flux(iSeg,3) = RCHFLX(iens,jSeg)%UPSTREAM_QI  ! Upstream accumulated flow
       flux(iSeg,4) = RCHFLX(iens,jSeg)%REACH_Q      ! KWT routed flow
       flux(iSeg,5) = RCHFLX(iens,jSeg)%REACH_Q_IRF  ! IRF routed flow
       flux(iSeg,6) = RCHFLX(iens,jSeg)%BASIN_QI     ! non-HRU routed flow (
      enddo

    end if ! end of root processor operation

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! Distribute global flux data to each process
    allocate(flux_local(nReach(pid),nFluxes), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux_local]'; return; endif
    do ix = 1,nFluxes
      associate(vec_in => reshape(flux(:,ix),[nSeg]))
      call shr_mpi_scatterV(vec_in, nReach, vec_out, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      flux_local(:,ix) = vec_out(:)
      end associate
    end do

    ! update RCHFLX_trib data structure
    do iSeg =1,nReach(pid) ! Loop through reaches per proc

      jSeg = rchIdxLocal(iSeg)

      RCHFLX_trib(iens,jSeg)%BASIN_QR(0) = flux(iSeg,1)  ! HRU routed flow (previous time step)
      RCHFLX_trib(iens,jSeg)%BASIN_QR(1) = flux(iSeg,2)  ! HRU routed flow (current time step)
      RCHFLX_trib(iens,jSeg)%UPSTREAM_QI = flux(iSeg,3)  ! Upstream accumulated flow
      RCHFLX_trib(iens,jSeg)%REACH_Q     = flux(iSeg,4)  ! KWT routed flow
      RCHFLX_trib(iens,jSeg)%REACH_Q_IRF = flux(iSeg,5)  ! IRF routed flow
      RCHFLX_trib(iens,jSeg)%BASIN_QI    = flux(iSeg,6)  ! non-HRU routed flow (

    end do

  elseif (commType == gather) then

    allocate(flux_local(nReach(pid),nFluxes), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux_local]'; return; endif

    do iSeg =1,nReach(pid)  ! Loop through (selected) tributary reaches

      jSeg = rchIdxLocal(iSeg)
      ! Transfer reach fluxes to 2D arrays
      flux_local(iSeg,1) = RCHFLX_trib(iens,jSeg)%BASIN_QR(0)  ! HRU routed flow (previous time step)
      flux_local(iSeg,2) = RCHFLX_trib(iens,jSeg)%BASIN_QR(1)  ! HRU routed flow (current time step)
      flux_local(iSeg,3) = RCHFLX_trib(iens,jSeg)%UPSTREAM_QI  ! Upstream accumulated flow
      flux_local(iSeg,4) = RCHFLX_trib(iens,jSeg)%REACH_Q      ! KWT routed flow
      flux_local(iSeg,5) = RCHFLX_trib(iens,jSeg)%REACH_Q_IRF  ! IRF routed flow
      flux_local(iSeg,6) = RCHFLX_trib(iens,jSeg)%BASIN_QI     ! non-HRU routed flow (

    end do

    allocate(flux(nSeg,nFluxes), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux]'; return; endif

    do ix = 1,nFluxes
      associate(vec_in => reshape(flux_local(:,ix),[nReach(pid)]))
      call shr_mpi_gatherV(vec_in, nReach, vec_out, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      flux(:,ix) = vec_out(:)
      end associate
    end do

    ! put it in global RCHFLX data structure
    if (pid==root) then
      do iSeg =1,nSeg ! Loop through all the reaches involved into communication

        jSeg = rchIdxGlobal(iSeg)

        RCHFLX(iens,jSeg)%BASIN_QR(0) = flux(iSeg,1)
        RCHFLX(iens,jSeg)%BASIN_QR(1) = flux(iSeg,2)
        RCHFLX(iens,jSeg)%UPSTREAM_QI = flux(iSeg,3)
        RCHFLX(iens,jSeg)%REACH_Q     = flux(iSeg,4)
        RCHFLX(iens,jSeg)%REACH_Q_IRF = flux(iSeg,5)
        RCHFLX(iens,jSeg)%BASIN_QI    = flux(iSeg,6)
      end do
    endif

  endif

 end subroutine mpi_comm_flux

 ! *********************************************************************
 ! subroutine: kinematic wave state communication
 ! *********************************************************************
 subroutine mpi_comm_kwt_state(pid,          &
                               nNodes,       &
                               iens,         &
                               nReach,       &
                               rchIdxGlobal, &
                               rchIdxLocal,  &
                               commType,     &
                               ierr, message)

  USE dataTypes,        ONLY: KREACH                             ! derived data type
  USE public_var,       ONLY: root
  USE globalData,       ONLY: KROUTE                             ! entire river reach kwt sate structure
  USE globalData,       ONLY: KROUTE_trib                        ! Reach k-wave data structures (entire river network and tributary only)

  ! input variables
  integer(i4b),             intent(in)  :: pid                   ! process id (MPI)
  integer(i4b),             intent(in)  :: nNodes                ! number of processes (MPI)
  integer(i4b),             intent(in)  :: iens                  ! ensemble index
  integer(i4b),             intent(in)  :: nReach(0:nNodes-1)    ! number of reaches communicate per node (dimension size == number of proc)
  integer(i4b),             intent(in)  :: rchIdxGlobal(:)       ! reach indices (w.r.t. global) to be transfer (dimension size == sum of nRearch)
  integer(i4b),             intent(in)  :: rchIdxLocal(:)        ! reach indices (w.r.t. local) (dimension size depends on procs )
  integer(i4b),             intent(in)  :: commType              ! communication type 1->scatter, 2->gather otherwise error
  ! output variables
  integer(i4b),             intent(out) :: ierr                  ! error code
  character(len=strLen),    intent(out) :: message               ! error message
  ! local variables
  character(len=strLen)                 :: cmessage              ! error message from a subroutine
  type(KREACH), allocatable             :: KROUTE0(:,:)          ! temp KROUTE data structure to hold updated states
  real(dp),     allocatable             :: QF(:),QF_trib(:)
  real(dp),     allocatable             :: QM(:),QM_trib(:)
  real(dp),     allocatable             :: TI(:),TI_trib(:)
  real(dp),     allocatable             :: TR(:),TR_trib(:)
  logical(lgt), allocatable             :: RF(:),RF_trib(:)
  integer(i4b)                          :: ix1, ix2
  integer(i4b)                          :: myid
  integer(i4b)                          :: nSeg                  ! number of reaches
  integer(i4b)                          :: iSeg, jSeg
  integer(i4b)                          :: ixWave
  integer(i4b), allocatable             :: nWave(:)
  integer(i4b), allocatable             :: nWave_trib(:)
  integer(i4b)                          :: totWave(0:nNodes-1)

  ierr=0; message='mpi_comm_kwt_state/'

  ! Number of total reaches to be communicated
  nSeg = sum(nReach)

  if (commType == scatter) then

    ! allocate nWave (number the same at all procs) at each proc
    allocate(nWave(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [nWave]'; return; endif

    if (pid==root) then

     ! extract only tributary reaches
     allocate(KROUTE0(1,nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE0]'; return; endif
     do iSeg =1,nSeg ! Loop through tributary reaches
      jSeg = rchIdxGlobal(iSeg)
      KROUTE0(1, iSeg) = KROUTE(iens,jSeg)
     enddo

     ! convert KROUTE data strucutre to state arrays
     call kwt_struc2array(iens, KROUTE0,  & !input: input state data structure
                          QF,QM,TI,TR,RF, & !output: states array
                          nWave,          & !output: number of waves per reach
                          ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif ! end of root process

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! will have to broadcast updated nWave to all proc
    call MPI_BCAST(nWave, nSeg, MPI_INT, root, MPI_COMM_WORLD, ierr)

    ! total waves from all the tributary reaches in each proc
    ix2=0
    do myid = 0, nNodes-1
      ix1=ix2+1
      ix2=ix1+nReach(myid)-1
      totWave(myid) = sum(nWave(ix1:ix2))
    enddo

    call shr_mpi_scatterV(nWave, nReach(0:nNodes-1), nWave_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Distribute modified KROUTE data to each process
    call shr_mpi_scatterV(QF, totWave(0:nNodes-1), QF_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_scatterV(QM, totWave(0:nNodes-1), QM_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_scatterV(TI, totWave(0:nNodes-1), TI_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_scatterV(TR, totWave(0:nNodes-1), TR_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_scatterV(RF, totWave(0:nNodes-1), RF_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! update KROUTE_trib data structure
    ixWave=1
    do iSeg =1,nReach(pid) ! Loop through reaches per proc

     jSeg = rchIdxLocal(iSeg)

     if (allocated(KROUTE_trib(iens,jSeg)%KWAVE)) then
      deallocate(KROUTE_trib(iens,jSeg)%KWAVE, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [KROUTE_trib(iens,jSeg)%KWAVE]'; return; endif
     endif

     allocate(KROUTE_trib(iens,jSeg)%KWAVE(0:nWave_trib(iSeg)-1),stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE_out(iens,iRch)%KWAVE]'; return; endif

     KROUTE_trib(iens,jSeg)%KWAVE(0:nWave_trib(iSeg)-1)%QF = QF_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
     KROUTE_trib(iens,jSeg)%KWAVE(0:nWave_trib(iSeg)-1)%QM = QM_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
     KROUTE_trib(iens,jSeg)%KWAVE(0:nWave_trib(iSeg)-1)%TI = TI_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
     KROUTE_trib(iens,jSeg)%KWAVE(0:nWave_trib(iSeg)-1)%TR = TR_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
     KROUTE_trib(iens,jSeg)%KWAVE(0:nWave_trib(iSeg)-1)%RF = RF_trib(ixWave:ixWave+nWave_trib(iSeg)-1)

     ixWave=ixWave+nWave_trib(iSeg) !update 1st idex of array

    end do

  elseif (commType == gather) then

    ! allocate nWave (number the same at all procs) and nWave_trib (number dependent on proc) at each proc
    allocate(nWave_trib(nReach(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [nWave_trib]'; return; endif

    ! extract only tributary reaches
    allocate(KROUTE0(1,nReach(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE0]'; return; endif
    do iSeg =1,nReach(pid)  ! Loop through tributary reaches
      jSeg = rchIdxLocal(iSeg)
      KROUTE0(1, iSeg) = KROUTE_trib(iens,jSeg)
    enddo

    ! Transfer KWT state data structure to flat arrays
    call kwt_struc2array(iens,KROUTE0,                            &
                         QF_trib,QM_trib,TI_trib,TR_trib,RF_trib, &
                         nWave_trib,                              &
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_gatherV(nWave_trib, nReach, nWave, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call MPI_BCAST(nWave, nSeg, MPI_INT, root, MPI_COMM_WORLD, ierr)

    ! total waves in reaches in each proc
    ix2=0
    do myid = 0, nNodes-1
      ix1=ix2+1
      ix2=ix1+nReach(myid)-1
      totWave(myid) = sum(nWave(ix1:ix2))
    enddo

    call shr_mpi_gatherV(QF_trib, totWave, QF, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_gatherV(QM_trib, totWave, QM, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_gatherV(TI_trib, totWave, TI, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_gatherV(TR_trib, totWave, TR, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_gatherV(RF_trib, totWave, RF, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! put it in global RCHFLX data structure
    if (pid==root) then
      ixWave=1
      do iSeg =1,nSeg ! Loop through all the reaches involved into communication

        jSeg = rchIdxGlobal(iSeg)

        if (allocated(KROUTE(iens,jSeg)%KWAVE)) then
          deallocate(KROUTE(iens,jSeg)%KWAVE, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [KROUTE(iens,jSeg)%KWAVE]'; return; endif
        endif

        allocate(KROUTE(iens,jSeg)%KWAVE(0:nWave(iSeg)-1),stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE_out(iens,iRch)%KWAVE]'; return; endif

        KROUTE(iens,jSeg)%KWAVE(0:nWave(iSeg)-1)%QF = QF(ixWave:ixWave+nWave(iSeg)-1)
        KROUTE(iens,jSeg)%KWAVE(0:nWave(iSeg)-1)%QM = QM(ixWave:ixWave+nWave(iSeg)-1)
        KROUTE(iens,jSeg)%KWAVE(0:nWave(iSeg)-1)%TI = TI(ixWave:ixWave+nWave(iSeg)-1)
        KROUTE(iens,jSeg)%KWAVE(0:nWave(iSeg)-1)%TR = TR(ixWave:ixWave+nWave(iSeg)-1)
        KROUTE(iens,jSeg)%KWAVE(0:nWave(iSeg)-1)%RF = RF(ixWave:ixWave+nWave(iSeg)-1)
        ixWave=ixWave+nWave(iSeg) !update 1st idex of array
      end do
    endif

  endif

 end subroutine mpi_comm_kwt_state


 ! *********************************************************************
 ! private subroutine
 ! *********************************************************************
 subroutine kwt_struc2array(iens, KROUTE_in,     &  ! input:
                            QF,QM,TI,TR,RF,      &  ! output:
                            nWave,               &
                            ierr, message)
  USE dataTypes,  only : KREACH             ! collection of particles in a given reach
  implicit none
  ! Input
  integer(i4b),          intent(in)              :: iens           ! ensemble index
  type(KREACH),          intent(in), allocatable :: KROUTE_in(:,:) ! reach state data
  ! Output error handling variables
  real(dp),              intent(out),allocatable :: QF(:)          ! flat array for wave Q
  real(dp),              intent(out),allocatable :: QM(:)          ! Flat array for modified Q
  real(dp),              intent(out),allocatable :: TI(:)          ! flat array for entiry time
  real(dp),              intent(out),allocatable :: TR(:)          ! flat array for exit time
  logical(lgt),          intent(out),allocatable :: RF(:)          ! flat array for exiting wave logical
  integer(i4b),          intent(out),allocatable :: nWave(:)       ! number of waves at each reach
  integer(i4b),          intent(out)             :: ierr           ! error code
  character(len=strLen), intent(out)             :: message        ! error message
  ! local variables
  integer(i4b)                                   :: ixWave         ! 1st indix of each reach
  integer(i4b)                                   :: iSeg           ! loop indix
  integer(i4b)                                   :: nSeg           ! number of reaches
  integer(i4b)                                   :: totWave        ! total number of waves from all the reaches

  ierr=0; message='kwt_struc2array/'

  nSeg = size(KROUTE_in(iens,:))
  if (.not.allocated(nWave)) then
    allocate(nWave(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [nWave]'; return; endif
  end if

  do iSeg = 1,nSeg
   nWave(iSeg) = size(KROUTE_in(iens,iSeg)%KWAVE)
  enddo

  totWave=sum(nWave)
  allocate(QF(totWave),QM(totWave),TI(totWave),TR(totWave),RF(totWave), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [QF,QM,TI,TF,RF]'; return; endif

  ixWave = 1
  do iSeg=1,nSeg
   QF(ixWave:ixWave+nWave(iSeg)-1) = KROUTE_in(iens,iSeg)%KWAVE(0:nWave(iSeg)-1)%QF
   QM(ixWave:ixWave+nWave(iSeg)-1) = KROUTE_in(iens,iSeg)%KWAVE(0:nWave(iSeg)-1)%QM
   TI(ixWave:ixWave+nWave(iSeg)-1) = KROUTE_in(iens,iSeg)%KWAVE(0:nWave(iSeg)-1)%TI
   TR(ixWave:ixWave+nWave(iSeg)-1) = KROUTE_in(iens,iSeg)%KWAVE(0:nWave(iSeg)-1)%TR
   RF(ixWave:ixWave+nWave(iSeg)-1) = KROUTE_in(iens,iSeg)%KWAVE(0:nWave(iSeg)-1)%RF
   ixWave = ixWave+nWave(iSeg)
  end do

 end subroutine kwt_struc2array

 ! *********************************************************************
 ! public subroutine: send global data
 ! *********************************************************************
 ! send all the necessary public variables to slave procs
 subroutine pass_global_data(ierr,message)   ! output: error control
  USE public_var, only : root
  USE globalData, only : timeVar           ! time variable
  USE globalData, only : iTime             ! time index
  USE globalData, only : refJulday         ! julian day: reference
  USE globalData, only : startJulday       ! julian day: start
  USE globalData, only : endJulday         ! julian day: end
  USE globalData, only : modJulday         ! julian day: at simulation time step
  USE globalData, only : TSEC              ! beginning/ending of simulation time step [sec]
  USE globalData, only : length_conv
  USE globalData, only : time_conv
  implicit none
  ! Input variables
  ! None
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                   ! error message
  ! Local variables
  ! None

  ierr=0; message='pass_global_data/'

  ! send scalars
  call MPI_BCAST(iTime,       1,     MPI_INT,              root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(refJulday,   1,     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(startJulday, 1,     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(endJulday,   1,     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(modJulday,   1,     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(TSEC,        2,     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(length_conv, 1,     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(time_conv,   1,     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

  call shr_mpi_bcast(timeVar,ierr, message)

 end subroutine pass_global_data

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

