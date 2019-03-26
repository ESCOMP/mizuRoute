MODULE mpi_routine

USE mpi

USE nrtype
USE dataTypes,         ONLY: KREACH
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
USE dataTypes,         ONLY: var_dlength   ! double precision type: var(:)%dat
USE dataTypes,         ONLY: var_clength   ! character type:        var(:)%dat

! named variables
USE var_lookup,        ONLY:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup,        ONLY:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup,        ONLY:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup,        ONLY:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology
USE var_lookup,        ONLY:ixPFAF,   nVarsPFAF    ! index of variables for the pfafstetter code

implicit none

private

public :: comm_ntopo_data
public :: mpi_route
public :: pass_global_data
public :: pass_public_var

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
                            ierr,message)   ! output: error control

  USE public_var
  USE globalData,        ONLY: ixDesire             ! desired reach index
  USE globalData,        ONLY: domains              ! domain data structure - for each domain, pfaf codes and list of segment indices
  USE globalData,        ONLY: nDomain              ! count of decomposed domains (tributaries + mainstems)
  USE globalData,        ONLY: RCHFLX, RCHFLX_trib  ! Reach flux data structures (entire river network and tributary only)
  USE globalData,        ONLY: KROUTE, KROUTE_trib  ! Reach k-wave data structures (entire river network and tributary only)
  USE globalData,        ONLY: NETOPO_trib, NETOPO_main
  USE globalData,        ONLY: RPARAM_trib, RPARAM_main
  USE globalData,        ONLY: fshape, tscale      ! parameters used for basin UH
  USE globalData,        ONLY: velo, diff          ! parameters used for UH
  USE globalData,        ONLY: mann_n, wscale      ! parameters used for KWT
  USE globalData,        ONLY: nEns
  USE globalData,        ONLY: hru_per_proc
  USE globalData,        ONLY: rch_per_proc
  USE alloc_data,        ONLY: alloc_struct
  USE process_ntopo,     ONLY: augment_ntopo        ! compute all the additional network topology (only compute option = on)
  USE process_ntopo,     ONLY: put_data_struct      !
  USE nr_utility_module, ONLY: indexx               ! sorted index array
  USE nr_utility_module, ONLY: arth                 !
  USE nr_utility_module, ONLY: findIndex            ! find index within a vector

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
  type(var_dlength), allocatable              :: structHRU_main(:)         ! ancillary data for HRUs
  type(var_dlength), allocatable              :: structSEG_main(:)         ! ancillary data for stream segments
  type(var_ilength), allocatable              :: structNTOPO_main(:)       ! network topology
  type(var_ilength), allocatable              :: structHRU2seg_main(:)     ! ancillary data for mapping hru2basin
  type(var_clength), allocatable              :: structPFAF_main(:)        ! ancillary data for pfafstetter code
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
  logical(lgt)                                :: isTribSeg(nRch_in)        ! logical to indicate tributary seg or not
  ! flat array for decomposed river network per domain (sub-basin)
  integer(i4b)                                :: idNode(nDomain)           ! node id array for each domain
  integer(i4b)                                :: rnkIdNode(nDomain)        ! ranked node id array for each domain
  ! mpi related variables
  integer(i4b)                                :: displs_hru(0:nNodes-1)       ! entry indices in receiving buffer (routedRunoff) at which to place the array from each proc
  integer(i4b)                                :: displs_rch(0:nNodes-1)       ! entry indices in receiving buffer (routedRunoff) at which to place the array from each proc
  integer(i4b)                                :: iSeg,iHru                 ! reach and hru loop indices
  integer(i4b)                                :: ix,ixx                    ! loop indices
  integer(i4b)                                :: myid                      ! process id indices
  integer(i4b)                                :: ixSeg1,ixSeg2             ! starting index and ending index, respectively, for reach array
  integer(i4b)                                :: ixHru1,ixHru2             ! starting index and ending index, respectively, for HRU array
  integer(i4b)                                :: idx                       ! node indix (1, ... , nNodes)
  character(len=strLen)                       :: cmessage                  ! error message from subroutine

  ierr=0; message='comm_ntopo_data/'

  ! spatial constant routing parameters
  call MPI_BCAST(fshape,  1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tscale,  1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(velo,    1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(diff,    1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(mann_n,  1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(wscale,  1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

  if (pid == root) then ! this is a root process

    ! allocate for the entire river network
    allocate(RCHFLX(nEns,nRch_in), KROUTE(nEns,nRch_in), stat=ierr)
    allocate(ixSubHRU(nHRU_in),ixSubSEG(nRch_in), stat=ierr)

    ! Create segIndex array from domains derived type. The array is sorted from node 0 through nNodes-1
    ! SegIndex Array needs to be contiguous when a chunk is sent to computing node (use sort function...)
    ! start with mainstem domain assigned to root node

    forall(ix=1:nDomain) idNode(ix) = domains(ix)%idNode
    call indexx(idNode,rnkIdNode)

    ixSeg2=0; ixHru2=0 ! last indices of domain chunks
    domain:do ix = 1, nDomain

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

     isTribSeg(ixSeg1:ixSeg2) = domains(ixx)%isTrib                 ! if domain is tributary, T otherwise, F
     ixNode(ixSeg1:ixSeg2)    = domains(ixx)%idNode                 ! node id
     pfaf(ixSeg1:ixSeg2)      = adjustl(trim(domains(ixx)%pfaf))    ! basin pfaf code

     end associate
    end do domain

    ! Count the number of reaches and hrus in each node
    ! index of seg_per_proc and hru_per_proc: -1 -> mainstem, 0 -> small tributaries, 1 through nNodes-1 -> large tributaries
    allocate(rch_per_proc(-1:nNodes-1), hru_per_proc(-1:nNodes-1), stat=ierr)
    rch_per_proc = 0
    hru_per_proc = 0
    do ix = 1,nDomain
     idx = domains(ix)%idNode
     rch_per_proc(idx) = rch_per_proc(idx) + size(domains(ix)%segIndex)
     hru_per_proc(idx) = hru_per_proc(idx) + size(domains(ix)%hruIndex)
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

  endif

  if (pid/=root) then
    allocate(rch_per_proc(-1:nNodes-1), hru_per_proc(-1:nNodes-1), stat=ierr)
  endif
  call MPI_BCAST(rch_per_proc, nNodes+1, MPI_INT, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(hru_per_proc, nNodes+1, MPI_INT, root, MPI_COMM_WORLD, ierr)

  ! Need to compute displacements
  displs_hru(0) = 0
  do myid = 1, nNodes-1
   displs_hru(myid) = sum(hru_per_proc(0:myid-1))
  end do
  displs_rch(0) = 0
  do myid = 1, nNodes-1
   displs_rch(myid) = sum(rch_per_proc(0:myid-1))
  end do

  allocate(segId_local    (rch_per_proc(pid)), &
           downSegId_local(rch_per_proc(pid)), &
           slope_local    (rch_per_proc(pid)), &
           length_local   (rch_per_proc(pid)), &
           hruId_local    (hru_per_proc(pid)), &
           hruSegId_local (hru_per_proc(pid)), &
           area_local     (hru_per_proc(pid)), &
           stat=ierr)

  ! Distribute tributary river data to each process
  call MPI_SCATTERV(segId(rch_per_proc(-1)+1:nRch_in),     rch_per_proc(0:nNodes-1), displs_rch, MPI_INT,      & ! flows from proc
                    segId_local,                           rch_per_proc(pid),                    MPI_INT, root,& ! gathered flows at root node
                    MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(downSegId(rch_per_proc(-1)+1:nRch_in), rch_per_proc(0:nNodes-1), displs_rch, MPI_INT,      & ! flows from proc
                    downSegId_local,                       rch_per_proc(pid),                    MPI_INT, root,& ! gathered flows at root node
                    MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(slope(rch_per_proc(-1)+1:nRch_in),     rch_per_proc(0:nNodes-1), displs_rch, MPI_DOUBLE_PRECISION,      & ! flows from proc
                    slope_local,                           rch_per_proc(pid),                    MPI_DOUBLE_PRECISION, root,& ! gathered flows at root node
                    MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(length(rch_per_proc(-1)+1:nRch_in),    rch_per_proc(0:nNodes-1), displs_rch, MPI_DOUBLE_PRECISION,      & ! flows from proc
                    length_local,                          rch_per_proc(pid),                    MPI_DOUBLE_PRECISION, root,& ! gathered flows at root node
                    MPI_COMM_WORLD, ierr)

  ! Distribute tributary hru data to each process
  call MPI_SCATTERV(hruId(hru_per_proc(-1)+1:nHRU_in),    hru_per_proc(0:nNodes-1), displs_hru, MPI_INT,      & ! flows from proc
                    hruId_local,                          hru_per_proc(pid),                    MPI_INT, root,& ! gathered flows at root node
                    MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(hruSegId(hru_per_proc(-1)+1:nHRU_in), hru_per_proc(0:nNodes-1), displs_hru, MPI_INT,      & ! flows from proc
                    hruSegId_local,                       hru_per_proc(pid),                    MPI_INT, root,& ! gathered flows at root node
                    MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(area(hru_per_proc(-1)+1:nHRU_in),     hru_per_proc(0:nNodes-1), displs_hru, MPI_DOUBLE_PRECISION,      & ! flows from proc
                    area_local,                           hru_per_proc(pid),                    MPI_DOUBLE_PRECISION, root,& ! gathered flows at root node
                    MPI_COMM_WORLD, ierr)

  allocate(RCHFLX_trib(nEns,rch_per_proc(pid)), KROUTE_trib(nEns,rch_per_proc(pid)), stat=ierr)

  call alloc_struct(hru_per_proc(pid),     & ! output: number of HRUs
                    rch_per_proc(pid),     & ! output: number of stream segments
                    structHRU_local,       & ! inout: ancillary data for HRUs
                    structSEG_local,       & ! inout: ancillary data for stream segments
                    structHRU2seg_local,   & ! inout: ancillary data for mapping hru2basin
                    structNTOPO_local,     & ! inout: ancillary data for network toopology
                    structPFAF_local,      & ! inout: ancillary data for pfafstetter code
                    ierr,cmessage)           ! output: error control

  ! Populate data structure
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
  if (ixDesire/=integerMissing) then
   ixDesire = findIndex(segId_local, desireId, integerMissing)
  endif

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

  call put_data_struct(rch_per_proc(pid), structSEG_local, structNTOPO_local, & ! input
                       RPARAM_trib, NETOPO_trib,                             & ! output:
                       ierr, cmessage)

  ! mainstem segments
  if (pid == root) then ! this is a root process

    associate (nSeg_main => rch_per_proc(-1), &
               nHRU_main => hru_per_proc(-1) )

    if (ixDesire/=integerMissing) then
     ixDesire = findIndex(segId(1:nSeg_main), desireId, integerMissing)
    endif

    call alloc_struct(nHRU_main,            & ! output: number of HRUs
                      nSeg_main,            & ! output: number of stream segments
                      structHRU_main,       & ! inout: ancillary data for HRUs
                      structSEG_main,       & ! inout: ancillary data for stream segments
                      structHRU2seg_main,   & ! inout: ancillary data for mapping hru2basin
                      structNTOPO_main,     & ! inout: ancillary data for network toopology
                      structPFAF_main,      & ! inout: ancillary data for pfafstetter code
                      ierr,cmessage)           ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Populate data structure
    ! reach
    do ix = 1, nSeg_main
      structNTOPO_main(ix)%var(ixNTOPO%segId)%dat(1)     = segId(ix)
      structNTOPO_main(ix)%var(ixNTOPO%downSegId)%dat(1) = downSegId(ix)
      structSEG_main  (ix)%var(ixSEG%length)%dat(1)      = length(ix)
      structSEG_main  (ix)%var(ixSEG%slope)%dat(1)       = slope(ix)
    end do
    ! hru
    do ix = 1, nHRU_main
      structHRU2SEG_main(ix)%var(ixHRU2SEG%HRUid)%dat(1)    = hruId(ix)
      structHRU2SEG_main(ix)%var(ixHRU2SEG%hruSegId)%dat(1) = hruSegId(ix)
      structHRU_main    (ix)%var(ixHRU%area)%dat(1)         = area(ix)
    end do

    call augment_ntopo(&
                       ! input: model control
                       nHRU_main,                   & ! number of stream segments
                       nSeg_main,                   & ! number of HRUs
                       ! inout: populate data structures
                       structHRU_main,              & ! ancillary data for HRUs
                       structSEG_main,              & ! ancillary data for stream segments
                       structHRU2seg_main,          & ! ancillary data for mapping hru2basin
                       structNTOPO_main,            & ! ancillary data for network toopology
                       ! output: error control
                       ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call put_data_struct(nSeg_main, structSEG_main, structNTOPO_main, & ! input
                         RPARAM_main, NETOPO_main,                    & ! output
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end associate

!    print*, 'ix, hruId, ixSubHRU, iySubHRU, hruSegId'
!    do ix = 1,nHRU_in
!      print*, ix, hruId(ix), ixSubHRU(ix), iySubHRU(ix), hruSegId(ix)
!    end do
!    print*, 'ix, segId, ixSubSEG, iySubSEG, ixNode, pfaf'
!    do ix = 1,nRch_in
!      print*, ix, segId(ix), ixSubSEG(ix), iySubSEG(ix), ixNode(ix), pfaf(ix)
!    enddo

  ! for other procs
  endif

!deleteme
!   if (pid==7) then
!   do ix =1, size(NETOPO_trib)
!   print*, NETOPO_trib(ix)%REACHID, RPARAM_trib(ix)%BASAREA,  NETOPO_trib(ix)%HRUID
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
  USE globalData, only : NETOPO_trib, NETOPO_main  ! tributary and mainstem reach netowrk topology structure
  USE globalData, only : NETOPO                    ! entire river reach netowrk topology structure
  USE globalData, only : RPARAM_trib, RPARAM_main  ! tributary and mainstem reach parameter structure
  USE globalData, only : RPARAM                    ! entire river reach parameter structure
  USE globalData, only : RCHFLX_trib               ! tributary reach flux structure
  USE globalData, only : RCHFLX                    ! entire reach flux structure
  USE globalData, only : KROUTE_trib               ! tributary reach kwt data structure
  USE globalData, only : KROUTE                    ! entire river reach kwt sate structure
  USE globalData, only : river_basin               ! OMP domain decomposition
  USE globalData, only : ixDesire                  ! desired reach index
  USE globalData, only : runoff_data               ! runoff data structure
  USE globalData, only : nHRU                      ! number of HRUs in the whoel river network
  USE globalData, only : ixHRU_order               ! global HRU index in the order of proc assignment
  USE globalData, only : ixRch_order               ! global reach index in the order of proc assignment
  USE globalData, only : hru_per_proc              ! number of hrus assigned to each proc (i.e., node)
  USE globalData, only : rch_per_proc              ! number of reaches assigned to each proc (i.e., node)
  USE globalData, only : TSEC                      ! beginning/ending of simulation time step [sec]
  ! external subroutines:
  ! mapping basin runoff to reach runoff
  USE remapping,  only : basin2reach               ! remap runoff from routing basins to routing reaches
  ! basin routing
  USE basinUH_module, only : IRF_route_basin       ! perform UH convolution for basin routing
  ! river routing
  USE accum_runoff_module, only : accum_runoff     ! upstream flow accumulation
  !USE kwt_route_module,    only : kwt_route     ! kinematic wave routing method
  !USE irf_route_module,    only : irf_route     ! unit hydrograph (impulse response function) routing method
  USE kwt_route_module,    only : kwt_route => kwt_route_orig   ! kinematic wave routing method
  USE irf_route_module,    only : irf_route => irf_route_orig    ! river network unit hydrograph method

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
  real(dp),     allocatable             :: reachRunoff_local(:)     ! reach runoff (m/s) for tributaries
  real(dp),     allocatable             :: reachRunoff_local1(:)     ! reach runoff (m/s) for tributaries
  real(dp),     allocatable             :: basinRunoff_main(:)      ! basin runoff (m/s) for mainstems (processed at root)
  real(dp),     allocatable             :: reachRunoff_main(:)      ! reach runoff (m/s) for mainstems (processed at root)
  real(dp),     allocatable             :: routedRunoff_local(:,:)  ! tributary routed runoff (m/s) for each proc
  real(dp),     allocatable             :: routedRunoff(:,:)        ! tributary routed runoff (m/s) gathered from each proc
  real(dp)                              :: T0,T1                    ! beginning/ending of simulation time step [sec]
  real(dp),     allocatable             :: QF(:), QF_trib(:)
  real(dp),     allocatable             :: QM(:), QM_trib(:)
  real(dp),     allocatable             :: TI(:), TI_trib(:)
  real(dp),     allocatable             :: TR(:), TR_trib(:)
  logical(lgt), allocatable             :: RF(:), RF_trib(:)
  type(KREACH), allocatable             :: KROUTE0(:,:)             !temp KROUTE data structure to hold updated states
  integer(i4b), allocatable             :: nWave(:), nWave_trib(:)
  integer(i4b), allocatable             :: ixWave
  integer(i4b)                          :: displs(0:nNodes-1)       ! entry indices in receiving buffer (routedRunoff) at which to place the array from each proc
  integer(i4b)                          :: displs_kw(0:nNodes-1)    ! entry indices in receiving buffer (state arrays) at which to place the array from each proc
  integer(i4b)                          :: totWave(0:nNodes-1)
  integer(i4b)                          :: ix1, ix2
  integer(i4b)                          :: iHru,iSeg, myid          ! loop indices
  integer(i4b)                          :: nRchTrib                 ! number of tributary reaches
  character(len=strLen)                 :: cmessage                 ! error message from subroutine

  ierr=0; message='mpi_route/'

  T0=TSEC(0); T1=TSEC(1)

  if (pid == root) then ! this is a root process
    ! Reaches/HRU assigned to root node include BOTH small tributaries and mainstem
    ! First, route "small tributaries" while routing over other bigger tributaries (at slave nodes).
    do iHru = 1,nHRU
      basinRunoff_sorted(iHru) = runoff_data%basinRunoff(ixHRU_order(iHru))
    enddo
  end if

  ! Need to compute displacements
  displs(0) = 0
  do myid = 1, nNodes-1
   displs(myid) = sum(hru_per_proc(0:myid-1))
  end do
  allocate(basinRunoff_local(hru_per_proc(pid)),    &
           reachRunoff_local(rch_per_proc(pid)),    &
           reachRunoff_local1(rch_per_proc(pid)),    &
           routedRunoff_local(rch_per_proc(pid),6), &
           stat=ierr)
  ! Distribute modified KROUTE data to each process
  call MPI_SCATTERV(basinRunoff_sorted(hru_per_proc(-1)+1:nHRU), hru_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION,      & ! flows from proc
                    basinRunoff_local,                           hru_per_proc(pid),                MPI_DOUBLE_PRECISION, root,& ! gathered flows at root node
                    MPI_COMM_WORLD, ierr)

! --------------------------------
! BEGIN: Put this in seprate routine
! --------------------------------
  ! 1. subroutine: map basin runoff to river network HRUs
  ! map the basin runoff to the stream network...
  call basin2reach(&
                  ! input
                  basinRunoff_local,       & ! basin runoff (m/s)
                  NETOPO_trib,             & ! reach topology
                  RPARAM_trib,             & ! reach parameter
                  ! output
                  reachRunoff_local,       & ! intent(out): reach runoff (m3/s)
                  ierr, cmessage)            ! intent(out): error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! 2. subroutine: basin route
  if (doesBasinRoute == 1) then
    ! instantaneous runoff volume (m3/s) to data structure
    RCHFLX_trib(iens,:)%BASIN_QI = reachRunoff_local(:)
    routedRunoff_local(:,6) = RCHFLX_trib(iens,:)%BASIN_QI
    ! perform Basin routing
    call IRF_route_basin(iens, RCHFLX_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  else
    ! no basin routing required (handled outside mizuRoute))
    RCHFLX_trib(iens,:)%BASIN_QR(0) = RCHFLX_trib(iens,:)%BASIN_QR(1)   ! streamflow from previous step
    RCHFLX_trib(iens,:)%BASIN_QR(1) = reachRunoff_local(:)                     ! streamflow (m3/s)
  end if
  ! populate reach fluxes in 2D array
  routedRunoff_local(:,1) = RCHFLX_trib(iens,:)%BASIN_QR(0)
  routedRunoff_local(:,2) = RCHFLX_trib(iens,:)%BASIN_QR(1)

  ! 3. subroutine: river reach routing
  ! perform upstream flow accumulation
  call accum_runoff(iens,              &  ! input: ensemble index
                    ixDesire,          &  ! input: index of verbose reach
                    NETOPO_trib,       &  ! input: reach topology data structure
                    RCHFLX_trib,       &  ! inout: reach flux data structure
                    ierr, cmessage)       ! output: error controls
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! populate reach fluxes in 2D array
  routedRunoff_local(:,3) = RCHFLX_trib(iens,:)%UPSTREAM_QI

  ! perform KWT routing
  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
   call kwt_route(iens,                 & ! input: ensemble index
                  river_basin,          & ! input: river basin data type
                  T0,T1,                & ! input: start and end of the time step
                  ixDesire,             & ! input: index of the desired reach
                  NETOPO_trib,          & ! input: reach topology data structure
                  RPARAM_trib,          & ! input: reach parameter data structure
                  KROUTE_trib,          & ! inout: reach state data structure
                  RCHFLX_trib,          & ! inout: reach flux data structure
                  ierr,cmessage)          ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! populate reach fluxes in 2D array
   routedRunoff_local(:,4) = RCHFLX_trib(iens,:)%REACH_Q
   call kwt_struc2array(iens,KROUTE_trib,                        &
                        QF_trib,QM_trib,TI_trib,TR_trib,RF_trib, &
                        nWave_trib,                              &
                        ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   totWave(pid) = sum(nWave_trib)
  endif

  ! perform IRF routing
  if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
   call irf_route(iens,                & ! input: ensemble index
                  river_basin,         & ! input: river basin data type
                  ixDesire,            & ! input: index of the desired reach
                  NETOPO_trib,         & ! input: reach topology data structure
                  RCHFLX_trib,         & ! inout: reach flux data structure
                  ierr,cmessage)         ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! populate reach fluxes in 2D array
   routedRunoff_local(:,5) = RCHFLX_trib(iens,:)%REACH_Q_IRF
  endif

  ! make sure that routing at all the procs finished
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! --------------------------------
  ! Subroutine: Collect all the tributary flows and states
  ! --------------------------------
  ! Need to compute displacements
  displs(0) = 0
  do myid = 1, nNodes-1
   displs(myid) = sum(rch_per_proc(0:myid-1))
  end do

  ! GATHER tributary routed flow and states from slave procs
  ! Can GATHER all together
  nRchTrib = sum(rch_per_proc(0:nNodes-1))
  allocate(routedRunoff(nRchTrib,6), stat=ierr)
  call MPI_GATHERV(routedRunoff_local(:,1), rch_per_proc(pid),                MPI_DOUBLE_PRECISION,       & ! flows from proc
                   routedRunoff(:,1),       rch_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                   MPI_COMM_WORLD, ierr)
  call MPI_GATHERV(routedRunoff_local(:,2), rch_per_proc(pid),                MPI_DOUBLE_PRECISION,       & ! flows from proc
                   routedRunoff(:,2),       rch_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                   MPI_COMM_WORLD, ierr)
  call MPI_GATHERV(routedRunoff_local(:,3), rch_per_proc(pid),                MPI_DOUBLE_PRECISION,       & ! flows from proc
                   routedRunoff(:,3),       rch_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                   MPI_COMM_WORLD, ierr)
  call MPI_GATHERV(routedRunoff_local(:,4), rch_per_proc(pid),                MPI_DOUBLE_PRECISION,       & ! flows from proc
                   routedRunoff(:,4),       rch_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                   MPI_COMM_WORLD, ierr)
  call MPI_GATHERV(routedRunoff_local(:,5), rch_per_proc(pid),                MPI_DOUBLE_PRECISION,       & ! flows from proc
                   routedRunoff(:,5),       rch_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                   MPI_COMM_WORLD, ierr)
  call MPI_GATHERV(routedRunoff_local(:,6), rch_per_proc(pid),                MPI_DOUBLE_PRECISION,       & ! flows from proc
                   routedRunoff(:,6),       rch_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                   MPI_COMM_WORLD, ierr)

  if (pid==root) then
    do iSeg =1,nRchTrib ! Loop through tributary reaches
      associate(iRch => ixRch_order(rch_per_proc(-1)+iSeg)) ! the first "rch_per_proc(-1)" reaches are mainstems
      ! flux
      RCHFLX(iens,iRch)%BASIN_QI    = routedRunoff(iSeg,6)
      RCHFLX(iens,iRch)%BASIN_QR(0) = routedRunoff(iSeg,1)
      RCHFLX(iens,iRch)%BASIN_QR(1) = routedRunoff(iSeg,2)
      RCHFLX(iens,iRch)%UPSTREAM_QI = routedRunoff(iSeg,3)
      RCHFLX(iens,iRch)%REACH_Q     = routedRunoff(iSeg,4)
      RCHFLX(iens,iRch)%REACH_Q_IRF = routedRunoff(iSeg,5)
      end associate
    end do
  endif

  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
    !KW states
    ! collect arrays storing number of waves for each reach from each proc
    allocate(nWave(nRchTrib), stat=ierr)
    call MPI_GATHERV(nWave_trib, rch_per_proc(pid),                MPI_INT,       & ! number of wave from proc
                     nWave,      rch_per_proc(0:nNodes-1), displs, MPI_INT, root, & ! gathered number of wave at root node
                     MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nWave, nRchTrib, MPI_INT, root, MPI_COMM_WORLD, ierr)

    ! total waves in reaches in each proc
    ix2=0
    do myid = 0, nNodes-1
      ix1=ix2+1
      ix2=ix1+rch_per_proc(myid)-1
      totWave(myid) = sum(nWave(ix1:ix2))
    enddo

    ! displacement of wave array
    displs_kw(0) = 0
    do myid = 1, nNodes-1
     displs_kw(myid) = sum(totWave(0:myid-1))
    end do

    ! collect state arrays from each proc
    allocate(QF(sum(totWave)), QM(sum(totWave)), TI(sum(totWave)), TR(sum(totWave)), RF(sum(totWave)), stat=ierr)
    ! Can GATHER all together for QF, QM, TI, TR
    call MPI_GATHERV(QF_trib, totWave(pid),                   MPI_DOUBLE_PRECISION,       & ! flows from proc
                     QF,      totWave(0:nNodes-1), displs_kw, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(QM_trib, totWave(pid),                   MPI_DOUBLE_PRECISION,       & ! flows from proc
                     QM,      totWave(0:nNodes-1), displs_kw, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(TI_trib, totWave(pid),                   MPI_DOUBLE_PRECISION,       & ! flows from proc
                     TI,      totWave(0:nNodes-1), displs_kw, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(TR_trib, totWave(pid),                   MPI_DOUBLE_PRECISION,       & ! flows from proc
                     TR,      totWave(0:nNodes-1), displs_kw, MPI_DOUBLE_PRECISION, root, & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(RF_trib, totWave(pid),                   MPI_LOGICAL,       & ! flows from proc
                     RF,      totWave(0:nNodes-1), displs_kw, MPI_LOGICAL, root, & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)

    ! clear tribuary state arrays for all procs
    deallocate(QF_trib, QM_trib, TI_trib, TR_trib, RF_trib, stat=ierr)

    ! put it in global RCHFLX data structure
    if (pid==root) then
      ixWave=1
      do iSeg =1,nRchTrib ! Loop through tributary reaches
        associate(iRch => ixRch_order(rch_per_proc(-1)+iSeg)) ! the first "rch_per_proc(-1)" reaches are mainstems
        ! states
        if (allocated(KROUTE(iens,iRch)%KWAVE)) then
          deallocate(KROUTE(iens,iRch)%KWAVE, stat=ierr)
        endif
        allocate(KROUTE(iens,iRch)%KWAVE(0:nWave(iSeg)-1),stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE_out(iens,iRch)%KWAVE]'; return; endif
        KROUTE(iens,iRch)%KWAVE(0:nWave(iSeg)-1)%QF = QF(ixWave:ixWave+nWave(iSeg)-1)
        KROUTE(iens,iRch)%KWAVE(0:nWave(iSeg)-1)%QM = QM(ixWave:ixWave+nWave(iSeg)-1)
        KROUTE(iens,iRch)%KWAVE(0:nWave(iSeg)-1)%TI = TI(ixWave:ixWave+nWave(iSeg)-1)
        KROUTE(iens,iRch)%KWAVE(0:nWave(iSeg)-1)%TR = TR(ixWave:ixWave+nWave(iSeg)-1)
        KROUTE(iens,iRch)%KWAVE(0:nWave(iSeg)-1)%RF = RF(ixWave:ixWave+nWave(iSeg)-1)
        end associate
        ixWave=ixWave+nWave(iSeg) !update 1st idex of array
      end do
      deallocate(QF, QM, TI, TR, RF, stat=ierr)
    endif
  endif

  ! make sure that routing at all the procs finished
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! --------------------------------
! mainstem routing
! --------------------------------
  if (pid==root) then
    allocate(basinRunoff_main(hru_per_proc(root-1)), stat=ierr)
    allocate(reachRunoff_main(rch_per_proc(root-1)), stat=ierr)

    basinRunoff_main(1:hru_per_proc(root-1)) = basinRunoff_sorted(1:hru_per_proc(root-1))

    ! 1. subroutine: map basin runoff to river network HRUs
    ! map the basin runoff to the stream network...
    call basin2reach(&
                    ! input
                    basinRunoff_main,       & ! basin runoff (m/s)
                    NETOPO_main,            & ! reach topology
                    RPARAM_main,            & ! reach parameter
                    ! output
                    reachRunoff_main,       & ! intent(out): reach runoff (m3/s)
                    ierr, cmessage)           ! intent(out): error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! 2. subroutine: basin route
    if (doesBasinRoute == 1) then
      ! instantaneous runoff volume (m3/s) to data structure
      do iSeg = 1,rch_per_proc(root-1)
       RCHFLX(iens,ixRch_order(iSeg))%BASIN_QI = reachRunoff_main(iSeg)
      end do
      ! perform Basin routing
      call IRF_route_basin(iens,                 &              ! input:  ensemble index
                           RCHFLX,               &              ! inout:  reach flux data structure
                           ierr, cmessage,       &              ! output: error controls
                           ixRch_order(1:rch_per_proc(root-1))) ! optional input: indices of reach to be routed
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    else
      ! no basin routing required (handled outside mizuRoute))
      do iSeg = 1,rch_per_proc(root-1)
      RCHFLX(iens,ixRch_order(iSeg))%BASIN_QR(0) = RCHFLX(iens,iSeg)%BASIN_QR(1)   ! streamflow from previous step
      RCHFLX(iens,ixRch_order(iSeg))%BASIN_QR(1) = reachRunoff_main(iSeg)          ! streamflow (m3/s)
      end do
    end if

    ! 3. subroutine: river reach routing
    ! perform upstream flow accumulation
    call accum_runoff(iens,                &  ! input: ensemble index
                      ixDesire,            &  ! input: index of verbose reach
                      NETOPO,              &  ! input: reach topology data structure
                      RCHFLX,              &  ! inout: reach flux data structure
                      ierr, cmessage,      &  ! output: error controls
                      ixRch_order(1:rch_per_proc(root-1)))      ! optional input: indices of reach to be routed
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! perform KWT routing
!    if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
!     call kwt_route(iens,                 & ! input: ensemble index
!                    river_basin,          & ! input: river basin data type
!                    T0,T1,                & ! input: start and end of the time step
!                    ixDesire,             & ! input: index of the desired reach
!                    NETOPO,               & ! input: reach topology data structure
!                    RPARAM,               & ! input: reach parameter data structure
!                    KROUTE,               & ! inout: reach state data structure
!                    RCHFLX,               & ! inout: reach flux data structure
!                    ierr,cmessage,        & ! output: error control
!                    ixRch_order(1:rch_per_proc(root-1))) ! optional input: indices of reach to be routed
!     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    endif

    ! perform IRF routing
    if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
     call irf_route(iens,                 & ! input: ensemble index
                    river_basin,          & ! input: river basin data type
                    ixDesire,             & ! input: index of the desired reach
                    NETOPO,               & ! input: reach topology data structure
                    RCHFLX,               & ! inout: reach flux data structure
                    ierr,cmessage,        & ! output: error control
                    ixRch_order(1:rch_per_proc(root-1)))      ! optional input: indices of reach to be routed
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

  endif ! end of root proc

  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
   ! get all tributary kwt states (excluding mainstems)
   if (pid==root) then
    ! extract only tributary reaches
    allocate(KROUTE0(1,nRchTrib), stat=ierr)
    do iSeg =1,nRchTrib ! Loop through tributary reaches
     associate(iRch => ixRch_order(rch_per_proc(-1)+iSeg)) ! the first "rch_per_proc(-1)" reaches are mainstems
     KROUTE0(1,iSeg) = KROUTE(iens,iRch)
     end associate
    enddo
    ! convert KROUTE data strucutre to state arrays
    call kwt_struc2array(iens, KROUTE0,  & !input: input state data structure
                         QF,QM,TI,TR,RF, & !output: states array
                         nWave,          & !output: number of waves per reach
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    deallocate(KROUTE0, stat=ierr)
   endif

   ! total waves from all the tributary reaches in each proc
   ix2=0
   do myid = 0, nNodes-1
     ix1=ix2+1
     ix2=ix1+rch_per_proc(myid)-1
     totWave(myid) = sum(nWave(ix1:ix2))
   enddo

   ! displacement of wave array
   displs_kw(0) = 0
   do myid = 1, nNodes-1
    displs_kw(myid) = sum(totWave(0:myid-1))
   end do

   allocate(QF_trib(totWave(pid)),QM_trib(totWave(pid)),TI_trib(totWave(pid)),TR_trib(totWave(pid)),RF_trib(totWave(pid)), stat=ierr)

   call MPI_SCATTERV(nWave,      rch_per_proc(0:nNodes-1), displs, MPI_INT,       &   ! number of wave from proc
                     nWave_trib, rch_per_proc(pid),                MPI_INT, root, &   ! gathered number of wave at root node
                     MPI_COMM_WORLD, ierr)

   ! Distribute modified KROUTE data to each process
   call MPI_SCATTERV(QF,      totWave(0:nNodes-1), displs_kw, MPI_DOUBLE_PRECISION, & ! flows from proc
                     QF_trib, totWave(pid),        MPI_DOUBLE_PRECISION, root,      & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)
   call MPI_SCATTERV(QM,      totWave(0:nNodes-1), displs_kw, MPI_DOUBLE_PRECISION, & ! flows from proc
                     QM_trib, totWave(pid),        MPI_DOUBLE_PRECISION, root,      & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)
   call MPI_SCATTERV(TI,      totWave(0:nNodes-1), displs_kw, MPI_DOUBLE_PRECISION, & ! flows from proc
                     TI_trib, totWave(pid),        MPI_DOUBLE_PRECISION, root,      & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)
   call MPI_SCATTERV(TR,      totWave(0:nNodes-1), displs_kw, MPI_DOUBLE_PRECISION, & ! flows from proc
                     TR_trib, totWave(pid),        MPI_DOUBLE_PRECISION, root,      & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)
   call MPI_SCATTERV(RF,      totWave(0:nNodes-1), displs_kw, MPI_LOGICAL,          & ! flows from proc
                     RF_trib, totWave(pid),        MPI_LOGICAL, root,               & ! gathered flows at root node
                     MPI_COMM_WORLD, ierr)

   ! update KROUTE_trib data structure
   ixWave=1
   do iSeg =1,rch_per_proc(pid) ! Loop through reaches per proc

    if (allocated(KROUTE_trib(iens,iSeg)%KWAVE)) then
     deallocate(KROUTE_trib(iens,iSeg)%KWAVE, stat=ierr)
    endif

    allocate(KROUTE_trib(iens,iSeg)%KWAVE(0:nWave_trib(iSeg)-1),stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE_out(iens,iRch)%KWAVE]'; return; endif

    KROUTE_trib(iens,iSeg)%KWAVE(0:nWave_trib(iSeg)-1)%QF = QF_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
    KROUTE_trib(iens,iSeg)%KWAVE(0:nWave_trib(iSeg)-1)%QM = QM_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
    KROUTE_trib(iens,iSeg)%KWAVE(0:nWave_trib(iSeg)-1)%TI = TI_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
    KROUTE_trib(iens,iSeg)%KWAVE(0:nWave_trib(iSeg)-1)%TR = TR_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
    KROUTE_trib(iens,iSeg)%KWAVE(0:nWave_trib(iSeg)-1)%RF = RF_trib(ixWave:ixWave+nWave_trib(iSeg)-1)

    ixWave=ixWave+nWave_trib(iSeg) !update 1st idex of array

   end do
 endif

 end subroutine mpi_route

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
  allocate(nWave(nSeg), stat=ierr)
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
 subroutine pass_global_data(pid, nNodes, ierr,message)   ! output: error control
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
  integer(i4b),                   intent(in)  :: pid                      ! process id (MPI)
  integer(i4b),                   intent(in)  :: nNodes                   ! number of processes (MPI)
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                   ! error message
  ! Local variables
  integer(i4b)                                :: nTime                     ! number of time step in runoff data
  integer(i4b)                                :: nTime_recv                ! number of time step in runoff data in slave proc
  integer(i4b)                                :: myid                      ! process id indices
  integer(i4b), parameter                     :: send_data_tag=2001
  integer(i4b), parameter                     :: return_data_tag=2002
  integer(i4b)                                :: status(MPI_STATUS_SIZE)

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

  ! send allocatable arrays
  if (pid == root) then ! this is a root process
    ! number of nTime
    nTime = size(timeVar)

    do myid = 1, nNodes-1
     ! number of nTime
     call MPI_SEND(nTime,                   1,  MPI_INT,              myid, send_data_tag, MPI_COMM_WORLD, ierr)
     call MPI_SEND(timeVar(1),          nTime,  MPI_DOUBLE_PRECISION, myid, send_data_tag, MPI_COMM_WORLD, ierr)
    end do
  else
     ! number of nTime
     call MPI_RECV(nTime_recv, 1, MPI_INT, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

     allocate(timeVar(nTime_recv), stat=ierr)
     call MPI_RECV(timeVar, nTime_recv, MPI_DOUBLE_PRECISION, root, send_data_tag, MPI_COMM_WORLD, status, ierr)

  endif

 end subroutine pass_global_data


 ! *********************************************************************
 ! public subroutine: send public var to tasks
 ! *********************************************************************
 ! (temporarily)
 ! send all the necessary public variables updated from control file to slave procs
 subroutine pass_public_var(ierr,message)   ! output: error control

  USE public_var, only : root
  USE public_var, only : hydGeometryOption       !
  USE public_var, only : topoNetworkOption       !
  USE public_var, only : computeReachList        !
  USE public_var, only : routOpt                 ! reach routing options
  USE public_var, only : desireId                ! reach id desired for print out
  USE public_var, only : doesBasinRoute          ! logical whether basin routing is performed or not
  USE public_var, only : dt                      ! runoff time step [sec]

  implicit none
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                   ! error message

  ierr=0; message='pass_public_var/'

  call MPI_BCAST(hydGeometryOption, 1, MPI_LOGICAL,          root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(topoNetworkOption, 1, MPI_LOGICAL,          root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(computeReachList,  1, MPI_LOGICAL,          root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(routOpt,           1, MPI_INT,              root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(desireId,          1, MPI_INT,              root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(doesBasinRoute,    1, MPI_INT,              root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(dt,                1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

 end subroutine pass_public_var

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

