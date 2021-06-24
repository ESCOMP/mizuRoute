MODULE mpi_routine

USE mpi

! general public variable
USE public_var                   ! iulog, integerMissing, realMissing, root, calendar, time_units
USE globalData, ONLY: masterproc

! numeric definition
USE nrtype

USE dataTypes, ONLY: var_ilength   ! integer type:     var(:)%dat
USE dataTypes, ONLY: var_dlength   ! double precision type: var(:)%dat
USE dataTypes, ONLY: var_clength   ! character type:        var(:)%dat
USE dataTypes, ONLY: subbasin_omp  ! openMP domain data structure

! named variables
USE var_lookup, ONLY:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup, ONLY:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup, ONLY:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup, ONLY:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology
USE var_lookup, ONLY:ixPFAF,   nVarsPFAF    ! index of variables for the pfafstetter code

! general utility
USE nr_utility_module, ONLY: indexx         ! create sorted index array
USE nr_utility_module, ONLY: arth           ! build a vector of regularly spaced numbers
USE nr_utility_module, ONLY: sizeo          ! get size of allocatable array (if not allocated, zero)
USE nr_utility_module, ONLY: findIndex      ! find index within a vector
USE perf_mod,          ONLY: t_startf       ! timing start
USE perf_mod,          ONLY: t_stopf        ! timing stop

! MPI utility
USE mpi_mod, ONLY: shr_mpi_bcast
USE mpi_mod, ONLY: shr_mpi_gatherV
USE mpi_mod, ONLY: shr_mpi_scatterV
USE mpi_mod, ONLY: shr_mpi_allgather
USE mpi_mod, ONLY: shr_mpi_barrier

implicit none

! common parameters within this module
! commputation types
integer(i4b),parameter  :: scatter=1
integer(i4b),parameter  :: gather=2
! omp domain decomposition methods
integer(i4b),parameter  :: upstream_size=1
integer(i4b),parameter  :: stream_order=2
! sub domain types
integer(i4b),parameter  :: mainstem=1
integer(i4b),parameter  :: tributary=2
integer(i4b),parameter  :: endorheic=3

logical(lgt), parameter :: debug_mpi=.false.
logical(lgt), parameter :: debug_trib_omp=.false.
logical(lgt), parameter :: debug_main_omp=.false.

private

public :: comm_ntopo_data
public :: mpi_restart
public :: mpi_route
public :: pass_global_data

contains

 ! *********************************************************************
 ! public subroutine: send reach/hru information to tasks and populate data structures
 ! *********************************************************************
 subroutine comm_ntopo_data(pid,                & ! input: proc id
                            nNodes,             & ! input: number of procs
                            comm,               & ! input: communicator
                            nRch_in,            & ! input: number of stream segments in whole domain
                            nHRU_in,            & ! input: number of HRUs that are connected to reaches
                            structHRU,          & ! input: data structure for HRUs
                            structSEG,          & ! input: data structure for stream segments
                            structHRU2seg,      & ! input: data structure for mapping hru2basin
                            structNTOPO,        & ! input: data structure for network toopology
                            ierr,message)         ! output: error control

  USE globalData,          ONLY: ixPrint                  ! reach index to be checked by on-screen pringing
  USE globalData,          ONLY: domains                  ! domain data structure - for each domain, pfaf codes and list of segment indices
  USE globalData,          ONLY: nDomain                  ! count of decomposed domains (tributaries + mainstems)
  USE globalData,          ONLY: river_basin_main         ! OMP domain data structure for mainstem
  USE globalData,          ONLY: river_basin_trib         ! OMP domain data structure for tributaries
  USE globalData,          ONLY: RCHFLX_trib              ! Reach flux data structures (per proc, tributary)
  USE globalData,          ONLY: RCHSTA_trib              ! Reach state data structures (per proc, tributary)
  USE globalData,          ONLY: RCHFLX_main              ! Reach flux data structures (master proc, mainstem)
  USE globalData,          ONLY: RCHSTA_main              ! Reach state data structures (master proc, mainstem)
  USE globalData,          ONLY: NETOPO_trib              ! network topology data structure (per proc, tributary)
  USE globalData,          ONLY: RPARAM_trib              ! reach physical parameter data structure (per proc, tributary)
  USE globalData,          ONLY: NETOPO_main              ! network topology data structure (per proc, tributary)
  USE globalData,          ONLY: RPARAM_main              ! reach physical parameter data structure (per proc, tributary)
  USE globalData,          ONLY: nEns                     ! ensemble numbers (currently only 1)
  USE globalData,          ONLY: nRch_mainstem            ! number of mainstem reaches
  USE globalData,          ONLY: nHRU_mainstem            ! number of mainstem HRUs
  USE globalData,          ONLY: hru_per_proc             ! number of hrus assigned to each proc (size = num of procs+1)
  USE globalData,          ONLY: rch_per_proc             ! number of reaches assigned to each proc (size = num of procs+1)
  USE globalData,          ONLY: ixHRU_order              ! global HRU index in the order of proc assignment (size = total number of HRUs contributing to any reaches, nContribHRU)
  USE globalData,          ONLY: ixRch_order              ! global reach index in the order of proc assignment (size = total number of reaches in the entire network)
  USE globalData,          ONLY: tribOutlet_per_proc      ! number of tributary outlets per proc (size = nNodes)
  USE globalData,          ONLY: global_ix_main           ! reach index in mainstem array linked to reach outlets in tributary array (size = sum of tributary outlets)
  USE globalData,          ONLY: global_ix_comm           ! global reach index at tributary reach outlets to mainstem (size = sum of tributary outlets within entire network)
  USE globalData,          ONLY: local_ix_comm            ! local reach index at tributary reach outlets to mainstem (size = sum of tributary outlets within entire network)
  USE globalData,          ONLY: reachID
  USE globalData,          ONLY: basinID
  USE alloc_data,          ONLY: alloc_struct
  USE process_ntopo,       ONLY: augment_ntopo            ! compute all the additional network topology (only compute option = on)
  USE process_ntopo,       ONLY: put_data_struct          !
  USE domain_decomposition,ONLY: omp_domain_decomposition ! domain decomposition for omp

  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: pid                      ! process id (MPI)
  integer(i4b),                   intent(in)  :: nNodes                   ! number of processes (MPI)
  integer(i4b),                   intent(in)  :: comm                     ! communicator
  integer(i4b),                   intent(in)  :: nRch_in                  ! number of total reaches
  integer(i4b),                   intent(in)  :: nHRU_in                  ! number of total HRUs that are connected to reaches
  type(var_dlength), allocatable, intent(in)  :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(in)  :: structSEG(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(in)  :: structHRU2SEG(:)         ! HRU to SEG mapping
  type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)           ! network topology
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                   ! error message
  ! Local variables
  ! data structure for tributary river network data
  type(var_dlength), allocatable              :: structHRU_local(:)        ! tributary ancillary data for HRUs
  type(var_dlength), allocatable              :: structSEG_local(:)        ! tributary ancillary data for stream segments
  type(var_ilength), allocatable              :: structNTOPO_local(:)      ! tributary network topology
  type(var_ilength), allocatable              :: structHRU2seg_local(:)    ! tributary ancillary data for mapping hru2basin
  type(var_clength), allocatable              :: structPFAF_local(:)       ! tributary ancillary data for pfafstetter code
  ! data structure for mainstem river network data (root node)
  type(var_dlength), allocatable              :: structHRU_main(:)         ! mainstem ancillary data for HRUs
  type(var_dlength), allocatable              :: structSEG_main(:)         ! mainstem ancillary data for stream segments
  type(var_ilength), allocatable              :: structNTOPO_main(:)       ! mainstem network topology
  type(var_ilength), allocatable              :: structHRU2seg_main(:)     ! mainstem ancillary data for mapping hru2basin
  type(var_clength), allocatable              :: structPFAF_main(:)        ! mainstem ancillary data for pfafstetter code
  ! 1D array for decomposed river network data
  integer(i4b),      allocatable              :: segId_local(:)            ! reach id for decomposed network
  integer(i4b),      allocatable              :: downSegId_local(:)        ! downstream reach id for decomposed network
  integer(i4b),      allocatable              :: hruId_local(:)            ! hru id array in decomposed network
  integer(i4b),      allocatable              :: hruSegId_local(:)         ! downstream reach id array in decomposed network
  integer(i4b),      allocatable              :: islake_local(:)           ! islake flag local
  integer(i4b),      allocatable              :: LakeTargVol_local(:)      ! LakeTargVol flag local
  integer(i4b),      allocatable              :: LakeModelType_local(:)    ! LakeModelType flag local
  real(dp),          allocatable              :: slope_local(:)            ! reach slope array in decomposed network
  real(dp),          allocatable              :: length_local(:)           ! reach length array in decomposed network
  real(dp),          allocatable              :: area_local(:)             ! hru area in decomposed network
  real(dp),          allocatable              :: D03_MaxStorage_local(:)   ! Doll 2003; maximum active storage for Doll 2003
  real(dp),          allocatable              :: D03_Coefficient_local(:)  ! Doll 2003; coefficient for Doll 2003
  real(dp),          allocatable              :: D03_Power_local(:)        ! Doll 2003; power for Doll 2003
  real(dp),          allocatable              :: H06_Smax_local(:)         ! Hanasaki 2006; maximume reservoir storage [m3]
  real(dp),          allocatable              :: H06_alpha_local(:)        ! Hanasaki 2006; fraction of active storage compared to total storage [-]
  real(dp),          allocatable              :: H06_envfact_local(:)      ! Hanasaki 2006; fraction of inflow that can be used to meet demand [-]
  real(dp),          allocatable              :: H06_S_ini_local(:)        ! Hanasaki 2006; initial storage used for initial estimation of release coefficient [m3]
  real(dp),          allocatable              :: H06_c1_local(:)           ! Hanasaki 2006; coefficient 1 for target release for irrigation reseroir [-]
  real(dp),          allocatable              :: H06_c2_local(:)           ! Hanasaki 2006; coefficient 2 for target release for irrigation reseroir [-]
  real(dp),          allocatable              :: H06_exponent_local(:)     ! Hanasaki 2006; Exponenet of actual release for "within-a-year" reservoir [-]
  real(dp),          allocatable              :: H06_denominator_local(:)  ! Hanasaki 2006; Denominator of actual release for "within-a-year" reservoir [-]
  real(dp),          allocatable              :: H06_c_compare_local(:)    ! Hanasaki 2006; Criterion for distinguish of "within-a-year" or "multi-year" reservoir [-]
  real(dp),          allocatable              :: H06_frac_Sdead_local(:)   ! Hanasaki 2006; Fraction of dead storage to maximume storage [-]
  real(dp),          allocatable              :: H06_E_rel_ini_local(:)    ! Hanasaki 2006; Initial release coefficient [-]
  real(dp),          allocatable              :: H06_I_Jan_local(:)        ! Hanasaki 2006; Average January   inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Feb_local(:)        ! Hanasaki 2006; Average Februrary inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Mar_local(:)        ! Hanasaki 2006; Average March     inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Apr_local(:)        ! Hanasaki 2006; Average April     inflow [m3/s]
  real(dp),          allocatable              :: H06_I_May_local(:)        ! Hanasaki 2006; Average May       inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Jun_local(:)        ! Hanasaki 2006; Average June      inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Jul_local(:)        ! Hanasaki 2006; Average July      inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Aug_local(:)        ! Hanasaki 2006; Average August    inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Sep_local(:)        ! Hanasaki 2006; Average September inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Oct_local(:)        ! Hanasaki 2006; Average October   inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Nov_local(:)        ! Hanasaki 2006; Average November  inflow [m3/s]
  real(dp),          allocatable              :: H06_I_Dec_local(:)        ! Hanasaki 2006; Average December  inflow [m3/s]
  real(dp),          allocatable              :: H06_D_Jan_local(:)        ! Hanasaki 2006; Average January   demand [m3/s]
  real(dp),          allocatable              :: H06_D_Feb_local(:)        ! Hanasaki 2006; Average Februrary demand [m3/s]
  real(dp),          allocatable              :: H06_D_Mar_local(:)        ! Hanasaki 2006; Average March     demand [m3/s]
  real(dp),          allocatable              :: H06_D_Apr_local(:)        ! Hanasaki 2006; Average April     demand [m3/s]
  real(dp),          allocatable              :: H06_D_May_local(:)        ! Hanasaki 2006; Average May       demand [m3/s]
  real(dp),          allocatable              :: H06_D_Jun_local(:)        ! Hanasaki 2006; Average June      demand [m3/s]
  real(dp),          allocatable              :: H06_D_Jul_local(:)        ! Hanasaki 2006; Average July      demand [m3/s]
  real(dp),          allocatable              :: H06_D_Aug_local(:)        ! Hanasaki 2006; Average Agust     demand [m3/s]
  real(dp),          allocatable              :: H06_D_Sep_local(:)        ! Hanasaki 2006; Average September demand [m3/s]
  real(dp),          allocatable              :: H06_D_Oct_local(:)        ! Hanasaki 2006; Average October   demand [m3/s]
  real(dp),          allocatable              :: H06_D_Nov_local(:)        ! Hanasaki 2006; Average Novermber demand [m3/s]
  real(dp),          allocatable              :: H06_D_Dec_local(:)        ! Hanasaki 2006; Average December  demand [m3/s]
  integer(i4b),      allocatable              :: H06_purpose_local(:)      ! Hanasaki 2006; reservoir purpose; (0= non-irrigation, 1=irrigation) [-]
  integer(i4b),      allocatable              :: H06_I_mem_F_local(:)      ! Hanasaki 2006; Flag to transition to modelled inflow [-]
  integer(i4b),      allocatable              :: H06_D_mem_F_local(:)      ! Hanasaki 2006; Flag to transition to modelled/provided demand [-]
  integer(i4b),      allocatable              :: H06_I_mem_L_local(:)      ! Hanasaki 2006; Memory length in years for inflow [year]
  integer(i4b),      allocatable              :: H06_D_mem_L_local(:)      ! Hanasaki 2006; Memory length in years for demand [year]



  logical(lgt),      allocatable              :: tribOutlet_local(:)       ! logical to indicate tributary outlet to mainstems
  ! 1D array for the entire river network
  integer(i4b)                                :: hruId(nHRU_in)            ! hru id for all the HRUs
  integer(i4b)                                :: hruSegId(nHRU_in)         ! hru-to-seg mapping for each hru
  integer(i4b)                                :: segId(nRch_in)            ! reach id for all the segments
  integer(i4b)                                :: downSegId(nRch_in)        ! downstream reach ID for each reach
  integer(i4b)                                :: islake(nRch_in)           ! islake flag
  integer(i4b)                                :: LakeTargVol(nRch_in)      ! LakeTargVol flag
  integer(i4b)                                :: LakeModelType(nRch_in)    ! LakeModelType flag local
  real(dp)                                    :: slope(nRch_in)            ! reach slope array for each reach
  real(dp)                                    :: length(nRch_in)           ! reach length array for each reach
  real(dp)                                    :: area(nHRU_in)             ! hru area for each hru
  real(dp)                                    :: D03_MaxStorage(nRch_in)   ! Doll 2003; maximum active storage for Doll 2003
  real(dp)                                    :: D03_Coefficient(nRch_in)  ! Doll 2003; coefficient for Doll 2003
  real(dp)                                    :: D03_Power(nRch_in)        ! Doll 2003; power for Doll 2003
  real(dp)                                    :: H06_Smax(nRch_in)         ! Hanasaki 2006; maximume reservoir storage [m3]
  real(dp)                                    :: H06_alpha(nRch_in)        ! Hanasaki 2006; fraction of active storage compared to total storage [-]
  real(dp)                                    :: H06_envfact(nRch_in)      ! Hanasaki 2006; fraction of inflow that can be used to meet demand [-]
  real(dp)                                    :: H06_S_ini(nRch_in)        ! Hanasaki 2006; initial storage used for initial estimation of release coefficient [m3]
  real(dp)                                    :: H06_c1(nRch_in)           ! Hanasaki 2006; coefficient 1 for target release for irrigation reseroir [-]
  real(dp)                                    :: H06_c2(nRch_in)           ! Hanasaki 2006; coefficient 2 for target release for irrigation reseroir [-]
  real(dp)                                    :: H06_exponent(nRch_in)     ! Hanasaki 2006; Exponenet of actual release for "within-a-year" reservoir [-]
  real(dp)                                    :: H06_denominator(nRch_in)  ! Hanasaki 2006; Denominator of actual release for "within-a-year" reservoir [-]
  real(dp)                                    :: H06_c_compare(nRch_in)    ! Hanasaki 2006; Criterion for distinguish of "within-a-year" or "multi-year" reservoir [-]
  real(dp)                                    :: H06_frac_Sdead(nRch_in)   ! Hanasaki 2006; Fraction of dead storage to maximume storage [-]
  real(dp)                                    :: H06_E_rel_ini(nRch_in)    ! Hanasaki 2006; Initial release coefficient [-]
  real(dp)                                    :: H06_I_Jan(nRch_in)        ! Hanasaki 2006; Average January   inflow [m3/s]
  real(dp)                                    :: H06_I_Feb(nRch_in)        ! Hanasaki 2006; Average Februrary inflow [m3/s]
  real(dp)                                    :: H06_I_Mar(nRch_in)        ! Hanasaki 2006; Average March     inflow [m3/s]
  real(dp)                                    :: H06_I_Apr(nRch_in)        ! Hanasaki 2006; Average April     inflow [m3/s]
  real(dp)                                    :: H06_I_May(nRch_in)        ! Hanasaki 2006; Average May       inflow [m3/s]
  real(dp)                                    :: H06_I_Jun(nRch_in)        ! Hanasaki 2006; Average June      inflow [m3/s]
  real(dp)                                    :: H06_I_Jul(nRch_in)        ! Hanasaki 2006; Average July      inflow [m3/s]
  real(dp)                                    :: H06_I_Aug(nRch_in)        ! Hanasaki 2006; Average August    inflow [m3/s]
  real(dp)                                    :: H06_I_Sep(nRch_in)        ! Hanasaki 2006; Average September inflow [m3/s]
  real(dp)                                    :: H06_I_Oct(nRch_in)        ! Hanasaki 2006; Average October   inflow [m3/s]
  real(dp)                                    :: H06_I_Nov(nRch_in)        ! Hanasaki 2006; Average November  inflow [m3/s]
  real(dp)                                    :: H06_I_Dec(nRch_in)        ! Hanasaki 2006; Average December  inflow [m3/s]
  real(dp)                                    :: H06_D_Jan(nRch_in)        ! Hanasaki 2006; Average January   demand [m3/s]
  real(dp)                                    :: H06_D_Feb(nRch_in)        ! Hanasaki 2006; Average Februrary demand [m3/s]
  real(dp)                                    :: H06_D_Mar(nRch_in)        ! Hanasaki 2006; Average March     demand [m3/s]
  real(dp)                                    :: H06_D_Apr(nRch_in)        ! Hanasaki 2006; Average April     demand [m3/s]
  real(dp)                                    :: H06_D_May(nRch_in)        ! Hanasaki 2006; Average May       demand [m3/s]
  real(dp)                                    :: H06_D_Jun(nRch_in)        ! Hanasaki 2006; Average June      demand [m3/s]
  real(dp)                                    :: H06_D_Jul(nRch_in)        ! Hanasaki 2006; Average July      demand [m3/s]
  real(dp)                                    :: H06_D_Aug(nRch_in)        ! Hanasaki 2006; Average Agust     demand [m3/s]
  real(dp)                                    :: H06_D_Sep(nRch_in)        ! Hanasaki 2006; Average September demand [m3/s]
  real(dp)                                    :: H06_D_Oct(nRch_in)        ! Hanasaki 2006; Average October   demand [m3/s]
  real(dp)                                    :: H06_D_Nov(nRch_in)        ! Hanasaki 2006; Average Novermber demand [m3/s]
  real(dp)                                    :: H06_D_Dec(nRch_in)        ! Hanasaki 2006; Average December  demand [m3/s]
  integer(i4b)                                :: H06_purpose(nRch_in)      ! Hanasaki 2006; reservoir purpose; (0= non-irrigation, 1=irrigation) [-]
  integer(i4b)                                :: H06_I_mem_F(nRch_in)      ! Hanasaki 2006; Flag to transition to modelled inflow [-]
  integer(i4b)                                :: H06_D_mem_F(nRch_in)      ! Hanasaki 2006; Flag to transition to modelled/provided demand [-]
  integer(i4b)                                :: H06_I_mem_L(nRch_in)      ! Hanasaki 2006; Memory length in years for inflow [year]
  integer(i4b)                                :: H06_D_mem_L(nRch_in)      ! Hanasaki 2006; Memory length in years for demand [year]
  integer(i4b)                                :: ixNode(nRch_in)           ! node assignment for each reach
  integer(i4b)                                :: ixDomain(nRch_in)         ! domain index for each reach
  logical(lgt),      allocatable              :: tribOutlet(:)             ! logical to indicate tributary outlet to mainstems over entire network
  integer(i4b)                                :: nTribOutlet               ! number of tributary outlet connected to mainstem
  integer(i4b)                                :: nRch_trib                 ! number of tributary outlets for each proc (scalar)
  integer(i4b),     allocatable               :: nRch_trib_array(:)        ! number of tributary outlets for each proc (array)
  ! flat array for decomposed river network per domain (sub-basin)
  integer(i4b)                                :: idNode(nDomain)           ! node id array for each domain
  integer(i4b)                                :: rnkIdNode(nDomain)        ! ranked node id array for each domain
  integer(i4b)                                :: jHRU,jSeg                 ! ranked indices
  integer(i4b)                                :: iUps                      ! immediate upstream loop indices
  integer(i4b),     allocatable               :: jUps(:)                   ! immediate upstream loop indices
  ! miscellaneous
  integer(i4b), allocatable                   :: seq_array(:)
  integer(i4b)                                :: nUps
  integer(i4b)                                :: iSeg,iHru                 ! reach and hru loop indices
  integer(i4b)                                :: ix1, ix2
  integer(i4b)                                :: ix,ixx                    ! loop indices
  integer(i4b)                                :: ixSeg1,ixSeg2             ! starting index and ending index, respectively, for reach array
  integer(i4b)                                :: ixHru1,ixHru2             ! starting index and ending index, respectively, for HRU array
  integer(i4b)                                :: idx                       ! node indix (-1, 0, 1, ... , nNodes-1)
  character(len=strLen)                       :: cmessage                  ! error message from subroutine

  ierr=0; message='comm_ntopo_data/'

  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! Part 1: define routing vectors ordered by domain/node
  !  - define the global indices ordered by domain/node
  !  - define the number of reaches/hrus on each processor
  !  - copy the data from the data structures to the ordered routing vectors
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************

  if (masterproc) then ! this is a root process

    ! allocate local and global indices
    allocate(rch_per_proc(-1:nNodes-1), hru_per_proc(-1:nNodes-1), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [rch_per_proc, hru_per_proc]'; return; endif

    allocate(ixHRU_order(nHRU_in),ixRch_order(nRch_in), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ixHRU_order, ixRch_order]'; return; endif

    ! Create segIndex array from domains derived type. The array is sorted from node 0 through nNodes-1
    ! SegIndex Array needs to be contiguous when a chunk is sent to computing node (use sort function...)
    ! start with mainstem domain assigned to root node

    ! domain is a contiguous collection of reaches/HRUs -- multiple domains may be on a single processor

    do ix=1,nDomain
      idNode(ix) = domains(ix)%idNode ! extracts the processing node from the "domain" data structire
    enddo
    call indexx(idNode,rnkIdNode) ! rank the processor nodes

    ! loop through the domains
    ixSeg2=0; ixHru2=0 ! last indices of domain chunks
    domain:do ix = 1, nDomain

     ! get the number of stream segments and HRUs in each domain
     ixx = rnkIdNode(ix)
     associate (nSubSeg => size(domains(ixx)%segIndex), nSubHru => size(domains(ixx)%hruIndex) )

     ! define reach global index array in order of node assignment
     if (domains(ixx)%basinType /= endorheic ) then   ! endorheic domain does not have reaches
       ixSeg1 = ixSeg2+1
       ixSeg2 = ixSeg1+nSubSeg-1
       ixRch_order(ixSeg1:ixSeg2) = domains(ixx)%segIndex(1:nSubSeg)
       ! extra information (debugging)
       ixNode(ixSeg1:ixSeg2)      = domains(ixx)%idNode
       ixDomain(ixSeg1:ixSeg2)    = ixx
     end if

     ! define hru index array in order of node assignment
     if (nSubHru>0) then
       ixHru1 = ixHru2+1
       ixHru2 = ixHru1+nSubHru-1
       ixHRU_order(ixHru1:ixHru2) = domains(ixx)%hruIndex(1:nSubHru) ! global hru index per node
     end if

     end associate
    end do domain

    ! Count the number of reaches and hrus in each node
    ! index of seg_per_proc and hru_per_proc: -1 -> mainstem, 0 -> small tributaries, 1 through nNodes-1 -> large tributaries
    rch_per_proc = 0
    hru_per_proc = 0
    do ix = 1,nDomain
     idx = domains(ix)%idNode
     rch_per_proc(idx) = rch_per_proc(idx) + sizeo(domains(ix)%segIndex)
     hru_per_proc(idx) = hru_per_proc(idx) + sizeo(domains(ix)%hruIndex)
    end do

    ! define routing vectors ordered by domain/node
    ! reach array
    do iSeg = 1,nRch_in
     jSeg = ixRch_order(iSeg) ! global index, ordered by domain/node
     segId(iSeg)           = structNTOPO(jSeg)%var(ixNTOPO%segId)%dat(1)
     downSegId(iSeg)       = structNTOPO(jSeg)%var(ixNTOPO%downSegId)%dat(1)
     slope(iSeg)           = structSEG(  jSeg)%var(ixSEG%slope)%dat(1)
     length(iSeg)          = structSEG(  jSeg)%var(ixSEG%length)%dat(1)
     if (is_lake_sim) then
       islake(iSeg)          = structNTOPO(jSeg)%var(ixNTOPO%islake)%dat(1)
       LakeTargVol(iSeg)     = structNTOPO(jSeg)%var(ixNTOPO%LakeTargVol)%dat(1)
       LakeModelType(iSeg)   = structNTOPO(jSeg)%var(ixNTOPO%LakeModelType)%dat(1)
       D03_MaxStorage(iSeg)  = structSEG(  jSeg)%var(ixSEG%D03_MaxStorage)%dat(1)
       D03_Coefficient(iSeg) = structSEG(  jSeg)%var(ixSEG%D03_Coefficient)%dat(1)
       D03_Power(iSeg)       = structSEG(  jSeg)%var(ixSEG%D03_Power)%dat(1)
       H06_Smax(iSeg)        = structSEG(  jSeg)%var(ixSEG%H06_Smax)%dat(1)
       H06_alpha(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_alpha)%dat(1)
       H06_envfact(iSeg)     = structSEG(  jSeg)%var(ixSEG%H06_envfact)%dat(1)
       H06_S_ini(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_S_ini)%dat(1)
       H06_c1(iSeg)          = structSEG(  jSeg)%var(ixSEG%H06_c1)%dat(1)
       H06_c2(iSeg)          = structSEG(  jSeg)%var(ixSEG%H06_c2)%dat(1)
       H06_exponent(iSeg)    = structSEG(  jSeg)%var(ixSEG%H06_exponent)%dat(1)
       H06_denominator(iSeg) = structSEG(  jSeg)%var(ixSEG%H06_denominator)%dat(1)
       H06_c_compare(iSeg)   = structSEG(  jSeg)%var(ixSEG%H06_c_compare)%dat(1)
       H06_frac_Sdead(iSeg)  = structSEG(  jSeg)%var(ixSEG%H06_frac_Sdead)%dat(1)
       H06_E_rel_ini(iSeg)   = structSEG(  jSeg)%var(ixSEG%H06_E_rel_ini)%dat(1)
       H06_I_Jan(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Jan)%dat(1)
       H06_I_Feb(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Feb)%dat(1)
       H06_I_Mar(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Mar)%dat(1)
       H06_I_Apr(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Apr)%dat(1)
       H06_I_May(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_May)%dat(1)
       H06_I_Jun(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Jun)%dat(1)
       H06_I_Jul(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Jul)%dat(1)
       H06_I_Aug(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Aug)%dat(1)
       H06_I_Sep(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Sep)%dat(1)
       H06_I_Oct(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Oct)%dat(1)
       H06_I_Nov(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Nov)%dat(1)
       H06_I_Dec(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_I_Dec)%dat(1)
       H06_D_Jan(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Jan)%dat(1)
       H06_D_Feb(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Feb)%dat(1)
       H06_D_Mar(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Mar)%dat(1)
       H06_D_Apr(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Apr)%dat(1)
       H06_D_May(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_May)%dat(1)
       H06_D_Jun(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Jun)%dat(1)
       H06_D_Jul(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Jul)%dat(1)
       H06_D_Aug(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Aug)%dat(1)
       H06_D_Sep(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Sep)%dat(1)
       H06_D_Oct(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Oct)%dat(1)
       H06_D_Nov(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Nov)%dat(1)
       H06_D_Dec(iSeg)       = structSEG(  jSeg)%var(ixSEG%H06_D_Dec)%dat(1)
       H06_purpose(iSeg)     = structSEG(  jSeg)%var(ixSEG%H06_purpose)%dat(1)
       H06_I_mem_F(iSeg)     = structSEG(  jSeg)%var(ixSEG%H06_I_mem_F)%dat(1)
       H06_D_mem_F(iSeg)     = structSEG(  jSeg)%var(ixSEG%H06_D_mem_F)%dat(1)
       H06_I_mem_L(iSeg)     = structSEG(  jSeg)%var(ixSEG%H06_I_mem_L)%dat(1)
       H06_D_mem_L(iSeg)     = structSEG(  jSeg)%var(ixSEG%H06_D_mem_L)%dat(1)
     end if
    end do

    ! hru array
    do iHru = 1,nHRU_in
      jHRU = ixHRU_order(iHru)  ! global index, ordered by domain/node
      hruId(iHru)    = structHRU2SEG(jHRU)%var(ixHRU2SEG%HRUid)%dat(1)
      hruSegId(iHru) = structHRU2SEG(jHRU)%var(ixHRU2SEG%hruSegId)%dat(1)
      area(iHru)     = structHRU(    jHRU)%var(ixHRU%area)%dat(1)
    enddo

    ! Reorder reachID and basinID for output to match up with order of RCHFLX/RCHSTA reach order and basinRunoff hru order
    reachID(1:nRch_in) = reachID( ixRch_order )
    basinID(1:nHRU_in) = basinID( ixHRU_order )

    if (debug_mpi) then
      write(iulog,'(a)') 'ix, segId, ixRch_order, domain-index, proc-id'
      do ix = 1,nRch_in
        write(iulog,*) ix, segId(ix), ixRch_order(ix), ixDomain(ix), ixNode(ix)
      enddo
    endif

  endif  ! ( masterproc )

  call shr_mpi_barrier(comm, message)

  ! sends the number of reaches/hrus per proc to all processors
  call shr_mpi_bcast(rch_per_proc, ierr, cmessage)
  call shr_mpi_bcast(hru_per_proc, ierr, cmessage)

  call shr_mpi_bcast(ixRch_order, ierr, cmessage)
  call shr_mpi_bcast(ixHRU_order, ierr, cmessage)

  ! define the number of reaches/hrus on the mainstem
  nRch_mainstem = rch_per_proc(-1)
  nHRU_mainstem = hru_per_proc(-1)

  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! Optional procedures: Multiple MPI tasks requested
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  if (nNodes>1) then

    ! -----------------------------------------------------------------------------
    !  Send the information for tributaries to individual processors
    ! -----------------------------------------------------------------------------
    call shr_mpi_scatterV(segId                 (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), segId_local,          ierr, cmessage)
    call shr_mpi_scatterV(downSegId             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), downSegId_local,      ierr, cmessage)
    call shr_mpi_scatterV(slope                 (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), slope_local,          ierr, cmessage)
    call shr_mpi_scatterV(length                (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), length_local,         ierr, cmessage)
    if (is_lake_sim) then
      call shr_mpi_scatterV(islake                (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), islake_local,         ierr, cmessage)
      call shr_mpi_scatterV(LakeTargVol           (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), LakeTargVol_local,    ierr, cmessage)
      call shr_mpi_scatterV(LakeModelType         (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), LakeModelType_local,  ierr, cmessage)
      call shr_mpi_scatterV(D03_MaxStorage        (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), D03_MaxStorage_local, ierr, cmessage)
      call shr_mpi_scatterV(D03_Coefficient       (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), D03_Coefficient_local,ierr, cmessage)
      call shr_mpi_scatterV(D03_Power             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), D03_Power_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_Smax              (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_Smax_local,       ierr, cmessage)
      call shr_mpi_scatterV(H06_alpha             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_alpha_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_envfact           (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_envfact_local,    ierr, cmessage)
      call shr_mpi_scatterV(H06_S_ini             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_S_ini_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_c1                (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_c1_local,         ierr, cmessage)
      call shr_mpi_scatterV(H06_c2                (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_c2_local,         ierr, cmessage)
      call shr_mpi_scatterV(H06_exponent          (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_exponent_local,   ierr, cmessage)
      call shr_mpi_scatterV(H06_denominator       (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_denominator_local,ierr, cmessage)
      call shr_mpi_scatterV(H06_c_compare         (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_c_compare_local,  ierr, cmessage)
      call shr_mpi_scatterV(H06_frac_Sdead        (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_frac_Sdead_local, ierr, cmessage)
      call shr_mpi_scatterV(H06_E_rel_ini         (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_E_rel_ini_local,  ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Jan             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Jan_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Feb             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Feb_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Mar             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Mar_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Apr             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Apr_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_May             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_May_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Jun             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Jun_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Jul             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Jul_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Aug             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Aug_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Sep             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Sep_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Oct             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Oct_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Nov             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Nov_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_Dec             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_Dec_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Jan             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Jan_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Feb             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Feb_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Mar             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Mar_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Apr             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Apr_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_May             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_May_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Jun             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Jun_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Jul             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Jul_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Aug             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Aug_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Sep             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Sep_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Oct             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Oct_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Nov             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Nov_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_Dec             (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_Dec_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_purpose           (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_purpose_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_mem_F           (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_mem_F_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_mem_F           (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_mem_F_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_I_mem_L           (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_I_mem_L_local,      ierr, cmessage)
      call shr_mpi_scatterV(H06_D_mem_L           (nRch_mainstem+1:nRch_in), rch_per_proc(0:nNodes-1), H06_D_mem_L_local,      ierr, cmessage)
    end if
    call shr_mpi_scatterV(hruId         (nHRU_mainstem+1:nHRU_in), hru_per_proc(0:nNodes-1), hruId_local,         ierr, cmessage)
    call shr_mpi_scatterV(hruSegId      (nHRU_mainstem+1:nHRU_in), hru_per_proc(0:nNodes-1), hruSegId_local,      ierr, cmessage)
    call shr_mpi_scatterV(area          (nHRU_mainstem+1:nHRU_in), hru_per_proc(0:nNodes-1), area_local,          ierr, cmessage)

    ! -----------------------------------------------------------------------------
    !  populate local (tributary) data structures and compute additional ancillary information
    ! -----------------------------------------------------------------------------

    allocate(RCHFLX_trib(nEns,rch_per_proc(pid)), RCHSTA_trib(nEns,rch_per_proc(pid)), stat=ierr)

    call alloc_struct(hru_per_proc(pid),     & ! input: number of HRUs
                      rch_per_proc(pid),     & ! input: number of stream segments
                      structHRU_local,       & ! inout: ancillary data for HRUs
                      structSEG_local,       & ! inout: ancillary data for stream segments
                      structHRU2seg_local,   & ! inout: ancillary data for mapping hru2basin
                      structNTOPO_local,     & ! inout: ancillary data for network toopology
                      structPFAF_local,      & ! inout: ancillary data for pfafstetter code
                      ierr,cmessage)           ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Populate local data structures needed for tributary network augumentation
    reach: do ix = 1,rch_per_proc(pid)
     structNTOPO_local(ix)%var(ixNTOPO%segId)%dat(1)          = segId_local(ix)
     structNTOPO_local(ix)%var(ixNTOPO%downSegId)%dat(1)      = downSegId_local(ix)
     structSEG_local  (ix)%var(ixSEG%length)%dat(1)           = length_local(ix)
     structSEG_local  (ix)%var(ixSEG%slope)%dat(1)            = slope_local(ix)
     if (is_lake_sim) then
       structNTOPO_local(ix)%var(ixNTOPO%islake)%dat(1)         = islake_local(ix)
       structNTOPO_local(ix)%var(ixNTOPO%LakeTargVol)%dat(1)    = LakeTargVol_local(ix)
       structNTOPO_local(ix)%var(ixNTOPO%LakeModelType)%dat(1)  = LakeModelType_local(ix)
       structSEG_local  (ix)%var(ixSEG%D03_MaxStorage)%dat(1)   = D03_MaxStorage_local(ix)
       structSEG_local  (ix)%var(ixSEG%D03_Coefficient)%dat(1)  = D03_Coefficient_local(ix)
       structSEG_local  (ix)%var(ixSEG%D03_Power)%dat(1)        = D03_Power_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_Smax)%dat(1)         = H06_Smax_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_alpha)%dat(1)        = H06_alpha_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_envfact)%dat(1)      = H06_envfact_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_S_ini)%dat(1)        = H06_S_ini_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_c1)%dat(1)           = H06_c1_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_c2)%dat(1)           = H06_c2_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_exponent)%dat(1)     = H06_exponent_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_denominator)%dat(1)  = H06_denominator_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_c_compare)%dat(1)    = H06_c_compare_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_frac_Sdead)%dat(1)   = H06_frac_Sdead_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_E_rel_ini)%dat(1)    = H06_E_rel_ini_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Jan)%dat(1)        = H06_I_Jan_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Feb)%dat(1)        = H06_I_Feb_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Mar)%dat(1)        = H06_I_Mar_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Apr)%dat(1)        = H06_I_Apr_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_May)%dat(1)        = H06_I_May_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Jun)%dat(1)        = H06_I_Jun_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Jul)%dat(1)        = H06_I_Jul_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Aug)%dat(1)        = H06_I_Aug_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Sep)%dat(1)        = H06_I_Sep_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Oct)%dat(1)        = H06_I_Oct_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Nov)%dat(1)        = H06_I_Nov_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_Dec)%dat(1)        = H06_I_Dec_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Jan)%dat(1)        = H06_D_Jan_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Feb)%dat(1)        = H06_D_Feb_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Mar)%dat(1)        = H06_D_Mar_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Apr)%dat(1)        = H06_D_Apr_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_May)%dat(1)        = H06_D_May_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Jun)%dat(1)        = H06_D_Jun_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Jul)%dat(1)        = H06_D_Jul_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Aug)%dat(1)        = H06_D_Aug_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Sep)%dat(1)        = H06_D_Sep_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Oct)%dat(1)        = H06_D_Oct_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Nov)%dat(1)        = H06_D_Nov_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_Dec)%dat(1)        = H06_D_Dec_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_purpose)%dat(1)      = H06_purpose_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_mem_F)%dat(1)      = H06_I_mem_F_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_mem_F)%dat(1)      = H06_D_mem_F_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_I_mem_L)%dat(1)      = H06_I_mem_L_local(ix)
       structSEG_local  (ix)%var(ixSEG%H06_D_mem_L)%dat(1)      = H06_D_mem_L_local(ix)
     end if
    end do reach

    hru: do ix=1,hru_per_proc(pid)
     structHRU2SEG_local(ix)%var(ixHRU2SEG%HRUid)%dat(1)    = hruId_local(ix)
     structHRU2SEG_local(ix)%var(ixHRU2SEG%hruSegId)%dat(1) = hruSegId_local(ix)
     structHRU_local    (ix)%var(ixHRU%area)%dat(1)         = area_local(ix)
    end do hru

    ! find index of desired reach
    if (desireId/=integerMissing) ixPrint(2) = findIndex(segId_local, desireId, integerMissing)

    ! compute additional ancillary infomration
    call augment_ntopo(hru_per_proc(pid),            & ! input: number of HRUs
                       rch_per_proc(pid),            & ! input: number of stream segments
                       structHRU_local,              & ! inout: ancillary data for HRUs
                       structSEG_local,              & ! inout: ancillary data for stream segments
                       structHRU2seg_local,          & ! inout: ancillary data for mapping hru2basin
                       structNTOPO_local,            & ! inout: ancillary data for network toopology
                       ierr, cmessage)                 ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! copy data to routing structres RPARAM_trib and NETOPO_trib
    call put_data_struct(rch_per_proc(pid), structSEG_local, structNTOPO_local, & ! input
                         RPARAM_trib, NETOPO_trib,                              & ! output:
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! OMP domain decomposition
    call omp_domain_decomposition(upstream_size, rch_per_proc(pid), structNTOPO_local, river_basin_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    if (debug_trib_omp .and. pid==3) then
      write(iulog,'(a)') 'segid, branch, order'
      do ix =1,size(river_basin_trib)
        do ixx = 1, size(river_basin_trib(ix)%branch)
          do iSeg = 1, river_basin_trib(ix)%branch(ixx)%nRch
            associate (idx_tmp => river_basin_trib(ix)%branch(ixx)%segIndex(iSeg))
            write(iulog,"(I15,A,I9,A,I9)") structNTOPO_local(idx_tmp)%var(ixNTOPO%segId)%dat(1),',',ixx,',',ix
            end associate
          enddo
        enddo
      enddo
    endif

    ! -----------------------------------------------------------------------------
    ! Find "dangling reach/hru", or tributary outlet reaches/hrus that link to mainstems
    ! -----------------------------------------------------------------------------
    allocate(tribOutlet_local(rch_per_proc(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [tribOutlet_local]'; return; endif
    tribOutlet_local = .false.
    do ix=1,rch_per_proc(pid)
      if (structNTOPO_local(ix)%var(ixNTOPO%downSegIndex)%dat(1) == -1 .and. &
          structNTOPO_local(ix)%var(ixNTOPO%downSegId)%dat(1) > 0) then  !index == -1 but there is donwstream id actually
        tribOutlet_local(ix) = .true.
      endif
    enddo

    ! gather array for tributary outlet reaches per each proc
    call shr_mpi_allgather(tribOutlet_local, rch_per_proc(root:nNodes-1), tribOutlet, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    nTribOutlet = count(tribOutlet)

    ! gather array for number of tributary outlet reaches from each proc
    nRch_trib = count(tribOutlet_local) ! number of tributary outlets for each proc
    call shr_mpi_allgather(nRch_trib, 1, nRch_trib_array, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! need to use 0:nNodes-1 bound instead of 1:nNodes
    allocate(tribOutlet_per_proc(0:nNodes-1), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [tribOutlet_per_proc]'; return; endif
    tribOutlet_per_proc(0:nNodes-1) = nRch_trib_array

    allocate(global_ix_comm(nTribOutlet), global_ix_main(nTribOutlet), &
             seq_array(rch_per_proc(pid)), local_ix_comm(nRch_trib), stat=ierr,  errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! mask non-tributary outlet reaches
    ! mainstem part (1:nRch_mainstem) has to be removed from ixRch_order
    ix1 = nRch_mainstem+1; ix2 = sum(rch_per_proc)
    global_ix_comm = pack(ixRch_order(ix1:ix2), tribOutlet)   ! size = number of tributary outlet reaches from all the procs

    seq_array = arth(1,1,rch_per_proc(pid))
    local_ix_comm = pack(seq_array, tribOutlet_local) ! size = number of tributary outlet reaches depending on proc

!       if (pid==0) then
!       print*, 'ix, local_ix_comm, NETOPO_trib%REACHIX, reachId, downstream ID, downstream local ix'
!       do ix =1, nRch_trib
!        print*,  ix, local_ix_comm(ix), NETOPO_trib(local_ix_comm(ix))%REACHIX, NETOPO_trib(local_ix_comm(ix))%REACHID, NETOPO_trib(local_ix_comm(ix))%DREACHK, NETOPO_trib(local_ix_comm(ix))%DREACHI
!       enddo
!       do ix =1, size(global_ix_comm)
!        print*,  ix, global_ix_comm(ix), NETOPO(global_ix_comm(ix))%REACHIX, NETOPO(global_ix_comm(ix))%REACHID, NETOPO(global_ix_comm(ix))%DREACHK, NETOPO(global_ix_comm(ix))%DREACHI
!       enddo
!       do ix =1, size(NETOPO_trib)
!       print*, NETOPO_trib(ix)%REACHID, NETOPO_trib(ix)%DREACHK, NETOPO_trib(ix)%DREACHI, tribOutlet(ix)
!       enddo
!       endif

  endif ! (nNodes>1)

  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! Optional procedures:  Mainstem decomposition for OMP if mainstem exists as a result of MPI domain docomposition
  ! ********************************************************************************************************************
  ! ********************************************************************************************************************
  ! Specifics:
  ! Build NTOPO, RPARAM, RCHFLX, RCHSTA (used for routing) data structures (only master proc)
  ! NEED TO TAKE INTO ACCOUNT FOR upstream Area, width, goodBasin

  if (nRch_mainstem > 0) then
   if (masterproc) then

     allocate(RCHFLX_main(nEns, nRch_mainstem+nTribOutlet), RCHSTA_main(nEns, nRch_mainstem+nTribOutlet), stat=ierr, errmsg=cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     call alloc_struct(nHRU_mainstem,              & ! input: number of HRUs
                       nRch_mainstem+nTribOutlet,  & ! input: number of stream segments
                       structHRU_main,             & ! inout: ancillary data for HRUs
                       structSEG_main,             & ! inout: ancillary data for stream segments
                       structHRU2seg_main,         & ! inout: ancillary data for mapping hru2basin
                       structNTOPO_main,           & ! inout: ancillary data for network toopology
                       structPFAF_main,            & ! inout: ancillary data for pfafstetter code
                       ierr,cmessage)                ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! Populate mainstem data structures
     main_rch: do ix = 1, nRch_mainstem
       structNTOPO_main(ix)%var(ixNTOPO%segId)%dat(1)         = segId(ix)
       structNTOPO_main(ix)%var(ixNTOPO%downSegId)%dat(1)     = downSegId(ix)
       structSEG_main  (ix)%var(ixSEG%length)%dat(1)          = length(ix)
       structSEG_main  (ix)%var(ixSEG%slope)%dat(1)           = slope(ix)
       if (is_lake_sim) then
         structNTOPO_main(ix)%var(ixNTOPO%islake)%dat(1)        = islake(ix)
         structNTOPO_main(ix)%var(ixNTOPO%LakeTargVol)%dat(1)   = LakeTargVol(ix)
         structNTOPO_main(ix)%var(ixNTOPO%LakeModelType)%dat(1) = LakeModelType(ix)
         structSEG_main  (ix)%var(ixSEG%D03_MaxStorage)%dat(1)  = D03_MaxStorage(ix)
         structSEG_main  (ix)%var(ixSEG%D03_Coefficient)%dat(1) = D03_Coefficient(ix)
         structSEG_main  (ix)%var(ixSEG%D03_Power)%dat(1)       = D03_Power(ix)
         structSEG_main  (ix)%var(ixSEG%H06_Smax)%dat(1)        = H06_Smax(ix)
         structSEG_main  (ix)%var(ixSEG%H06_alpha)%dat(1)       = H06_alpha(ix)
         structSEG_main  (ix)%var(ixSEG%H06_envfact)%dat(1)     = H06_envfact(ix)
         structSEG_main  (ix)%var(ixSEG%H06_S_ini)%dat(1)       = H06_S_ini(ix)
         structSEG_main  (ix)%var(ixSEG%H06_c1)%dat(1)          = H06_c1(ix)
         structSEG_main  (ix)%var(ixSEG%H06_c2)%dat(1)          = H06_c2(ix)
         structSEG_main  (ix)%var(ixSEG%H06_exponent)%dat(1)    = H06_exponent(ix)
         structSEG_main  (ix)%var(ixSEG%H06_denominator)%dat(1) = H06_denominator(ix)
         structSEG_main  (ix)%var(ixSEG%H06_c_compare)%dat(1)   = H06_c_compare(ix)
         structSEG_main  (ix)%var(ixSEG%H06_frac_Sdead)%dat(1)  = H06_frac_Sdead(ix)
         structSEG_main  (ix)%var(ixSEG%H06_E_rel_ini)%dat(1)   = H06_E_rel_ini(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Jan)%dat(1)       = H06_I_Jan(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Feb)%dat(1)       = H06_I_Feb(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Mar)%dat(1)       = H06_I_Mar(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Apr)%dat(1)       = H06_I_Apr(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_May)%dat(1)       = H06_I_May(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Jun)%dat(1)       = H06_I_Jun(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Jul)%dat(1)       = H06_I_Jul(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Aug)%dat(1)       = H06_I_Aug(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Sep)%dat(1)       = H06_I_Sep(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Oct)%dat(1)       = H06_I_Oct(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Nov)%dat(1)       = H06_I_Nov(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_Dec)%dat(1)       = H06_I_Dec(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Jan)%dat(1)       = H06_D_Jan(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Feb)%dat(1)       = H06_D_Feb(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Mar)%dat(1)       = H06_D_Mar(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Apr)%dat(1)       = H06_D_Apr(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_May)%dat(1)       = H06_D_May(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Jun)%dat(1)       = H06_D_Jun(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Jul)%dat(1)       = H06_D_Jul(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Aug)%dat(1)       = H06_D_Aug(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Sep)%dat(1)       = H06_D_Sep(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Oct)%dat(1)       = H06_D_Oct(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Nov)%dat(1)       = H06_D_Nov(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_Dec)%dat(1)       = H06_D_Dec(ix)
         structSEG_main  (ix)%var(ixSEG%H06_purpose)%dat(1)     = H06_purpose(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_mem_F)%dat(1)     = H06_I_mem_F(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_mem_F)%dat(1)     = H06_D_mem_F(ix)
         structSEG_main  (ix)%var(ixSEG%H06_I_mem_L)%dat(1)     = H06_I_mem_L(ix)
         structSEG_main  (ix)%var(ixSEG%H06_D_mem_L)%dat(1)     = H06_D_mem_L(ix)
       end if
     end do main_rch

     ups_trib: do ix = 1, nTribOutlet
       ixx = global_ix_comm(ix)
       structNTOPO_main(nRch_mainstem+ix)%var(ixNTOPO%segId)%dat(1)         = structNTOPO(ixx)%var(ixNTOPO%segId)%dat(1)
       structNTOPO_main(nRch_mainstem+ix)%var(ixNTOPO%downSegId)%dat(1)     = structNTOPO(ixx)%var(ixNTOPO%downSegId)%dat(1)
       structSEG_main  (nRch_mainstem+ix)%var(ixSEG%length)%dat(1)          = structSEG(ixx)%var(ixSEG%length)%dat(1)
       structSEG_main  (nRch_mainstem+ix)%var(ixSEG%slope)%dat(1)           = structSEG(ixx)%var(ixSEG%slope)%dat(1)
       if (is_lake_sim) then
         structNTOPO_main(nRch_mainstem+ix)%var(ixNTOPO%islake)%dat(1)        = structNTOPO(ixx)%var(ixNTOPO%islake)%dat(1)
         structNTOPO_main(nRch_mainstem+ix)%var(ixNTOPO%LakeTargVol)%dat(1)   = structNTOPO(ixx)%var(ixNTOPO%LakeTargVol)%dat(1)
         structNTOPO_main(nRch_mainstem+ix)%var(ixNTOPO%LakeModelType)%dat(1) = structNTOPO(ixx)%var(ixNTOPO%LakeModelType)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%D03_MaxStorage)%dat(1)  = structSEG(ixx)%var(ixSEG%D03_MaxStorage)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%D03_Coefficient)%dat(1) = structSEG(ixx)%var(ixSEG%D03_Coefficient)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%D03_Power)%dat(1)       = structSEG(ixx)%var(ixSEG%D03_Power)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_Smax)%dat(1)        = structSEG(ixx)%var(ixSEG%H06_Smax)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_alpha)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_alpha)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_envfact)%dat(1)     = structSEG(ixx)%var(ixSEG%H06_envfact)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_S_ini)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_S_ini)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_c1)%dat(1)          = structSEG(ixx)%var(ixSEG%H06_c1)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_c2)%dat(1)          = structSEG(ixx)%var(ixSEG%H06_c1)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_exponent)%dat(1)    = structSEG(ixx)%var(ixSEG%H06_exponent)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_denominator)%dat(1) = structSEG(ixx)%var(ixSEG%H06_denominator)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_c_compare)%dat(1)   = structSEG(ixx)%var(ixSEG%H06_c_compare)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_frac_Sdead)%dat(1)  = structSEG(ixx)%var(ixSEG%H06_frac_Sdead)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_E_rel_ini)%dat(1)   = structSEG(ixx)%var(ixSEG%H06_E_rel_ini)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Jan)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Jan)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Feb)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Feb)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Mar)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Mar)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Apr)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Apr)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_May)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_May)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Jun)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Jun)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Jul)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Jul)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Aug)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Aug)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Sep)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Sep)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Oct)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Oct)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Nov)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Nov)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_Dec)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Dec)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Jan)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_I_Jan)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Feb)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Feb)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Mar)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Mar)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Apr)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Apr)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_May)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_May)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Jun)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Jun)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Jul)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Jul)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Aug)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Aug)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Sep)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Sep)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Oct)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Oct)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Nov)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Nov)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_Dec)%dat(1)       = structSEG(ixx)%var(ixSEG%H06_D_Dec)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_purpose)%dat(1)     = structSEG(ixx)%var(ixSEG%H06_purpose)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_mem_F)%dat(1)     = structSEG(ixx)%var(ixSEG%H06_I_mem_F)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_D_mem_F)%dat(1)     = structSEG(ixx)%var(ixSEG%H06_D_mem_F)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_mem_L)%dat(1)     = structSEG(ixx)%var(ixSEG%H06_I_mem_L)%dat(1)
         structSEG_main  (nRch_mainstem+ix)%var(ixSEG%H06_I_mem_L)%dat(1)     = structSEG(ixx)%var(ixSEG%H06_I_mem_L)%dat(1)
       end if
       global_ix_main(ix) = nRch_mainstem+ix   ! index in mainstem array that is link to tributary outlet
     end do ups_trib

     main_hru: do ix=1, nHRU_mainstem
       structHRU2SEG_main(ix)%var(ixHRU2SEG%HRUid)%dat(1)    = hruId(ix)
       structHRU2SEG_main(ix)%var(ixHRU2SEG%hruSegId)%dat(1) = hruSegId(ix)
       structHRU_main    (ix)%var(ixHRU%area)%dat(1)         = area(ix)
     end do main_hru

     ! compute additional ancillary infomration
     call augment_ntopo(nHRU_mainstem,             & ! input: number of HRUs
                        nRch_mainstem+nTribOutlet, & ! input: number of stream segments
                        structHRU_main,            & ! inout: ancillary data for HRUs
                        structSEG_main,            & ! inout: ancillary data for stream segments
                        structHRU2seg_main,        & ! inout: ancillary data for mapping hru2basin
                        structNTOPO_main,          & ! inout: ancillary data for network toopology
                        ierr, cmessage)              ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     do iSeg = 1, nRch_mainstem
       jSeg = ixRch_order(iSeg) ! global index, ordered by domain/node
       structSEG_main(iSeg) = structSEG(jSeg)
       nUps = size(structNTOPO(jSeg)%var(ixNTOPO%upSegIds)%dat)
       deallocate(structNTOPO_main(iSeg)%var(ixNTOPO%goodBasin)%dat)
       allocate(structNTOPO_main(iSeg)%var(ixNTOPO%goodBasin)%dat(nUps))
       structNTOPO_main(iSeg)%var(ixNTOPO%goodBasin)%dat(:) = structNTOPO(jSeg)%var(ixNTOPO%goodBasin)%dat(:)
       allocate(jUps(nUps))
       do iUps = 1,nUps
         do ix = 1,nUps
           if (structNTOPO_main(iSeg)%var(ixNTOPO%upSegIds)%dat(ix) == structNTOPO(jSeg)%var(ixNTOPO%upSegIds)%dat(iUps)) then
             jUps(iUps) = ix
           end if
         end do
       end do
       structNTOPO_main(iSeg)%var(ixNTOPO%upSegIds    )%dat(:) = structNTOPO_main(iSeg)%var(ixNTOPO%upSegIds)%dat(jUps)
       structNTOPO_main(iSeg)%var(ixNTOPO%upSegIndices)%dat(:) = structNTOPO_main(iSeg)%var(ixNTOPO%upSegIndices)%dat(jUps)
       deallocate(jUps)
     end do

     do iSeg = 1, nTribOutlet
       jSeg = global_ix_comm(iSeg)
       structSEG_main(global_ix_main(iSeg)) = structSEG(jSeg)
       nUps = size(structNTOPO(jSeg)%var(ixNTOPO%upSegIds)%dat)
       deallocate(structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%goodBasin   )%dat, &
                  structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%upSegIds    )%dat, &
                  structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%upSegIndices)%dat)
       allocate(structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%goodBasin   )%dat(nUps),&
                structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%upSegIds    )%dat(nUps),&
                structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%upSegIndices)%dat(nUps))
       structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%goodBasin   )%dat(:) = structNTOPO(jSeg)%var(ixNTOPO%goodBasin)%dat(:)
       structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%upSegIds    )%dat(:) = structNTOPO(jSeg)%var(ixNTOPO%upSegIds)%dat(:)
       structNTOPO_main(global_ix_main(iSeg))%var(ixNTOPO%upSegIndices)%dat(:) = integerMissing
     enddo

     ! find index of desired reach
     if (desireId/=integerMissing) ixPrint(1) = findIndex(segId, desireId, integerMissing)

     ! copy data to routing structres RPARAM_trib and NETOPO_trib
     call put_data_struct(nRch_mainstem+nTribOutlet, structSEG_main, structNTOPO_main, & ! input
                          RPARAM_main, NETOPO_main,                                    & ! output:
                          ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! OMP domain decomposition
     call omp_domain_decomposition(stream_order, nRch_mainstem+nTribOutlet, structNTOPO_main, river_basin_main, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    if (debug_main_omp) then
      write(iulog,'(a)') 'segid, branch, order'
      do ix = 1, size(river_basin_main)
        do ixx = 1, size(river_basin_main(ix)%branch)
          do iSeg = 1, river_basin_main(ix)%branch(ixx)%nRch
            associate (idx_tmp => river_basin_main(ix)%branch(ixx)%segIndex(iSeg))
            write(iulog,"(I15,A,I9,A,I9)") structNTOPO_main(idx_tmp)%var(ixNTOPO%segId)%dat(1),',',ixx,',',ix
            end associate
          end do
        end do
      enddo
    endif

   end if ! (masterproc)
  end if ! (nRch_mainstem > 0)

  call shr_mpi_barrier(comm, message)

 end subroutine comm_ntopo_data


 ! *********************************************************************
 ! public subroutine: restart flux/state communication
 ! *********************************************************************
 subroutine mpi_restart(pid,           & ! input: proc id
                        nNodes,        & ! input: number of procs
                        comm,          & ! input: communicator
                        iens,          & ! input: ensemble index
                        ierr,message)    ! output: error control
  ! shared data
  USE globalData, ONLY: nRch             ! number of reaches in the whoel river network
  USE globalData, ONLY: rch_per_proc     ! number of reaches assigned to each proc (i.e., node)
  USE globalData, ONLY: nRch_mainstem    ! number of mainstem reaches
  USE globalData, ONLY: ixRch_order      ! global reach index in the order of proc assignment (size = total number of reaches in the entire network)
  USE globalData, ONLY: global_ix_comm   ! global reach index at tributary reach outlets to mainstem (size = sum of tributary outlets within entire network)
  USE globalData, ONLY: RCHFLX_trib      ! tributary reach flux structure
  USE globalData, ONLY: RCHFLX_main      ! mainstem reach flux structure
  USE globalData, ONLY: RCHSTA_trib      ! tributary reach state structure
  USE globalData, ONLY: RCHSTA_main      ! mainstem reach state structure
  USE globalData, ONLY: RCHFLX           ! entire reach flux structure
  USE globalData, ONLY: RCHSTA           !
  USE globalData, ONLY: TSEC             ! beginning/ending of simulation time step [sec]

  implicit none

  ! input variables
  integer(i4b),             intent(in)  :: pid                      ! process id (MPI)
  integer(i4b),             intent(in)  :: nNodes                   ! number of processes (MPI)
  integer(i4b),             intent(in)  :: comm                     ! communicator
  integer(i4b),             intent(in)  :: iens                     ! ensemble index
  ! Output variables
  integer(i4b),             intent(out) :: ierr
  character(len=strLen),    intent(out) :: message                  ! error message
  ! local variables
  character(len=strLen)                 :: cmessage                 ! error message from subroutine
  integer(i4b)                          :: iSeg,jSeg,iTbound        ! loop indices
  integer(i4b)                          :: nTbound = 2
  real(dp),     allocatable             :: flux_global(:)           ! basin runoff (m/s) for entire reaches
  real(dp),     allocatable             :: vol_global(:,:)          ! reach/lake volume (m3) for entire network
  real(dp),     allocatable             :: vol_global_tmp(:)        ! temporary reach/lake volume (m3) for entire network
  real(dp),     allocatable             :: flux_local(:)            ! basin runoff (m/s) for tributaries
  real(dp),     allocatable             :: vol_local(:)             ! reach/lake volume (m3) for tributaries

  ierr=0; message='mpi_restart/'

  call MPI_BCAST(TSEC, nTbound, MPI_DOUBLE_PRECISION, root, comm, ierr)

  allocate(flux_global(nRch), vol_global(nRch,nTbound), vol_global_tmp(nRch))
  allocate(vol_local(rch_per_proc(pid)), flux_local(rch_per_proc(pid)))

  if (masterproc) then

    do iSeg = 1, nRch
      flux_global(iSeg) = RCHFLX(iens,iSeg)%BASIN_QR(1)
      vol_global(iSeg,1:2) = RCHFLX(iens,iSeg)%REACH_VOL(0:1)
    enddo

    ! Distribute global flux/state (RCHFLX & RCHSTA) to mainstem (RCHFLX_main & RCHSTA_main)
    do iSeg = 1, nRch_mainstem
      jSeg = ixRch_order(iSeg)
      RCHFLX_main(iens, iSeg) = RCHFLX(iens, jSeg)
      RCHSTA_main(iens, iSeg) = RCHSTA(iens, jSeg)
    end do
    do iSeg = 1,size(global_ix_comm)
      jSeg = global_ix_comm(iSeg)
      RCHFLX_main(iens, iSeg+nRch_mainstem) = RCHFLX(iens, jSeg)
      RCHSTA_main(iens, iSeg+nRch_mainstem) = RCHSTA(iens, jSeg)
    end do

  else

    flux_global(:)  = realMissing
    vol_global(:,:) = realMissing

  endif

  ! Distribute global flux/state (RCHFLX & RCHSTA) to tributary (RCHFLX_trib & RCHSTA_trib)
  ! flux communication (only basin delayed runoff flux)
  call mpi_comm_single_flux(pid, nNodes, comm,                        &
                            flux_global,                              &
                            flux_local,                               &
                            rch_per_proc(root:nNodes-1),              &
                            ixRch_order(rch_per_proc(root-1)+1:nRch), &
                            arth(1,1,rch_per_proc(pid)),              &
                            scatter,                                  &
                            ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  do iSeg = 1, rch_per_proc(pid)
    RCHFLX_trib(iens,iSeg)%BASIN_QR(1) = flux_local(iSeg)
  enddo

  ! basin IRF state communication
  if (doesBasinRoute == 1) then
    call mpi_comm_irf_bas_state(pid, nNodes, comm,                        &
                                iens,                                     &
                                rch_per_proc(root:nNodes-1),              &
                                RCHFLX,                                   &
                                RCHFLX_trib,                              &
                                ixRch_order(rch_per_proc(root-1)+1:nRch), &
                                arth(1,1,rch_per_proc(pid)),              &
                                scatter,                                  & ! communication type
                                ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! KWT state communication
  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
    call mpi_comm_kwt_state(pid, nNodes, comm,                        &
                            iens,                                     &
                            rch_per_proc(root:nNodes-1),              &
                            RCHSTA,                                   &
                            RCHSTA_trib,                              &
                            ixRch_order(rch_per_proc(root-1)+1:nRch), &
                            arth(1,1,rch_per_proc(pid)),              &
                            scatter,                                  & ! communication type
                            ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! KWE state communication
  if (routOpt==kinematicWaveEuler) then
    call mpi_comm_kwe_state(pid, nNodes, comm,                        &
                            iens,                                     &
                            rch_per_proc(root:nNodes-1),              &
                            RCHSTA,                                   &
                            RCHSTA_trib,                              &
                            ixRch_order(rch_per_proc(root-1)+1:nRch), &
                            arth(1,1,rch_per_proc(pid)),              &
                            scatter,                                  & ! communication type
                            ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! IRF state communication
  if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
    call mpi_comm_irf_state(pid, nNodes, comm,                        &
                            iens,                                     &
                            rch_per_proc(root:nNodes-1),              &
                            RCHFLX,                                   &
                            RCHFLX_trib,                              &
                            ixRch_order(rch_per_proc(root-1)+1:nRch), &
                            arth(1,1,rch_per_proc(pid)),              &
                            scatter,                                  & ! communication type
                            ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! volume communication
    do iTbound =1,nTbound
      vol_global_tmp(:) = vol_global(:,iTbound)
      call mpi_comm_single_flux(pid, nNodes, comm,                        &
                                vol_global_tmp,                           &
                                vol_local,                                &
                                rch_per_proc(root:nNodes-1),              &
                                ixRch_order(rch_per_proc(root-1)+1:nRch), &
                                arth(1,1,rch_per_proc(pid)),              &
                                scatter,                                  &
                                ierr, message)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      do iSeg = 1, rch_per_proc(pid)
        RCHFLX_trib(iens,iSeg)%REACH_VOL(iTbound-1) = vol_local(iSeg)
      enddo

    end do

  endif

  ! no need for the entire domain flux/state data strucure
  if (masterproc) then
    deallocate(RCHFLX, RCHSTA)
  end if

 end subroutine mpi_restart

 ! *********************************************************************
 ! public subroutine: routing routine
 ! *********************************************************************
 subroutine mpi_route(pid,           & ! input: proc id
                      nNodes,        & ! input: number of procs
                      comm,          & ! input: communicator
                      iens,          & ! input: ensemble index
                      ierr,message,  & ! output: error control
                      scatter_ro)      ! optional input: logical. .true. => scatter global hru runoff to mainstem and tributary
  ! shared data
  USE globalData, ONLY: ixPrint             ! desired reach index
  USE globalData, ONLY: NETOPO_trib         ! tributary river netowrk topology structure
  USE globalData, ONLY: NETOPO_main         ! mainstem river netowrk topology structure
  USE globalData, ONLY: RPARAM_trib         ! tributary reach parameter structure
  USE globalData, ONLY: RPARAM_main         ! mainstem reach parameter structure
  USE globalData, ONLY: RCHFLX_trib         ! tributary reach flux structure
  USE globalData, ONLY: RCHFLX_main         ! Reach flux data structures (master proc, mainstem)
  USE globalData, ONLY: RCHSTA_trib         ! tributary reach state data structure
  USE globalData, ONLY: RCHSTA_main         ! Reach state data structures (master proc, mainstem)
  USE globalData, ONLY: nHRU_mainstem       ! number of mainstem HRUs
  USE globalData, ONLY: basinRunoff_main    ! mainstem only HRU runoff
  USE globalData, ONLY: basinRunoff_trib    ! tributary only HRU runoff
  USE globalData, ONLY: basinEvapo_main     ! mainstem only HRU Evaporation
  USE globalData, ONLY: basinEvapo_trib     ! tributary only HRU Evaporation
  USE globalData, ONLY: basinPrecip_main    ! mainstem only HRU Precipitation
  USE globalData, ONLY: basinPrecip_trib    ! tributary only HRU Precipitation
  USE globalData, ONLY: river_basin_trib    ! tributary OMP domain data structure
  USE globalData, ONLY: river_basin_main    ! mainstem OMP domain data structure
  USE globalData, ONLY: nRch_mainstem       ! number of mainstem reaches
  USE globalData, ONLY: flux_wm_main        ! nRch flux holder for mainstem
  USE globalData, ONLY: flux_wm_trib        ! nRch flux holder for tributary
  USE globalData, ONLY: vol_wm_main         ! nRch target vol holder for mainstem
  USE globalData, ONLY: vol_wm_trib         ! nRch target vol holder for tributary
  USE globalData, ONLY: rch_per_proc        ! number of reaches assigned to each proc
  USE globalData, ONLY: tribOutlet_per_proc ! number of tributary outlets per proc (array size = nNodes)
  USE globalData, ONLY: global_ix_main      ! reach index at tributary reach outlets to mainstem (size = sum of tributary outlets in all the procs)
  USE globalData, ONLY: local_ix_comm       ! local reach index at tributary reach outlets to mainstem for each proc (size = sum of tributary outlets in proc)
  USE public_var, ONLY: is_lake_sim         ! logical whether or not lake should be simulated
  USE public_var, ONLY: is_flux_wm          ! logical whether or not fluxes should be passed
  USE public_var, ONLY: is_vol_wm           ! logical whether or not target volume should be passed
  ! routing driver
  USE main_route_module, ONLY: main_route   ! routing driver

  ! Details:
  ! Reaches/HRUs assigned to master proc include BOTH small tributaries and mainstem. Reaches/HRUs assigned to other procs includ ONLY tributaries
  ! 1. Route "small" tributaries" assigned to master proc and "bigger tributaries" on other procs the same time
  ! 2 (if mainstem exist). pass flow/state variables at "tributary-mainstem junction", which is outlet of tributaries into mainstem to master proc
  ! 3 (if mainstem exist). Route mainstem at master proc
  ! 4 (if mainstem exist & for a certain routing method). pass updated flow/state variables at "tributary-mainstem junction to other procs

  implicit none

  ! input variables
  integer(i4b),             intent(in)  :: pid                      ! process id (MPI)
  integer(i4b),             intent(in)  :: nNodes                   ! number of processes (MPI)
  integer(i4b),             intent(in)  :: comm                     ! communicator
  integer(i4b),             intent(in)  :: iens                     ! ensemble index
  ! Output variables
  integer(i4b),             intent(out) :: ierr
  character(len=strLen),    intent(out) :: message                  ! error message
  ! optional input variables
  logical(lgt), optional,   intent(in)  :: scatter_ro               ! logical to indicate if scattering global runoff is required
  ! local variables
  integer(i4b), allocatable             :: ixRchProcessed(:)        ! reach indice list to be processed
  integer(i4b)                          :: nRch_trib                ! number of reaches on one tributary (in a proc)
  logical(lgt)                          :: doesScatterRunoff        ! temporal logical to indicate if scattering global runoff is required
  character(len=strLen)                 :: cmessage                 ! error message from subroutine

  ierr=0; message='mpi_route/'

  if (present(scatter_ro)) then
    doesScatterRunoff = scatter_ro
  else
    doesScatterRunoff = .true.
  end if

  ! --------------------------------
  ! distribute (scatter) global runoff, evaporation
  ! and precipitation array to local runoff arrays
  !  - basinRunoff_main (only at master proc)
  !  - basinRunoff_trib (all procs)
  !  - basinEvapo_main  (only at master proc)
  !  - basinEvapo_trib  (all procs)
  !  - basinPrecip_main (only at master proc)
  !  - basinPrecip_trib (all procs)
  ! --------------------------------
  if (doesScatterRunoff) then

    call t_startf ('route/scatter-runoff')
    call scatter_runoff(nNodes, comm, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call t_stopf ('route/scatter-runoff')

  else

    if (nNodes==1) then

      if (.not.allocated(basinRunoff_main)) then
        ierr=10; message=trim(message)//'master proc: mainstem runoff array is not allocated'; return
      end if
      if (is_lake_sim) then
        if (.not.allocated(basinEvapo_main)) then
          ierr=10; message=trim(message)//'master proc: mainstem evaporation array is not allocated'; return
        end if
        if (.not.allocated(basinPrecip_main)) then
          ierr=10; message=trim(message)//'master proc: mainstem precipitation array is not allocated'; return
        end if
      end if

    else

      if (masterproc) then

        write(iulog,*) 'NOTE: HRU Runoff is already decomposed into mainstems and tributaries. No need for scatter_runoff'

        if (nHRU_mainstem > 0 .and. .not.allocated(basinRunoff_main)) then
          ierr=10; message=trim(message)//'mainstem runoff array is not allocated/populated'; return
        end if
        if (.not.allocated(basinRunoff_trib)) then
          ierr=10; message=trim(message)//'master proc: tributary runoff array is not allocated'; return
        end if

        if (is_lake_sim) then
          if (nHRU_mainstem > 0 .and. .not.allocated(basinEvapo_main)) then
            ierr=10; message=trim(message)//'mainstem evaporation array is not allocated/populated'; return
          end if
          if (.not.allocated(basinEvapo_trib)) then
            ierr=10; message=trim(message)//'master proc: tributary evaporation array is not allocated'; return
          end if
          if (.not.allocated(basinPrecip_trib)) then
            ierr=10; message=trim(message)//'master proc: tributary precipitation array is not allocated'; return
          end if
          if (nHRU_mainstem > 0 .and. .not.allocated(basinPrecip_main)) then
            ierr=10; message=trim(message)//'mainstem precipitation array is not allocated/populated'; return
          end if
        end if

      else

        if(.not.allocated(basinRunoff_trib)) then
          ierr=10; message=trim(message)//'tributary runoff array is not allocated/populated'; return
        end if
        if (is_lake_sim) then
          if (.not.allocated(basinEvapo_trib)) then
            ierr=10; message=trim(message)//'tributary evaporation array is not allocated/populated'; return
          end if
          if (.not.allocated(basinPrecip_trib)) then
            ierr=10; message=trim(message)//'tributary precipitation array is not allocated/populated'; return
          end if
        end if
      end if
    end if
  end if


  ! --------------------------------
  ! distribute (scatter) water managemnt flux and
  ! volume in the
  !  - flux_wm_main (only at master proc)
  !  - flux_wm_trib (all procs)
  !  - vol_wm_main  (only at master proc)
  !  - vol_wm_trib  (all procs)
  ! --------------------------------

  if (is_flux_wm.or.(is_vol_wm.and.is_lake_sim)) then

    if (doesScatterRunoff) then

      call t_startf ('route/scatter-wm')
      call scatter_wm(nNodes, comm, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      call t_stopf ('route/scatter-wm')

    else

      if (nNodes==1) then

        if (is_flux_wm) then
          if (.not.allocated(flux_wm_main)) then
            ierr=10; message=trim(message)//'master proc: mainstem water management flux array is not allocated'; return
          end if
        end if
        if (is_vol_wm.and.is_lake_sim) then
          if (.not.allocated(vol_wm_main)) then
            ierr=10; message=trim(message)//'master proc: mainstem target volume array (for lakes) is not allocated'; return
          end if
        end if

      else

        if (masterproc) then

          write(iulog,*) 'NOTE: reach flux and target volume are already decomposed into mainstems and tributaries. No need for scatter_wm'

          if (is_flux_wm) then
            if (nRch_mainstem > 0 .and. .not.allocated(flux_wm_main)) then
              ierr=10; message=trim(message)//'mainstem water management flux array is not allocated/populated'; return
            end if
            if (.not.allocated(flux_wm_trib)) then
              ierr=10; message=trim(message)//'master proc: tributary water management flux array is not allocated'; return
            end if
          end if

          if (is_vol_wm.and.is_lake_sim) then
            if (nRch_mainstem > 0 .and. .not.allocated(vol_wm_main)) then
              ierr=10; message=trim(message)//'mainstem target volume array (for lakes) is not allocated/populated'; return
            end if
            if (.not.allocated(vol_wm_trib)) then
              ierr=10; message=trim(message)//'master proc: tributary target volume array (for lakes) is not allocated'; return
            end if
          end if


        else

          if (is_flux_wm) then
            if(.not.allocated(flux_wm_trib)) then
              ierr=10; message=trim(message)//'tributary water management fluxes array is not allocated/populated'; return
            end if
          end if

          if (is_vol_wm.and.is_lake_sim) then
            if(.not.allocated(vol_wm_trib)) then
              ierr=10; message=trim(message)//'tributary target volume array (for lakes) is not allocated/populated'; return
            end if
          end if
        end if
      end if
    end if
  end if


  ! --------------------------------
  ! Perform tributary routing (for all procs)
  ! --------------------------------
  if (nNodes>1) then

  call t_startf ('route/tributary-route')

  !Idenfity number of tributary reaches for each procs
  nRch_trib = rch_per_proc(pid)
  allocate(ixRchProcessed(nRch_trib), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [ixRchProcessed]'; return; endif

  ! Define processing reach indices in terms of tributary data sets
  ixRchProcessed = arth(1,1,nRch_trib)

  ! Perform routing
  call main_route(iens,              &  ! input: ensemble index
                  basinRunoff_trib,  &  ! input: basin (i.e.,HRU) runoff (m/s)
                  basinEvapo_trib,   &  ! input: basin (i.e. HRU) Evapo  (m/s)
                  basinPrecip_trib,  &  ! input: basin (i.e. HRU) Precip (m/s)
                  flux_wm_trib,      &  ! reach (i.e.,reach) flux (m3/s)
                  vol_wm_trib,       &  ! reach (i.e.,reach) target volume for lakes (m3)
                  ixRchProcessed,    &  ! input: indices of reach to be routed
                  river_basin_trib,  &  ! input: OMP basin decomposition
                  NETOPO_trib,       &  ! input: reach topology data structure
                  RPARAM_trib,       &  ! input: reach parameter data structure
                  ixPrint(2),        &  ! input: reach index to be checked by on-screen pringing
                  RCHFLX_trib,       &  ! inout: reach flux data structure
                  RCHSTA_trib,       &  ! inout: reach state data structure
                  ierr, message)        ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call t_stopf ('route/tributary-route')

  ! make sure that routing at all the procs finished
  call shr_mpi_barrier(comm, message)

  if (nRch_mainstem==0) return   ! if there are no mainstem reaches, finished

  ! --------------------------------
  ! Collect all the tributary flows
  ! --------------------------------
  call t_startf('route/gather-state-flux')

  ! flux communication
  call mpi_comm_flux(pid, nNodes, comm,   & ! input: mpi rank, number of tasks, and communicator
                     iens,                & ! input: ensemble index (not used now)
                     tribOutlet_per_proc, & ! input: number of reaches communicate per node (dimension size == number of proc)
                     RCHFLX_main,         & ! input:
                     RCHFLX_trib,         & ! input:
                     global_ix_main,      & ! input:
                     local_ix_comm,       & ! input: local reach indices per proc (dimension size depends on procs )
                     gather,              & ! input: communication type
                     ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! KWT state communication
  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then

    call mpi_comm_kwt_state(pid, nNodes, comm,   & ! input: mpi rank, number of tasks, and communicator
                            iens,                & ! input:
                            tribOutlet_per_proc, & ! input: number of reaches communicate per node (dimension size == number of proc)
                            RCHSTA_main,         & ! input:
                            RCHSTA_trib,         & ! input:
                            global_ix_main,      & ! input:
                            local_ix_comm,       & ! input: local reach indices per proc (dimension size depends on procs )
                            gather,              & ! communication type
                            ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif

  call t_stopf('route/gather-state-flux')

  ! make sure that routing at all the procs finished
  call shr_mpi_barrier(comm, message)

  end if ! (nNodes>1)

  ! --------------------------------
  ! perform mainstem routing
  ! --------------------------------
  if (masterproc) then

    call t_startf ('route/mainstem_route')

    if (allocated(ixRchProcessed)) then
      deallocate(ixRchProcessed, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem deallocating array for [ixRchProcessed]'; return; endif
    end if
    allocate(ixRchProcessed(nRch_mainstem), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ixRchProcessed]'; return; endif

    ! Define processing reach indices
    ixRchProcessed = arth(1,1,nRch_mainstem)

    call main_route(iens,                    &  ! input: ensemble index
                    basinRunoff_main,        &  ! input: basin (i.e.,HRU) runoff (m/s)
                    basinEvapo_main,         &  ! input: basin (i.e. HRU) Evapo  (m/s)
                    basinPrecip_main,        &  ! input: basin (i.e. HRU) Precip (m/s)
                    flux_wm_main,            &  ! reach (i.e.,reach) flux (m3/s)
                    vol_wm_main,             &  ! reach (i.e.,reach) target volume for lakes (m3)
                    ixRchProcessed,          &  ! input: indices of reach to be routed
                    river_basin_main,        &  ! input: OMP basin decomposition
                    NETOPO_main,             &  ! input: reach topology data structure
                    RPARAM_main,             &  ! input: reach parameter data structure
                    ixPrint(1),              &  ! input: reach index to be checked by on-screen pringing
                    RCHFLX_main,             &  ! inout: reach flux data structure
                    RCHSTA_main,             &  ! inout: reach state data structure
                    ierr, message)              ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call t_stopf ('route/mainstem_route')

  endif ! end of root proc

  ! make sure that routing at all the procs finished
  call shr_mpi_barrier(comm, message)

  if (nNodes>1) then
  ! --------------------------------
  ! Distribute updated tributary states (only tributary reaches flowing into mainstem) to processors to update states upstream reaches
  ! --------------------------------
  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then

    call t_startf ('route/scatter-kwt-state')

    call mpi_comm_kwt_state(pid, nNodes, comm,   & ! input: mpi rank, number of tasks, and communicator
                            iens,                & ! input:
                            tribOutlet_per_proc, & ! input: number of reaches communicate per node (dimension size == number of proc)
                            RCHSTA_main,         & ! input:
                            RCHSTA_trib,         & ! input:
                            global_ix_main,      &
                            local_ix_comm,       & ! input: local reach indices per proc (dimension size depends on procs )
                            scatter,             & ! input: 1 = scatter
                            ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call t_stopf ('route/scatter-kwt-state')

  endif ! end of kwt option

  ! make sure that routing at all the procs finished
  call shr_mpi_barrier(comm, message)

  endif !(nNodes>1)

 end subroutine mpi_route

 ! *********************************************************************
 ! private subroutine: scatter global domain runoff data to local domain
 ! *********************************************************************
 SUBROUTINE scatter_runoff(nNodes, comm, & ! mpi variables: number nodes, communicator
                           ierr, message)  ! error controls

  USE globalData, ONLY: nHRU              ! number of all HRUs
  USE globalData, ONLY: nHRU_mainstem     ! number of mainstem HRUs
  USE globalData, ONLY: runoff_data       ! runoff data structure
  USE globalData, ONLY: hru_per_proc      ! number of hrus assigned to each proc (i.e., node)
  USE globalData, ONLY: basinRunoff_main  ! HRU runoff holder for mainstem
  USE globalData, ONLY: basinRunoff_trib  ! HRU runoff holder for tributary
  USE globalData, ONLY: basinEvapo_main   ! HRU evaporation holder for mainstem
  USE globalData, ONLY: basinEvapo_trib   ! HRU evaporation holder for tributary
  USE globalData, ONLY: basinPrecip_main  ! HRU precipitation holder for mainstem
  USE globalData, ONLY: basinPrecip_trib  ! HRU precipitation holder for tributary
  USE public_var, ONLY: is_lake_sim       ! logical whether or not lake should be simulated


  ! input variables
  integer(i4b),           intent(in)      :: nNodes                    ! number of processes (MPI)
  integer(i4b),           intent(in)      :: comm                      ! communicator
  ! output variables
  integer(i4b),           intent(out)     :: ierr                      ! error code
  character(len=strLen),  intent(out)     :: message                   ! error message
  ! local variables
  integer(i4b)                            :: iHru,jHru                 ! loop indices
  real(dp)                                :: basinRunoff_local(nHRU)   ! temporal basin runoff (m/s) for whole domain
  real(dp)                                :: basinEvapo_local(nHRU)    ! temporal basin evaporation (m/s) for whole domain
  real(dp)                                :: basinPrecip_local(nHRU)   ! temporal basin precipitation (m/s) for whole domain
  character(len=strLen)                   :: cmessage                  ! error message from a subroutine

  ierr=0; message='scatter_runoff/'

  if (nNodes==1) then

    ! if only single proc is used, all runoff is stored in mainstem runoff array
    if (.not. allocated(basinRunoff_main)) then
      allocate(basinRunoff_main(nHRU), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating array for [basinRunoff_main]'; return; endif
    end if
    basinRunoff_main(:) = runoff_data%basinRunoff(:)

    if (is_lake_sim) then

      ! if only single proc is used, all evaporation is stored in mainstem evaporation array
      if (.not. allocated(basinEvapo_main)) then
        allocate(basinEvapo_main(nHRU), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [basinEvapo_main]'; return; endif
      end if
      basinEvapo_main(:) = runoff_data%basinEvapo(:)

      ! if only single proc is used, all precipitation is stored in mainstem precipitation array
      if (.not. allocated(basinPrecip_main)) then
        allocate(basinPrecip_main(nHRU), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [basinPrecip_main]'; return; endif
      end if
      basinPrecip_main(:) = runoff_data%basinPrecip(:)

    end if

  else

    ! sort the basin runoff, precipitation and evaporation in terms of nodes/domains
    if (masterproc) then ! this is a root process

      if (.not. allocated(basinRunoff_main)) then
        allocate(basinRunoff_main(nHRU_mainstem), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [basinRunoff_main]'; return; endif
      endif

      if (is_lake_sim) then

        if (.not. allocated(basinEvapo_main)) then
          allocate(basinEvapo_main(nHRU_mainstem), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating array for [basinEvapo_main]'; return; endif

        endif

        if (.not. allocated(basinPrecip_main)) then
          allocate(basinPrecip_main(nHRU_mainstem), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating array for [basinPrecip_main]'; return; endif
        endif

      end if


      ! runoff at hru in mainstem and tributaries
      basinRunoff_local(1:nHRU) = runoff_data%basinRunoff(1:nHRU)
      basinRunoff_main(1:nHRU_mainstem) = basinRunoff_local(1:nHRU_mainstem)


      ! evaporation and precipitation at main channel and tributaries
      if (is_lake_sim) then
        basinEvapo_local (1:nHRU) = runoff_data%basinEvapo(1:nHRU)
        basinPrecip_local(1:nHRU) = runoff_data%basinPrecip(1:nHRU)
        basinEvapo_main (1:nHRU_mainstem) = basinEvapo_local (1:nHRU_mainstem)
        basinPrecip_main(1:nHRU_mainstem) = basinPrecip_local(1:nHRU_mainstem)
      end if

    end if

    call shr_mpi_barrier(comm, message)

    ! Distribute the basin runoff to each process
    call shr_mpi_scatterV(basinRunoff_local(nHRU_mainstem+1:nHRU), &
                          hru_per_proc(0:nNodes-1),                &
                          basinRunoff_trib,                        &
                          ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    if (is_lake_sim) then
      call shr_mpi_scatterV(basinEvapo_local(nHRU_mainstem+1:nHRU),  &
                            hru_per_proc(0:nNodes-1),                &
                            basinEvapo_trib,                         &
                            ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call shr_mpi_scatterV(basinPrecip_local(nHRU_mainstem+1:nHRU), &
                            hru_per_proc(0:nNodes-1),                &
                            basinPrecip_trib,                        &
                            ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end if

  end if

 END SUBROUTINE scatter_runoff


 ! *********************************************************************
 ! private subroutine: scatter global domain water management data to local domain
 ! *********************************************************************
 SUBROUTINE scatter_wm(nNodes, comm, &    ! mpi variables: number nodes, communicator
                       ierr, message)     ! error controls

  USE globalData, ONLY: nRch              ! number of all reach
  USE globalData, ONLY: nRch_mainstem     ! number of mainstem reach
  USE globalData, ONLY: wm_data           ! water management data structure
  USE globalData, ONLY: rch_per_proc      ! number of reach assigned to each proc (i.e., node)
  USE globalData, ONLY: flux_wm_main      ! nRch flux holder for mainstem
  USE globalData, ONLY: flux_wm_trib      ! nRch flux holder for tributary
  USE globalData, ONLY: vol_wm_main       ! nRch target vol holder for mainstem
  USE globalData, ONLY: vol_wm_trib       ! nRch target vol holder for tributary
  USE public_var, ONLY: is_lake_sim       ! logical whether or not fluxes should be passed
  USE public_var, ONLY: is_flux_wm        ! logical whether or not fluxes should be passed
  USE public_var, ONLY: is_vol_wm         ! logical whether or not target volume should be passed

  ! input variables
  integer(i4b),           intent(in)  :: nNodes                          ! number of processes (MPI)
  integer(i4b),           intent(in)  :: comm                            ! communicator
  ! output variables
  integer(i4b),           intent(out) :: ierr                            ! error code
  character(len=strLen),  intent(out) :: message                         ! error message
  ! local variables
  integer(i4b)                        :: iHru,jHru                       ! loop indices
  real(dp)                            :: Rch_flux_local(nRch)            ! temporal reach flux (m3/s) for whole domain
  real(dp)                            :: Rch_vol_local(nRch)             ! temporal reach (lake) volume (m3) for whole domain
  character(len=strLen)               :: cmessage                        ! error message from a subroutine

  ierr=0; message='scatter_runoff/'

  if (nNodes==1) then

    if (is_flux_wm) then

      ! if only single proc is used, all fluxes are stored in mainstem runoff array
      if (.not. allocated(flux_wm_main)) then
        allocate(flux_wm_main(nRch), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [flux_wm_main]'; return; endif
      end if
      flux_wm_main(:) = wm_data%flux_wm(:)

    endif

    if (is_vol_wm.and.is_lake_sim) then

      ! if only single proc is used, all target volumes are stored in mainstem runoff array
      if (.not. allocated(vol_wm_main)) then
        allocate(vol_wm_main(nRch), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [vol_wm_main]'; return; endif
      end if
      vol_wm_main(:) = wm_data%vol_wm(:)

    endif

  else

    ! sort the reach flux, target vol in terms of nodes/domains
    if (masterproc) then ! this is a root process

      if (is_flux_wm) then
        if (.not. allocated(flux_wm_main)) then
          allocate(flux_wm_main(nRch_mainstem), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating array for [flux_wm_main]'; return; endif
        endif
        ! flux or target vol at reach in mainstem and tributaries
        Rch_flux_local(1:nRch) = wm_data%flux_wm(1:nRch)
        flux_wm_main(1:nRch_mainstem) = Rch_flux_local(1:nRch_mainstem)
      endif

      if (is_vol_wm.and.is_lake_sim) then
        if (.not. allocated(vol_wm_main)) then
          allocate(vol_wm_main(nRch_mainstem), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating array for [vol_wm_main]'; return; endif
        endif
        ! target vol at reach in mainstem and tributaries
        Rch_vol_local(1:nRch) = wm_data%vol_wm(1:nRch)
        vol_wm_main(1:nRch_mainstem) = Rch_vol_local(1:nRch_mainstem)
      end if

    end if

    call shr_mpi_barrier(comm, message)

    ! Distribute the read flux to each process
    if (is_flux_wm) then
      call shr_mpi_scatterV(Rch_flux_local(nRch_mainstem+1:nRch),  &
                            rch_per_proc(0:nNodes-1),              &
                            flux_wm_trib,                          &
                            ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

    if (is_vol_wm.and.is_lake_sim) then
      call shr_mpi_scatterV(Rch_vol_local(nRch_mainstem+1:nRch),   &
                            rch_per_proc(0:nNodes-1),              &
                            vol_wm_trib,                           &
                            ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

  end if

 END SUBROUTINE scatter_wm

 ! *********************************************************************
 ! subroutine: single flux communication
 ! *********************************************************************
 subroutine mpi_comm_single_flux(pid,          &
                                 nNodes,       &
                                 comm,         & ! input: communicator
                                 flux_global,  &
                                 flux_local,   &
                                 nReach,       &
                                 rchIdxGlobal, &
                                 rchIdxLocal,  &
                                 commType,     &
                                 ierr, message)

  ! input variables
  integer(i4b),             intent(in)    :: pid                   ! process id (MPI)
  integer(i4b),             intent(in)    :: nNodes                ! number of processes (MPI)
  integer(i4b),             intent(in)    :: comm                  ! communicator
  real(dp),     allocatable,intent(inout) :: flux_global(:)        ! global flux vectors
  real(dp),     allocatable,intent(inout) :: flux_local(:)         ! local flux vectors
  integer(i4b),             intent(in)    :: nReach(0:nNodes-1)    ! number of reaches communicate per node (dimension size == number of proc)
  integer(i4b),             intent(in)    :: rchIdxGlobal(:)       ! reach indices (w.r.t. global) to be transfer (dimension size == sum of nRearch)
  integer(i4b),             intent(in)    :: rchIdxLocal(:)        ! reach indices (w.r.t. local) (dimension size depends on procs )
  integer(i4b),             intent(in)    :: commType              ! communication type 1->scatter, 2->gather otherwise error
  ! output variables
  integer(i4b),             intent(out)   :: ierr                  ! error code
  character(len=strLen),    intent(out)   :: message               ! error message
  ! local variables
  character(len=strLen)                   :: cmessage              ! error message from a subroutine
  real(dp),     allocatable               :: flux_global_tmp(:)    ! temporal global flux array
  real(dp),     allocatable               :: flux_local_tmp(:)     ! temporal local flux array
  integer(i4b)                            :: nSeg                  ! number of reaches involved into communication
  integer(i4b)                            :: iSeg, jSeg            ! reach looping index

  ierr=0; message='mpi_comm_single_flux/'

  ! Number of total reaches to be communicated
  nSeg = sum(nReach)

  if (commType == scatter) then

    if (masterproc) then

      allocate(flux_global_tmp(nSeg), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux_global_tmp]'; return; endif

      do iSeg =1,nSeg ! Loop through global reach
       jSeg = rchIdxGlobal(iSeg)
       flux_global_tmp(iSeg) = flux_global(jSeg)     !
      enddo

    end if ! end of root processor operation

    call shr_mpi_barrier(comm, message)

    ! Distribute global flux data to each process
    allocate(flux_local_tmp(nReach(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux_local_tmp]'; return; endif

    call shr_mpi_scatterV(flux_global_tmp, nReach, flux_local_tmp, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! put flux at correct location in local flux array
    do iSeg =1,nReach(pid) ! Loop through reaches per proc
      jSeg = rchIdxLocal(iSeg)
      flux_local(jSeg) = flux_local_tmp(iSeg)
    end do

  elseif (commType == gather) then

    allocate(flux_local_tmp(nReach(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux_local_tmp]'; return; endif

    do iSeg =1,nReach(pid)
      jSeg = rchIdxLocal(iSeg)
      flux_local_tmp(iSeg) = flux_local(jSeg)
    end do

    allocate(flux_global_tmp(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux_global_tmp]'; return; endif

    call shr_mpi_gatherV(flux_local_tmp, nReach, flux_global_tmp, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! put flux at correct location in global flux array
    if (masterproc) then
      do iSeg =1,nSeg ! Loop through all the reaches involved into communication
        jSeg = rchIdxGlobal(iSeg)
        flux_global(jSeg) = flux_global_tmp(iSeg)
      end do
    endif

  endif

 end subroutine mpi_comm_single_flux

 ! *********************************************************************
 ! subroutine: all fluxes communication
 ! *********************************************************************
 subroutine mpi_comm_flux(pid,          &
                          nNodes,       &
                          comm,         & ! input: communicator
                          iens,         &
                          nReach,       &
                          RCHFLX_global,&
                          RCHFLX_local, &
                          rchIdxGlobal, &
                          rchIdxLocal,  &
                          commType,     &
                          ierr, message)

  USE dataTypes,  ONLY: STRFLX              ! reach flux data structure

  ! input variables
  integer(i4b),             intent(in)    :: pid                   ! process id (MPI)
  integer(i4b),             intent(in)    :: nNodes                ! number of processes (MPI)
  integer(i4b),             intent(in)    :: comm                  ! communicator
  integer(i4b),             intent(in)    :: iens                  ! ensemble index
  integer(i4b),             intent(in)    :: nReach(0:nNodes-1)    ! number of reaches communicate per node (dimension size == number of proc)
  type(STRFLX),allocatable, intent(inout) :: RCHFLX_global(:,:)
  type(STRFLX),allocatable, intent(inout) :: RCHFLX_local(:,:)
  integer(i4b),             intent(in)    :: rchIdxGlobal(:)       ! reach indices (w.r.t. global) to be transfer (dimension size == sum of nRearch)
  integer(i4b),             intent(in)    :: rchIdxLocal(:)        ! reach indices (w.r.t. local) (dimension size depends on procs )
  integer(i4b),             intent(in)    :: commType              ! communication type 1->scatter, 2->gather otherwise error
  ! output variables
  integer(i4b),             intent(out)   :: ierr                  ! error code
  character(len=strLen),    intent(out)   :: message               ! error message
  ! local variables
  character(len=strLen)                   :: cmessage              ! error message from a subroutine
  real(dp),     allocatable               :: flux(:,:)             !
  real(dp),     allocatable               :: flux_local(:,:)       !
  real(dp),     allocatable               :: vec_out(:)            ! output vector from mpi gather/scatter routine
  integer(i4b)                            :: nSeg                  ! number of reaches
  integer(i4b)                            :: iSeg, jSeg            ! reach looping index
  integer(i4b)                            :: ix                    ! general looping index
  integer(i4b),parameter                  :: nFluxes=6             ! number of flux variables

  ierr=0; message='mpi_comm_flux/'

  ! Number of total reaches to be communicated
  nSeg = sum(nReach)

  if (commType == scatter) then

    if (masterproc) then

      allocate(flux(nSeg, nFluxes), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux]'; return; endif

      do iSeg =1,nSeg ! Loop through tributary reaches
       jSeg = rchIdxGlobal(iSeg)
       flux(iSeg,1) = RCHFLX_global(iens,jSeg)%BASIN_QR(0)  ! HRU routed flow (previous time step)
       flux(iSeg,2) = RCHFLX_global(iens,jSeg)%BASIN_QR(1)  ! HRU routed flow (current time step)
       flux(iSeg,3) = RCHFLX_global(iens,jSeg)%UPSTREAM_QI  ! Upstream accumulated flow
       flux(iSeg,4) = RCHFLX_global(iens,jSeg)%REACH_Q      ! KWT routed flow
       flux(iSeg,5) = RCHFLX_global(iens,jSeg)%REACH_Q_IRF  ! IRF routed flow
       flux(iSeg,6) = RCHFLX_global(iens,jSeg)%BASIN_QI     ! HRU non-routed flow
      enddo

    end if ! end of root processor operation

    call shr_mpi_barrier(comm, message)

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

      RCHFLX_local(iens,jSeg)%BASIN_QR(0) = flux_local(iSeg,1)  ! HRU routed flow (previous time step)
      RCHFLX_local(iens,jSeg)%BASIN_QR(1) = flux_local(iSeg,2)  ! HRU routed flow (current time step)
      RCHFLX_local(iens,jSeg)%UPSTREAM_QI = flux_local(iSeg,3)  ! Upstream accumulated flow
      RCHFLX_local(iens,jSeg)%REACH_Q     = flux_local(iSeg,4)  ! KWT routed flow
      RCHFLX_local(iens,jSeg)%REACH_Q_IRF = flux_local(iSeg,5)  ! IRF routed flow
      RCHFLX_local(iens,jSeg)%BASIN_QI    = flux_local(iSeg,6)  ! HRU non-routed flow

    end do

  elseif (commType == gather) then

    allocate(flux_local(nReach(pid),nFluxes), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [flux_local]'; return; endif

    do iSeg =1,nReach(pid)  ! Loop through (selected) tributary reaches

      jSeg = rchIdxLocal(iSeg)
      ! Transfer reach fluxes to 2D arrays
      flux_local(iSeg,1) = RCHFLX_local(iens,jSeg)%BASIN_QR(0)  ! HRU routed flow (previous time step)
      flux_local(iSeg,2) = RCHFLX_local(iens,jSeg)%BASIN_QR(1)  ! HRU routed flow (current time step)
      flux_local(iSeg,3) = RCHFLX_local(iens,jSeg)%UPSTREAM_QI  ! Upstream accumulated flow
      flux_local(iSeg,4) = RCHFLX_local(iens,jSeg)%REACH_Q      ! KWT routed flow
      flux_local(iSeg,5) = RCHFLX_local(iens,jSeg)%REACH_Q_IRF  ! IRF routed flow
      flux_local(iSeg,6) = RCHFLX_local(iens,jSeg)%BASIN_QI     ! non-HRU routed flow (

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
    if (masterproc) then
      do iSeg =1,nSeg ! Loop through all the reaches involved into communication

        jSeg = rchIdxGlobal(iSeg)
        RCHFLX_global(iens,jSeg)%BASIN_QR(0) = flux(iSeg,1)
        RCHFLX_global(iens,jSeg)%BASIN_QR(1) = flux(iSeg,2)
        RCHFLX_global(iens,jSeg)%UPSTREAM_QI = flux(iSeg,3)
        RCHFLX_global(iens,jSeg)%REACH_Q     = flux(iSeg,4)
        RCHFLX_global(iens,jSeg)%REACH_Q_IRF = flux(iSeg,5)
        RCHFLX_global(iens,jSeg)%BASIN_QI    = flux(iSeg,6)
      end do
    endif

  endif

 end subroutine mpi_comm_flux

 ! *********************************************************************
 ! subroutine: basin IRF state communication
 ! *********************************************************************
 subroutine mpi_comm_irf_bas_state(pid,          &
                                   nNodes,       &
                                   comm,         & ! input: communicator
                                   iens,         &
                                   nReach,       &
                                   RCHFLX_global,&
                                   RCHFLX_local, &
                                   rchIdxGlobal, &
                                   rchIdxLocal,  &
                                   commType,     &
                                   ierr, message)

  USE dataTypes,  ONLY: STRFLX              ! reach flux data structure

  ! input variables
  integer(i4b),             intent(in)    :: pid                   ! process id (MPI)
  integer(i4b),             intent(in)    :: nNodes                ! number of processes (MPI)
  integer(i4b),             intent(in)    :: comm                  ! communicator
  integer(i4b),             intent(in)    :: iens                  ! ensemble index
  integer(i4b),             intent(in)    :: nReach(0:nNodes-1)    ! number of reaches communicate per node (dimension size == number of proc)
  type(STRFLX),allocatable, intent(inout) :: RCHFLX_global(:,:)
  type(STRFLX),allocatable, intent(inout) :: RCHFLX_local(:,:)
  integer(i4b),             intent(in)    :: rchIdxGlobal(:)       ! reach indices (w.r.t. global) to be transfer (dimension size == sum of nRearch)
  integer(i4b),             intent(in)    :: rchIdxLocal(:)        ! reach indices (w.r.t. local) (dimension size depends on procs )
  integer(i4b),             intent(in)    :: commType              ! communication type 1->scatter, 2->gather otherwise error
  ! output variables
  integer(i4b),             intent(out)   :: ierr                  ! error code
  character(len=strLen),    intent(out)   :: message               ! error message
  ! local variables
  character(len=strLen)                   :: cmessage              ! error message from a subroutine
  type(STRFLX), allocatable               :: RCHFLX0(:,:)          ! temp RCHFLX data structure to hold updated states
  real(dp),     allocatable               :: qfuture(:),qfuture_trib(:)
  integer(i4b)                            :: ix1, ix2
  integer(i4b)                            :: myid
  integer(i4b)                            :: nSeg                  ! number of reaches
  integer(i4b)                            :: iSeg, jSeg
  integer(i4b)                            :: ixTdh
  integer(i4b), allocatable               :: ntdh(:)
  integer(i4b), allocatable               :: ntdh_trib(:)
  integer(i4b)                            :: totTdh(0:nNodes-1)

  ierr=0; message='mpi_comm_irf_bas_state/'

  ! Number of total reaches to be communicated
  nSeg = sum(nReach)

  if (commType == scatter) then

    ! allocate nWave (number the same at all procs) at each proc
    allocate(ntdh(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ntdh]'; return; endif

    if (masterproc) then

     ! extract only tributary reaches
     allocate(RCHFLX0(1,nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [RCHFLX0]'; return; endif
     do iSeg =1,nSeg ! Loop through tributary reaches
      jSeg = rchIdxGlobal(iSeg)
      RCHFLX0(1, iSeg) = RCHFLX_global(iens,jSeg)
     enddo

     ! convert RCHFLX data strucutre to state arrays
     call irf_bas_struc2array(iens, RCHFLX0,  & !input: input state data structure
                          qfuture,            & !output: states array
                          ntdh,               & !output: number of qfuture per reach
                          ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif ! end of root process

    call shr_mpi_barrier(comm, message)

    ! will have to broadcast updated ntdh to all proc
    call MPI_BCAST(ntdh, nSeg, MPI_INTEGER, root, comm, ierr)

    ! total waves from all the tributary reaches in each proc
    ix2=0
    do myid = 0, nNodes-1
      ix1=ix2+1
      ix2=ix1+nReach(myid)-1
      totTdh(myid) = sum(ntdh(ix1:ix2))
    enddo

    call shr_mpi_scatterV(ntdh, nReach(0:nNodes-1), ntdh_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Distribute modified RCHFLX data to each process
    call shr_mpi_scatterV(qfuture, totTdh(0:nNodes-1), qfuture_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! update RCHFLX_local data structure
    ixTdh=1
    do iSeg =1,nReach(pid) ! Loop through reaches per proc

     jSeg = rchIdxLocal(iSeg)

     if (allocated(RCHFLX_local(iens,jSeg)%QFUTURE)) then
      deallocate(RCHFLX_local(iens,jSeg)%QFUTURE, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [RCHFLX_local(iens,jSeg)%QFUTURE_IRF]'; return; endif
     endif

     allocate(RCHFLX_local(iens,jSeg)%QFUTURE(1:ntdh_trib(iSeg)),stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [RCHFLX_local(iens,iRch)%QFUTURE_IRF]'; return; endif

     RCHFLX_local(iens,jSeg)%QFUTURE(1:ntdh_trib(iSeg)) = qfuture_trib(ixTdh:ixTdh+ntdh_trib(iSeg)-1)

     ixTdh=ixTdh+ntdh_trib(iSeg) !update 1st idex of array

    end do

  elseif (commType == gather) then

    ! allocate ntdh (number the same at all procs) and ntdh_trib (number dependent on proc) at each proc
    allocate(ntdh_trib(nReach(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ntdh_trib]'; return; endif

    ! extract only tributary reaches
    allocate(RCHFLX0(1,nReach(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [RCHFLX0]'; return; endif
    do iSeg =1,nReach(pid)  ! Loop through tributary reaches
      jSeg = rchIdxLocal(iSeg)
      RCHFLX0(1, iSeg) = RCHFLX_local(iens,jSeg)
    enddo

    ! Transfer KWT state data structure to flat arrays
    call irf_bas_struc2array(iens,RCHFLX0,       &
                         qfuture_trib,       &
                         ntdh_trib,          &
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_gatherV(ntdh_trib, nReach, ntdh, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call MPI_BCAST(ntdh, nSeg, MPI_INTEGER, root, comm, ierr)

    ! total waves in reaches in each proc
    ix2=0
    do myid = 0, nNodes-1
      ix1=ix2+1
      ix2=ix1+nReach(myid)-1
      totTdh(myid) = sum(ntdh(ix1:ix2))
    enddo

    call shr_mpi_gatherV(qfuture_trib, totTdh, qfuture, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! put it in global RCHFLX data structure
    if (masterproc) then
      ixTdh=1
      do iSeg =1,nSeg ! Loop through all the reaches involved into communication

        jSeg = rchIdxGlobal(iSeg)

        if (allocated(RCHFLX_global(iens,jSeg)%QFUTURE)) then
          deallocate(RCHFLX_global(iens,jSeg)%QFUTURE, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [RCHFLX_global(iens,jSeg)%QFUTURE]'; return; endif
        endif

        allocate(RCHFLX_global(iens,jSeg)%QFUTURE(1:ntdh(iSeg)),stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [RCHFLX_global(iens,iRch)%QFUTURE]'; return; endif

        RCHFLX_global(iens,jSeg)%QFUTURE(1:ntdh(iSeg)) = qfuture(ixTdh:ixTdh+ntdh(iSeg)-1)
        ixTdh=ixTdh+ntdh(iSeg) !update 1st idex of array
      end do
    endif

  endif

 end subroutine mpi_comm_irf_bas_state

 ! *********************************************************************
 ! subroutine: IRF state communication
 ! *********************************************************************
 subroutine mpi_comm_irf_state(pid,          &
                               nNodes,       &
                               comm,         & ! input: communicator
                               iens,         &
                               nReach,       &
                               RCHFLX_global,&
                               RCHFLX_local, &
                               rchIdxGlobal, &
                               rchIdxLocal,  &
                               commType,     &
                               ierr, message)

  USE dataTypes,  ONLY: STRFLX              ! reach flux data structure

  ! input variables
  integer(i4b),             intent(in)    :: pid                   ! process id (MPI)
  integer(i4b),             intent(in)    :: nNodes                ! number of processes (MPI)
  integer(i4b),             intent(in)    :: comm                  ! communicator
  integer(i4b),             intent(in)    :: iens                  ! ensemble index
  integer(i4b),             intent(in)    :: nReach(0:nNodes-1)    ! number of reaches communicate per node (dimension size == number of proc)
  type(STRFLX),allocatable, intent(inout) :: RCHFLX_global(:,:)
  type(STRFLX),allocatable, intent(inout) :: RCHFLX_local(:,:)
  integer(i4b),             intent(in)    :: rchIdxGlobal(:)       ! reach indices (w.r.t. global) to be transfer (dimension size == sum of nRearch)
  integer(i4b),             intent(in)    :: rchIdxLocal(:)        ! reach indices (w.r.t. local) (dimension size depends on procs )
  integer(i4b),             intent(in)    :: commType              ! communication type 1->scatter, 2->gather otherwise error
  ! output variables
  integer(i4b),             intent(out)   :: ierr                  ! error code
  character(len=strLen),    intent(out)   :: message               ! error message
  ! local variables
  character(len=strLen)                   :: cmessage              ! error message from a subroutine
  type(STRFLX), allocatable               :: RCHFLX0(:,:)          ! temp RCHFLX data structure to hold updated states
  real(dp),     allocatable               :: qfuture(:),qfuture_trib(:)
  integer(i4b)                            :: ix1, ix2
  integer(i4b)                            :: myid
  integer(i4b)                            :: nSeg                  ! number of reaches
  integer(i4b)                            :: iSeg, jSeg
  integer(i4b)                            :: ixTdh
  integer(i4b), allocatable               :: ntdh(:)
  integer(i4b), allocatable               :: ntdh_trib(:)
  integer(i4b)                            :: totTdh(0:nNodes-1)

  ierr=0; message='mpi_comm_irf_state/'

  ! Number of total reaches to be communicated
  nSeg = sum(nReach)

  if (commType == scatter) then

    ! allocate nWave (number the same at all procs) at each proc
    allocate(ntdh(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ntdh]'; return; endif

    if (masterproc) then

     ! extract only tributary reaches
     allocate(RCHFLX0(1,nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [RCHFLX0]'; return; endif
     do iSeg =1,nSeg ! Loop through tributary reaches
      jSeg = rchIdxGlobal(iSeg)
      RCHFLX0(1, iSeg) = RCHFLX_global(iens,jSeg)
     enddo

     ! convert RCHFLX data strucutre to state arrays
     call irf_struc2array(iens, RCHFLX0,  & !input: input state data structure
                          qfuture,        & !output: states array
                          ntdh,           & !output: number of qfuture per reach
                          ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif ! end of root process

    call shr_mpi_barrier(comm, message)

    ! will have to broadcast updated ntdh to all proc
    call MPI_BCAST(ntdh, nSeg, MPI_INTEGER, root, comm, ierr)

    ! total waves from all the tributary reaches in each proc
    ix2=0
    do myid = 0, nNodes-1
      ix1=ix2+1
      ix2=ix1+nReach(myid)-1
      totTdh(myid) = sum(ntdh(ix1:ix2))
    enddo

    call shr_mpi_scatterV(ntdh, nReach(0:nNodes-1), ntdh_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Distribute modified RCHFLX data to each process
    call shr_mpi_scatterV(qfuture, totTdh(0:nNodes-1), qfuture_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! update RCHFLX_local data structure
    ixTdh=1
    do iSeg =1,nReach(pid) ! Loop through reaches per proc

     jSeg = rchIdxLocal(iSeg)

     if (allocated(RCHFLX_local(iens,jSeg)%QFUTURE_IRF)) then
      deallocate(RCHFLX_local(iens,jSeg)%QFUTURE_IRF, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [RCHFLX_local(iens,jSeg)%QFUTURE_IRF]'; return; endif
     endif

     allocate(RCHFLX_local(iens,jSeg)%QFUTURE_IRF(1:ntdh_trib(iSeg)),stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [RCHFLX_local(iens,iRch)%QFUTURE_IRF]'; return; endif

     RCHFLX_local(iens,jSeg)%QFUTURE_IRF(1:ntdh_trib(iSeg)) = qfuture_trib(ixTdh:ixTdh+ntdh_trib(iSeg)-1)

     ixTdh=ixTdh+ntdh_trib(iSeg) !update 1st idex of array

    end do

  elseif (commType == gather) then

    ! allocate ntdh (number the same at all procs) and ntdh_trib (number dependent on proc) at each proc
    allocate(ntdh_trib(nReach(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ntdh_trib]'; return; endif

    ! extract only tributary reaches
    allocate(RCHFLX0(1,nReach(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [RCHFLX0]'; return; endif
    do iSeg =1,nReach(pid)  ! Loop through tributary reaches
      jSeg = rchIdxLocal(iSeg)
      RCHFLX0(1, iSeg) = RCHFLX_local(iens,jSeg)
    enddo

    ! Transfer KWT state data structure to flat arrays
    call irf_struc2array(iens,RCHFLX0,       &
                         qfuture_trib,       &
                         ntdh_trib,          &
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_gatherV(ntdh_trib, nReach, ntdh, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call MPI_BCAST(ntdh, nSeg, MPI_INTEGER, root, comm, ierr)

    ! total waves in reaches in each proc
    ix2=0
    do myid = 0, nNodes-1
      ix1=ix2+1
      ix2=ix1+nReach(myid)-1
      totTdh(myid) = sum(ntdh(ix1:ix2))
    enddo

    call shr_mpi_gatherV(qfuture_trib, totTdh, qfuture, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! put it in global RCHFLX data structure
    if (masterproc) then
      ixTdh=1
      do iSeg =1,nSeg ! Loop through all the reaches involved into communication

        jSeg = rchIdxGlobal(iSeg)

        if (allocated(RCHFLX_global(iens,jSeg)%QFUTURE_IRF)) then
          deallocate(RCHFLX_global(iens,jSeg)%QFUTURE_IRF, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [RCHFLX_global(iens,jSeg)%QFUTURE_IRF]'; return; endif
        endif

        allocate(RCHFLX_global(iens,jSeg)%QFUTURE_IRF(1:ntdh(iSeg)),stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [RCHFLX_global(iens,iRch)%QFUTURE_IRF]'; return; endif

        RCHFLX_global(iens,jSeg)%QFUTURE_IRF(1:ntdh(iSeg)) = qfuture(ixTdh:ixTdh+ntdh(iSeg)-1)
        ixTdh=ixTdh+ntdh(iSeg) !update 1st idex of array
      end do
    endif

  endif

 end subroutine mpi_comm_irf_state


 ! *********************************************************************
 ! subroutine: KWE state communication
 ! *********************************************************************
 subroutine mpi_comm_kwe_state(pid,          &
                               nNodes,       &
                               comm,         & ! input: communicator
                               iens,         &
                               nReach,       &
                               RCHSTA_global,&
                               RCHSTA_local, &
                               rchIdxGlobal, &
                               rchIdxLocal,  &
                               commType,     &
                               ierr, message)

  USE dataTypes,  ONLY: EKWRCH
  USE dataTypes,  ONLY: STRSTA

  ! input variables
  integer(i4b),             intent(in)    :: pid                   ! process id (MPI)
  integer(i4b),             intent(in)    :: nNodes                ! number of processes (MPI)
  integer(i4b),             intent(in)    :: comm                  ! communicator
  integer(i4b),             intent(in)    :: iens                  ! ensemble index
  integer(i4b),             intent(in)    :: nReach(0:nNodes-1)    ! number of reaches communicate per node (dimension size == number of proc)
  type(STRSTA), allocatable,intent(inout) :: RCHSTA_global(:,:)
  type(STRSTA), allocatable,intent(inout) :: RCHSTA_local(:,:)
  integer(i4b),             intent(in)    :: rchIdxGlobal(:)       ! reach indices (w.r.t. global) to be transfer (dimension size == sum of nRearch)
  integer(i4b),             intent(in)    :: rchIdxLocal(:)        ! reach indices (w.r.t. local) (dimension size depends on procs )
  integer(i4b),             intent(in)    :: commType              ! communication type 1->scatter, 2->gather otherwise error
  ! output variables
  integer(i4b),             intent(out)   :: ierr                  ! error code
  character(len=strLen),    intent(out)   :: message               ! error message
  ! local variables
  character(len=strLen)                   :: cmessage              ! error message from a subroutine
  type(EKWRCH), allocatable               :: EKW0(:,:)             ! temp KEW data structure to hold updated states
  real(dp),     allocatable               :: Q(:),Q_trib(:)
  real(dp),     allocatable               :: A(:),A_trib(:)
  integer(i4b)                            :: myid
  integer(i4b)                            :: nSeg                  ! number of reaches
  integer(i4b)                            :: iSeg, jSeg
  integer(i4b)                            :: ixMesh
  integer(i4b)                            :: totMesh(0:nNodes-1)
  integer(i4b)                            :: totMeshAll

  ierr=0; message='mpi_comm_kwe_state/'

  ! Number of total reaches to be communicated
  nSeg = sum(nReach)

  if (commType == scatter) then

    if (masterproc) then

     ! extract only tributary reaches
     allocate(EKW0(1,nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [EKW0]'; return; endif
     do iSeg =1,nSeg ! Loop through tributary reaches
      jSeg = rchIdxGlobal(iSeg)
      EKW0(1, iSeg) = RCHSTA_global(iens,jSeg)%EKW_ROUTE
     enddo

     ! convert EKWRCH data strucutre to state arrays
     call kwe_struc2array(iens, EKW0,   & !input: input state data structure
                          Q,A,          & !output: states array
                          ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif ! end of root process

    call shr_mpi_barrier(comm, message)

    ! total waves from all the tributary reaches in each proc
    do myid = 0, nNodes-1
      totMesh(myid) = 4*nReach(myid)
    enddo

    ! need to allocate global array to be scattered at the other tasks
    totMeshAll = nSeg*4
    if (.not.masterproc) then
     allocate(Q(totMeshAll), A(totMeshAll), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [Q,A]'; return; endif
    endif
    call shr_mpi_barrier(comm, message)

    ! Distribute modified LKW_ROUTE%KWAVE data to each process
    call shr_mpi_scatterV(Q, totMesh(0:nNodes-1), Q_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_scatterV(A, totMesh(0:nNodes-1), A_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! update RCHSTA_local%EKW_ROUTE data structure
    ixMesh=1
    do iSeg =1,nReach(pid) ! Loop through reaches per proc

     jSeg = rchIdxLocal(iSeg)

     RCHSTA_local(iens,jSeg)%EKW_ROUTE%Q(1:4) = Q_trib(ixMesh:ixMesh+3)
     RCHSTA_local(iens,jSeg)%EKW_ROUTE%A(1:4) = A_trib(ixMesh:ixMesh+3)

     ixMesh=ixMesh+4 !update 1st idex of array

    end do

  elseif (commType == gather) then


  endif

 end subroutine mpi_comm_kwe_state

 ! *********************************************************************
 ! subroutine: kinematic wave state communication
 ! *********************************************************************
 subroutine mpi_comm_kwt_state(pid,          & ! input:
                               nNodes,       & ! input: number of node
                               comm,         & ! input: communicator
                               iens,         &
                               nReach,       &
                               RCHSTA_global,&
                               RCHSTA_local, &
                               rchIdxGlobal, &
                               rchIdxLocal,  &
                               commType,     &
                               ierr, message)

  USE dataTypes,  ONLY: LKWRCH
  USE dataTypes,  ONLY: STRSTA

  ! input variables
  integer(i4b),             intent(in)    :: pid                   ! process id (MPI)
  integer(i4b),             intent(in)    :: nNodes                ! number of processes (MPI)
  integer(i4b),             intent(in)    :: comm                  ! communicator
  integer(i4b),             intent(in)    :: iens                  ! ensemble index
  integer(i4b),             intent(in)    :: nReach(0:nNodes-1)    ! number of reaches communicate per node (dimension size == number of proc)
  type(STRSTA), allocatable,intent(inout) :: RCHSTA_global(:,:)
  type(STRSTA), allocatable,intent(inout) :: RCHSTA_local(:,:)
  integer(i4b),             intent(in)    :: rchIdxGlobal(:)       ! reach indices (w.r.t. global) to be transfer (dimension size == sum of nRearch)
  integer(i4b),             intent(in)    :: rchIdxLocal(:)        ! reach indices (w.r.t. local) (dimension size depends on procs )
  integer(i4b),             intent(in)    :: commType              ! communication type 1->scatter, 2->gather otherwise error
  ! output variables
  integer(i4b),             intent(out)   :: ierr                  ! error code
  character(len=strLen),    intent(out)   :: message               ! error message
  ! local variables
  character(len=strLen)                   :: cmessage              ! error message from a subroutine
  type(LKWRCH), allocatable               :: KROUTE0(:,:)          ! temp KROUTE data structure to hold updated states
  real(dp),     allocatable               :: QF(:),QF_trib(:)
  real(dp),     allocatable               :: QM(:),QM_trib(:)
  real(dp),     allocatable               :: TI(:),TI_trib(:)
  real(dp),     allocatable               :: TR(:),TR_trib(:)
  logical(lgt), allocatable               :: RF(:),RF_trib(:)
  integer(i4b)                            :: ix1, ix2
  integer(i4b)                            :: myid
  integer(i4b)                            :: nSeg                  ! number of reaches
  integer(i4b)                            :: iSeg, jSeg
  integer(i4b)                            :: ixWave
  integer(i4b), allocatable               :: nWave(:)
  integer(i4b), allocatable               :: nWave_trib(:)
  integer(i4b)                            :: totWave(0:nNodes-1)
  integer(i4b)                            :: totWaveAll

  ierr=0; message='mpi_comm_kwt_state/'

  ! Number of total reaches to be communicated
  nSeg = sum(nReach)

  if (commType == scatter) then

    ! allocate nWave (number the same at all procs) at each proc
    allocate(nWave(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [nWave]'; return; endif

    if (masterproc) then

     ! extract only tributary reaches
     allocate(KROUTE0(1,nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE0]'; return; endif
     do iSeg =1,nSeg ! Loop through tributary reaches
      jSeg = rchIdxGlobal(iSeg)
      KROUTE0(1, iSeg) = RCHSTA_global(iens,jSeg)%LKW_ROUTE
     enddo

     ! convert KROUTE data strucutre to state arrays
     call kwt_struc2array(iens, KROUTE0,  & !input: input state data structure
                          QF,QM,TI,TR,RF, & !output: states array
                          nWave,          & !output: number of waves per reach
                          ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif ! end of root process

    call shr_mpi_barrier(comm, message)

    ! will have to broadcast updated nWave to all proc
    call MPI_BCAST(nWave, nSeg, MPI_INTEGER, root, comm, ierr)

    ! total waves from all the tributary reaches in each proc
    ix2=0
    do myid = 0, nNodes-1
      ix1=ix2+1
      ix2=ix1+nReach(myid)-1
      totWave(myid) = sum(nWave(ix1:ix2))
    enddo

    totWaveAll = sum(nWave)
    ! need to allocate global array to be scattered at the other tasks
    if (.not.masterproc) then
     allocate(QF(totWaveAll), QM(totWaveAll), TI(totWaveAll), TR(totWaveAll), RF(totWaveAll), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [QF,..,RF]'; return; endif
    endif
    call shr_mpi_barrier(comm, message)

    call shr_mpi_scatterV(nWave, nReach(0:nNodes-1), nWave_trib, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Distribute modified LKW_ROUTE%KWAVE data to each process
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

    ! update RCHSTA_local%LKW_ROUTE data structure
    ixWave=1
    do iSeg =1,nReach(pid) ! Loop through reaches per proc

     jSeg = rchIdxLocal(iSeg)

     if (allocated(RCHSTA_local(iens,jSeg)%LKW_ROUTE%KWAVE)) then
      deallocate(RCHSTA_local(iens,jSeg)%LKW_ROUTE%KWAVE, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [KROUTE_trib(iens,jSeg)%KWAVE]'; return; endif
     endif

     allocate(RCHSTA_local(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave_trib(iSeg)-1),stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE_trib(iens,iRch)%KWAVE]'; return; endif

     RCHSTA_local(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave_trib(iSeg)-1)%QF = QF_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
     RCHSTA_local(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave_trib(iSeg)-1)%QM = QM_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
     RCHSTA_local(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave_trib(iSeg)-1)%TI = TI_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
     RCHSTA_local(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave_trib(iSeg)-1)%TR = TR_trib(ixWave:ixWave+nWave_trib(iSeg)-1)
     RCHSTA_local(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave_trib(iSeg)-1)%RF = RF_trib(ixWave:ixWave+nWave_trib(iSeg)-1)

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
      KROUTE0(1, iSeg) = RCHSTA_local(iens,jSeg)%LKW_ROUTE
    enddo

    ! Transfer KWT state data structure to flat arrays
    call kwt_struc2array(iens,KROUTE0,                            &
                         QF_trib,QM_trib,TI_trib,TR_trib,RF_trib, &
                         nWave_trib,                              &
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call shr_mpi_gatherV(nWave_trib, nReach, nWave, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    call MPI_BCAST(nWave, nSeg, MPI_INTEGER, root, comm, ierr)

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
    if (masterproc) then
      ixWave=1
      do iSeg =1,nSeg ! Loop through all the reaches involved into communication

        jSeg = rchIdxGlobal(iSeg)

        if (allocated(RCHSTA_global(iens,jSeg)%LKW_ROUTE%KWAVE)) then
          deallocate(RCHSTA_global(iens,jSeg)%LKW_ROUTE%KWAVE, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [KROUTE(iens,jSeg)%KWAVE]'; return; endif
        endif

        allocate(RCHSTA_global(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave(iSeg)-1),stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating array for [KROUTE_out(iens,iRch)%KWAVE]'; return; endif

        RCHSTA_global(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave(iSeg)-1)%QF = QF(ixWave:ixWave+nWave(iSeg)-1)
        RCHSTA_global(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave(iSeg)-1)%QM = QM(ixWave:ixWave+nWave(iSeg)-1)
        RCHSTA_global(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave(iSeg)-1)%TI = TI(ixWave:ixWave+nWave(iSeg)-1)
        RCHSTA_global(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave(iSeg)-1)%TR = TR(ixWave:ixWave+nWave(iSeg)-1)
        RCHSTA_global(iens,jSeg)%LKW_ROUTE%KWAVE(0:nWave(iSeg)-1)%RF = RF(ixWave:ixWave+nWave(iSeg)-1)
        ixWave=ixWave+nWave(iSeg) !update 1st idex of array
      end do
    endif

  endif

 end subroutine mpi_comm_kwt_state

 ! *********************************************************************
 ! private subroutine
 ! *********************************************************************
 subroutine irf_bas_struc2array(iens, RCHFLX_in,    &  ! input:
                                qfuture_bas,        &  ! output:
                                ntdh_bas,               &  ! output:
                                ierr, message)
  USE dataTypes, ONLY: STRFLX              ! reach flux data structure
  implicit none
  ! Input
  integer(i4b),          intent(in)              :: iens           ! ensemble index
  type(STRFLX),          intent(in), allocatable :: RCHFLX_in(:,:) ! reach state data
  ! Output error handling variables
  real(dp),              intent(out),allocatable :: qfuture_bas(:) ! flat array for wave Q
  integer(i4b),          intent(out),allocatable :: ntdh_bas(:)    ! number of waves at each reach
  integer(i4b),          intent(out)             :: ierr           ! error code
  character(len=strLen), intent(out)             :: message        ! error message
  ! local variables
  integer(i4b)                                   :: ixTdh          ! 1st indix of each reach
  integer(i4b)                                   :: iSeg           ! loop indix
  integer(i4b)                                   :: nSeg           ! number of reaches
  integer(i4b)                                   :: totTdh         ! total number of waves from all the reaches

  ierr=0; message='irf_bas_struc2array/'

  nSeg = size(RCHFLX_in(iens,:))
  if (.not.allocated(ntdh_bas)) then
    allocate(ntdh_bas(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ntdh_bas]'; return; endif
  end if

  do iSeg = 1,nSeg
   ntdh_bas(iSeg) = size(RCHFLX_in(iens,iSeg)%QFUTURE)
  enddo

  totTdh=sum(ntdh_bas)
  allocate(qfuture_bas(totTdh), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [qfuture_bas]'; return; endif

  ixTdh = 1
  do iSeg=1,nSeg
   qfuture_bas(ixTdh:ixTdh+ntdh_bas(iSeg)-1) = RCHFLX_in(iens,iSeg)%QFUTURE(1:ntdh_bas(iSeg))
   ixTdh = ixTdh+ntdh_bas(iSeg)
  end do

 end subroutine irf_bas_struc2array


 ! *********************************************************************
 ! private subroutine
 ! *********************************************************************
 subroutine irf_struc2array(iens, RCHFLX_in,    &  ! input:
                            qfuture,            &  ! output:
                            ntdh,               &  ! output:
                            ierr, message)
  USE dataTypes, ONLY: STRFLX              ! reach flux data structure
  implicit none
  ! Input
  integer(i4b),          intent(in)              :: iens           ! ensemble index
  type(STRFLX),          intent(in), allocatable :: RCHFLX_in(:,:) ! reach state data
  ! Output error handling variables
  real(dp),              intent(out),allocatable :: qfuture(:)     ! flat array for wave Q
  integer(i4b),          intent(out),allocatable :: ntdh(:)        ! number of waves at each reach
  integer(i4b),          intent(out)             :: ierr           ! error code
  character(len=strLen), intent(out)             :: message        ! error message
  ! local variables
  integer(i4b)                                   :: ixTdh          ! 1st indix of each reach
  integer(i4b)                                   :: iSeg           ! loop indix
  integer(i4b)                                   :: nSeg           ! number of reaches
  integer(i4b)                                   :: totTdh         ! total number of waves from all the reaches

  ierr=0; message='irf_struc2array/'

  nSeg = size(RCHFLX_in(iens,:))
  if (.not.allocated(ntdh)) then
    allocate(ntdh(nSeg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [ntdh]'; return; endif
  end if

  do iSeg = 1,nSeg
   ntdh(iSeg) = size(RCHFLX_in(iens,iSeg)%QFUTURE_IRF)
  enddo

  totTdh=sum(ntdh)
  allocate(qfuture(totTdh), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [qfuture]'; return; endif

  ixTdh = 1
  do iSeg=1,nSeg
   qfuture(ixTdh:ixTdh+ntdh(iSeg)-1) = RCHFLX_in(iens,iSeg)%QFUTURE_IRF(1:ntdh(iSeg))
   ixTdh = ixTdh+ntdh(iSeg)
  end do

 end subroutine irf_struc2array

 ! *********************************************************************
 ! private subroutine
 ! *********************************************************************
 subroutine kwt_struc2array(iens, KROUTE_in,     &  ! input:
                            QF,QM,TI,TR,RF,      &  ! output:
                            nWave,               &
                            ierr, message)
  USE dataTypes, ONLY: LKWRCH             ! collection of particles in a given reach
  implicit none
  ! Input
  integer(i4b),          intent(in)              :: iens           ! ensemble index
  type(LKWRCH),          intent(in), allocatable :: KROUTE_in(:,:) ! reach state data
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
 ! private subroutine
 ! *********************************************************************
 subroutine kwe_struc2array(iens, EKW_in,     &  ! input:
                            Q,A,              &  ! output:
                            ierr, message)
  USE dataTypes, ONLY: EKWRCH             ! collection of particles in a given reach
  implicit none
  ! Input
  integer(i4b),          intent(in)              :: iens           ! ensemble index
  type(EKWRCH),          intent(in), allocatable :: EKW_in(:,:)    ! reach state data
  ! Output error handling variables
  real(dp),              intent(out),allocatable :: Q(:)           ! flat array for wave Q
  real(dp),              intent(out),allocatable :: A(:)           ! Flat array for modified Q
  integer(i4b),          intent(out)             :: ierr           ! error code
  character(len=strLen), intent(out)             :: message        ! error message
  ! local variables
  integer(i4b)                                   :: ixMesh         ! 1st indix of each reach
  integer(i4b)                                   :: iSeg           ! loop indix
  integer(i4b)                                   :: nSeg           ! number of reaches
  integer(i4b)                                   :: totMesh        ! total number of waves from all the reaches

  ierr=0; message='kwe_struc2array/'

  nSeg = size(EKW_in(iens,:))

  totMesh = 4*nSeg

  allocate(Q(totMesh),A(totMesh), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [Q,A]'; return; endif

  ixMesh = 1
  do iSeg=1,nSeg
   Q(ixMesh:ixMesh+3) = EKW_in(iens,iSeg)%Q(1:4)
   A(ixMesh:ixMesh+3) = EKW_in(iens,iSeg)%A(1:4)
   ixMesh = ixMesh+4
  end do

 end subroutine kwe_struc2array

 ! *********************************************************************
 ! public subroutine: send global data
 ! *********************************************************************
 ! send all the necessary public/global variables neccesary in task
 subroutine pass_global_data(comm, ierr,message)   ! output: error control
  USE globalData, ONLY: nRch,nHRU         ! number of reaches and hrus in whole network
  USE globalData, ONLY: timeVar           ! time variable
  USE globalData, ONLY: iTime             ! time index
  USE globalData, ONLY: refJulday         ! julian day: reference
  USE globalData, ONLY: startJulday       ! julian day: start
  USE globalData, ONLY: endJulday         ! julian day: end
  USE globalData, ONLY: modJulday         ! julian day: at simulation time step
  USE globalData, ONLY: roJulday          ! julian day: runoff input time
  USE globalData, ONLY: reachID
  USE globalData, ONLY: basinID
  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: comm    ! communicator
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message ! error message
  ! Local variables
  ! None

  ierr=0; message='pass_global_data/'

  ! send scalars
  call MPI_BCAST(iTime,       1,     MPI_INTEGER,          root, comm, ierr)
  call MPI_BCAST(nRch,        1,     MPI_INTEGER,          root, comm, ierr)
  call MPI_BCAST(nHRU,        1,     MPI_INTEGER,          root, comm, ierr)
  call MPI_BCAST(calendar,  strLen,  MPI_CHARACTER,        root, comm, ierr)
  call MPI_BCAST(time_units,strLen,  MPI_CHARACTER,        root, comm, ierr)
  call MPI_BCAST(refJulday,     1,   MPI_DOUBLE_PRECISION, root, comm, ierr)
  call MPI_BCAST(startJulday,   1,   MPI_DOUBLE_PRECISION, root, comm, ierr)
  call MPI_BCAST(endJulday,     1,   MPI_DOUBLE_PRECISION, root, comm, ierr)
  call MPI_BCAST(modJulday,     1,   MPI_DOUBLE_PRECISION, root, comm, ierr)

  call shr_mpi_bcast(roJulday,ierr, message)
  call shr_mpi_bcast(timeVar,ierr, message)
  call shr_mpi_bcast(reachID,ierr, message)
  call shr_mpi_bcast(basinID,ierr, message)

 end subroutine pass_global_data

end module mpi_routine

