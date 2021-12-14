module globalData
  ! This module includes shared data
  USE pio

  USE public_var, ONLY: integerMissing

  ! data types
  USE nrtype

  ! metadata types
  USE dataTypes, ONLY: struct_info   ! metadata type
  USE dataTypes, ONLY: dim_info      ! metadata type
  USE dataTypes, ONLY: var_info      ! metadata type
  USE objTypes,  ONLY: var_info_new  ! metadata type

  ! input data file info
  USE dataTypes, ONLY: infileinfo    ! strture for getting the information of input files

  ! parameter structures
  USE dataTypes, ONLY: RCHPRP        ! Reach parameters (properties)
  USE dataTypes, ONLY: RCHTOPO       ! Network topology

  ! routing structures
  USE dataTypes, ONLY: STRFLX        ! fluxes in each reach
  USE dataTypes, ONLY: STRSTA        ! states in each reach

  ! lake structures
  USE dataTypes, ONLY: LAKPRP        ! lake properties
  USE dataTypes, ONLY: LAKTOPO       ! lake topology
  USE dataTypes, ONLY: LKFLX         ! lake fluxes

  ! remapping structures
  USE dataTypes, ONLY: remap         ! remapping data type
  USE dataTypes, ONLY: runoff        ! runoff data type
  USE dataTypes, ONLY: wm            ! water management (flux to/from segment, target volume) data type

  ! basin data structure
  USE dataTypes, ONLY: subbasin_omp  ! mainstem+tributary data structures
  USE dataTypes, ONLY: subbasin_mpi  ! mpi domain data structure (list of segments and hru and type of basins)

  ! time data structure
  USE datetime_data, ONLY: datetime  ! datetime data class

  ! number of variables for data structure
  USE var_lookup, ONLY: nStructures
  USE var_lookup, ONLY: nDimensions
  USE var_lookup, ONLY: nStateDims
  USE var_lookup, ONLY: nQdims
  USE var_lookup, ONLY: nVarsHRU
  USE var_lookup, ONLY: nVarsHRU2SEG
  USE var_lookup, ONLY: nVarsSEG
  USE var_lookup, ONLY: nVarsNTOPO
  USE var_lookup, ONLY: nVarsPFAF
  USE var_lookup, ONLY: nVarsRFLX
  USE var_lookup, ONLY: nVarsBasinQ
  USE var_lookup, ONLY: nVarsIRFbas
  USE var_lookup, ONLY: nVarsIRF
  USE var_lookup, ONLY: nVarsKWT
  USE var_lookup, ONLY: nVarsKWE

  implicit none

  save

  ! ---------- Catchment and reach IDs  -------------------------------------------------------------------------

  integer(i4b)                   , public :: nHRU                 ! number of HRUs in the whole river network
  integer(i4b)                   , public :: nContribHRU          ! number of HRUs that are connected to any reaches
  integer(i4b)                   , public :: nRch                 ! number of reaches in the whole river network
  integer(i4b)                   , public :: nRch_mainstem        ! number of reaches in mainstems
  integer(i4b)                   , public :: nHRU_mainstem        ! number of HRUs in mainstems
  integer(i4b)    , allocatable  , public :: basinID(:)           ! HRU id
  integer(i4b)    , allocatable  , public :: reachID(:)           ! reach id

  ! ---------- Date/Time data  -------------------------------------------------------------------------

  integer(i4b)                   , public :: iTime                ! time index at simulation time step
  real(dp)        , allocatable  , public :: timeVar(:)           ! time variables (unit given by runoff data)
  real(dp)                       , public :: TSEC(0:1)            ! begning and end of time step since simulation started (sec)
  type(datetime)                 , public :: simDatetime(0:1)     ! previous and current simulation time (yyyy:mm:dd:hh:mm:ss)
  type(datetime)                 , public :: begDatetime          ! simulation start date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime)                 , public :: endDatetime          ! simulation end date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime)                 , public :: restDatetime         ! desired restart date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime)                 , public :: dropDatetime         ! restart dropoff date/time (yyyy:mm:dd:hh:mm:ss)

  ! ---------- input file information -------------------------------------------------------------------

  type(infileinfo), allocatable,  public :: infileinfo_data(:)    ! input runoff file information
  type(infileinfo), allocatable,  public :: infileinfo_data_wm(:) ! input water management (abstaction/injection) file information

  ! ---------- Misc. data -------------------------------------------------------------------------
  ! standalone mode
  logical(lgt)                   , public :: isStandalone=.true.       ! flag to indicate model is running in standalone mode (True), otherwise coupled mode

  ! I/O stuff
  logical(lgt)                   , public :: isFileOpen                ! flag to indicate output netcdf is open
  integer(i4b)                   , public :: ixPrint(1:2)=integerMissing    ! index of desired reach to be on-screen print
  ! ennsemble number (maybe to be removed)
  integer(i4b)                   , public :: nEns=1                    ! number of ensemble

  ! ---------- MPI/OMP/PIO variables ----------------------------------------------------------------

  integer(i4b)                   , public :: mpicom_route              ! communicator for this program
  integer(i4b)                   , public :: pid                       ! process id
  integer(i4b)                   , public :: nNodes                    ! number of MPI processors
  integer(i4b)                   , public :: nThreads                  ! number of threads
  logical(lgt)                   , public :: masterproc                ! root logical. root processor => true, other => false
  logical(lgt)                   , public :: multiProcs                ! MPI multi-processors logical. => number of processors>1 true, other => false
  character(len=strLen)          , public :: pio_netcdf_format = "64bit_offset"
  character(len=strLen)          , public :: pio_typename      = "pnetcdf"
  integer(i4b)                   , public :: pio_numiotasks    = -99
  integer(i4b)                   , public :: pio_rearranger    = 2     ! 0=>PIO_rearr_none 1=> PIO_rearr_box 2=> PIO_rearr_subset
  integer(i4b)                   , public :: pio_root          = 1
  integer(i4b)                   , public :: pio_stride        = 1
  type(iosystem_desc_t)          , public :: pioSystem                  ! PIO I/O system data

  ! ---------- conversion factors -------------------------------------------------------------------

  real(dp)                       , public :: time_conv                  ! time conversion factor -- used to convert to mm/s
  real(dp)                       , public :: length_conv                ! length conversion factor -- used to convert to mm/s

  ! ---------- routing parameter names -------------------------------------------------------------------
  ! spatial constant ....
  real(dp)                       , public :: fshape                     ! shape parameter in time delay histogram (=gamma distribution) [-]
  real(dp)                       , public :: tscale                     ! scaling factor for the time delay histogram [sec]
  real(dp)                       , public :: velo                       ! velocity [m/s] for Saint-Venant equation   added by NM
  real(dp)                       , public :: diff                       ! diffusivity [m2/s] for Saint-Venant equation   added by NM
  real(dp)                       , public :: mann_n                     ! manning's roughness coefficient [unitless]  added by NM
  real(dp)                       , public :: wscale                     ! scaling factor for river width [-] added by NM

  ! ---------- general structure information --------------------------------------------------------

  type(struct_info)              , public :: meta_struct(nStructures)   ! metadata on the data structures
  type(dim_info)                 , public :: meta_dims(nDimensions)     ! metadata on the dimensions for network topology
  type(dim_info)                 , public :: meta_stateDims(nStateDims) ! metadata on the dimensions for state variables
  type(dim_info)                 , public :: meta_qDims(nQdims)         ! metadata on the dimensions for state variables

  ! ---------- metadata structures ------------------------------------------------------------------

  ! define vectors of metadata
  type(var_info)                 , public :: meta_HRU    (nVarsHRU    ) ! HRU properties
  type(var_info)                 , public :: meta_HRU2SEG(nVarsHRU2SEG) ! HRU-to-segment mapping
  type(var_info)                 , public :: meta_SEG    (nVarsSEG    ) ! stream segment properties
  type(var_info)                 , public :: meta_NTOPO  (nVarsNTOPO  ) ! network topology
  type(var_info)                 , public :: meta_PFAF   (nVarsPFAF   ) ! pfafstetter code
  type(var_info_new)             , public :: meta_rflx   (nVarsRFLX   ) ! reach flux variables
  type(var_info_new)             , public :: meta_basinQ (nVarsBasinQ ) ! reach inflow from basin
  type(var_info_new)             , public :: meta_irf_bas(nVarsIRFbas ) ! basin IRF routing fluxes/states
  type(var_info_new)             , public :: meta_kwt    (nVarsKWT    ) ! KWT routing fluxes/states
  type(var_info_new)             , public :: meta_irf    (nVarsIRF    ) ! IRF routing fluxes/states
  type(var_info_new)             , public :: meta_kwe    (nVarsKWE    ) ! KWE routing fluxes/states

  ! ---------- shared data structures ----------------------------------------------------------------------

  ! river topology and parameter structures
  type(RCHPRP)    , allocatable  , public :: RPARAM(:)            ! Reach Parameters for whole domain
  type(RCHTOPO)   , allocatable  , public :: NETOPO(:)            ! River Network topology for whole domain
  type(RCHPRP)    , allocatable  , public :: RPARAM_trib(:)       ! Reach Parameters for tributaries
  type(RCHTOPO)   , allocatable  , public :: NETOPO_trib(:)       ! River Network topology tributaries
  type(RCHPRP)    , allocatable  , public :: RPARAM_main(:)       ! Reach Parameters for mainstems
  type(RCHTOPO)   , allocatable  , public :: NETOPO_main(:)       ! River Network topology for mainstems

  ! time delay histogram
  REAL(DP)        , allocatable  , public :: FRAC_FUTURE(:)       ! fraction of runoff in future time steps

  ! routing data structures
  TYPE(STRSTA)    , allocatable  , public :: RCHSTA(:,:)          ! Routing state variables (ensembles, space [reaches]) for the entire river network
  TYPE(STRFLX)    , allocatable  , public :: RCHFLX(:,:)          ! Reach fluxes (ensembles, space [reaches]) for entire river network
  TYPE(STRSTA)    , allocatable  , public :: RCHSTA_trib(:,:)     ! Routing state variables (ensembles, space [reaches]) for tributary
  TYPE(STRFLX)    , allocatable  , public :: RCHFLX_trib(:,:)     ! Reach fluxes (ensembles, space [reaches]) for tributaries
  TYPE(STRSTA)    , allocatable  , public :: RCHSTA_main(:,:)     ! Routing state variables (ensembles, space [reaches]) for mainstem
  TYPE(STRFLX)    , allocatable  , public :: RCHFLX_main(:,:)     ! Reach fluxes (ensembles, space [reaches]) for mainstem

  ! lakes data structures
  TYPE(LAKPRP)    , allocatable  , public :: LPARAM(:)            ! Lake parameters
  TYPE(LAKTOPO)   , allocatable  , public :: LKTOPO(:)            ! Lake topology
  TYPE(LKFLX)     , allocatable  , public :: LAKFLX(:,:)          ! Lake fluxes
  TYPE(LKFLX)     , allocatable  , public :: LAKFLX_local(:,:)    ! Lake fluxes for decomposed domains

  ! mapping structures
  type(remap)                    , public :: remap_data           ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)

  ! hru runoff data
  type(runoff)                   , public :: runoff_data          ! HRU runoff data structure for one time step for LSM HRUs and River network HRUs
  real(dp)        , allocatable  , public :: basinRunoff_trib(:)  ! HRU runoff array (m/s) for tributaries
  real(dp)        , allocatable  , public :: basinRunoff_main(:)  ! HRU runoff array (m/s) for mainstem
  real(dp)        , allocatable  , public :: basinEvapo_trib(:)   ! HRU evaporation array (m/s) for tributaries
  real(dp)        , allocatable  , public :: basinEvapo_main(:)   ! HRU evaporation array (m/s) for mainstem
  real(dp)        , allocatable  , public :: basinPrecip_trib(:)  ! HRU precipitation array (m/s) for tributaries
  real(dp)        , allocatable  , public :: basinPrecip_main(:)  ! HRU precipitation array (m/s) for mainstem

  ! seg water management fluxes and target volume
  type(wm)                       , public :: wm_data              ! SEG flux and target vol data structure for one time step for river network
  real(dp)        , allocatable  , public :: flux_wm_trib(:)      ! SEG flux array (m3/s) for tributaries
  real(dp)        , allocatable  , public :: flux_wm_main(:)      ! SEG flux array (m3/s) for mainstem
  real(dp)        , allocatable  , public :: vol_wm_trib(:)       ! SEG target volume (for lakes) (m3) for tributaries
  real(dp)        , allocatable  , public :: vol_wm_main(:)       ! SEG target volume (for lakes) (m3) for mainstem

  ! domain data
  ! MPI
  type(subbasin_mpi),allocatable , public :: domains_mpi(:)       ! mpi domain decomposition data structure
  integer(i4b)                   , public :: nDomain_mpi          ! number of mpi domains
  ! OMP
  type(subbasin_omp), allocatable, public :: river_basin_main(:)   ! openMP domain decomposition for mainstem
  type(subbasin_omp), allocatable, public :: river_basin_trib(:)   ! openMP domain decomposition for tributary

  integer(i4b)    , allocatable  , public :: ixHRU_order(:)       ! global HRU index in the order of proc assignment (size = num of hrus contributing reach in entire network)
  integer(i4b)    , allocatable  , public :: ixRch_order(:)       ! global reach index in the order of proc assignment (size = num of reaches in entire network))
  integer(i4b)    , allocatable  , public :: hru_per_proc(:)      ! number of hrus assigned to each proc (size = num of procs
  integer(i4b)    , allocatable  , public :: rch_per_proc(:)      ! number of reaches assigned to each proc (size = num of procs)

  integer(i4b)    , allocatable  , public :: global_ix_main(:)    ! index array in mainstem array for tributary reach outlet (size = num of tributary outlets)
  integer(i4b)    , allocatable  , public :: global_ix_comm(:)    ! global index array for tributary reach outlet (size = num of tributary outlets)
  integer(i4b)    , allocatable  , public :: local_ix_comm(:)     ! local index array for tributary reach outlet (size = num of tributary outlets per proc)
  integer(i4b)    , allocatable  , public :: tribOutlet_per_proc(:)! number of tributary outlet reaches assigned to each proc (size = num of procs)

end module globalData
