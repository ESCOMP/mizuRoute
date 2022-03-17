MODULE globalData

  USE pio

  USE public_var

  USE nrtype

  USE dataTypes, ONLY: struct_info   ! metadata type - data structure
  USE dataTypes, ONLY: dim_info      ! metadata type - variable dimensions
  USE dataTypes, ONLY: var_info      ! metadata type - variable
  USE objTypes,  ONLY: var_info_new  ! metadata type - variable

  USE dataTypes, ONLY: infileinfo    ! data strture - information of input files

  USE dataTypes, ONLY: RCHPRP        ! data structure - Reach parameters (properties)
  USE dataTypes, ONLY: RCHTOPO       ! data structure - Network topology

  USE dataTypes, ONLY: STRFLX        ! data structure - fluxes in each reach
  USE dataTypes, ONLY: STRSTA        ! data structure - states in each reach

  USE dataTypes, ONLY: LAKPRP        ! data structure - lake properties
  USE dataTypes, ONLY: LAKTOPO       ! data structure - lake topology
  USE dataTypes, ONLY: LKFLX         ! data structure - lake fluxes

  USE dataTypes, ONLY: remap         ! data structure - remapping data type
  USE dataTypes, ONLY: runoff        ! data structure - runoff data type
  USE dataTypes, ONLY: wm            ! data structure - water management (flux to/from segment, target volume) data type

  USE dataTypes, ONLY: subbasin_omp  ! data structure - omp domain decomposition
  USE dataTypes, ONLY: subbasin_mpi  ! data structure - mpi domain decomposition

  USE dataTypes, ONLY : cMolecule    ! data structure - computational molecule number

  USE datetime_data, ONLY: datetime  ! datetime data class

  USE var_lookup, ONLY: nStructures  ! number of variables in data structure (struct_info)
  USE var_lookup, ONLY: nDimensions  ! number of variables in dimensions related to network topology
  USE var_lookup, ONLY: nStateDims   ! number of variables in dimensions related to restart variables
  USE var_lookup, ONLY: nQdims       ! number of variables in dimensions related to fluxes/states variables
  USE var_lookup, ONLY: nVarsHRU     ! number of variables in data structure (catchment propoerties)
  USE var_lookup, ONLY: nVarsHRU2SEG ! number of variables in data structure (river-catchment topology)
  USE var_lookup, ONLY: nVarsSEG     ! number of variables in data structure (river reach propeties)
  USE var_lookup, ONLY: nVarsNTOPO   ! number of variables in data structure (river network topology)
  USE var_lookup, ONLY: nVarsPFAF    ! number of variables in data structure (pfaffstetter related variable)
  USE var_lookup, ONLY: nVarsRFLX    ! number of variables in data structure (river flux/state)
  USE var_lookup, ONLY: nVarsBasinQ  ! number of variables in data structure (restart vars for
  USE var_lookup, ONLY: nVarsIRFbas  ! number of variables in data structure (restart vars for overland unit-hydrograph routing)
  USE var_lookup, ONLY: nVarsIRF     ! number of variables in data structure (restart vars for unit-hydrograph routing)
  USE var_lookup, ONLY: nVarsKWT     ! number of variables in data structure (restart vars for lagrangian kinematic wave)
  USE var_lookup, ONLY: nVarsKW      ! number of variables in data structure (restart vars for kinematic wave routing)
  USE var_lookup, ONLY: nVarsMC      ! number of variables in data structure (restart vars for muskingum-cunge routing)
  USE var_lookup, ONLY: nVarsDW      ! number of variables in data structure (restart vars for diffusive wave routing)

  implicit none

  save

  ! ---------- Catchment and reach IDs  -------------------------------------------------------------------------

  integer(i4b),                    public :: nHRU                 ! number of HRUs in the whole river network
  integer(i4b),                    public :: nContribHRU          ! number of HRUs that are connected to any reaches
  integer(i4b),                    public :: nRch                 ! number of reaches in the whole river network
  integer(i4b),                    public :: nRch_mainstem        ! number of reaches in mainstems
  integer(i4b),                    public :: nHRU_mainstem        ! number of HRUs in mainstems
  integer(i4b),     allocatable,   public :: basinID(:)           ! HRU id in the whole river network
  integer(i4b),     allocatable,   public :: reachID(:)           ! reach id in the whole river network

  ! ---------- routing methods  -------------------------------------------------------------------------
  integer(i4b)                   , public :: nRoutes                    ! number of active routing methods
  integer(i4b)    , allocatable  , public :: routeMethods(:)            ! active routing method id
  logical(lgt)                   , public :: onRoute(0:nRouteMethods-1) ! logical to indicate active routing method(s)
  integer(i4b)                   , public :: idxSUM                     ! index of SUM method
  integer(i4b)                   , public :: idxIRF                     ! index of IRF method
  integer(i4b)                   , public :: idxKWT                     ! index of KWT method
  integer(i4b)                   , public :: idxKW                      ! index of KW method
  integer(i4b)                   , public :: idxMC                      ! index of MC method
  integer(i4b)                   , public :: idxDW                      ! index of DW method

  ! ---------- Date/Time data  -------------------------------------------------------------------------

  integer(i4b),                    public :: iTime                ! time index at simulation time step
  real(dp),         allocatable,   public :: timeVar(:)           ! time variables (unit given by runoff data)
  real(dp),                        public :: TSEC(0:1)            ! begning and end of time step since simulation started (sec)
  type(datetime),                  public :: simDatetime(0:1)     ! previous and current simulation time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: begDatetime          ! simulation start date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: endDatetime          ! simulation end date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: restDatetime         ! desired restart date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: dropDatetime         ! restart dropoff date/time (yyyy:mm:dd:hh:mm:ss)

  ! ---------- input file information -------------------------------------------------------------------

  type(infileinfo), allocatable,   public :: infileinfo_data(:)    ! input runoff file information
  type(infileinfo), allocatable,   public :: infileinfo_data_wm(:) ! input water management (abstaction/injection) file information

  ! ---------- Misc. data -------------------------------------------------------------------------
  character(len=strLen),           public :: runMode='standalone'        ! run options: standalone or cesm-coupling
  logical(lgt),                    public :: isHistFileOpen=.false.      ! flag to indicate history output netcdf is open
  integer(i4b),                    public :: ixPrint(1:2)=integerMissing ! index of desired reach to be on-screen print
  integer(i4b),                    public :: nEns=1                      ! number of ensemble
  type(cMolecule),                 public :: nMolecule                   ! number of computational molecule (used for KW, MC, DW)

  ! ---------- MPI/OMP/PIO variables ----------------------------------------------------------------

  integer(i4b),                    public :: mpicom_route                       ! communicator for this program
  integer(i4b),                    public :: pid                                ! process id
  integer(i4b),                    public :: nNodes                             ! number of MPI processors
  integer(i4b),                    public :: nThreads                           ! number of threads
  logical(lgt),                    public :: masterproc                         ! root logical. root processor => true, other => false
  logical(lgt),                    public :: multiProcs                         ! MPI multi-processors logical. => number of processors>1 true, other => false
  character(len=strLen),           public :: pio_netcdf_format = "64bit_offset" !
  character(len=strLen),           public :: pio_typename      = "pnetcdf"      !
  integer(i4b),                    public :: pio_numiotasks    = -99            !
  integer(i4b),                    public :: pio_rearranger    = 2              ! 0=>PIO_rearr_none 1=> PIO_rearr_box 2=> PIO_rearr_subset
  integer(i4b),                    public :: pio_root          = 1              !
  integer(i4b),                    public :: pio_stride        = 1              !
  type(iosystem_desc_t),           public :: pioSystem                          ! PIO I/O system data

  ! ---------- conversion factors -------------------------------------------------------------------

  real(dp),                        public :: time_conv                  ! time conversion factor -- used to convert to mm/s
  real(dp),                        public :: length_conv                ! length conversion factor -- used to convert to mm/s

  ! ---------- routing parameter names -------------------------------------------------------------------
  ! spatial constant ....
  real(dp),                        public :: fshape                     ! shape parameter in time delay histogram (=gamma distribution) [-]
  real(dp),                        public :: tscale                     ! scaling factor for the time delay histogram [sec]
  real(dp),                        public :: velo                       ! velocity [m/s] for Saint-Venant equation   added by NM
  real(dp),                        public :: diff                       ! diffusivity [m2/s] for Saint-Venant equation   added by NM
  real(dp),                        public :: mann_n                     ! manning's roughness coefficient [unitless]  added by NM
  real(dp),                        public :: wscale                     ! scaling factor for river width [-] added by NM

  ! ---------- general structure information --------------------------------------------------------

  type(struct_info),               public :: meta_struct(nStructures)   ! metadata on the data structures
  type(dim_info),                  public :: meta_dims(nDimensions)     ! metadata on the dimensions for network topology
  type(dim_info),                  public :: meta_stateDims(nStateDims) ! metadata on the dimensions for state variables
  type(dim_info),                  public :: meta_qDims(nQdims)         ! metadata on the dimensions for state variables

  ! ---------- metadata structures ------------------------------------------------------------------

  ! define vectors of metadata
  type(var_info),                  public :: meta_HRU    (nVarsHRU    ) ! HRU properties
  type(var_info),                  public :: meta_HRU2SEG(nVarsHRU2SEG) ! HRU-to-segment mapping
  type(var_info),                  public :: meta_SEG    (nVarsSEG    ) ! stream segment properties
  type(var_info),                  public :: meta_NTOPO  (nVarsNTOPO  ) ! network topology
  type(var_info),                  public :: meta_PFAF   (nVarsPFAF   ) ! pfafstetter code
  type(var_info_new),              public :: meta_rflx   (nVarsRFLX   ) ! reach flux variables
  type(var_info_new),              public :: meta_basinQ (nVarsBasinQ ) ! reach inflow from basin
  type(var_info_new),              public :: meta_irf_bas(nVarsIRFbas ) ! basin IRF routing fluxes/states
  type(var_info_new),              public :: meta_irf    (nVarsIRF    ) ! IRF routing fluxes/states
  type(var_info_new),              public :: meta_kwt    (nVarsKWT    ) ! KWT routing fluxes/states
  type(var_info_new),              public :: meta_kw     (nVarsKW     ) ! KW routing fluxes/states
  type(var_info_new),              public :: meta_mc     (nVarsMC     ) ! MC routing restart fluxes/states
  type(var_info_new),              public :: meta_dw     (nVarsDW     ) ! DW routing restart fluxes/states

  ! ---------- shared data structures ----------------------------------------------------------------------

  ! river topology and parameter structures
  type(RCHPRP),       allocatable, public :: RPARAM(:)              ! Reach Parameters for whole domain
  type(RCHTOPO),      allocatable, public :: NETOPO(:)              ! River Network topology for whole domain
  type(RCHPRP),       allocatable, public :: RPARAM_trib(:)         ! Reach Parameters for tributaries
  type(RCHTOPO),      allocatable, public :: NETOPO_trib(:)         ! River Network topology tributaries
  type(RCHPRP),       allocatable, public :: RPARAM_main(:)         ! Reach Parameters for mainstems
  type(RCHTOPO),      allocatable, public :: NETOPO_main(:)         ! River Network topology for mainstems

  ! time delay histogram
  real(dp),           allocatable, public :: FRAC_FUTURE(:)         ! fraction of runoff in future time steps

  ! routing data structures
  type(STRSTA),       allocatable, public :: RCHSTA(:,:)            ! restart variables (ensembles, space [reaches]) for the entire river network
  type(STRFLX),       allocatable, public :: RCHFLX(:,:)            ! Reach fluxes (ensembles, space [reaches]) for entire river network
  type(STRSTA),       allocatable, public :: RCHSTA_trib(:,:)       ! restart variables (ensembles, space [reaches]) for tributary
  type(STRFLX),       allocatable, public :: RCHFLX_trib(:,:)       ! Reach fluxes (ensembles, space [reaches]) for tributaries
  type(STRSTA),       allocatable, public :: RCHSTA_main(:,:)       ! restart variables (ensembles, space [reaches]) for mainstem
  type(STRFLX),       allocatable, public :: RCHFLX_main(:,:)       ! Reach fluxes (ensembles, space [reaches]) for mainstem

  ! lakes data structures
  type(LAKPRP),       allocatable, public :: LPARAM(:)              ! Lake parameters
  type(LAKTOPO),      allocatable, public :: LKTOPO(:)              ! Lake topology
  type(LKFLX),        allocatable, public :: LAKFLX(:,:)            ! Lake fluxes
  type(LKFLX),        allocatable, public :: LAKFLX_local(:,:)      ! Lake fluxes for decomposed domains

  ! mapping structures
  type(remap),                     public :: remap_data             ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)

  ! hru runoff data
  type(runoff),                    public :: runoff_data            ! HRU runoff data structure for one time step for LSM HRUs and River network HRUs
  real(dp),           allocatable, public :: basinRunoff_trib(:)    ! HRU runoff array (m/s) for tributaries
  real(dp),           allocatable, public :: basinRunoff_main(:)    ! HRU runoff array (m/s) for mainstem
  real(dp),           allocatable, public :: basinEvapo_trib(:)     ! HRU evaporation array (m/s) for tributaries
  real(dp),           allocatable, public :: basinEvapo_main(:)     ! HRU evaporation array (m/s) for mainstem
  real(dp),           allocatable, public :: basinPrecip_trib(:)    ! HRU precipitation array (m/s) for tributaries
  real(dp),           allocatable, public :: basinPrecip_main(:)    ! HRU precipitation array (m/s) for mainstem

  ! seg water management fluxes and target volume
  type(wm),                        public :: wm_data                ! SEG flux and target vol data structure for one time step for river network
  real(dp),           allocatable, public :: flux_wm_trib(:)        ! SEG flux array (m3/s) for tributaries
  real(dp),           allocatable, public :: flux_wm_main(:)        ! SEG flux array (m3/s) for mainstem
  real(dp),           allocatable, public :: vol_wm_trib(:)         ! SEG target volume (for lakes) (m3) for tributaries
  real(dp),           allocatable, public :: vol_wm_main(:)         ! SEG target volume (for lakes) (m3) for mainstem

  ! domain data
  ! MPI
  type(subbasin_mpi), allocatable, public :: domains_mpi(:)         ! mpi domain decomposition data structure
  integer(i4b),                    public :: nDomain_mpi            ! number of mpi domains
  ! OMP
  type(subbasin_omp), allocatable, public :: river_basin_main(:)    ! openMP domain decomposition for mainstem
  type(subbasin_omp), allocatable, public :: river_basin_trib(:)    ! openMP domain decomposition for tributary

  integer(i4b),       allocatable, public :: ixHRU_order(:)         ! global HRU index in the order of proc assignment (size = num of hrus contributing reach in entire network)
  integer(i4b),       allocatable, public :: ixRch_order(:)         ! global reach index in the order of proc assignment (size = num of reaches in entire network))
  integer(i4b),       allocatable, public :: hru_per_proc(:)        ! number of hrus assigned to each proc (size = num of procs
  integer(i4b),       allocatable, public :: rch_per_proc(:)        ! number of reaches assigned to each proc (size = num of procs)

  integer(i4b),       allocatable, public :: global_ix_main(:)      ! index array in mainstem array for tributary reach outlet (size = num of tributary outlets)
  integer(i4b),       allocatable, public :: global_ix_comm(:)      ! global index array for tributary reach outlet (size = num of tributary outlets)
  integer(i4b),       allocatable, public :: local_ix_comm(:)       ! local index array for tributary reach outlet (size = num of tributary outlets per proc)
  integer(i4b),       allocatable, public :: tribOutlet_per_proc(:) ! number of tributary outlet reaches assigned to each proc (size = num of procs)

END MODULE globalData
