MODULE globalData

  ! General rules
  ! This module include data/variables that can be accessed from any other modules
  ! User do not directly interface the data/variables
  ! Data/values can be altered throughout the runtime
  ! See public_var.f90 for difference

  USE pio

  USE public_var

  USE nrtype

  USE dataTypes, ONLY: struct_info   ! metadata type - data structure
  USE dataTypes, ONLY: dim_info      ! metadata type - variable dimensions
  USE dataTypes, ONLY: var_info      ! metadata type - variable
  USE objTypes,  ONLY: var_info_new  ! metadata type - variable

  USE dataTypes, ONLY: inFileInfo    ! data strture - information of input files
  USE dataTypes, ONLY: gage          ! data structure - gauge metadata

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

  USE dataTypes, ONLY: commLink      ! data structure - arbitrary element-to-element link

  USE dataTypes, ONLY: subbasin_omp  ! data structure - omp domain decomposition
  USE dataTypes, ONLY: subbasin_mpi  ! data structure - mpi domain decomposition

  USE dataTypes, ONLY : cMolecule    ! data structure - computational molecule number

  USE datetime_data, ONLY: datetime  ! datetime data class

  USE base_route, ONLY: routeContainer ! a container of instantiated routing methods

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
  USE var_lookup, ONLY: nVarsHFLX    ! number of variables in data structure (HRU flux/state)
  USE var_lookup, ONLY: nVarsBasinQ  ! number of variables in data structure (restart vars for
  USE var_lookup, ONLY: nVarsIRFbas  ! number of variables in data structure (restart vars for overland unit-hydrograph routing)
  USE var_lookup, ONLY: nVarsIRF     ! number of variables in data structure (restart vars for unit-hydrograph routing)
  USE var_lookup, ONLY: nVarsKWT     ! number of variables in data structure (restart vars for lagrangian kinematic wave)
  USE var_lookup, ONLY: nVarsKW      ! number of variables in data structure (restart vars for kinematic wave routing)
  USE var_lookup, ONLY: nVarsMC      ! number of variables in data structure (restart vars for muskingum-cunge routing)
  USE var_lookup, ONLY: nVarsDW      ! number of variables in data structure (restart vars for diffusive wave routing)

  implicit none

  save

  ! ---------- mizuRoute version -------------------------------------------------------------------

  character(len=strLen)          , public :: version='undefined'        ! mizuRoute version
  character(len=strLen)          , public :: gitBranch='undefined'      ! github branch
  character(len=strLen)          , public :: gitHash='undefined'        ! github commit hash

  ! ---------- Catchment and reach IDs  -------------------------------------------------------------------------

  integer(i4b),                    public :: nHRU                 ! number of HRUs in the whole river network
  integer(i4b),                    public :: nContribHRU          ! number of HRUs that are connected to any reaches
  integer(i4b),                    public :: nRch                 ! number of reaches in the whole river network
  integer(i4b),                    public :: nRch_mainstem        ! number of reaches in mainstems
  integer(i4b),                    public :: nHRU_mainstem        ! number of HRUs in mainstems
  integer(i4b),                    public :: nRch_trib            ! number of reaches in tributaries
  integer(i4b),                    public :: nHRU_trib            ! number of HRUs in tributaries
  integer(i4b),     allocatable,   public :: basinID(:)           ! HRU id in the whole river network
  integer(i4b),     allocatable,   public :: reachID(:)           ! reach id in the whole river network

  ! ---------- routing methods  -------------------------------------------------------------------------
  type(routeContainer), allocatable , public :: rch_routes(:)           ! a collection of routing method objects
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
  real(dp),                        public :: sec2tunit            ! time conversion- seconds per time unit used in simulation
  real(dp),                        public :: timeVar(1:2)         ! time variables [sec] at endpoints of a simulation time step
  real(dp),                        public :: TSEC(1:2)            ! seconds at endpoints of a simulation time step since simulation started
  type(datetime),                  public :: simDatetime(0:2)     ! previous, current and next simulation time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: begDatetime          ! simulation start date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: endDatetime          ! simulation end date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: restDatetime         ! desired restart date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: dropDatetime         ! restart dropoff date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: roBegDatetime        ! forcing data start date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime),                  public :: wmBegDatetime        ! water management data start date/time (yyyy:mm:dd:hh:mm:ss)

  ! ---------- input file information -------------------------------------------------------------------

  type(infileinfo), allocatable,   public :: inFileInfo_ro(:)    ! input runoff/evapo/precipi file information
  type(infileinfo), allocatable,   public :: inFileInfo_wm(:)    ! input water management (abstaction/injection) file information

  ! ---------- Misc. data -------------------------------------------------------------------------
  character(len=strLen),           public :: runMode='standalone'        ! run options: standalone or cesm-coupling
  character(len=12),               public :: hfile_dayStamp='period-start' ! day stamp in history file for daily file. period-start or period-end
  integer(i4b),                    public :: ixPrint(1:2)=integerMissing ! index of desired reach to be on-screen print
  integer(i4b),                    public :: nEns=1                      ! number of ensemble
  integer(i4b),                    public :: maxtdh                      ! maximum unit-hydrograph future time steps
  type(cMolecule),                 public :: nMolecule                   ! number of computational molecule (used for KW, MC, DW)
  character(FileStrLen),           public :: hfileout=charMissing        ! history output file name
  character(FileStrLen),           public :: hfileout_gage=charMissing   ! gage-only history output file name
  character(FileStrLen),           public :: rfileout=charMissing        ! restart output file name
  logical(lgt),                    public :: initHvars=.false.           ! status of history variable data initialization
  logical(lgt),                    public :: isColdStart=.true.          ! initial river state - cold start (T) or from restart file (F)

  ! ---------- MPI/OMP/PIO variables ----------------------------------------------------------------

  integer(i4b),                    public :: mpicom_route                       ! communicator for this program
  integer(i4b),                    public :: pid                                ! process id
  integer(i4b),                    public :: nNodes                             ! number of MPI processors
  integer(i4b),                    public :: nThreads                           ! number of threads
  logical(lgt),                    public :: masterproc                         ! root logical. root processor => true, other => false
  logical(lgt),                    public :: multiProcs                         ! MPI multi-processors logical. => number of processors>1 true, other => false
  integer(i4b),                    public :: pio_numiotasks    = -99            ! Number of iotasks (ntasks/stride) - see PIO documentation for more information
  integer(i4b),                    public :: pio_rearranger    = 2              ! 0=>PIO_rearr_none 1=> PIO_rearr_box 2=> PIO_rearr_subset
  integer(i4b),                    public :: pio_root          = 1              ! PIO root
  integer(i4b),                    public :: pio_stride        = 1              ! PIO stride - see PIO documentation for more information
  type(iosystem_desc_t), pointer,  public :: pioSystem                          ! PIO I/O system data
  ! pio decomposition used for history file variables
  type(io_desc_t),                 public :: ioDesc_hru_float                   ! [hru] (float)
  type(io_desc_t),                 public :: ioDesc_rch_float                   ! [reach] (float)
  type(io_desc_t),                 public :: ioDesc_gauge_float                 ! [gauge points] (float)
  ! decomposition used for restart variables
  type(io_desc_t),                 public :: ioDesc_rch_int                     ! [reach x ensemble] (integer)
  type(io_desc_t),                 public :: ioDesc_rch_double                  ! [reach x ensemble] (double precision)
  type(io_desc_t),                 public :: ioDesc_hist_rch_double             ! [reach] (double precision)
  type(io_desc_t),                 public :: ioDesc_hru_double                  ! [hru] (double precision)
  type(io_desc_t),                 public :: ioDesc_wave_int                    ! [reach x max. number of waves x ensemble] (integer)
  type(io_desc_t),                 public :: ioDesc_wave_double                 ! [reach x max. number of waves x ensemble] (double precision)
  type(io_desc_t),                 public :: ioDesc_mesh_kw_double              ! [reach x Euler kinematic wave computational meshes x ensemble] (double precision)
  type(io_desc_t),                 public :: ioDesc_mesh_mc_double              ! [reach x Muskingum-Cunge computational meshes x ensemble] (double precision)
  type(io_desc_t),                 public :: ioDesc_mesh_dw_double              ! [reach x Duffusive wave computational meshes x ensemble] (double precision)
  type(io_desc_t),                 public :: ioDesc_irf_double                  ! [reach x IRF future timer-steps x ensemble] (double precision)
  type(io_desc_t),                 public :: ioDesc_vol_double                  ! [reach x time-bounds x ensemble] (double precision)
  type(io_desc_t),                 public :: ioDesc_irf_bas_double              ! [reach x basin UH future time-steps x ensemble] (double precision)
  integer(i4b),   allocatable,     public :: index_write_gage(:)                ! reach indices to gauge points w.r.t. distributed domains

  ! ---------- conversion factors -------------------------------------------------------------------

  real(dp),                        public :: time_conv                  ! time conversion factor -- used to convert to mm/s
  real(dp),                        public :: length_conv                ! length conversion factor -- used to convert to mm/s

  ! ---------- routing parameter names -------------------------------------------------------------------
  ! spatial constant ....
  real(dp),                        public :: fshape                     ! shape parameter in time delay histogram (=gamma distribution) [-]
  real(dp),                        public :: tscale                     ! scaling factor for the time delay histogram [sec]
  real(dp),                        public :: velo                       ! velocity [m/s] for Saint-Venant equation
  real(dp),                        public :: diff                       ! diffusivity [m2/s] for Saint-Venant equation
  real(dp),                        public :: mann_n                     ! manning's roughness coefficient [-]
  real(dp),                        public :: wscale                     ! scaling factor for river width [-]
  real(dp),                        public :: channelDepth=10.           ! channel bankfull depth [m]
  real(dp),                        public :: floodplainSlope=0.001      ! floodplain down slope [-]

  ! ---------- general structure information --------------------------------------------------------

  type(struct_info),               public :: meta_struct(nStructures)   ! metadata on the data structures
  type(dim_info),                  public :: meta_dims(nDimensions)     ! metadata on the dimensions for network topology
  type(dim_info),                  public :: meta_stateDims(nStateDims) ! metadata on the dimensions for state variables
  type(dim_info),                  public :: meta_qDims(nQdims)         ! metadata on the dimensions for flux variables
  type(dim_info),                  public :: meta_qDims_gage(nQdims)    ! metadata on the dimensions for flux variables

  ! ---------- metadata structures ------------------------------------------------------------------

  ! define vectors of metadata
  type(var_info),                  public :: meta_HRU    (nVarsHRU    ) ! HRU properties
  type(var_info),                  public :: meta_HRU2SEG(nVarsHRU2SEG) ! HRU-to-segment mapping
  type(var_info),                  public :: meta_SEG    (nVarsSEG    ) ! stream segment properties
  type(var_info),                  public :: meta_NTOPO  (nVarsNTOPO  ) ! network topology
  type(var_info),                  public :: meta_PFAF   (nVarsPFAF   ) ! pfafstetter code
  type(var_info_new),              public :: meta_rflx   (nVarsRFLX   ) ! reach flux variables
  type(var_info_new),              public :: meta_hflx   (nVarsHFLX   ) ! hru flux variables
  type(var_info_new),              public :: meta_basinQ (nVarsBasinQ ) ! reach inflow from basin
  type(var_info_new),              public :: meta_irf_bas(nVarsIRFbas ) ! basin IRF routing fluxes/states
  type(var_info_new),              public :: meta_irf    (nVarsIRF    ) ! IRF routing fluxes/states
  type(var_info_new),              public :: meta_kwt    (nVarsKWT    ) ! KWT routing fluxes/states
  type(var_info_new),              public :: meta_kw     (nVarsKW     ) ! KW routing fluxes/states
  type(var_info_new),              public :: meta_mc     (nVarsMC     ) ! MC routing restart fluxes/states
  type(var_info_new),              public :: meta_dw     (nVarsDW     ) ! DW routing restart fluxes/states

  ! ---------- shared data structures ----------------------------------------------------------------------
  type(gage),                      public :: gage_data              ! gauge metadata

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
  type(STRSTA),       allocatable, public :: RCHSTA(:,:)            ! restart variables (ensembles, reaches) for the entire domain
  type(STRFLX),       allocatable, public :: RCHFLX(:,:)            ! Reach fluxes (ensembles, reaches) for entire domain
  type(STRSTA),       allocatable, public :: RCHSTA_trib(:,:)       ! restart variables (ensembles, reaches) for decomposed domains
  type(STRFLX),       allocatable, public :: RCHFLX_trib(:,:)       ! Reach fluxes (ensembles, reaches) for decomposed domains

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
  ! -- MPI
  type(subbasin_mpi), allocatable, public :: domains_mpi(:)         ! mpi domain decomposition data structure
  integer(i4b),                    public :: nDomain_mpi            ! number of mpi domains

  ! -- OMP
  type(subbasin_omp), allocatable, public :: river_basin_main(:)    ! openMP domain decomposition for mainstem
  type(subbasin_omp), allocatable, public :: river_basin_trib(:)    ! openMP domain decomposition for tributary

  ! -- reach and hru ordered index
  integer(i4b),       allocatable, public :: ixHRU_order(:)         ! global HRU index in the order of proc assignment (size = num of hrus contributing reach in entire network)
  integer(i4b),       allocatable, public :: ixRch_order(:)         ! global reach index in the order of proc assignment (size = num of reaches in entire network))
  integer(i4b),       allocatable, public :: hru_per_proc(:)        ! number of hrus assigned to each proc (size = num of procs
  integer(i4b),       allocatable, public :: rch_per_proc(:)        ! number of reaches assigned to each proc (size = num of procs)
  integer(i4b),       allocatable, public :: nTribOutlet            ! number of tributary reaches flowing into mainstem
  integer(i4b),       allocatable, public :: global_ix_main(:)      ! index array in mainstem array for tributary reach outlet (size = num of tributary outlets)
  integer(i4b),       allocatable, public :: global_ix_comm(:)      ! global index array for tributary reach outlet (size = num of tributary outlets)
  integer(i4b),       allocatable, public :: local_ix_comm(:)       ! local index array for tributary reach outlet (size = num of tributary outlets per proc)
  integer(i4b),       allocatable, public :: tribOutlet_per_proc(:) ! number of tributary outlet reaches assigned to each proc (size = num of procs)

  ! -- reach to reach connection
  type(commLink),     allocatable, public :: commRch(:)             ! reach-reach connections for reach flux transfer

END MODULE globalData
