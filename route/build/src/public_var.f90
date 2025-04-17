MODULE public_var

  ! General rules
  ! This module include variables that can be accessed from any other modules
  ! User can set values of the variables through control file (read_control.f90).
  ! Values should not be altered during the runtime.
  ! See globalData.f90 for difference

  USE nrtype

  implicit none

  save

  ! physical constants
  real(dp),    parameter,public   :: pi=3.14159265359_dp    ! pi
  real(dp),    parameter,public   :: Cw=4190.0_dp           ! heat capacity of water [J/kg/K]
  real(dp),    parameter,public   :: RoW=0.99975_dp         ! density of water [kg/m³] at 10 C-degree

  ! some common constant variables (not likely to change value)
  real(dp),    parameter,public   :: secprmin=60._dp        ! number of seconds in a minute
  real(dp),    parameter,public   :: secprhour=3600._dp     ! number of seconds in an hour
  real(dp),    parameter,public   :: secprday=86400._dp     ! number of seconds in a day
  integer(i4b),parameter,public   :: months_per_yr=12       ! number of months in a year
  integer(i4b),parameter,public   :: days_per_yr=365        ! number of days in a year
  real(dp),    parameter,public   :: hr_per_day = 24.0_dp   ! hours per days
  real(dp),    parameter,public   :: min_per_hour = 60.0_dp ! minutes per hour
  real(dp),    parameter,public   :: maxTimeDiff=1/secprday ! time difference tolerance for input checks
  real(dp),    parameter,public   :: verySmall=tiny(1.0_dp) ! a very small number
  real(dp),    parameter,public   :: min_slope=1.e-6_dp     ! minimum slope
  real(dp),    parameter,public   :: runoffMin=1.e-15_dp    ! minimum runoff from each basin
  real(dp),    parameter,public   :: negRunoffTol=-1.e-3_dp ! nagative runoff tolerance
  real(dp),    parameter,public   :: lakeWBtol=1.e-3_dp     ! lake water balance tolerance

  ! routing related constants
  integer(i4b),parameter,public   :: MAXQPAR=20             ! maximum number of particles

  ! openMP domain decompostion parameters
  integer(i4b),parameter,public   :: maxSegs=100            ! maximum reach numbers within tributaries
  integer(i4b),parameter,public   :: maxLevel=20            ! maximum mainstem levels used for OMP domain decomposition

  ! constants for general use
  real(dp),    parameter,public   :: MaxPosVal=1.e36_dp     ! maximum value for positive value
  real(dp),    parameter,public   :: MinPosVal=1.e-10_dp    ! minimum value for positive value
  integer(i4b),parameter,public   :: integerMissing=-9999   ! missing value for integers
  real(sp),    parameter,public   :: floatMissing=-9999._sp ! missing value for real32 numbers
  real(dp),    parameter,public   :: realMissing=-9999._dp  ! missing value for real64 numbers
  character(5),parameter,public   :: charMissing='empty'    ! missing value for character

  ! mpi related parameters
  integer(i4b),parameter,public   :: root=0                 ! root node id

  ! I/O related parameters
  integer(i4b),          public   :: iulog=6                            ! logical unit identifier
  character(len=strLen), public   :: rpntfil='rpointer.rof'             ! file name for local restart pointer file
  character(len=strLen), public   :: pio_netcdf_format="64bit_offset"   ! netCDF format - use '64bit_offset' for PIO use or 'netCDF-4'
  character(len=strLen), public   :: pio_typename="pnetcdf"             ! netcdf, pnetcdf, netcdf4c, or netcdf4p

  ! ---------- named variables ----------------------------------------------------------------------

  ! true/false
  integer(i4b), parameter, public :: true=1                     ! true
  integer(i4b), parameter, public :: false=0                    ! false

  ! variable types
  integer(i4b), parameter, public :: varType_integer   = 1001   ! named variable for an integer
  integer(i4b), parameter, public :: varType_double    = 1002   ! named variable for a double precision
  integer(i4b), parameter, public :: varType_character = 1003   ! named variable for a double precision

  ! compute versus read from file
  integer(i4b), parameter,public  :: compute=1              ! compute given variable
  integer(i4b), parameter,public  :: doNotCompute=0         ! do not compute given variable
  integer(i4b), parameter,public  :: readFromFile=0         ! read given variable from a file

  ! routing methods
  integer(i4b), parameter,public  :: nRouteMethods=6         ! number of routing methods available
  integer(i4b), parameter,public  :: accumRunoff=0           ! runoff accumulation over all the upstream reaches
  integer(i4b), parameter,public  :: impulseResponseFunc=1   ! impulse response function
  integer(i4b), parameter,public  :: kinematicWaveTracking=2 ! Lagrangian kinematic wave
  integer(i4b), parameter,public  :: kinematicWave=3         ! kinematic wave
  integer(i4b), parameter,public  :: muskingumCunge=4        ! muskingum-cunge
  integer(i4b), parameter,public  :: diffusiveWave=5         ! diffusiveWave

  ! ---------- variables in the control file --------------------------------------------------------

  ! Control file variables
  ! DIRECTORIES
  character(len=strLen),public    :: ancil_dir            = charMissing     ! directory containing ancillary data (network, mapping, namelist)
  character(len=strLen),public    :: input_dir            = charMissing     ! directory containing input forcing data
  character(len=strLen),public    :: output_dir           = charMissing     ! directory for routed flow output (netCDF)
  character(len=strLen),public    :: restart_dir          = charMissing     ! directory for restart output (netCDF)
  ! RUN CONTROL
  character(len=strLen),public    :: case_name            = ''              ! name of simulation
  logical(lgt),public             :: continue_run         = .false.         ! T-> append output in existing history files. F-> write output in new history file
  character(len=strLen),public    :: simStart             = ''              ! date string defining the start of the simulation
  character(len=strLen),public    :: simEnd               = ''              ! date string defining the end of the simulation
  character(len=strLen),public    :: newFileFrequency     = 'yearly'        ! frequency for new output files (daily, monthly, yearly, single)
  character(len=strLen),public    :: outputFrequency      = '1'             ! output frequency (integer for multiple of simulation time step or daily, monthly or yearly)
  character(len=strLen),public    :: outputNameOption     = 'specific'      ! option for routing method dependent output names (e.g., routedRunoff) - generic or specific (default)
  integer(i4b)         ,public    :: nOutFreq             = integerMissing  ! integer output frequency
  character(len=10)    ,public    :: routOpt              = '0'             ! routing scheme options  0: accum runoff, 1:IRF, 2:KWT, 3:KW, 4:MC, 5:DW
  integer(i4b)         ,public    :: doesBasinRoute       = 1               ! basin routing options   0-> no, 1->IRF, otherwise error
  logical(lgt)         ,public    :: floodplain           = .false.         ! logical if flood water is computed or not (floodplain is added)
  integer(i4b)         ,public    :: hw_drain_point       = 2               ! how to add inst. runoff in reach for headwater HRUs. 1->top of reach, 2->bottom of reach (default)
  logical(lgt)         ,public    :: is_lake_sim          = .false.         ! logical if lakes are activated in simulation
  logical(lgt)         ,public    :: lakeRegulate         = .true.          ! logical: F -> turn all the lakes into natural (lakeType=1) regardless of lakeModelType defined individually
  integer(i4b)         ,public    :: LakeInputOption      = 0               ! fluxes for lake simulation; 0->evaporation+precipitation, 1->runoff, 3->evaporation+precipitation+runoff
  logical(lgt)         ,public    :: is_flux_wm           = .false.         ! logical if flow is added or removed from a reach
  logical(lgt)         ,public    :: is_vol_wm            = .false.         ! logical if target volume is considered for a lake
  logical(lgt)         ,public    :: is_vol_wm_jumpstart  = .false.         ! logical if true the volume is reset to target volume for the first time step of modeling
  real(dp)             ,public    :: scale_factor_runoff  = realMissing     ! float scale to scale the runoff
  real(dp)             ,public    :: offset_value_runoff  = realMissing     ! float offset for runoff
  real(dp)             ,public    :: scale_factor_Ep      = realMissing     ! float scale to scale the potential evaporation
  real(dp)             ,public    :: offset_value_Ep      = realMissing     ! float offset for evaporation
  logical(lgt)         ,public    :: is_Ep_upward_negative= .false.         ! set to true when convention of upward evaporation is negative
  real(dp)             ,public    :: scale_factor_prec    = realMissing     ! float scale to scale the precipitation
  real(dp)             ,public    :: offset_value_prec    = realMissing     ! float offset for precipitation
  real(dp)             ,public    :: min_length_route     = 0.0_dp          ! float; minimum reach length for routing to be performed. pass-through is performed for length less than this threshold
  logical(lgt)         ,public    :: compWB               = .false.         ! logical if entire domain water balance is computed
  real(dp)             ,public    :: dt                   = realMissing     ! simulation time step (seconds)
  ! RIVER NETWORK TOPOLOGY
  character(len=strLen),public    :: fname_ntopOld        = ''              ! old filename containing stream network topology information
  logical(lgt)         ,public    :: ntopAugmentMode      = .false.         ! option for river network augmentation mode. terminate the program after writing augmented ntopo.
  character(len=strLen),public    :: fname_ntopNew        = ''              ! new filename containing stream network topology information
  character(len=strLen),public    :: dname_sseg           = ''              ! dimension name of segment in river network data
  character(len=strLen),public    :: dname_nhru           = ''              ! dimension name of hru in river network data
  ! RUNOFF, EVAPORATION AND PRECIPITATION FILE
  character(len=strLen),public    :: fname_qsim           = charMissing     ! runoff netCDF name
  character(len=strLen),public    :: vname_qsim           = charMissing     ! variable name for runoff
  character(len=strLen),public    :: vname_evapo          = charMissing     ! variable name for actual evapoartion
  character(len=strLen),public    :: vname_precip         = charMissing     ! variable name for precipitation
  character(len=strLen),public    :: vname_time           = charMissing     ! variable name for time
  character(len=strLen),public    :: vname_hruid          = charMissing     ! variable name for runoff hru id
  character(len=strLen),public    :: dname_time           = charMissing     ! dimension name for time
  character(len=strLen),public    :: dname_hruid          = charMissing     ! dimension name for hru in runoff data
  character(len=strLen),public    :: dname_xlon           = charMissing     ! dimension name for x (j, longitude) dimension
  character(len=strLen),public    :: dname_ylat           = charMissing     ! dimension name for y (i, latitude) dimension
  character(len=strLen),public    :: units_qsim           = charMissing     ! units of runoff data
  real(dp)             ,public    :: dt_ro                = realMissing     ! runoff time step (seconds)
  real(dp)             ,public    :: input_fillvalue      = realMissing     ! fillvalue used for input variables (runoff, precipitation, evaporation)
  character(len=strLen),public    :: ro_time_units        = charMissing     ! time units used in ro netcdf. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted
  character(len=strLen),public    :: ro_calendar          = charMissing     ! calendar used in ro netcdf
  character(len=strLen),public    :: ro_time_stamp        = 'start'         ! time stamp used in runoff input - start (default), middle, or end, otherwise error
  ! Water-management input netCDF - water abstraction/infjection or lake target volume
  character(len=strLen),public    :: fname_wm             = ''              ! the txt file name that includes nc files holesing the abstraction, injection, target volume values
  character(len=strLen),public    :: vname_flux_wm        = ''              ! variable name for abstraction or injection from or to a river segment
  character(len=strLen),public    :: vname_vol_wm         = ''              ! variable name for target volume when lake is_lake_sim is on
  character(len=strLen),public    :: vname_time_wm        = ''              ! variable name for time
  character(len=strLen),public    :: vname_segid_wm       = ''              ! variable name for runoff hru id
  character(len=strLen),public    :: dname_time_wm        = ''              ! dimension name for time
  character(len=strLen),public    :: dname_segid_wm       = ''              ! dimension name for hru in runoff data
  real(dp)             ,public    :: dt_wm                = realMissing     ! water-management time step (seconds)
  ! RUNOFF REMAPPING
  logical(lgt),         public    :: is_remap             = .false.         ! logical whether or not runnoff needs to be mapped to river network HRU
  character(len=strLen),public    :: fname_remap          = charMissing     ! runoff mapping netCDF name
  character(len=strLen),public    :: vname_hruid_in_remap = charMissing     ! variable name for river network hru id
  character(len=strLen),public    :: vname_weight         = charMissing     ! variable name for areal weights of runoff HRUs within each river network
  character(len=strLen),public    :: vname_qhruid         = charMissing     ! variable name for runoff HRU ID
  character(len=strLen),public    :: vname_num_qhru       = charMissing     ! variable for numbers of runoff HRUs within each river network HRU
  character(len=strLen),public    :: vname_i_index        = charMissing     ! variable for numbers of y (latitude) index if runoff file is grid
  character(len=strLen),public    :: vname_j_index        = charMissing     ! variable for numbers of x (longitude) index if runoff file is grid
  character(len=strLen),public    :: dname_hru_remap      = charMissing     ! dimension name for river network HRU
  character(len=strLen),public    :: dname_data_remap     = charMissing     ! dimension name for runoff HRU ID
  ! RESTART OPTION
  character(len=strLen),public    :: restart_write        = 'never'         ! restart write option (case-insensitive): never, last, specified, yearly, monthly, daily
  character(len=strLen),public    :: restart_date         = charMissing     ! specifed restart date
  integer(i4b)         ,public    :: restart_month        = 1               ! restart periodic month. Default Jan (write every January of year)
  integer(i4b)         ,public    :: restart_day          = 1               ! restart periodic day.   Default 1st (write every 1st of month)
  integer(i4b)         ,public    :: restart_hour         = 0               ! restart periodic hour.  Default 0hr (write every 00 hr of day)
  character(len=strLen),public    :: fname_state_in       = charMissing     ! name of state file
  ! SPATIAL CONSTANT PARAMETERS
  character(len=strLen),public    :: param_nml            = ''              ! name of the namelist file
  ! GAUGE DATA
  character(len=strLen),public    :: gageMetaFile         = charMissing     ! name of the gauge metadata csv
  logical(lgt),public             :: outputAtGage         = .false.         ! logical; T-> history file output at only gauge points
  character(len=strLen),public    :: fname_gageObs        = charMissing     ! gauge data netcdf name
  character(len=strLen),public    :: vname_gageFlow       = charMissing     ! variable name for gauge flow data
  character(len=strLen),public    :: vname_gageSite       = charMissing     ! variable name for site name data
  character(len=strLen),public    :: vname_gageTime       = charMissing     ! variable name for time data
  character(len=strLen),public    :: dname_gageSite       = charMissing     ! dimension name for gauge site
  character(len=strLen),public    :: dname_gageTime       = charMissing     ! dimension name for time
  integer(i4b)         ,public    :: strlen_gageSite      = 30              ! maximum character length for site name
  ! OUTPUT OPTIONS
  real(dp)             ,public    :: histTimeStamp_offset = 0._dp           ! time stamp offset [second] from a start of time step
  logical(lgt)         ,public    :: outputInflow         = .false.         ! logical; T-> write upstream inflow in history file output
  ! USER OPTIONS
  integer(i4b)         ,public    :: qmodOption           = 0               ! options for streamflow modification (DA): 0-> no DA, 1->direct insertion
  integer(i4b)         ,public    :: QerrTrend            = 1               ! temporal discharge error decreasing trend: 1->constant, 2->linear, 3->logistic, 4->exponential
  integer(i4b)         ,public    :: qBlendPeriod         = 10              ! number of time steps for which streamflow modification is performed through blending observation
  integer(i4b)         ,public    :: hydGeometryOption    = readFromFile    ! option for hydraulic geometry calculations (0=read from file, 1=compute)
  integer(i4b)         ,public    :: topoNetworkOption    = compute         ! option for network topology calculations (0=read from file, 1=compute)
  integer(i4b)         ,public    :: computeReachList     = compute         ! option to compute list of upstream reaches (0=do not compute, 1=compute)
  ! TIME
  character(len=strLen),public    :: time_units           = charMissing     ! time units. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted
  character(len=strLen),public    :: calendar             = charMissing     ! calendar name
  ! MISCELLANEOUS
  logical(lgt)         ,public    :: debug                = .false.         ! print out detaled information
  integer(i4b)         ,public    :: idSegOut             = integerMissing  ! id of outlet stream segment
  integer(i4b)         ,public    :: desireId             = integerMissing  ! turn off checks or speficy reach ID if necessary to print on screen
  ! PFAFCODE
  integer(i4b)         ,public    :: maxPfafLen           = 32              ! maximum digit of pfafstetter code (default 32).
  character(len=1)     ,public    :: pfafMissing          = '0'             ! missing pfafcode (e.g., reach without any upstream area)

  ! CESM Coupling variables
  character(len=32)    ,public    :: bypass_routing_option = 'direct_in_place' ! bypass routing model method: direct_in_place or direct_to_outlet
  character(len=32)    ,public    :: qgwl_runoff_option    = 'threshold'       ! method for handling qgwl runoff: all, negative, or threshold

END MODULE public_var
