module public_var
  ! This module include variables that can be accessed from any other modules and values not altered
  ! except that variables read from control file are populated.

  use nrtype, only: i4b,dp,lgt
  use nrtype, only: strLen  ! string length
  implicit none

  save

  ! ---------- mizuRoute version -------------------------------------------------------------------

  character(len=strLen), parameter, public    :: mizuRouteVersion='v1.2.2'

  ! ---------- common constants ---------------------------------------------------------------------

  ! physical constants
  real(dp),    parameter,public   :: pi=3.14159265359_dp   ! pi

  ! some common constant variables (not likely to change value)
  real(dp),    parameter,public   :: secprmin=60._dp        ! number of seconds in a minute
  real(dp),    parameter,public   :: secprhour=3600._dp     ! number of seconds in an hour
  real(dp),    parameter,public   :: secprday=86400._dp     ! number of seconds in a day
  integer(i4b),parameter,public   :: months_per_yr=12       ! number of months in a year
  integer(i4b),parameter,public   :: days_per_yr=365        ! number of days in a year
  real(dp),    parameter,public   :: hr_per_day = 24.0_dp   ! hours per days
  real(dp),    parameter,public   :: min_per_hour = 60.0_dp ! minutes per hour
  real(dp),    parameter,public   :: verySmall=tiny(1.0_dp) ! a very small number
  real(dp),    parameter,public   :: min_slope=1.e-6_dp     ! minimum slope
  real(dp),    parameter,public   :: runoffMin=1.e-15_dp    ! minimum runoff from each basin
  real(dp),    parameter,public   :: negRunoffTol=-1.e-3_dp ! nagative runoff tolerance

  ! routing related constants
  integer(i4b),parameter,public   :: MAXQPAR=20             ! maximum number of particles

  ! domain decomposition parameters
  integer(i4b),parameter,public   :: maxDomain=150000       ! maximum sub-domains

  ! constants for general use
  real(dp),    parameter,public   :: MinPosVal=1.e-10_dp    ! minimum value for positive value
  integer(i4b),parameter,public   :: integerMissing=-9999   ! missing value for integers
  real(dp),    parameter,public   :: realMissing=-9999._dp  ! missing value for real numbers
  character(5),parameter,public   :: charMissing='empty'    ! missing value for character

  ! I/O related parameters
  integer(i4b),          public   :: iulog=6                ! logical unit identifier

  ! ---------- named variables ----------------------------------------------------------------------

  ! true/false
  integer(i4b), parameter, public :: true=1001                  ! true
  integer(i4b), parameter, public :: false=1002                 ! false

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
  character(len=strLen),public    :: ancil_dir            = ''              ! directory containing ancillary data (network, mapping, namelist)
  character(len=strLen),public    :: input_dir            = ''              ! directory containing input runoff netCDF
  character(len=strLen),public    :: output_dir           = ''              ! directory for routed flow output (netCDF)
  character(len=strLen),public    :: restart_dir          = charMissing     ! directory for restart output (netCDF)
  ! RIVER NETWORK TOPOLOGY
  character(len=strLen),public    :: fname_ntopOld        = ''              ! old filename containing stream network topology information
  logical(lgt)         ,public    :: ntopAugmentMode      = .false.         ! option for river network augmentation mode. terminate the program after writing augmented ntopo.
  character(len=strLen),public    :: fname_ntopNew        = ''              ! new filename containing stream network topology information
  character(len=strLen),public    :: dname_sseg           = ''              ! dimension name of segment in river network data
  character(len=strLen),public    :: dname_nhru           = ''              ! dimension name of hru in river network data
  ! RUNOFF FILE
  character(len=strLen),public    :: fname_qsim           = ''              ! simulated runoff netCDF name
  character(len=strLen),public    :: vname_qsim           = ''              ! variable name for simulated runoff
  character(len=strLen),public    :: vname_time           = ''              ! variable name for time
  character(len=strLen),public    :: vname_hruid          = ''              ! variable name for runoff hru id
  character(len=strLen),public    :: dname_time           = ''              ! dimension name for time
  character(len=strLen),public    :: dname_hruid          = ''              ! dimension name for hru in runoff data
  character(len=strLen),public    :: dname_xlon           = ''              ! dimension name for x (j, longitude) dimension
  character(len=strLen),public    :: dname_ylat           = ''              ! dimension name for y (i, latitude) dimension
  character(len=strLen),public    :: units_qsim           = ''              ! units of simulated runoff data
  real(dp)             ,public    :: dt                   = realMissing     ! time step (seconds)
  real(dp)             ,public    :: ro_fillvalue         = realMissing     ! fillvalue used for runoff depth variable
  logical(lgt)         ,public    :: userRunoffFillvalue  = .false.         ! true -> runoff depth fillvalue used in netcdf is specified here, otherwise -> false
  ! RUNOFF REMAPPING
  logical(lgt),public             :: is_remap             = .false.         ! logical whether or not runnoff needs to be mapped to river network HRU
  character(len=strLen),public    :: fname_remap          = ''              ! runoff mapping netCDF name
  character(len=strLen),public    :: vname_hruid_in_remap = ''              ! variable name for river network hru id
  character(len=strLen),public    :: vname_weight         = ''              ! variable name for areal weights of runoff HRUs within each river network
  character(len=strLen),public    :: vname_qhruid         = ''              ! variable name for runoff HRU ID
  character(len=strLen),public    :: vname_num_qhru       = ''              ! variable for numbers of runoff HRUs within each river network HRU
  character(len=strLen),public    :: vname_i_index        = ''              ! variable for numbers of y (latitude) index if runoff file is grid
  character(len=strLen),public    :: vname_j_index        = ''              ! variable for numbers of x (longitude) index if runoff file is grid
  character(len=strLen),public    :: dname_hru_remap      = ''              ! dimension name for river network HRU
  character(len=strLen),public    :: dname_data_remap     = ''              ! dimension name for runoff HRU ID
  ! RUN CONTROL
  character(len=strLen),public    :: case_name            = ''              ! name of simulation. used as head of model output and restart file
  character(len=strLen),public    :: simStart             = ''              ! date string defining the start of the simulation
  character(len=strLen),public    :: simEnd               = ''              ! date string defining the end of the simulation
  character(len=10)    ,public    :: routOpt              = '0'             ! routing scheme options  0: accum runoff, 1:IRF, 2:KWT, 3:KW, 4:MC, 5:DW
  integer(i4b)         ,public    :: doesBasinRoute       = 1               ! basin routing options   0-> no, 1->IRF, otherwise error
  character(len=strLen),public    :: newFileFrequency     = 'annual'        ! frequency for new output files (day, month, annual, single)
  ! STATES
  character(len=strLen),public    :: restart_write        = 'never'         ! restart write option: N[n]ever-> never write, L[l]ast -> write at last time step, S[s]pecified, Monthly, Daily
  character(len=strLen),public    :: restart_date         = charMissing     ! specifed restart date
  integer(i4b)         ,public    :: restart_month        = 1               ! restart periodic month. Default Jan (write every January of year)
  integer(i4b)         ,public    :: restart_day          = 1               ! restart periodic day.   Default 1st (write every 1st of month)
  integer(i4b)         ,public    :: restart_hour         = 0               ! restart periodic hour.  Default 0hr (write every 00 hr of day)
  character(len=strLen),public    :: fname_state_in       = charMissing     ! name of state file
  ! SPATIAL CONSTANT PARAMETERS
  character(len=strLen),public    :: param_nml            = ''              ! name of the namelist file
  ! USER OPTIONS
  logical(lgt)         ,public    :: qtakeOption          = .false.         ! option for abstraction/injection
  integer(i4b)         ,public    :: hydGeometryOption    = compute         ! option for hydraulic geometry calculations (0=read from file, 1=compute)
  integer(i4b)         ,public    :: topoNetworkOption    = compute         ! option for network topology calculations (0=read from file, 1=compute)
  integer(i4b)         ,public    :: computeReachList     = compute         ! option to compute list of upstream reaches (0=do not compute, 1=compute)
  ! TIME
  character(len=strLen),public    :: time_units           = charMissing     ! time units time units. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted
  character(len=strLen),public    :: calendar             = charMissing     ! calendar name
  ! MISCELLANEOUS
  logical(lgt)         ,public    :: debug                = .false.         ! print out detaled information
  integer(i4b)         ,public    :: idSegOut             = integerMissing  ! id of outlet stream segment
  integer(i4b)         ,public    :: desireId             = integerMissing  ! turn off checks or speficy reach ID if necessary to print on screen
  character(len=strLen),public    :: netcdf_format        = 'netcdf4'       ! netcdf format for output
  ! PFAFCODE
  integer(i4b)         ,public    :: maxPfafLen           = 32              ! maximum digit of pfafstetter code (default 32).
  character(len=1)     ,public    :: pfafMissing          = '0'             ! missing pfafcode (e.g., reach without any upstream area)

end module public_var
