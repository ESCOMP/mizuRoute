module public_var
  ! This module include variables that can be accessed from any other modules and values not altered
  ! Examples of variables are physical parameteres, namelist variable, variables in control etc.
  use nrtype, only: i4b,dp,lgt
  use nrtype, only: strLen  ! string length
  implicit none

  save

  ! ---------- common constants ---------------------------------------------------------------------

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
  real(dp),    parameter,public   :: MinPosVal=1.e-10_dp    ! minimum value for positive value
  integer(i4b),parameter,public   :: MAXQPAR=20             ! maximum number of particles
  integer(i4b),          public   :: maxPfafLen=32          ! maximum digit of pfafstetter code (default 32)
  integer(i4b),          public   :: nThresh=10000          ! maximum river segments in each river basin
  integer(i4b),parameter,public   :: integerMissing=-9999   ! missing value for integers
  real(dp),    parameter,public   :: realMissing=-9999._dp  ! missing value for real numbers

  ! ---------- named variables ----------------------------------------------------------------------

  ! output file frequency
  integer(i4b) , parameter        :: annual=1001            ! named variable for yearly output files
  integer(i4b) , parameter        :: month=1002             ! named variable for monthly output files
  integer(i4b) , parameter        :: day=1003               ! named variable for daily output files

  ! compute versus read from file
  integer(i4b), parameter,public  :: compute=1              ! compute given variable
  integer(i4b), parameter,public  :: doNotCompute=0         ! do not compute given variable
  integer(i4b), parameter,public  :: readFromFile=0         ! read given variable from a file

  ! routing methods
  integer(i4b), parameter,public  :: allRoutingMethods=0    ! all routing methods
  integer(i4b), parameter,public  :: impulseResponseFunc=1  ! impulse response function
  integer(i4b), parameter,public  :: kinematicWave=2        ! kinematic wave

  ! ---------- variables in the control file --------------------------------------------------------

  ! Control file variables
  ! DIRECTORIES
  character(len=strLen),public    :: ancil_dir            = ''              ! directory containing ancillary data
  character(len=strLen),public    :: input_dir            = ''              ! directory containing input data
  character(len=strLen),public    :: output_dir           = ''              ! directory containing output data
  ! SIMULATION TIME
  character(len=strLen),public    :: simStart             = ''              ! date string defining the start of the simulation
  character(len=strLen),public    :: simEnd               = ''              ! date string defining the end of the simulation
  ! RUNOFF FILE
  character(len=strLen),public    :: fname_qsim           = ''              ! simulated runoff netCDF name
  character(len=strLen),public    :: vname_qsim           = ''              ! variable name for simulated runoff
  character(len=strLen),public    :: vname_time           = ''              ! variable name for time
  character(len=strLen),public    :: vname_hruid          = ''              ! variable name for runoff hru id
  character(len=strLen),public    :: units_qsim           = ''              ! units of simulated runoff data
  character(len=strLen),public    :: dname_time           = ''              ! dimension name for time
  character(len=strLen),public    :: dname_hruid          = ''              ! dimension name for hru in runoff data
  character(len=strLen),public    :: dname_xlon           = ''              ! dimension name for x (j, longitude) dimension
  character(len=strLen),public    :: dname_ylat           = ''              ! dimension name for y (i, latitude) dimension
  ! RIVER NETWORK TOPOLOGY
  character(len=strLen),public    :: fname_ntopOld        = ''              ! old filename containing stream network topology information
  logical(lgt)         ,public    :: ntopWriteOption      = .false.         ! option for writing augmented network topology (false=no,true=yes)
  character(len=strLen),public    :: fname_ntopNew        = ''              ! new filename containing stream network topology information
  character(len=strLen),public    :: dname_sseg           = ''              ! dimension name of segment in river network data
  character(len=strLen),public    :: dname_nhru           = ''              ! dimension name of hru in river network data
  ! ROUTED FLOW OUTPUT
  character(len=strLen),public    :: fname_output         = ''              ! name of output file
  integer(i4b)         ,public    :: newFileFrequency     = annual          ! frequency for new output files (day, month, annual)
  ! USER OPTIONS
  integer(i4b)         ,public    :: hydGeometryOption    = compute         ! option for hydraulic geometry calculations (0=read from file, 1=compute)
  integer(i4b)         ,public    :: topoNetworkOption    = compute         ! option for network topology calculations (0=read from file, 1=compute)
  integer(i4b)         ,public    :: computeReachList     = doNotCompute    ! option to compute list of upstream reaches (0=do not compute, 1=compute)
  ! STATES
  logical(lgt)         ,public    :: isRestart            = .false.         ! restart option: True-> model run with restart, F -> model run with empty channels
  character(len=strLen),public    :: fname_state_in       = ''              ! name of state file
  character(len=strLen),public    :: fname_state_out      = ''              ! name of state file
  ! PARAMETER
  character(len=strLen),public    :: param_nml            = ''              ! name of the namelist file
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
  ! TIME
  real(dp)             ,public    :: dt                   = realMissing     ! time step (seconds)
  character(len=strLen),public    :: time_units           = ''              ! time units (seconds, hours, or days)
  character(len=strLen),public    :: calendar             = ''              ! calendar name
  real(dp)             ,public    :: startJulday          = realMissing     ! julian day: start
  real(dp)             ,public    :: endJulday            = realMissing     ! julian day: end
  real(dp)             ,public    :: refJulday            = realMissing     ! julian day: reference
  ! MISCELLANEOUS
  integer(i4b)         ,public    :: idSegOut             = integerMissing  ! id of outlet stream segment
  integer(i4b)         ,public    :: routOpt              = integerMissing  ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error
  integer(i4b)         ,public    :: desireId             = integerMissing  ! turn off checks or speficy reach ID if necessary to print on screen
  integer(i4b)         ,public    :: doesBasinRoute       = 1               ! basin routing options   0-> no, 1->IRF, otherwise error

end module public_var
