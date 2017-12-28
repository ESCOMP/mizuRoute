module public_var
  ! This module include variables that can be accessed from any other modules and values not altered
  ! Examples of variables are physical parameteres, namelist variable, variables in control etc.
  use nrtype

  implicit none

  save

  ! some common constant variables (not likely to change value)
  integer(i4b),parameter,public   :: strLen=256             ! length of character string
  integer(i4b),parameter,public   :: imiss=-999             ! missing value for integer value
  real(dp),    parameter,public   :: dmiss=-999.0_dp        ! missing value for floating value
  real(dp),    parameter,public   :: secprhour=3600._dp     ! number of seconds in an hour
  real(dp),    parameter,public   :: secprday=86400._dp     ! number of seconds in a day
  real(dp),    parameter,public   :: verySmall=tiny(1.0_dp) ! a very small number
  real(dp),    parameter,public   :: runoffMin=1.e-15_dp    ! minimum runoff from each basin
  real(dp),    parameter,public   :: MinPosVal=1.e-10_dp    ! minimum value for positive value
  integer(i4b),parameter,public   :: MAXQPAR=20             ! maximum number of particles

  ! Control file variables
  ! DIRECTORIES
  character(len=strLen),public    :: ancil_dir              ! directory containing ancillary data
  character(len=strLen),public    :: input_dir              ! directory containing input data
  character(len=strLen),public    :: output_dir             ! directory containing output data
  ! RUNOFF FILE
  character(len=strLen),public    :: fname_qsim             ! simulated runoff netCDF name
  character(len=strLen),public    :: vname_qsim             ! variable name for simulated runoff
  character(len=strLen),public    :: vname_time             ! variable name for time
  character(len=strLen),public    :: vname_hruid            ! variable name for runoff hru id
  character(len=strLen),public    :: units_qsim             ! units of simulated runoff data
  character(len=strLen),public    :: dname_time             ! dimension name for time
  character(len=strLen),public    :: dname_hruid            ! dimension name for hru in runoff data
  character(len=strLen),public    :: units_time             ! time units
  ! RIVER NETWORK TOPOLOGY
  character(len=strLen),public    :: fname_ntop             ! filename containing stream network topology information
  ! ROUTED FLOW OUTPUT
  character(len=strLen),public    :: fname_output           ! name of output file
  ! STATES
  logical(lgt)         ,public    :: isRestart              ! restart option: True-> model run with restart, F -> model run with empty channels
  character(len=strLen),public    :: fname_state_in         ! name of state file
  character(len=strLen),public    :: fname_state_out        ! name of state file
  ! PARAMETER
  character(len=strLen),public    :: param_nml              ! name of the namelist file
  ! RUNOFF REMAPPING
  logical(lgt),public             :: is_remap               ! logical whether or not runnoff needs to be mapped to river network HRU
  character(len=strLen),public    :: fname_remap            ! runoff mapping netCDF name
  character(len=strLen),public    :: vname_hruid_in_remap   ! variable name for river network hru id
  character(len=strLen),public    :: vname_weight           ! variable name for areal weights of runoff HRUs within each river network
  character(len=strLen),public    :: vname_qhruid           ! variable name for runoff HRU ID
  character(len=strLen),public    :: vname_num_qhru         ! variable for numbers of runoff HRUs within each river network HRU
  character(len=strLen),public    :: dname_hru_remap        ! dimension name for river network HRU
  character(len=strLen),public    :: dname_data_remap       ! dimension name for runoff HRU ID
  ! MISCELLANEOUS
  real(dp)             ,public    :: dt                     ! time step (seconds)
  integer(i4b)         ,public    :: iSegOut                ! index of outlet stream segment
  integer(i4b)         ,public    :: routOpt                ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error

end module public_var
