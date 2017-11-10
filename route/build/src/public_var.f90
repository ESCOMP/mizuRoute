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
  character(len=strLen),public    :: ancil_dir              ! directory containing ancillary data
  character(len=strLen),public    :: input_dir              ! directory containing input data
  character(len=strLen),public    :: output_dir             ! directory containing output data
  character(len=strLen),public    :: fname_qsim             ! filename containing simulated runoff
  character(len=strLen),public    :: vname_qsim             ! variable name for simulated runoff
  character(len=strLen),public    :: vname_hruid            ! coordinate name for the HRUid
  character(len=strLen),public    :: vname_time             ! coordinate name for time
  character(len=strLen),public    :: units_qsim             ! units of simulated runoff data
  character(len=strLen),public    :: units_time             ! time units
  character(len=strLen),public    :: fname_ntop             ! filename containing stream network topology information
  character(len=strLen),public    :: fname_output           ! name of output file
  character(len=strLen),public    :: fname_state_in         ! name of state file
  character(len=strLen),public    :: fname_state_out        ! name of state file
  character(len=strLen),public    :: param_nml              ! name of the namelist file
  real(dp)             ,public    :: dt                     ! time step (seconds)
  integer(i4b)         ,public    :: iSegOut                ! index of outlet stream segment
  integer(i4b)         ,public    :: routOpt                ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error
  logical(lgt)         ,public    :: isRestart              ! restart option: True-> model run with restart, F -> model run with empty channels 

end module public_var
