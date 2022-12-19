MODULE model_setup

USE nrtype
USE public_var,        ONLY: iulog
USE public_var,        ONLY: debug
USE public_var,        ONLY: integerMissing
USE public_var,        ONLY: realMissing
USE public_var,        ONLY: charMissing
USE nr_utils,          ONLY: arth
USE nr_utils,          ONLY: unique         ! get unique element array
USE nr_utils,          ONLY: indexx         ! get rank of data value

implicit none

private
public :: init_mpi
public :: init_data

CONTAINS

 ! *********************************************************************
 ! public subroutine: initialize MPI for stand-alone program
 ! *********************************************************************
 SUBROUTINE init_mpi()

  ! Initialize MPI and get OMP thread

  USE globalData,      ONLY: mpicom_route ! communicator id
  USE init_model_data, ONLY: get_mpi_omp
  USE mpi_utils,       ONLY: shr_mpi_init

  implicit none
  ! input:  None
  ! output: None
  ! local variables
  character(len=strLen)       :: message             ! error message

  ! initialize error control
  message='init_mpi/'

  call shr_mpi_init(mpicom_route, message)

  call get_mpi_omp(mpicom_route)

 END SUBROUTINE init_mpi


 ! *********************************************************************
 ! public subroutine: initialize all the data - For stand-alone
 ! *********************************************************************
 SUBROUTINE init_data(pid,           & ! input: proc id
                      nNodes,        & ! input: number of procs
                      comm,          & ! input: communicator
                      ierr, message)   ! output: error control

   USE public_var,          ONLY: continue_run        ! T-> append output in existing history files. F-> write output in new history file
   USE mpi_process,         ONLY: pass_global_data    ! mpi globaldata copy to slave proc
   USE init_model_data,     ONLY: init_ntopo_data
   USE init_model_data,     ONLY: init_state_data
   USE write_simoutput_pio, ONLY: init_histFile       ! open existing history file to append (only continue_run is true)

   implicit none
   ! Argument variables
   integer(i4b),              intent(in)    :: pid              ! proc id
   integer(i4b),              intent(in)    :: nNodes           ! number of procs
   integer(i4b),              intent(in)    :: comm             ! communicator
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variables
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ierr=0; message='init_data/'

   ! runoff input files initialization
   call init_inFile_pop(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! time initialization
   call init_time(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! network topology data initialization
   call init_ntopo_data(pid, nNodes, comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! Initialize forcing (ro, pr, evap), water-management data and mapping data(only at main core)
   if (pid==0) then
     call init_forc_data(ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   ! broadcast public and some global variables
   call pass_global_data(comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! channel state initialization
   call init_state_data(pid, nNodes, comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (continue_run) then
     call init_histFile(ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

 END SUBROUTINE init_data


 ! *********************************************************************
 ! private subroutine: initiate the reading of the netcdf files for runoff
 ! or abstraction or injection
 ! *********************************************************************

 SUBROUTINE init_inFile_pop (ierr, message)  ! output

  USE public_var, ONLY: input_dir               ! directory containing the text files of fname_qsim and fname_wm
  USE public_var, ONLY: fname_qsim              ! simulated runoff txt file that includes the NetCDF file names
  USE public_var, ONLY: vname_time              ! variable name for time
  USE public_var, ONLY: dname_time              ! dimension name for time
  USE public_var, ONLY: fname_wm                ! simulated runoff txt file that includes the NetCDF file names
  USE public_var, ONLY: vname_time_wm           ! variable name for time
  USE public_var, ONLY: dname_time_wm           ! dimension name for time
  USE globalData, ONLY: inFileInfo_ro           ! metadata of the ro/evap/p input files
  USE globalData, ONLY: inFileInfo_wm           ! metadata of the input files for abstration, injection and target volume
  USE public_var, ONLY: is_lake_sim             ! logical whether or not lake should be simulated
  USE public_var, ONLY: is_flux_wm              ! logical whether or not abstraction and injection should be read from the file
  USE public_var, ONLY: is_vol_wm               ! logical whether or not target volume for lakes should be read

  ! Argument variables
  integer(i4b),         intent(out)    :: ierr             ! error code
  character(*),         intent(out)    :: message          ! error message
  ! local variables
  character(len=strLen)                :: cmessage         ! error message of downwind routine

  ierr=0; message='init_inFile_pop/'

  call inFile_pop(input_dir,         & ! input: name of the directory of the txt file
                  fname_qsim,        & ! input: name of the txt file hold the nc file names
                  vname_time,        & ! input: name of variable time in the nc files
                  dname_time,        & ! input: name of dimention time in the nc files
                  inFileInfo_ro,     & ! output: input file information
                  ierr, cmessage)      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

  if ((is_flux_wm).or.(is_vol_wm.and.is_lake_sim)) then     ! if either of abstraction injection or target volume is activated
    call inFile_pop(input_dir,            & ! input: name of the directory of the txt file
                    fname_wm,             & ! input: name of the txt file hold the nc file names
                    vname_time_wm,        & ! input: name of variable time in the nc files
                    dname_time_wm,        & ! input: name of dimention time in the nc files
                    inFileInfo_wm,        & ! output: input file information
                    ierr, cmessage)         ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

    call inFile_sync_time(inFileInfo_ro,      & ! input: the structure of simulated runoff, evapo and
                          inFileInfo_wm,      & ! inout: input file information
                          ierr, cmessage)       ! output: error control

  endif

  END SUBROUTINE init_inFile_pop

 ! *********************************************************************
 ! private subroutine: read the name of the netcdf that is specified
 ! in a text file, populates the filed of inFiledata dataType
 ! *********************************************************************
 SUBROUTINE inFile_pop(dir_name,         & ! input: name of the directory of the txt file
                       file_name,        & ! input: name of the txt file hold the nc file names
                       time_var_name,    & ! input: name of variable time in the nc files
                       time_dim_name,    & ! input: name of dimention time in the nc files
                       inputFileInfo,    & ! output: input file information
                       ierr, message)      ! output: error control

  USE dataTypes,           ONLY: inFileInfo     ! the data type for storing the infromation of the nc files and its attributes
  USE datetime_data,       ONLY: datetime       ! datetime data
  USE ascii_utils,         ONLY: file_open      ! open file (performs a few checks as well)
  USE ascii_utils,         ONLY: get_vlines     ! get a list of character strings from non-comment lines
  USE ncio_utils,          ONLY: get_nc         ! Read netCDF variable data
  USE ncio_utils,          ONLY: get_var_attr   ! Read attributes variables
  USE ncio_utils,          ONLY: get_nc_dim_len ! get the nc dimension length
  USE public_var,          ONLY: time_units     ! get the time units from control file and replace if not provided in nc files
  USE public_var,          ONLY: calendar       ! get the calendar from control file and replace if not provided in nc files

  ! Argument variables
  character(len=strLen), intent(in)                 :: dir_name         ! the name of the directory that the txt file located
  character(len=strLen), intent(in)                 :: file_name        ! the name of the file that include the nc file names
  character(len=strLen), intent(in)                 :: time_var_name    ! the name of the time variable
  character(len=strLen), intent(in)                 :: time_dim_name    ! the name of dimension time
  type(infileinfo),      intent(inout), allocatable :: inputFileInfo(:) ! the name of structure that hold the infile information
  integer(i4b),          intent(out)                :: ierr             ! error code
  character(*),          intent(out)                :: message          ! error message
  ! local varibales
  integer(i4b)                                      :: unit             ! file unit (free unit output from file_open)
  character(len=7)                                  :: t_unit           ! time units. "<time_step> since yyyy-MM-dd hh:mm:ss"
  integer(i4b)                                      :: iFile            ! counter for forcing files
  integer(i4b)                                      :: nFile            ! number of nc files identified in the text file
  integer(i4b)                                      :: nTime            ! hard coded for now
  type(datetime)                                    :: refDatetime      ! reference datetime for each file
  real(dp)                                          :: convTime2Days    ! conversion of the day to the local time
  character(len=strLen)                             :: infilename       ! input filename
  character(len=strLen),allocatable                 :: dataLines(:)     ! vector of lines of information (non-comment lines)
  character(len=strLen)                             :: filenameData     ! name of forcing datafile
  character(len=strLen)                             :: cmessage         ! error message of downwind routine

  ierr=0; message='inFile_pop/'

  ! build filename and its path containing list of NetCDF files
  infilename = trim(dir_name)//trim(file_name)

  ! open the text file
  call file_open(trim(infilename),unit,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! get a list of character strings from non-commented lines
  call get_vlines(unit,dataLines,ierr,cmessage)
  if(ierr/=0)then; ierr=20; message=trim(message)//trim(cmessage); return; end if
  nFile = size(dataLines) ! get the name of the lines in the file

  ! allocate space for forcing information
  allocate(inputFileInfo(nFile), stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for forcFileInfo'; return; end if

  ! poputate the forcingInfo structure with filenames, and time variables/attributes
  do iFile=1,nFile

   ! split the line into "words" (expect one word: the file describing forcing data for that index)
   read(dataLines(iFile),*,iostat=ierr) filenameData
   if(ierr/=0)then; message=trim(message)//'problem reading a line of data from file ['//trim(infilename)//']'; return; end if

   ! set forcing file name
   inputFileInfo(iFile)%infilename = trim(filenameData)

   ! get the time units, assuming the water managment nc files has the same calendar as the first
   if (trim(time_units) == charMissing) then
     call get_var_attr(trim(dir_name)//trim(inputFileInfo(iFile)%infilename), &
                       trim(time_var_name), 'units', inputFileInfo(iFile)%unit, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else
     inputFileInfo(iFile)%unit = trim(time_units)
   end if

   ! get the calendar, assuming the water managment nc files has the same calendar as the first
   if (trim(calendar) == charMissing) then
     call get_var_attr(trim(dir_name)//trim(inputFileInfo(iFile)%infilename), &
                       trim(time_var_name), 'calendar', inputFileInfo(iFile)%calendar, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else
     inputFileInfo(iFile)%calendar = trim(calendar)
   end if

   ! get the dimension of the time to populate nTime and pass it to the get_nc file
   call get_nc_dim_len(trim(dir_name)//trim(inputFileInfo(iFile)%infilename), &
                       trim(time_dim_name), inputFileInfo(iFile)%nTime, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   nTime = inputFileInfo(iFile)%nTime ! the length of time varibale for each nc file

   ! allocate space for time varibale of each file
   allocate(inputFileInfo(iFile)%timeVar(nTime))
   if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for inputFileInfo(:)%timeVar'; return; end if

   ! get the time varibale
   call get_nc(trim(dir_name)//trim(inputFileInfo(iFile)%infilename), &
               time_var_name, inputFileInfo(iFile)%timeVar, 1, nTime, ierr, cmessage) ! does it needs timeVar(:)

   ! get the time multiplier needed to convert time to units of days for each nc file
   t_unit = trim( inputFileInfo(iFile)%unit(1:index(inputFileInfo(iFile)%unit,' ')) )
   select case( trim(t_unit) )
    case('seconds','second','sec','s'); convTime2Days=86400._dp
    case('minutes','minute','min','m'); convTime2Days=1440._dp
    case('hours'  ,'hour'  ,'hr' ,'h'); convTime2Days=24._dp
    case('days'   ,'day'   ,'d');       convTime2Days=1._dp
    case default
      ierr=20; message=trim(message)//'<time_units>= '//trim(t_unit)//': <time_units> must be seconds, minutes, hours or days.'; return
   end select
   ! convert timeValue unit to day
   inputFileInfo(iFile)%timeVar(:) = inputFileInfo(iFile)%timeVar(:)/convTime2Days

   ! get the reference julian day from the nc file
   call refDatetime%str2datetime(trim(inputFileInfo(iFile)%unit), ierr, message)
   call refDatetime%julianday(trim(inputFileInfo(iFile)%calendar), inputFileInfo(iFile)%ncrefjulday, ierr, cmessage)
   if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [ncrefjulday]'; return; endif

   ! populated the index of the iTimebound for each nc file
   if (iFile==1) then
    inputFileInfo(iFile)%iTimebound(1) = 1
    inputFileInfo(iFile)%iTimebound(2) = nTime
   else ! if multiple files specfied in the txt file
    inputFileInfo(iFile)%iTimebound(1) = inputFileInfo(iFile-1)%iTimebound(2) + 1 ! the last index from the perivous nc file + 1
    inputFileInfo(iFile)%iTimebound(2) = inputFileInfo(iFile-1)%iTimebound(2) + nTime ! the last index from the perivous nc file + 1
   endif

  end do

  ! close ascii file
  close(unit=unit,iostat=ierr); if(ierr/=0)then;message=trim(message)//'problem closing forcing file list'; return; end if

 END SUBROUTINE inFile_pop


 ! *********************************************************************
 ! private subroutine: to synchronize the iTimebound of
 ! the inputFileInfo_wm to match the inputFileInfo
 ! *********************************************************************
 SUBROUTINE inFile_sync_time(inputFileInfo_ro,   & ! input: the structure of simulated runoff, evapo and
                             inputFileInfo_wm,   & ! inout: input file information
                             ierr, message)        ! output: error control

  USE dataTypes,  ONLY: infileinfo    ! the data type for storing the infromation of the nc files and its attributes
  USE public_var, ONLY: dt_ro         ! forcing time step in seconds
  USE public_var, ONLY: secprday      ! conversion of steps in days to seconds

  ! Argument variables
  type(infileinfo),    intent(in)       :: inputFileInfo_ro(:)   ! the name of structure that hold the infile information
  type(infileinfo),    intent(inout)    :: inputFileInfo_wm(:)   ! the name of structure that hold the infile information
  integer(i4b),        intent(out)      :: ierr                  ! error code
  character(*),        intent(out)      :: message               ! error message
  ! local variables
  integer(i4b)                          :: nt
  integer(i4b)                          :: nFile                 ! number of nc files for the simulated runoff
  integer(i4b)                          :: nFile_wm              ! number of nc files for the water managent
  integer(i4b)                          :: iFile                 ! for loop over the nc files
  real(dp)                              :: day_runoff_start      ! the Julian day that runoff starts
  real(dp)                              :: day_start_diff        ! conversion of the day to the local time
  real(dp)                              :: day_end_diff          ! conversion of the day to the local time

  ierr=0; message='inFile_sync_time/'

  ! set the reference julday based on the first nc file of simulation
  nFile               = size(inputFileInfo_ro)
  nFile_wm            = size(inputFileInfo_wm)
  day_runoff_start    = inputFileInfo_ro(1)%timeVar(1)+inputFileInfo_ro(1)%ncrefjulday

  do iFile=1,nFile_wm

    nt = inputFileInfo_wm(iFile)%nTime ! get the number of time

    day_start_diff = inputFileInfo_wm(iFile)%timeVar(1) +inputFileInfo_wm(iFile)%ncrefjulday - day_runoff_start
    day_end_diff   = inputFileInfo_wm(iFile)%timeVar(nt)+inputFileInfo_wm(iFile)%ncrefjulday - day_runoff_start

    inputFileInfo_wm(iFile)%iTimebound(1) = int(day_start_diff * secprday/dt_ro, kind=i4b) + 1 ! to convert the day difference into time step difference
    inputFileInfo_wm(iFile)%iTimebound(2) = int(day_end_diff   * secprday/dt_ro, kind=i4b) + 1 ! to convert the day difference into time step difference

  end do

  ! checks if the staring and ending iTime of the inputFileInfo_wm overlap with the inputFileInfo of simulated runoff, evapo and precip
  if (inputFileInfo_wm(1)%iTimebound(1) > inputFileInfo_ro(1)%iTimebound(1)) then
    print*, "The first water managment netCDF starts later than the first runoff, evapo and precip netCDF and may cause crash"
  endif
  if (inputFileInfo_wm(nFile_wm)%iTimebound(2) < inputFileInfo_ro(nFile)%iTimebound(2)) then
    print*, "The last water managment netCDf ends earlier than the last runoff, evapo and precip netCDF and may cause crash"
  endif
  if (inputFileInfo_wm(1)%iTimebound(1) < inputFileInfo_ro(1)%iTimebound(1)) then
    print*, "The water managment netCDF starts earlier than the last runoff, evapo and precip netCDF"
  endif
  if (inputFileInfo_wm(nFile_wm)%iTimebound(2) > inputFileInfo_ro(nFile)%iTimebound(2)) then
    print*, "The water managment netCDf ends later than the last runoff, evapo and precip netCDF"
  endif

 END SUBROUTINE inFile_sync_time

 ! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 SUBROUTINE init_time(ierr, message)
  ! purpose: Save the following time related global variables
  ! - time_units
  ! - calendar
  ! - timeVar
  ! - iTime
  ! - begDatetime, endDatetime:   simulationg start and end datetime
  ! - restDatetime, dropDatetime

  USE ascii_utils,   ONLY: lower                    ! convert string to lower case
  USE datetime_data, ONLY: datetime                 ! datetime data
  USE public_var,    ONLY: time_units               ! time units (seconds, hours, or days)
  USE public_var,    ONLY: simStart                 ! date string defining the start of the simulation
  USE public_var,    ONLY: simEnd                   ! date string defining the end of the simulation
  USE public_var,    ONLY: calendar                 ! calendar used for simulation
  USE public_var,    ONLY: dt                       ! simulation time step in second
  USE public_var,    ONLY: continue_run             ! logical to indicate sppend output in existing history file
  USE public_var,    ONLY: secprday                 ! unit conversion from day to sec
  USE public_var,    ONLY: restart_write            ! restart write option
  USE public_var,    ONLY: restart_date             ! restart datetime
  USE public_var,    ONLY: restart_month            ! periodic restart month
  USE public_var,    ONLY: restart_day              ! periodic restart day
  USE public_var,    ONLY: restart_hour             ! periodic restart hr
  USE public_var,    ONLY: maxTimeDiff              ! time difference tolerance for input checks
  USE globalData,    ONLY: timeVar                  ! time variables (unit given by runoff data)
  USE globalData,    ONLY: iTime                    ! time index at simulation time step
  USE globalData,    ONLY: simDatetime              ! model time data (yyyy:mm:dd:hh:mm:ss)
  USE globalData,    ONLY: begDatetime              ! simulation begin datetime data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,    ONLY: endDatetime              ! simulation end datetime data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,    ONLY: restDatetime             ! restart time data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,    ONLY: dropDatetime             ! restart dropoff calendar date/time
  USE globalData,    ONLY: roBegDatetime            ! forcing data start datetime data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,    ONLY: wmBegDatetime            ! water-managment data start datetime data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,    ONLY: infileinfo_ro            ! the information of the input files
  USE globalData,    ONLY: infileinfo_wm            ! the information of the input files
  USE public_var,    ONLY: is_lake_sim              ! logical whether or not lake simulations are activated
  USE public_var,    ONLY: is_flux_wm               ! logical whether or not abstraction and injection should be read from the file
  USE public_var,    ONLY: is_vol_wm                ! logical whether or not target volume for lakes should be read

  implicit none

  ! Argument variables
  integer(i4b),              intent(out)   :: ierr                ! error code
  character(*),              intent(out)   :: message             ! error message
  ! local variables
  integer(i4b)                             :: ix
  integer(i4b)                             :: counter
  integer(i4b)                             :: nTime
  integer(i4b)                             :: nTime_wm
  integer(i4b)                             :: nt
  integer(i4b)                             :: nFile               ! number of nc files
  integer(i4b)                             :: nFile_wm            ! number of nc files
  character(len=strLen)                    :: calendar_wm         ! calendar used in water management input file
  character(len=strLen)                    :: t_unit              ! time units. "<time_step> since yyyy-MM-dd hh:mm:ss"
  integer(i4b)                             :: iFile               ! for loop over the nc files
  type(datetime), allocatable              :: roCal(:)
  type(datetime), allocatable              :: wmCal(:)
  type(datetime)                           :: dummyDatetime       ! temp datetime
  integer(i4b)                             :: nDays               ! number of days in a month
  real(dp)                                 :: juldaySim           ! julian day of first simulation time
  real(dp), allocatable                    :: roJulday(:)         ! julian day in runoff data
  real(dp), allocatable                    :: roJulday_diff(:)    ! the difference of two concequative elements in roJulyday
  real(dp)                                 :: refJulday           ! reference julian day in runoff input file
  real(dp), allocatable                    :: roJulday_wm(:)      ! Julian day of concatenated netCDF for water management
  real(dp), allocatable                    :: roJulday_diff_wm(:) ! the difference of two concequative elements in roJulyday_wm
  character(len=strLen)                    :: cmessage            ! error message of downwind routine
  character(len=50)                        :: fmt1='(a,I4,a,I2.2,a,I2.2,x,I2.2,a,I2.2,a,F5.2)'

  ierr=0; message='init_time/'

  ! Set time attributes for continuous time variables (saved in globalData to use for output)
  ! Use from 1st file
  ! These are used for history output
  calendar   = inFileInfo_ro(1)%calendar
  time_units = inFileInfo_ro(1)%unit

  ! get reference julianday in the 1st file
  refJulday  = inFileInfo_ro(1)%ncrefjulday

  ! get the number of the total time length of all the nc files
  nFile = size(inFileInfo_ro)
  nTime = sum(inFileInfo_ro(:)%nTime)

  ! Define time varialbes: timeVar and roJulday
  allocate(roJulday(nTime), roCal(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! Get roJulday: Julian day series of concatenated netCDF
  counter = 1;
  do iFile=1,nFile
    nt = inFileInfo_ro(iFile)%nTime
    roJulday(counter:counter+nt-1) = &
    inFileInfo_ro(iFile)%timeVar(1:nt)+inFileInfo_ro(iFile)%ncrefjulday
    counter = counter + inFileInfo_ro(iFile)%nTime
  end do

  do ix=1,nTime
    call roCal(ix)%jul2datetime(roJulday(ix), calendar, ierr, cmessage)
  end do

  ! save water management data starting datetime
  roBegDatetime = roCal(1)

  call begDatetime%str2datetime(simStart, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [begDatetime]'; return; endif

  call endDatetime%str2datetime(simEnd, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endDatetime]'; return; endif

  ! check that the dates are aligned
  if(endDatetime < begDatetime) then
    write(cmessage,'(7a)') 'simulation end is before simulation start:', new_line('a'), '<sim_start>= ', trim(simStart), new_line('a'), '<sim_end>= ', trim(simEnd)
    ierr=20; message=trim(message)//trim(cmessage); return
  endif

  ! check sim_start is after the last time step in runoff data
  if (begDatetime > roCal(nTime)) then
    write(iulog,'(2a)') new_line('a'),'ERROR: <sim_start> is after the last time step in input runoff'
    write(iulog,fmt1)  ' runoff_end  : ', roCal(nTime)%year(),'-',roCal(nTime)%month(),'-',roCal(nTime)%day(), roCal(nTime)%hour(),':', roCal(nTime)%minute(),':',roCal(nTime)%sec()
    write(iulog,fmt1)  ' <sim_start> : ', begDatetime%year(),'-',begDatetime%month(),'-',begDatetime%day(), begDatetime%hour(),':', begDatetime%minute(),':',begDatetime%sec()
    ierr=20; message=trim(message)//'check <sim_start> against runoff input time'; return
  endif

  ! Compare sim_start vs. time at first time step in runoff data
  if (begDatetime < roCal(1)) then
    write(iulog,'(2a)') new_line('a'),'WARNING: <sim_start> is before the first time step in input runoff'
    write(iulog,fmt1)  ' runoff_start: ', roCal(1)%year(),'-',roCal(1)%month(),'-',roCal(1)%day(), roCal(1)%hour(),':', roCal(1)%minute(),':',roCal(1)%sec()
    write(iulog,fmt1)  ' <sim_start> : ', begDatetime%year(),'-',begDatetime%month(),'-',begDatetime%day(), begDatetime%hour(),':', begDatetime%minute(),':',begDatetime%sec()
    write(iulog,'(a)') ' Reset <sim_start> to runoff_start'
    begDatetime = roCal(1)
  endif

  ! Compare sim_end vs. time at last time step in runoff data
  if (endDatetime > roCal(nTime)) then
    write(iulog,'(2a)')  new_line('a'),'WARNING: <sim_end> is after the last time step in input runoff'
    write(iulog,fmt1)   ' runoff_end: ', roCal(nTime)%year(),'-',roCal(nTime)%month(),'-',roCal(nTime)%day(), roCal(nTime)%hour(),':', roCal(nTime)%minute(),':',roCal(nTime)%sec()
    write(iulog,fmt1)   ' <sim_end> : ', endDatetime%year(),'-',endDatetime%month(),'-',endDatetime%day(), endDatetime%hour(),':', endDatetime%minute(),':',endDatetime%sec()
    write(iulog,'(a)')  ' Reset <sim_end> to runoff_end'
    endDatetime = roCal(nTime)
  endif

  ! check if the julian day of contacenated files do not have overlap or gap if nTime is larger than 1
  if (nTime > 1) then
    allocate(roJulday_diff(nTime-1), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    ! calculate the difference of consequative time in julian day
    roJulday_diff = roJulday(1:nTime-1) - roJulday(2:nTime)
    ! check if the difference are identical otherwise error and terminate
    if ( any(abs(roJulday_diff-roJulday_diff(1)) > maxTimeDiff) ) then
      write(iulog,'(2a)') new_line('a'),'ERROR: time spacing in netCDF input(s) is not consistent within tolerance maxTimeDiff = ',maxTimeDiff
      ierr=20; message=trim(message)//'make sure the input netCDF files do not have time overlaps or gaps'; return
    end if
  endif

  ! water management options on
  if ((is_flux_wm).or.(is_vol_wm.and.is_lake_sim)) then

    calendar_wm = inFileInfo_wm(1)%calendar

    ! get the number of the total time length of all the water management nc files
    nFile_wm = size(inFileInfo_wm)
    nTime_wm = sum(inFileInfo_wm(:)%nTime)

    ! Define time varialbes: roJulday_wm
    allocate(roJulday_wm(nTime_wm), wmCal(nTime), stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! roJulday_wm: Julian day of concatenated netCDF for water management
    counter = 1;
    do iFile=1,nFile_wm
      nt = inFileInfo_wm(iFile)%nTime
      roJulday_wm(counter:counter+nt-1) = &
      inFileInfo_wm(iFile)%timeVar(1:nt)+inFileInfo_wm(iFile)%ncrefjulday
      counter = counter + inFileInfo_wm(iFile)%nTime
    enddo

    do ix=1,nTime
      call wmCal(ix)%jul2datetime(roJulday_wm(ix), calendar_wm, ierr, cmessage)
    end do

    ! save water management data starting datetime
    wmBegDatetime = wmCal(1)

    ! check sim_start is after the last time step in water management data
    if(begDatetime > wmCal(nTime_wm)) then
      write(iulog,'(2a)') new_line('a'),'ERROR: <sim_start> is after the last time step in input runoff'
      ierr=20; message=trim(message)//'check <sim_start> against water management input time'; return
    endif

    ! check sim_end is before the first time step in water management data
    if(endDatetime < wmCal(1)) then
      write(iulog,'(2a)') new_line('a'),'ERROR: <sim_end> is before the last time step in input runoff'
      ierr=20; message=trim(message)//'check <sim_start> against water management input time'; return
    endif

    ! check if the julian day of contacenated files do not have overlap or gap if nTime_wm is larger than 1
    if (nTime_wm > 1) then
      allocate(roJulday_diff_wm(nTime_wm-1), stat=ierr)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      ! calculate the difference of consequative time in julian day
      roJulday_diff_wm = roJulday_wm (1:nTime_wm-1) - roJulday_wm (2:nTime_wm)
      ! check if the difference are identical otherwise error and terminate
      if ( any(abs(roJulday_diff_wm-roJulday_diff_wm(1)) > maxTimeDiff) ) then
          write(iulog,'(2a)') new_line('a'),'ERROR: time spacing in water management netCDF input(s) is not consistent within tolerance maxTimeDiff = ',maxTimeDiff
          ierr=20; message=trim(message)//'make sure the water management input netCDF files do not have time overlaps or gaps'; return
      end if
    endif

  endif

  ! set initial model simulation time (beginning of simulation time step)
  simDatetime(1) = begDatetime

  ! set simulation time step index (should be one to start)
  iTime = 1

  ! set time variable first simulation time step
  call begDatetime%julianday(calendar, juldaySim,  ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  t_unit = trim( time_units(1:index(time_units,' ')) )
  select case( trim(t_unit) )
    case('seconds','second','sec','s'); timeVar = (juldaySim - refJulday)*86400._dp
    case('minutes','minute','min');     timeVar = (juldaySim - refJulday)*1440._dp
    case('hours','hour','hr','h');      timeVar = (juldaySim - refJulday)*24._dp
    case('days','day','d');             timeVar = (juldaySim - refJulday)
    case default
      ierr=20; message=trim(message)//'<tunit>= '//trim(t_unit)//': <tunit> must be seconds, minutes, hours or days.'; return
  end select

  if (continue_run) then
    simDatetime(0) = simDatetime(1)%add_sec(-dt, calendar, ierr, cmessage)
  end if

  ! Set restart calendar date/time and dropoff calendar date/time and
  ! -- For periodic restart options  ---------------------------------------------------------------------
  ! Ensure that user-input restart month, day are valid.
  ! "Annual" option:  if user input day exceed number of days given user input month, set to last day
  ! "Monthly" option: use 2000-01 as template calendar yr/month
  ! "Daily" option:   use 2000-01-01 as template calendar yr/month/day
  select case(lower(trim(restart_write)))
    case('annual')
      dummyDatetime = datetime(2000, restart_month, 1, 0, 0, 0.0_dp)
      nDays = dummyDatetime%ndays_month(calendar, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage); return; endif
      if (restart_day > nDays) restart_day=nDays
    case('monthly'); restart_month = 1
    case('daily');   restart_month = 1; restart_day = 1
  end select

  select case(lower(trim(restart_write)))
    case('last')
      dropDatetime = endDatetime
      restart_month = dropDatetime%month(); restart_day = dropDatetime%day(); restart_hour = dropDatetime%hour()
    case('specified')
      if (trim(restart_date) == charMissing) then
        ierr=20; message=trim(message)//'<restart_date> must be provided when <restart_write> option is "specified"'; return
      end if
      call restDatetime%str2datetime(restart_date, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [restart_date]'; return; endif
      dropDatetime = restDatetime%add_sec(-dt, calendar, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [restDatetime->dropDatetime]'; return; endif
      restart_month = dropDatetime%month(); restart_day = dropDatetime%day(); restart_hour = dropDatetime%hour()
    case('yearly','monthly','daily')
      restDatetime = datetime(2000, restart_month, restart_day, restart_hour, 0, 0._dp)
      dropDatetime = restDatetime%add_sec(-dt, calendar, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [ dropDatetime for periodical restart]'; return; endif
      restart_month = dropDatetime%month(); restart_day = dropDatetime%day(); restart_hour = dropDatetime%hour()
    case('never')
      dropDatetime = datetime(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
    case default
      ierr=20; message=trim(message)//'Accepted <restart_write> options (case insensitive): last, never, specified, yearly, monthly, or daily '; return
  end select

 END SUBROUTINE init_time

 ! ********************************************************************************
 ! private subroutine: initialize runoff, and runoff-mapping data - For stand-alone
 ! ********************************************************************************
 SUBROUTINE init_forc_data(ierr, message)   ! output: error control

   USE public_var,  ONLY: ancil_dir            ! name of the ancillary directory
   USE public_var,  ONLY: input_dir            ! name of the runoff input directory
   USE public_var,  ONLY: vname_qsim           ! name of simulated runoff varibale
   USE public_var,  ONLY: vname_evapo          ! name of simulated evaporation varibale
   USE public_var,  ONLY: vname_precip         ! name of simulated precipitation varibale
   USE public_var,  ONLY: vname_hruid          ! name of name of varibale hruid
   USE public_var,  ONLY: dname_hruid          ! name of dimension for varibale hruid
   USE public_var,  ONLY: dname_xlon           ! name of dimension for lon
   USE public_var,  ONLY: dname_ylat           ! name of dimension for lat
   USE public_var,  ONLY: vname_flux_wm        ! name of varibale abstraction/injection
   USE public_var,  ONLY: vname_vol_wm         ! name of varibale target volume
   USE public_var,  ONLY: vname_segid_wm       ! name of varibale river network hruid for abs/inj
   USE public_var,  ONLY: dname_segid_wm       ! name of dimension hruid
   USE public_var,  ONLY: fname_remap          ! name of runoff mapping netCDF name
   USE public_var,  ONLY: is_remap             ! logical whether or not runnoff needs to be mapped to river network HRU
   USE public_var,  ONLY: is_lake_sim          ! logical if lakes simulations are activated
   USE public_var,  ONLY: is_flux_wm           ! logical whether or not abstraction or injection should be read
   USE public_var,  ONLY: is_vol_wm            ! logical whether or not target volume should be read
   USE globalData,  ONLY: nHRU                 ! number of HRUs over the entire domain
   USE globalData,  ONLY: nRch                 ! number of reaches over the entire domain
   USE globalData,  ONLY: basinID              ! basin ID
   USE globalData,  ONLY: reachID              ! reach ID
   USE globalData,  ONLY: inFileInfo_ro        ! metadata of the ro/evap/p input files
   USE globalData,  ONLY: inFileInfo_wm        ! metadata of the input files for abstration, injection and target volume
   USE globalData,  ONLY: remap_data           ! runoff mapping data structure
   USE globalData,  ONLY: runoff_data          ! runoff data structure
   USE globalData,  ONLY: wm_data              ! abstraction injection data structure
   USE read_runoff, ONLY: read_runoff_metadata ! read meta data from runoff data
   USE read_remap,  ONLY: get_remap_data       ! read remap data

   implicit none
   ! Argument variables
   integer(i4b), intent(out)          :: ierr             ! error code
   character(*), intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)              :: fname            ! input file name
   integer(i4b), allocatable          :: unq_qhru_id(:)
   integer(i4b), allocatable          :: unq_idx(:)
   character(len=strLen)              :: cmessage         ! error message from subroutine

   ierr=0; message='init_forc_data/'

   ! passing the first nc file as global file name to read
   fname = trim(input_dir)//trim(inFileInfo_ro(1)%infilename)

   ! get runoff metadata for simulated runoff, evaporation and precipitation
   call read_runoff_metadata(fname,                           & ! input: filename
                             vname_qsim,                      & ! input: varibale name for simulated runoff
                             vname_hruid,                     & ! input: varibale hruid
                             dname_hruid,                     & ! input: dimension of varibale hru
                             dname_ylat,                      & ! input: dimension of lat
                             dname_xlon,                      & ! input: dimension of lon
                             runoff_data%nSpace,              & ! nSpace of the input in runoff or wm strcuture
                             runoff_data%sim,                 & ! 1D simulation
                             runoff_data%sim2D,               & ! 2D simulation
                             runoff_data%hru_id,              & ! ID of seg or hru in data
                             runoff_data%fillvalue,           & ! fillvalue for data
                             ierr, cmessage)                    ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! initialize routing catchment array (runoff_data%basinRunoff)
   allocate(runoff_data%basinRunoff(nHRU), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   runoff_data%basinRunoff(:) = realMissing

   if (is_lake_sim) then
     allocate(runoff_data%basinEvapo(nHRU), stat=ierr, errmsg=cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     runoff_data%basinEvapo(:) = realMissing

     allocate(runoff_data%basinPrecip(nHRU), stat=ierr, errmsg=cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     runoff_data%basinPrecip(:) = realMissing
   end if

   ! need to remap runoff to HRUs
   if (is_remap) then

     ! get runoff mapping file
     call get_remap_data(trim(ancil_dir)//trim(fname_remap),  & ! input: file name
                         runoff_data%nSpace,                  & ! input: number of spatial elements
                         remap_data,                          & ! output: data structure to remap data from a polygon
                         ierr, cmessage)                        ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! get indices of the HRU ids in the mapping file in the routing layer
     call get_qix(remap_data%hru_id, &  ! input: vector of ids in mapping file
                  basinID,           &  ! input: vector of ids in the routing layer
                  remap_data%hru_ix, &  ! output: indices of hru ids in routing layer
                  ierr, cmessage)       ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     if (debug) then
       write(iulog,'(2a)') new_line('a'), 'DEBUG: Corresponding between River-Network(RN) hru in mapping data and RN hru in river network data'
       write(iulog,'(2x,a,I15)') '(1) number of RN hru in river-network = ', size(basinID)
       write(iulog,'(2x,a,I15)') '(2) number of RN hru in mapping       = ', size(remap_data%hru_id)
       write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two  = ', count(remap_data%hru_ix/=integerMissing)
       if(count(remap_data%hru_ix/=integerMissing)/=size(basinID))then
         message=trim(message)//'(1) not equal (2)'
         ierr=20; return
       endif
     end if

     if ( runoff_data%nSpace(2) == integerMissing ) then
       ! get indices of the "overlap HRUs" (the runoff input) in the runoff vector
       call get_qix(remap_data%qhru_id, &  ! input: vector of ids in mapping file
                    runoff_data%hru_id, &  ! input: vector of ids in runoff file
                    remap_data%qhru_ix, &  ! output: indices of mapping ids in runoff file
                    ierr, cmessage)           ! output: error control
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

       if (debug) then
         call unique(remap_data%qhru_id, unq_qhru_id, unq_idx)
         write(iulog,'(2a)') new_line('a'),'DEBUG: corresponding between Hydro-Model (HM) hru in mapping data and HM hru in runoff data'
         write(iulog,'(2x,a,I15)') '(1) number of HM hru in hyrdo-model  = ', size(runoff_data%hru_id)
         write(iulog,'(2x,a,I15)') '(2) number of HM hru in mapping      = ', size(unq_qhru_id)
         write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two = ', count(remap_data%qhru_ix(unq_idx)/=integerMissing)
       end if
     end if

   else ! if runoff given in RN_HRU

     allocate(runoff_data%hru_ix(size(runoff_data%hru_id)), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data%hru_ix'; return; endif

     ! get indices of the HRU ids in the runoff file in the routing layer
     call get_qix(runoff_data%hru_id,  &    ! input: vector of ids in mapping file
                  basinID,             &    ! input: vector of ids in the routing layer
                  runoff_data%hru_ix,  &    ! output: indices of hru ids in routing layer
                  ierr, cmessage)           ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     if (debug) then
       write(iulog,'(2a)') new_line('a'), 'DEBUG: corresponding between River-Network (RN) hru in runoff data and RN hru in river network data'
       write(iulog,'(2x,a,I15)') '(1) number of RN hru in river-network = ', size(basinID)
       write(iulog,'(2x,a,I15)') '(2) number of RN hru in hyrdo-model   = ', size(runoff_data%hru_id)
       write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two  = ', count(runoff_data%hru_ix/=integerMissing)
       if(count(runoff_data%hru_ix/=integerMissing)/=size(basinID))then
         message=trim(message)//'(1) not equal (2)'
         ierr=20; return
       endif
     end if

   endif

   ! Optionals: lake simulation with target volume or water-management
   if ((is_flux_wm).or.((is_vol_wm).and.(is_lake_sim))) then

     ! passing the first netCDF as global file meta to read
     fname = trim(input_dir)//trim(inFileInfo_wm(1)%infilename)

     call read_runoff_metadata(fname,                         & ! input: filename
                               vname_flux_wm,                 & ! input: varibale name for simulated runoff
                               vname_segid_wm,                & ! input: varibale hruid
                               dname_segid_wm,                & ! input: dimension of varibale hru
                               dname_ylat,                    & ! input: dimension of lat
                               dname_xlon,                    & ! input: dimension of lon
                               wm_data%nSpace,                & ! nSpace of the input in runoff or wm strcuture
                               wm_data%sim,                   & ! 1D simulation
                               wm_data%sim2D,                 & ! 2D simulation
                               wm_data%seg_id,                & ! ID of seg in data
                               wm_data%fillvalue,             & ! fillvalue for data
                               ierr, cmessage)                  ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     if (is_flux_wm) then
       allocate(wm_data%flux_wm(nRch), stat=ierr, errmsg=cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
       wm_data%flux_wm(:) = realMissing
     end if

     if (is_vol_wm) then
       allocate(wm_data%vol_wm(nRch), stat=ierr, errmsg=cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
       wm_data%vol_wm(:) = realMissing
     end if

     ! allocate the hru_ix based on number of hru_id presented in the
     allocate(wm_data%seg_ix(size(wm_data%seg_id)), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating wm_data%hru_ix'; return; endif

     ! get indices of the seg ids in the input file in the routing layer
     call get_qix(wm_data%seg_id,  &    ! input: vector of ids in mapping file
                  reachID,         &    ! input: vector of ids in the routing layer
                  wm_data%seg_ix,  &    ! output: indices of hru ids in routing layer
                  ierr, cmessage)       ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   endif

 END SUBROUTINE init_forc_data

 ! *****
 ! private subroutine: get indices of mapping points within runoff file... (To be removed since the same routine in nr_utils)
 ! ***********************************************************************
 SUBROUTINE get_qix(qid,qidMaster,qix,ierr,message)

 implicit none
 ! input
 integer(i4b), intent(in)  :: qid(:)                       ! ID of input vector
 integer(i4b), intent(in)  :: qidMaster(:)                 ! ID of master vector
 ! output
 integer(i4b), intent(out) :: qix(:)                       ! index within master vector
 integer(i4b), intent(out) :: ierr                         ! error code
 character(*), intent(out) :: message                      ! error message
 ! local
 integer(i4b)             :: rankID( size(qid) )           ! rank of input vector
 integer(i4b)             :: rankMaster( size(qidMaster) ) ! rank of master vector
 integer(i4b)             :: ix,jx,ixMaster                ! array indices
 integer(i4b)             :: nx                            ! counter

 ierr=0; message='get_qix/'

 ! sort the data vector from smallest to largest
 call indexx(qid,       rankID)
 call indexx(qidMaster, rankMaster)

 qix(1:size(qid)) = integerMissing
 nx=0
 jx=1
 ! loop through id vector
 do ix=1,size(qid)

  ! find match
  do ixMaster=jx,size(qidMaster) ! normally a very short loop

   ! keep track of trials
   nx=nx+1
   !write(*,*) 'qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) ) = ', qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) )

   ! find match
   if( qid( rankId(ix) ) == qidMaster( rankMaster(ixMaster) ) )then
    qix( rankId(ix) ) = rankMaster(ixMaster)
    jx = ixMaster
    exit
   endif

   ! unable to find match
   if( qidMaster( rankMaster(ixMaster) ) > qid( rankId(ix) ) )then
    qix( rankId(ix) ) = integerMissing
    jx = ixMaster
    exit
   endif

  end do  ! ixMaster

  ! print progress
  if(qix( rankId(ix) )/=integerMissing .and. mod(ix,1000000)==0)then
   write(iulog,*) trim(message)//'matching ids: ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) ) = ', &
                                                ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) )
  endif

 end do  ! looping through the vector

 ! check again
 do ix=1,size(qid)
  if(qix(ix) /= integerMissing)then
   if(qid(ix) /= qidMaster( qix(ix) ) )then
    write(iulog,'(a,2(x,I10,x,I15))') 'ERROR Mapping: ix, qid(ix), qix(ix), qidMaster(qix(ix))=', ix, qid(ix), qix(ix), qidMaster(qix(ix))
    message=trim(message)//'unable to find the match'
    ierr=20; return
   endif
  endif
 end do

 END SUBROUTINE get_qix

END MODULE model_setup
