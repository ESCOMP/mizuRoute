MODULE model_setup
USE nrtype
USE public_var,        ONLY: iulog
USE public_var,        ONLY: debug
USE public_var,        ONLY: integerMissing
USE public_var,        ONLY: realMissing
USE public_var,        ONLY: charMissing
USE nr_utils,          ONLY: match_index
USE nr_utils,          ONLY: unique         ! get unique element array

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
  ! Argument variables:  None
  ! local variables
  character(len=strLen)       :: message             ! error message

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
   USE globalData,          ONLY: version             ! mizuRoute version
   USE globalData,          ONLY: gitBranch           ! git branch
   USE globalData,          ONLY: gitHash             ! git commit hash
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

  ! Get mizuRoute model information
#if defined(VERSION)
  version=VERSION
#endif
#if defined(HASH)
  gitHash=HASH
#endif
#if defined(BRANCH)
  gitBranch=BRANCH
#endif

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

   ! restart initialization
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
 SUBROUTINE init_inFile_pop(ierr, message)

  USE public_var, ONLY: input_dir               ! directory containing the text files of fname_qsim and fname_wm
  USE public_var, ONLY: fname_qsim              ! runoff netcdf file name. multiple files can be expressed with wildcard
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
  USE public_var, ONLY: ro_time_units           ! time units used in forcing input netcdfs
  USE public_var, ONLY: ro_calendar             ! calendar used in forcing input netcdfs

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
                  ro_calendar,       & ! input: name of calendar used in runoff netcdf
                  ro_time_units,     & ! input: time units used in runoff netcdf
                  inFileInfo_ro,     & ! output: input file information
                  ierr, cmessage)      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

  if ((is_flux_wm).or.(is_vol_wm.and.is_lake_sim)) then     ! if either of abstraction injection or target volume is activated
    call inFile_pop(input_dir,            & ! input: name of the directory of the txt file
                    fname_wm,             & ! input: name of the txt file hold the nc file names
                    vname_time_wm,        & ! input: name of variable time in the nc files
                    dname_time_wm,        & ! input: name of dimention time in the nc files
                    ro_calendar,          & ! input: name of calendar used in water-management netcdf (assume the same as runoff)
                    ro_time_units,        & ! input: time units used in water-management netcdf (assume the same as runoff)
                    inFileInfo_wm,        & ! output: input file information
                    ierr, cmessage)         ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if
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
                       calendar_in,      & ! input: name of calendar used in input file
                       time_units_in,    & ! input: time units used in input file
                       inputFileInfo,    & ! output: input file information
                       ierr, message)      ! output: error control

  USE globalData,          ONLY: masterproc
  USE dataTypes,           ONLY: inFileInfo     ! the data type for storing the infromation of the nc files and its attributes
  USE datetime_data,       ONLY: datetime       ! datetime data
  USE ascii_utils,         ONLY: file_open      ! open file (performs a few checks as well)
  USE ascii_utils,         ONLY: get_vlines     ! get a list of character strings from non-comment lines
  USE ncio_utils,          ONLY: is_netcdf_file ! check if a file is netcdf
  USE ncio_utils,          ONLY: get_nc         ! Read netCDF variable data
  USE ncio_utils,          ONLY: check_attr     ! check if the attribute exist for a variable
  USE ncio_utils,          ONLY: get_var_attr   ! Read attributes variables
  USE ncio_utils,          ONLY: get_nc_dim_len ! get the nc dimension length
  USE public_var,          ONLY: secprmin,  &   ! time conversion factor (min->sec)
                                 secprhour, &   ! time conversion factor (hour->sec)
                                 secprday       ! time conversion factor (day->sec)
  USE mpi_utils,           ONLY: shr_mpi_bcast

  ! Argument variables
  character(len=strLen), intent(in)                 :: dir_name         ! the name of the directory that the txt file located
  character(len=strLen), intent(in)                 :: file_name        ! the name of the file that include the nc file names
  character(len=strLen), intent(in)                 :: time_var_name    ! the name of the time variable
  character(len=strLen), intent(in)                 :: time_dim_name    ! the name of dimension time
  character(len=strLen), intent(in)                 :: calendar_in      ! name of calendar used in input file
  character(len=strLen), intent(in)                 :: time_units_in    ! time units used in input file
  type(infileinfo),      intent(inout), allocatable :: inputFileInfo(:) ! the name of structure that hold the infile information
  integer(i4b),          intent(out)                :: ierr             ! error code
  character(*),          intent(out)                :: message          ! error message
  ! local varibales
  integer(i4b)                                      :: funit            ! file unit (free unit output from file_open)
  character(len=strLen)                             :: tmp_file_list    ! temporal text listing all the input netcdf(s)
  character(len=7)                                  :: t_unit           ! time units. "<time_step> since yyyy-MM-dd hh:mm:ss"
  integer(i4b)                                      :: iFile            ! counter for forcing files
  integer(i4b)                                      :: nFile            ! number of nc files identified in the text file
  integer(i4b)                                      :: nTime            ! hard coded for now
  logical(lgt)                                      :: existAttr        ! attribute exists or not
  logical(lgt)                                      :: tmp_file_exists  ! tmp file exists or not
  logical(lgt)                                      :: is_nc            ! input file is netcdf and not ascii
  real(dp)                                          :: convTime2sec     ! time conversion to second
  character(len=strLen)                             :: infilename       ! input filename
  character(len=strLen),allocatable                 :: dataLines(:)     ! vector of lines of information (non-comment lines)
  character(len=strLen)                             :: cmessage         ! error message of downwind routine
  character(len=strLen)                             :: trim_file_name   ! temporal text keeping the trimmed file name

  ierr=0; message='inFile_pop/'

  ! build filename and its path containing list of NetCDF files
  ! then construct a character array including the file paths
  trim_file_name = trim(file_name)
  infilename = trim(dir_name)//trim_file_name
  tmp_file_list = trim(dir_name)//'tmp'
  is_nc = .fasle.
  ! remove possible tmp file
  call execute_command_line("rm -f "//trim(tmp_file_list))
  
  if (masterproc) then
    
    ! check if the infile is a wild card with * or ? such as file*.nc4 or file?.nc
    DO i = 1, len(trim_file_name)
      IF (trim_file_name(i:i) == '*' .OR. trim_file_name(i:i) == '?') THEN
        ! create the tmp_file_list on disk
        call execute_command_line("ls "//infilename//" > "//trim(tmp_file_list))
        EXIT
      END IF
    END DO

    ! check if tmp is created
    INQUIRE(FILE=tmp_file_list, EXIST=tmp_file_exists)
    
    ! is tmp file is not created then it should be file.nc or file_name.txt 
    if (.NOT. tmp_file_exists) then
      ! check the file is netcdf
      call is_netcdf_file (infilename, is_nc, ierr, message)
      ! check opening is successful
      if (is_nc) then
        call execute_command_line("ls "//infilename//" > "//trim(tmp_file_list))
      else
        call execute_command_line("cp "//infilename//" > "//trim(tmp_file_list))
      end if
    end if
    
    ! open the tmp file
    call file_open(tmp_file_list,funit,ierr,cmessage) ! open the text file
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! get a list of character strings from non-commented lines
    call get_vlines(funit,dataLines,ierr,cmessage)
    if(ierr/=0)then; ierr=20; message=trim(message)//trim(cmessage); return; end if
    ! remove the tmp file
    call execute_command_line("rm -f "//trim(tmp_file_list))

  end if

  call shr_mpi_bcast(dataLines, ierr, cmessage)
  if(ierr/=0)then; ierr=20; message=trim(message)//trim(cmessage); return; end if

  nFile = size(dataLines) ! get the name of the lines in the file

  ! allocate space for forcing information
  allocate(inputFileInfo(nFile), stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for forcFileInfo'; return; end if

  ! poputate the forcingInfo structure with filenames, and time variables/attributes
  do iFile=1,nFile

    ! set forcing file name
    inputFileInfo(iFile)%infilename = dataLines(iFile)

    ! get the time units. if not exsit in netcdfs, provided from the control file
    existAttr = check_attr(trim(inputFileInfo(iFile)%infilename), time_var_name, 'units')
    if (existAttr) then
      call get_var_attr(trim(inputFileInfo(iFile)%infilename), &
                        time_var_name, 'units', inputFileInfo(iFile)%unit, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    else
      if (trim(time_units_in)==charMissing) then
        write(cmessage, '(2A)')  trim(time_var_name), '. No units attribute exist nor provided by user in ro_time_units in control file'
        ierr=10; message=trim(message)//trim(cmessage); return
      end if
      inputFileInfo(iFile)%unit = time_units_in
    end if

    ! get the calendar. if not exsit in netcdfs, provided from the control file
    existAttr = check_attr(trim(inputFileInfo(iFile)%infilename), time_var_name, 'calendar')
    if (existAttr) then
      call get_var_attr(trim(inputFileInfo(iFile)%infilename), &
                        time_var_name, 'calendar', inputFileInfo(iFile)%calendar, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    else
      if (trim(calendar_in)==charMissing) then
        write(cmessage, '(2A)')  trim(time_var_name), '. No calendar attribute exist nor provided by user in ro_calendar in control file'
        ierr=10; message=trim(message)//trim(cmessage); return
      end if
      inputFileInfo(iFile)%calendar = calendar_in
    end if

    ! get the dimension of the time to populate nTime and pass it to the get_nc file
    call get_nc_dim_len(trim(inputFileInfo(iFile)%infilename), &
                        time_dim_name, inputFileInfo(iFile)%nTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    nTime = inputFileInfo(iFile)%nTime ! the length of time varibale for each nc file

    ! allocate space for time varibale of each file
    allocate(inputFileInfo(iFile)%timeVar(nTime))
    if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for inputFileInfo(:)%timeVar'; return; end if

    ! get the time varibale
    call get_nc(trim(inputFileInfo(iFile)%infilename), &
                time_var_name, inputFileInfo(iFile)%timeVar, 1, nTime, ierr, cmessage) ! does it needs timeVar(:)

    ! get the time multiplier needed to convert time to units of days for each nc file
    t_unit = trim( inputFileInfo(iFile)%unit(1:index(inputFileInfo(iFile)%unit,' ')) )
    select case( trim(t_unit) )
      case('seconds','second','sec','s'); convTime2sec=1._dp
      case('minutes','minute','min','m'); convTime2sec=secprmin
      case('hours'  ,'hour'  ,'hr' ,'h'); convTime2sec=secprhour
      case('days'   ,'day'   ,'d');       convTime2sec=secprday
      case default
        ierr=20; message=trim(message)//'<time_units>= '//trim(t_unit)//': <time_units> must be seconds, minutes, hours or days.'; return
    end select
    ! convert timeValue unit to second
    inputFileInfo(iFile)%timeVar(:) = inputFileInfo(iFile)%timeVar(:)*convTime2sec

    ! get the reference datetime from the nc file
    call inputFileInfo(iFile)%refDatetime%str2datetime(trim(inputFileInfo(iFile)%unit), inputFileInfo(iFile)%calendar, ierr, message)
    if(ierr/=0) then; message=trim(message)//trim(cmessage)//' inputFileInfo%refDatetime%str2datetime'; return; endif

    ! populated the index of the iTimebound for each nc file
    if (iFile==1) then
      inputFileInfo(iFile)%iTimebound(1) = 1
      inputFileInfo(iFile)%iTimebound(2) = nTime
    else ! if multiple files specfied in the txt file
      inputFileInfo(iFile)%iTimebound(1) = inputFileInfo(iFile-1)%iTimebound(2) + 1 ! the last index from the perivous nc file + 1
      inputFileInfo(iFile)%iTimebound(2) = inputFileInfo(iFile-1)%iTimebound(2) + nTime ! the last index from the perivous nc file + 1
    endif

  end do

  close(unit=funit,iostat=ierr) ! close ascii file
  if(ierr/=0)then;message=trim(message)//'problem closing forcing file list'; return; end if

 END SUBROUTINE inFile_pop

 ! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 SUBROUTINE init_time(ierr, message)

  ! purpose: Save the following time related global variables
  ! - dt_ro: forcing time step [sec]
  ! - dt_wm: water-management data time step [sec]
  ! - time_units: time units used for simulation datetime. format="<units> since yyyy-mm-dd hh:mm:ss"
  ! - calendar: calendar used for simulation datetime
  ! - timeVar: time variable [sec] for simulation (since referece datetime
  ! - iTime: time index of simulation time step
  ! - begDatetime: datetime at front of 1st simulation time step period
  ! - endDatetime: datetime at front of last simulation time step period
  ! - restDatetime: datetime at front of restart simulation time step period
  ! - dropDatetime: datetime at front of simulation time step period when restart file is written

  USE ascii_utils,   ONLY: lower                    ! convert string to lower case
  USE datetime_data, ONLY: datetime                 ! datetime data
  USE public_var,    ONLY: time_units               ! time units (seconds, hours, or days)
  USE public_var,    ONLY: simStart                 ! date string defining the start of the simulation
  USE public_var,    ONLY: simEnd                   ! date string defining the end of the simulation
  USE public_var,    ONLY: calendar                 ! calendar used for simulation
  USE public_var,    ONLY: dt                       ! simulation time step in second
  USE public_var,    ONLY: dt_ro                    ! forcing time step in second
  USE public_var,    ONLY: dt_wm                    ! water-management time step in second
  USE public_var,    ONLY: continue_run             ! logical to indicate sppend output in existing history file
  USE public_var,    ONLY: secprday                 ! unit conversion from day to sec
  USE public_var,    ONLY: secprhour                ! unit conversion from hour to sec
  USE public_var,    ONLY: secprmin                 ! unit conversion from minute to sec
  USE public_var,    ONLY: restart_write            ! restart write option
  USE public_var,    ONLY: restart_date             ! restart datetime
  USE public_var,    ONLY: restart_month            ! periodic restart month
  USE public_var,    ONLY: restart_day              ! periodic restart day
  USE public_var,    ONLY: restart_hour             ! periodic restart hr
  USE public_var,    ONLY: maxTimeDiff              ! time difference tolerance for input checks
  USE public_var,    ONLY: ro_time_stamp            ! time stamp for runoff input: front, end or middle of time step
  USE globalData,    ONLY: timeVar                  ! time variables at time step endpoints (unit given by runoff data)
  USE globalData,    ONLY: iTime                    ! time index at simulation time step
  USE globalData,    ONLY: sec2tunit                ! seconds per time unit
  USE globalData,    ONLY: simDatetime              ! model time data (yyyy:mm:dd:hh:mm:ss)
  USE globalData,    ONLY: begDatetime              ! simulation begin datetime data (yyyy:mm:dd:hh:mm:sec). front of time step
  USE globalData,    ONLY: endDatetime              ! simulation end datetime data (yyyy:mm:dd:hh:mm:sec). front of time step
  USE globalData,    ONLY: restDatetime             ! restart time data (yyyy:mm:dd:hh:mm:sec). front of time step
  USE globalData,    ONLY: dropDatetime             ! restart dropoff calendar date/time. front of time step
  USE globalData,    ONLY: roBegDatetime            ! forcing data start datetime data (yyyy:mm:dd:hh:mm:sec) front of time step
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
  character(len=strLen)                    :: t_unit              ! time units. "<time_step> since yyyy-MM-dd hh:mm:ss"
  integer(i4b)                             :: iFile               ! for loop over the nc files
  type(datetime), allocatable              :: roCal(:)            ! datetime in runoff data
  type(datetime), allocatable              :: wmCal(:)            ! datetime in water-management data
  type(datetime)                           :: roDatetime_end      ! end of datetime in runoff data
  type(datetime)                           :: refDatetime         ! reference datetime from unit time based on "time_units"
  type(datetime)                           :: dummyDatetime       ! temp datetime
  integer(i4b)                             :: nDays               ! number of days in a month
  real(dp), allocatable                    :: roTimeVar(:)        ! elapsed seconds from reference datetime in runoff data
  real(dp), allocatable                    :: roTimeVar_diff(:)   ! time difference in second between two concequative runoff time step
  real(dp), allocatable                    :: wmTimeVar(:)        ! elapsed seconds from reference datetime in water management
  real(dp), allocatable                    :: wmTimeVar_diff(:)   ! time difference in second between two concequative water-management time step
  character(len=strLen)                    :: cmessage            ! error message of downwind routine
  character(len=50)                        :: fmt1='(a,I4,a,I2.2,a,I2.2,x,I2.2,a,I2.2,a,F5.2)'

  ierr=0; message='init_time/'

  ! Set simulation time attributes-calendar and time units (used for history file)
  ! For time units, if they are not provided in control files, Use from 1st runoff input netCDF
  ! For calendar, always use the same as forcing input
  calendar = inFileInfo_ro(1)%calendar
  if (trim(time_units)==charMissing) then
    time_units = inFileInfo_ro(1)%unit
  end if
  call refDatetime%str2datetime(time_units, calendar, ierr, message)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refDatetime]'; return; endif

  ! get the number of the total time length of all the nc files
  nFile = size(inFileInfo_ro)
  nTime = sum(inFileInfo_ro(:)%nTime)

  ! runoff data: time variable-roTimeVar [sec] from reference datetime (refDatetime) and datetime (roCal)
  allocate(roTimeVar(nTime), roCal(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  counter = 1;
  do iFile=1,nFile
    nt = inFileInfo_ro(iFile)%nTime
    do ix=1, nt
      roCal(counter) = inFileInfo_ro(iFile)%refDatetime%add_sec(inFileInfo_ro(iFile)%timeVar(ix), ierr, cmessage)
      roTimeVar(counter) = roCal(counter) - refDatetime ! sec
      counter = counter + 1
    end do
  end do

  ! check if time step in runoff data is consistent
  if (nTime > 1) then
    allocate(roTimeVar_diff(nTime-1), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    ! calculate the difference of consequative time in julian day
    roTimeVar_diff = roTimeVar(2:nTime)-roTimeVar(1:nTime-1)
    ! check if the difference are identical otherwise error and terminate
    if ( any(abs(roTimeVar_diff-roTimeVar_diff(1)) > maxTimeDiff) ) then
      write(iulog,'(2a)') new_line('a'),'ERROR: time spacing in netCDF input(s) is not consistent within tolerance maxTimeDiff = ',maxTimeDiff
      ierr=20; message=trim(message)//'make sure the input netCDF files do not have time overlaps or gaps'; return
    end if
  endif

  ! runoff data time step [sec]- dt_ro is saved
  dt_ro = roTimeVar_diff(1)

  ! datetime of runoff time at front of the 1st time step:  roBegDatetime (global data)
  ! datetime of runoff time at end of last time step: roDatetime_end (local data)
  select case(trim(ro_time_stamp))
    case('start')
      roBegDatetime  = roCal(1)
      roDatetime_end = roCal(nTime)%add_sec(dt_ro, ierr, cmessage)
    case('end')
      roBegDatetime  = roCal(1)%add_sec(-dt_ro, ierr, cmessage)
      roDatetime_end = roCal(nTime)
    case('middle')
      roBegDatetime  = roCal(1)%add_sec(-dt_ro/2.0, ierr, cmessage)
      roDatetime_end = roCal(nTime)%add_sec(dt_ro/2.0, ierr, cmessage)
  end select

  call begDatetime%str2datetime(simStart, calendar, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [begDatetime]'; return; endif

  call endDatetime%str2datetime(simEnd, calendar, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endDatetime]'; return; endif

  ! check that the dates are aligned
  if(endDatetime < begDatetime) then
    write(cmessage,'(7a)') 'simulation end is before simulation start:', new_line('a'), '<sim_start>= ', trim(simStart), new_line('a'), '<sim_end>= ', trim(simEnd)
    ierr=20; message=trim(message)//trim(cmessage); return
  endif

  ! check sim_start is after the last time step in runoff data
  if (begDatetime > roDatetime_end) then
    write(iulog,'(2a)') new_line('a'),'ERROR: <sim_start> is after the last time step in input runoff'
    write(iulog,fmt1)  ' runoff_end  : ', roDatetime_end%year(),'-',roDatetime_end%month(),'-',roDatetime_end%day(), roDatetime_end%hour(),':', roDatetime_end%minute(),':',roDatetime_end%sec()
    write(iulog,fmt1)  ' <sim_start> : ', begDatetime%year(),'-',begDatetime%month(),'-',begDatetime%day(), begDatetime%hour(),':', begDatetime%minute(),':',begDatetime%sec()
    ierr=20; message=trim(message)//'check <sim_start> against runoff input time'; return
  endif

  ! Compare sim_start vs. time at first time step in runoff data
  if (begDatetime < roBegDatetime) then
    write(iulog,'(2a)') new_line('a'),'WARNING: <sim_start> is before the first time step in input runoff'
    write(iulog,fmt1)  ' runoff_start: ', roCal(1)%year(),'-',roCal(1)%month(),'-',roCal(1)%day(), roCal(1)%hour(),':', roCal(1)%minute(),':',roCal(1)%sec()
    write(iulog,fmt1)  ' <sim_start> : ', begDatetime%year(),'-',begDatetime%month(),'-',begDatetime%day(), begDatetime%hour(),':', begDatetime%minute(),':',begDatetime%sec()
    write(iulog,'(a)') ' Reset <sim_start> to runoff_start'
    begDatetime = roCal(1)
  endif

  ! Compare sim_end vs. time at last time step in runoff data
  if (endDatetime > roDatetime_end) then
    write(iulog,'(2a)')  new_line('a'),'WARNING: <sim_end> is after the last time step in input runoff'
    write(iulog,fmt1)   ' runoff_end: ', roDatetime_end%year(),'-',roDatetime_end%month(),'-',roDatetime_end%day(), roDatetime_end%hour(),':', roDatetime_end%minute(),':',roDatetime_end%sec()
    write(iulog,fmt1)   ' <sim_end> : ', endDatetime%year(),'-',endDatetime%month(),'-',endDatetime%day(), endDatetime%hour(),':', endDatetime%minute(),':',endDatetime%sec()
    write(iulog,'(a)')  ' Reset <sim_end> to runoff_end'
    endDatetime = roCal(nTime)
  endif

  ! water management options on
  if ((is_flux_wm).or.(is_vol_wm.and.is_lake_sim)) then

    ! get the number of the total time length of all the water management nc files
    nFile_wm = size(inFileInfo_wm)
    nTime_wm = sum(inFileInfo_wm(:)%nTime)

    ! elapsed time-wmTimeVar [sec] from reference datetime (refDatetime) and datetime (wmCal)
    allocate(wmTimeVar(nTime_wm), wmCal(nTime_wm), stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    counter = 1;
    do iFile=1,nFile_wm
      nt = inFileInfo_wm(iFile)%nTime
      do ix=1, nt
        wmCal(counter) = inFileInfo_wm(iFile)%refDatetime%add_sec(inFileInfo_wm(iFile)%timeVar(ix), ierr, cmessage)
        wmTimeVar(counter) = wmCal(counter) - refDatetime ! sec
        counter = counter + 1
      end do
    enddo

    ! check if time step in water-management data is consistent
    if (nTime_wm > 1) then
      allocate(wmTimeVar_diff(nTime_wm-1), stat=ierr)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      ! calculate the difference of consequative time in julian day
      wmTimeVar_diff = wmTimeVar(2:nTime_wm) - wmTimeVar(1:nTime_wm-1)
      ! check if the difference are identical otherwise error and terminate
      if ( any(abs(wmTimeVar_diff-wmTimeVar_diff(1)) > maxTimeDiff) ) then
        write(iulog,'(2a)') new_line('a'),'ERROR: time spacing in water management netCDF input(s) is not consistent within tolerance maxTimeDiff = ',maxTimeDiff
        ierr=20; message=trim(message)//'make sure the water management input netCDF files do not have time overlaps or gaps'; return
      end if
    endif

    ! water-management data time step [sec]
    dt_wm = wmTimeVar_diff(1)

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

  endif ! end of water manamgemnt option

  ! set initial model simulation time (beginning of simulation time step and next time step)
  simDatetime(1) = begDatetime
  simDatetime(2) = simDatetime(1)%add_sec(dt, ierr, cmessage)
  if (continue_run) then
    simDatetime(0) = simDatetime(1)%add_sec(-dt, ierr, cmessage)
  end if

  ! set simulation time step index (should be one to start)
  iTime = 1

  ! unit conversion from second to time unit used in the simulation
  t_unit = trim( time_units(1:index(time_units,' ')) )
  select case( trim(t_unit) )
    case('seconds','second','sec','s'); sec2tunit=1._dp
    case('minutes','minute','min');     sec2tunit=secprmin
    case('hours','hour','hr','h');      sec2tunit=secprhour
    case('days','day','d');             sec2tunit=secprday
    case default
      ierr=20; message=trim(message)//'<tunit>= '//trim(t_unit)//': <tunit> must be seconds, minutes, hours or days.'; return
  end select

  ! Initialize timeVar : model time (time step endpoints) in model time unit (t_unit), used for time output
  timeVar(1) = begDatetime - refDatetime  ! second since reference datetime
  timeVar(2) = timeVar(1) + dt

  ! Set restart calendar date/time and dropoff calendar date/time and
  ! -- For periodic restart options  ---------------------------------------------------------------------
  ! Ensure that user-input restart month, day are valid.
  ! "yearly" option:  if user input day exceed number of days given user input month, set to last day
  ! "monthly" option: use 2000-01 as template calendar yr/month
  ! "daily" option:   use 2000-01-01 as template calendar yr/month/day
  select case(lower(trim(restart_write)))
    case('yearly')
      dummyDatetime = datetime(2000, restart_month, 1, 0, 0, 0.0_dp, calendar=calendar)
      nDays = dummyDatetime%ndays_month()
      if(nDays==integerMissing) then; ierr=10; message=trim(message)//'dummyDatetime calendar may be invalid'; return; endif
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
      call restDatetime%str2datetime(restart_date, calendar, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [restart_date]'; return; endif
      dropDatetime = restDatetime%add_sec(-dt, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [restDatetime->dropDatetime]'; return; endif
      restart_month = dropDatetime%month(); restart_day = dropDatetime%day(); restart_hour = dropDatetime%hour()
    case('yearly','monthly','daily')
      restDatetime = datetime(2000, restart_month, restart_day, restart_hour, 0, 0._dp, calendar=calendar)
      dropDatetime = restDatetime%add_sec(-dt, ierr, cmessage)
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
   USE read_runoff, ONLY: read_forcing_metadata ! read meta data from runoff data
   USE read_remap,  ONLY: get_remap_data       ! read remap data

   implicit none
   ! Argument variables
   integer(i4b), intent(out)          :: ierr             ! error code
   character(*), intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)              :: fname            ! input file name
   integer(i8b), allocatable          :: unq_qhru_id(:)
   integer(i4b), allocatable          :: unq_idx(:)
   character(len=strLen)              :: cmessage         ! error message from subroutine

   ierr=0; message='init_forc_data/'

   ! passing the first nc file as global file name to read
   fname = trim(inFileInfo_ro(1)%infilename)

   if (trim(vname_hruid)==charMissing) then
     ! get runoff metadata for simulated runoff, evaporation and precipitation
     call read_forcing_metadata(fname,                           & ! input: filename
                                vname_qsim,                      & ! input: varibale name for simulated runoff
                                vname_hruid,                     & ! input: varibale hruid
                                dname_hruid,                     & ! input: dimension of varibale hru
                                dname_ylat,                      & ! input: dimension of lat
                                dname_xlon,                      & ! input: dimension of lon
                                runoff_data%nSpace,              & ! nSpace of the input in runoff or wm strcuture
                                runoff_data%sim2D,               & ! 2D simulation
                                runoff_data%hru_id,              & ! ID of seg or hru in data
                                runoff_data%fillvalue,           & ! fillvalue for data
                                ierr, cmessage)                    ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else if (trim(dname_ylat)==charMissing .and. trim(dname_xlon)==charMissing) then
     ! get runoff metadata for simulated runoff, evaporation and precipitation
     call read_forcing_metadata(fname,                           & ! input: filename
                                vname_qsim,                      & ! input: varibale name for simulated runoff
                                vname_hruid,                     & ! input: varibale hruid
                                dname_hruid,                     & ! input: dimension of varibale hru
                                dname_ylat,                      & ! input: dimension of lat
                                dname_xlon,                      & ! input: dimension of lon
                                runoff_data%nSpace,              & ! nSpace of the input in runoff or wm strcuture
                                runoff_data%sim,                 & ! 1D simulation
                                runoff_data%hru_id,              & ! ID of seg or hru in data
                                runoff_data%fillvalue,           & ! fillvalue for data
                                ierr, cmessage)                    ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

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
     remap_data%hru_ix = match_index(basinID, remap_data%hru_id, ierr, cmessage)
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
       remap_data%qhru_ix = match_index(runoff_data%hru_id, remap_data%qhru_id, ierr, cmessage)
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
     runoff_data%hru_ix = match_index(int(basinID,kind=i8b), runoff_data%hru_id, ierr, cmessage)
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
     fname = trim(inFileInfo_wm(1)%infilename)

     call read_forcing_metadata(fname,                         & ! input: filename
                               vname_flux_wm,                 & ! input: varibale name for simulated runoff
                               vname_segid_wm,                & ! input: varibale hruid
                               dname_segid_wm,                & ! input: dimension of varibale hru
                               dname_ylat,                    & ! input: dimension of lat
                               dname_xlon,                    & ! input: dimension of lon
                               wm_data%nSpace,                & ! nSpace of the input in runoff or wm strcuture
                               wm_data%sim,                   & ! 1D simulation
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
     wm_data%seg_ix = match_index(reachID, wm_data%seg_id, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   endif

 END SUBROUTINE init_forc_data

END MODULE model_setup
