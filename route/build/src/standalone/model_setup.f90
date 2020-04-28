MODULE model_setup

! data types
USE nrtype,    ONLY: i4b,dp,lgt          ! variable types, etc.
USE nrtype,    ONLY: strLen              ! length of characters

! Shared data
USE public_var, ONLY: iulog              ! i/o logical unit number
USE public_var, ONLY: integerMissing
USE public_var, ONLY: realMissing

! external routines
USE init_model_data,  ONLY: init_ntopo_data
USE init_model_data,  ONLY: init_state_data
USE init_model_data,  ONLY: get_mpi_omp
USE init_model_data,  ONLY: update_time

implicit none

! privacy -- everything private unless declared explicitly
private
public :: init_mpi
public :: init_data
public :: infile_name
public :: inFile_pop

CONTAINS

 ! *********************************************************************
 ! public subroutine: initialize MPI for stand-alone program
 ! *********************************************************************
 SUBROUTINE init_mpi()

  ! Initialize MPI and get OMP thread

  ! shared data used
  USE globalData, ONLY: mpicom_route ! communicator id
  ! subroutines: populate metadata
  USE mpi_mod, ONLY: shr_mpi_init

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

   USE globalData,  ONLY: isFileOpen             ! file open/close status
   ! external subroutines
   USE mpi_routine, ONLY: pass_global_data       ! mpi globaldata copy to slave proc

   IMPLICIT NONE
   integer(i4b),              intent(in)    :: pid              ! proc id
   integer(i4b),              intent(in)    :: nNodes           ! number of procs
   integer(i4b),              intent(in)    :: comm             ! communicator
   ! output: error control
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variable
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ! initialize error control
   ierr=0; message='init_data/'

   ! network topology data initialization
   call init_ntopo_data(pid, nNodes, comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! runoff and mapping data initialization
   call init_runoff_data(pid, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! broadcast public and some global variables
   call pass_global_data(comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! channel state initialization
   call init_state_data(pid, nNodes, comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   isFileOpen=.false.

 END SUBROUTINE init_data


 ! ********************************************************************************
 ! public subroutine: initialize runoff, and runoff-mapping data - For stand-alone
 ! ********************************************************************************
 SUBROUTINE init_runoff_data(pid,           & ! input: proc id
                             ierr, message)   ! output: error control

   USE public_var, ONLY: is_remap             ! logical whether or not runnoff needs to be mapped to river network HRU
   USE globalData, ONLY: remap_data           ! runoff mapping data structure
   USE globalData, ONLY: runoff_data          ! runoff data structure

   implicit none
   ! input:
   integer(i4b),              intent(in)    :: pid              ! proc id
   ! output: error control
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local:
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   if (pid==0) then
     ! runoff and remap data initialization (TO DO: split runoff and remap initialization)
     call init_runoff(is_remap,        & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                      remap_data,      & ! output: data structure to remap data
                      runoff_data,     & ! output: data structure for runoff
                      ierr, cmessage)    ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! DateTime initialization
     ! call init_time(runoff_data%ntime, ierr, cmessage)
     ! if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end if  ! if processor=0 (root)

 END SUBROUTINE init_runoff_data

 ! *********************************************************************
 ! private subroutine: read the name of the netcdf files, populates the
 ! filed of inFiledata dataType assuming the calendar, time units, e.g. days
 ! since yyyy-mm-dd hh:mm:ss, of all the files are identical
 ! *********************************************************************
 SUBROUTINE inFile_pop(ierr, message)  ! output


  ! data types
  USE dataTypes, ONLY: infileinfo
  USE io_netcdf, only:get_nc
  USE io_netcdf, only:get_var_attr_real
  USE io_netcdf, only:get_var_attr_char
  USE io_netcdf, only:get_nc_dim_len

  ! subroutines
  USE ascii_util_module, ONLY:file_open        ! open file (performs a few checks as well)
  USE ascii_util_module, ONLY:get_vlines       ! get a list of character strings from non-comment lines

  ! Shared data
  USE public_var, ONLY: input_dir       ! directory containing input data
  USE public_var, ONLY: fname_qsim      ! simulated runoff txt file that includes the NetCDF file names
  USE public_var, ONLY: vname_time      ! variable name for time
  USE public_var, ONLY: time_units      ! time units (seconds, hours, or days)
  USE public_var, ONLY: dname_time      !
  !USE public_var, ONLY: simStart        ! date string defining the start of the simulation
  !USE public_var, ONLY: simEnd          ! date string defining the end of the simulation
  USE public_var, ONLY: calendar        ! calendar name
  USE globalData, ONLY: timeVar         ! time variables (unit given by runoff data)
  !USE globalData, ONLY: iTime           ! time index at simulation time step
  USE globalData, ONLY: convTime2Days   ! conversion multipliers for time unit of runoff input to day
  !USE globalData, ONLY: refJulday       ! julian day: reference
  !USE globalData, ONLY: startJulday     ! julian day: start of routing simulation
  !USE globalData, ONLY: endJulday       ! julian day: end of routing simulation
  !USE globalData, ONLY: modJulday       ! julian day: at model time step
  !USE globalData, ONLY: modTime         ! model time data (yyyy:mm:dd:hh:mm:ss)
  USE globalData, ONLY: infileinfo_data ! the information of the input files

  ! output: error control
  integer(i4b),         intent(out)    :: ierr             ! error code
  character(*),         intent(out)    :: message          ! error message

  ! local varibales
  !integer(i4b)                        :: size_fname_qsim       ! error code
  character(len=strLen)                :: infilename            ! input filename
  integer(i4b)                         :: unt                   ! file unit (free unit output from file_open)
  character(len=strLen),allocatable    :: dataLines(:)          ! vector of lines of information (non-comment lines)
  integer(i4b)                         :: iFile                 ! counter for forcing files
  integer(i4b)                         :: nFile                 ! number of forcing files in forcing file list
  character(len=strLen)                :: filenameData          ! name of forcing datafile
  !type(infileinfo), allocatable       :: infileinfo_data(:)    ! the file data
  integer(i4b)                         :: nTime                 ! hard coded for now
  character(len=strLen)                :: cmessage              ! error message of downwind routine
  integer(i4b)                         :: counter               ! counter
  integer(i4b)                         :: i                     ! counter

  ! initialize error control
  ierr=0; message='inFile_pop/'

  ! build filename and its path containing list of NetCDF files
  infilename = trim(input_dir)//trim(fname_qsim)

  ! open file
  call file_open(trim(infilename),unt,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! get a list of character strings from non-comment lines
  call get_vlines(unt,dataLines,ierr,cmessage)
  if(ierr/=0)then; ierr=20; message=trim(message)//trim(cmessage); return; end if
  nFile = size(dataLines) ! get the name of the lines in the file
  print*, 'number of lines in the text file', nFile

  ! allocate space for forcing information
  !if(allocated(infileinfo_data)) deallocate(infileinfo_data)
  allocate(infileinfo_data(nFile)) ! allocate(infileinfo_data(nFile), stat=ierr)
  !if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for forcFileInfo'; return; end if
  print*, "infileinfo is allocated"

  ! poputate the forcingInfo structure with filenames, julian day of sart and end of the simulation
  do iFile=1,nFile

   ! split the line into "words" (expect one word: the file describing forcing data for that index)
   read(dataLines(iFile),*,iostat=ierr) filenameData
   if(ierr/=0)then; message=trim(message)//'problem reading a line of data from file ['//trim(infilename)//']'; return; end if

   ! set forcing file name
   infileinfo_data(iFile)%infilename = trim(filenameData)
   print*, infileinfo_data(iFile)%infilename

   ! get the time units
   call get_var_attr_char(trim(input_dir)//trim(infileinfo_data(iFile)%infilename), &
                          trim(vname_time), 'units', infileinfo_data(iFile)%unit, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   print*, "unit is", infileinfo_data(iFile)%unit

   ! get the calendar
   call get_var_attr_char(trim(input_dir)//trim(infileinfo_data(iFile)%infilename), &
                          trim(vname_time), 'calendar', infileinfo_data(iFile)%calendar, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   print*, "calendar is", infileinfo_data(iFile)%calendar

   ! how to get the dimension of the time to populate nTime and pass it to the get_nc file
   call get_nc_dim_len(trim(input_dir)//trim(infileinfo_data(iFile)%infilename), &
                       trim(dname_time), infileinfo_data(iFile)%nTime, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   print*, "time dimension is ", infileinfo_data(iFile)%nTime
   nTime = infileinfo_data(iFile)%nTime
   print*, "nTime ", nTime

   ! allocate space for time varibale of each file
   if(allocated(infileinfo_data(iFile)%timeVar)) deallocate(infileinfo_data(iFile)%timeVar)
   allocate(infileinfo_data(iFile)%timeVar(nTime))
   !allocated(infileinfo_data(iFile)%timeVar(infileinfo_data(iFile)%nTime))
   print*, "time var is allocated"

   ! get the time varibale
   call get_nc(trim(input_dir)//trim(infileinfo_data(iFile)%infilename), &
               vname_time, infileinfo_data(iFile)%timeVar, 1, nTime, ierr, cmessage) ! does it needs timeVar(:)
   print*, "time var is", infileinfo_data(iFile)%timeVar

   ! get the time multiplier needed to convert time to units of days for each file
   select case( trim( infileinfo_data(iFile)%unit(1:index(infileinfo_data(iFile)%unit,' ')) ) )
    case('seconds'); convTime2Days=86400._dp
    case('hours');   convTime2Days=24._dp
    case('days');    convTime2Days=1._dp
    case default;    ierr=20; message=trim(message)//'unable to identify time units'; return
   end select
   infileinfo_data(iFile)%convTime2Days = convTime2Days
   print*, "conversion", convTime2Days

   ! populated the index of the iTime for each nc file
   if (iFile.eq.1) then
    infileinfo_data(iFile)%iTimebound(1) = 1
    infileinfo_data(iFile)%iTimebound(2) = size(infileinfo_data(iFile)%timeVar)
   endif
   if (iFile.gt.1) then
    infileinfo_data(iFile)%iTimebound(1) = infileinfo_data(iFile-1)%iTimebound(2) + 1 ! the last index from the perivous nc file + 1
    infileinfo_data(iFile)%iTimebound(2) = infileinfo_data(iFile-1)%iTimebound(2) + nTime ! the last index from the perivous nc file + 1
   endif

   print*, "time bound = ",infileinfo_data(iFile)%iTimebound

  enddo

  ! find the total length of the timeVar
  counter = 0; ! counter
  do iFile=1,nFile
   do i = 1, infileinfo_data(iFile)%nTime
    counter = counter + 1
   enddo
  enddo
  print*, "size of the time var: ", counter ! the total number of time steps in input files from first to last

  ! allocate the timeVar
  allocate(timeVar(counter)) ! what is the allocate stat?

  ! pass the time var from each file into the global time var
  counter = 1; ! counter
  do iFile=1,nFile
   timeVar(counter:counter+infileinfo_data(iFile)%nTime-1) = infileinfo_data(iFile)%timeVar
   counter = counter + +infileinfo_data(iFile)%nTime
  enddo

  print*, "global time var : ", timeVar

  ! passing the first nc file as global netCDF file
  fname_qsim = trim(infileinfo_data(1)%infilename)
  calendar = infileinfo_data(1)%calendar

  print*, "name of the nc file from the pop in file", fname_qsim


  ! call init_time_new to get the first iTime
  call init_time_new(ierr, cmessage)

 END SUBROUTINE inFile_pop


 ! *********************************************************************
 ! private subroutine: initialize time data modifeid based on infile pop
 ! *********************************************************************
 SUBROUTINE init_time_new(ierr, message)  ! output

  ! subroutines:
  USE process_time_module, ONLY: process_time  ! process time information
  USE io_netcdf,           ONLY: get_nc        ! netcdf input
  ! derived datatype
  USE dataTypes, ONLY: time           ! time data type
  ! Shared data
  !USE public_var, ONLY: input_dir     ! directory containing input data
  !USE public_var, ONLY: fname_qsim    ! simulated runoff netCDF name
  !USE public_var, ONLY: vname_time    ! variable name for time
  !USE public_var, ONLY: time_units    ! time units (seconds, hours, or days)
  USE public_var, ONLY: simStart      ! date string defining the start of the simulation
  USE public_var, ONLY: simEnd        ! date string defining the end of the simulation
  USE public_var, ONLY: calendar      ! calendar name
  USE globalData, ONLY: timeVar       ! time variables (unit given by runoff data)
  USE globalData, ONLY: iTime         ! time index at simulation time step
  USE globalData, ONLY: convTime2Days ! conversion multipliers for time unit of runoff input to day
  USE globalData, ONLY: refJulday     ! julian day: reference
  USE globalData, ONLY: startJulday   ! julian day: start of routing simulation
  USE globalData, ONLY: endJulday     ! julian day: end of routing simulation
  USE globalData, ONLY: modJulday     ! julian day: at model time step
  USE globalData, ONLY: modTime       ! model time data (yyyy:mm:dd:hh:mm:ss)
  USE globalData, ONLY: infileinfo_data ! the information of the input files

  implicit none

  ! output: error control
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer(i4b)                             :: ix
  character(len=strLen)                    :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_time_new/'

  ! extract time information from the control information
  call process_time(trim(infileinfo_data(1)%unit),  calendar, refJulday,   ierr, cmessage)
  print*, refJulday
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refJulday]'; return; endif
  call process_time(trim(simStart),calendar, startJulday, ierr, cmessage)
  print*, startJulday
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startJulday]'; return; endif
  call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage)
  print*, endJulday
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endJulday]'; return; endif

  ! check that the dates are aligned
  if(endJulday<startJulday) then; ierr=20; message=trim(message)//'simulation end is before simulation start'; return; endif

  ! fast forward time to time index at simStart and save iTime and modJulday
  ! need to convert time unit in timeVar to day
  do ix = 1, size(timeVar)
    modJulday = refJulday + timeVar(ix)/convTime2Days
    if( modJulday < startJulday ) cycle
    exit
  enddo
  iTime = ix
  print*,  iTime

  ! initialize previous model time
  !modTime(0:1) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
  modTime(0) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

 END SUBROUTINE init_time_new


 ! *********************************************************************
 ! private subroutine: get the name of input file based on iTime, will be called
 ! in get_hru_runoff to ajust for file name given iTime
 ! *********************************************************************
 SUBROUTINE infile_name(ierr, message)  ! output

  ! subroutines:
  USE process_time_module, ONLY: process_time  ! process time information
  USE io_netcdf,           ONLY: get_nc        ! netcdf input
  ! derived datatype
  USE dataTypes, ONLY: time           ! time data type
  ! Shared data
  !USE public_var, ONLY: input_dir      ! directory containing input data
  USE public_var, ONLY: fname_qsim     ! simulated runoff netCDF name
  !USE public_var, ONLY: vname_time     ! variable name for time
  !USE public_var, ONLY: time_units     ! time units (seconds, hours, or days)
  !USE public_var, ONLY: simStart       ! date string defining the start of the simulation
  !USE public_var, ONLY: simEnd         ! date string defining the end of the simulation
  !USE public_var, ONLY: calendar       ! calendar name
  USE globalData, ONLY: timeVar         ! time variables (unit given by runoff data)
  USE globalData, ONLY: iTime           ! time index at simulation time step
  !USE globalData, ONLY: convTime2Days  ! conversion multipliers for time unit of runoff input to day
  !USE globalData, ONLY: refJulday      ! julian day: reference
  !USE globalData, ONLY: startJulday    ! julian day: start of routing simulation
  !USE globalData, ONLY: endJulday      ! julian day: end of routing simulation
  !USE globalData, ONLY: modJulday      ! julian day: at model time step
  USE globalData, ONLY: modTime         ! model time data (yyyy:mm:dd:hh:mm:ss)
  USE globalData, ONLY: infileinfo_data ! the information of the input files
  USE globalData, ONLY: iTime_local     ! iTime index for the given netcdf file

  implicit none

  ! output:
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer(i4b)                             :: ix
  character(len=strLen)                    :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_time/'

  ! initialize error control
  print*, "inside infile_name", fname_qsim

  ! fast forward time to time index at simStart and save iTime and modJulday
  ! need to convert time unit in timeVar to day
  ixloop: do ix = 1, size(infileinfo_data) !loop over number of file
   if (iTime.ge.infileinfo_data(ix)%iTimebound(1).and.iTime.le.infileinfo_data(ix)%iTimebound(2)) then
    iTime_local = iTime - infileinfo_data(ix)%iTimebound(1) + 1
    fname_qsim = trim(infileinfo_data(ix)%infilename)
    exit ixloop
   endif
  enddo ixloop

  print*, "inside infile_name iTime", iTime
  print*, "inside infile_name iTime_local", iTime_local
  print*, "inside infile_name file name", fname_qsim

  ! initialize previous model time
  !modTime(0:1) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
  modTime(0) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

 END SUBROUTINE infile_name

 ! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 SUBROUTINE init_time(nTime,     &    ! input: number of time steps
                      ierr, message)  ! output

  ! subroutines:
  USE process_time_module, ONLY: process_time  ! process time information
  USE io_netcdf,           ONLY: get_nc        ! netcdf input
  ! derived datatype
  USE dataTypes, ONLY: time           ! time data type
  ! Shared data
  USE public_var, ONLY: input_dir     ! directory containing input data
  USE public_var, ONLY: fname_qsim    ! simulated runoff netCDF name
  USE public_var, ONLY: vname_time    ! variable name for time
  USE public_var, ONLY: time_units    ! time units (seconds, hours, or days)
  USE public_var, ONLY: simStart      ! date string defining the start of the simulation
  USE public_var, ONLY: simEnd        ! date string defining the end of the simulation
  USE public_var, ONLY: calendar      ! calendar name
  USE globalData, ONLY: timeVar       ! time variables (unit given by runoff data)
  USE globalData, ONLY: iTime         ! time index at simulation time step
  USE globalData, ONLY: convTime2Days ! conversion multipliers for time unit of runoff input to day
  USE globalData, ONLY: refJulday     ! julian day: reference
  USE globalData, ONLY: startJulday   ! julian day: start of routing simulation
  USE globalData, ONLY: endJulday     ! julian day: end of routing simulation
  USE globalData, ONLY: modJulday     ! julian day: at model time step
  USE globalData, ONLY: modTime       ! model time data (yyyy:mm:dd:hh:mm:ss)

  implicit none

  ! input:
  integer(i4b),              intent(in)    :: nTime
  ! output: error control
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer(i4b)                             :: ix
  character(len=strLen)                    :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_time/'

  ! time initialization
  allocate(timeVar(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time data
  call get_nc(trim(input_dir)//trim(fname_qsim), vname_time, timeVar, 1, nTime, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time multiplier needed to convert time to units of days
  select case( trim( time_units(1:index(time_units,' ')) ) )
   case('seconds'); convTime2Days=86400._dp
   case('hours');   convTime2Days=24._dp
   case('days');    convTime2Days=1._dp
   case default;    ierr=20; message=trim(message)//'unable to identify time units'; return
  end select

  ! extract time information from the control information
  call process_time(time_units,    calendar, refJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refJulday]'; return; endif
  call process_time(trim(simStart),calendar, startJulday, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startJulday]'; return; endif
  call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endJulday]'; return; endif

  ! check that the dates are aligned
  if(endJulday<startJulday) then; ierr=20; message=trim(message)//'simulation end is before simulation start'; return; endif

  ! fast forward time to time index at simStart and save iTime and modJulday
  ! need to convert time unit in timeVar to day
  do ix = 1, nTime
    modJulday = refJulday + timeVar(ix)/convTime2Days
    if( modJulday < startJulday ) cycle
    exit
  enddo
  iTime = ix

  ! initialize previous model time
  !modTime(0:1) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
  modTime(0) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

 END SUBROUTINE init_time


 ! *****
 ! private subroutine: get mapping data between runoff hru and river network hru
 ! *********************************************************************
 SUBROUTINE init_runoff(remap_flag,      & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                        remap_data_in,   & ! output: data structure to remap data
                        runoff_data_in,  & ! output: data structure for runoff
                        ierr, message)     ! output: error control

 USE public_var,  ONLY: ancil_dir              ! name of the ancillary directory
 USE public_var,  ONLY: input_dir              ! name of the runoff input directory
 USE public_var,  ONLY: fname_qsim             ! name of simulated runoff netCDF
 USE public_var,  ONLY: fname_remap            ! name of runoff mapping netCDF name
 USE public_var,  ONLY: calendar               ! name of calendar
 USE public_var,  ONLY: time_units             ! time units
 USE globalData,  ONLY: basinID                ! basin ID
 USE dataTypes,   ONLY: remap                  ! remapping data type
 USE dataTypes,   ONLY: runoff                 ! runoff data type
 USE read_runoff, ONLY: read_runoff_metadata   ! read meta data from runoff data
 USE read_remap,  ONLY: get_remap_data         ! read remap data

 implicit none
 ! data structures
 logical(lgt),      intent(in)      :: remap_flag       ! logical whether or not runnoff needs to be mapped to river network HRU
 type(remap)  , intent(out)         :: remap_data_in    ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
 type(runoff) , intent(out)         :: runoff_data_in   ! runoff for one time step for all HRUs
 ! error control
 integer(i4b), intent(out)          :: ierr             ! error code
 character(*), intent(out)          :: message          ! error message
 ! local variables
 integer(i4b)                       :: iHRU, jHRU       ! hru loop indices
 character(len=strLen)              :: cmessage         ! error message from subroutine

 ! initialize error control
 ierr=0; message='init_runoff/'

 ! get runoff metadata
 call read_runoff_metadata(trim(input_dir)//trim(fname_qsim), & ! input: filename
                          runoff_data_in,                     & ! output: runoff data structure
                          time_units, calendar,               & ! output: number of time steps, time units, calendar
                          ierr, cmessage)                       ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 !write(*,*) 'runoff_data_in%nSpace, nTime, trim(time_units) = ', runoff_data_in%nSpace(:), runoff_data_in%nTime, trim(time_units)

 ! need to remap runoff to HRUs
 if (remap_flag) then

   ! get runoff mapping file
   call get_remap_data(trim(ancil_dir)//trim(fname_remap),     & ! input: file name
                       runoff_data_in%nSpace,                  & ! input: number of spatial elements
                       remap_data_in,                          & ! output: data structure to remap data from a polygon
                       ierr, cmessage)                           ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! get indices of the HRU ids in the mapping file in the routing layer
   call get_qix(remap_data_in%hru_id, &  ! input: vector of ids in mapping file
                basinID,              &  ! input: vector of ids in the routing layer
                remap_data_in%hru_ix, &  ! output: indices of hru ids in routing layer
                ierr, cmessage)          ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if ( runoff_data_in%nSpace(2) == integerMissing ) then
     ! get indices of the "overlap HRUs" (the runoff input) in the runoff vector
     call get_qix(remap_data_in%qhru_id, &  ! input: vector of ids in mapping file
                  runoff_data_in%hru_id, &  ! input: vector of ids in runoff file
                  remap_data_in%qhru_ix, &  ! output: indices of mapping ids in runoff file
                  ierr, cmessage)           ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   ! check
   if(count(remap_data_in%hru_ix/=integerMissing)/=size(basinID))then
    message=trim(message)//'unable to identify all polygons in the mapping file'
    ierr=20; return
   endif

   ! check that the basins match
   do iHRU = 1, size(remap_data_in%hru_ix)
     jHRU = remap_data_in%hru_ix(iHRU)
     if (jHRU == integerMissing) cycle
     if( remap_data_in%hru_id(iHRU) /= basinID(jHRU) )then
      message=trim(message)//'mismatch in HRU ids for basins in the routing layer'
      ierr=20; return
     endif
   end do

 else ! if runoff given in RN_HRU

   allocate(runoff_data_in%hru_ix(size(runoff_data_in%hru_id)), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data_in%hru_ix'; return; endif

   ! get indices of the HRU ids in the runoff file in the routing layer
   call get_qix(runoff_data_in%hru_id,  &    ! input: vector of ids in mapping file
                basinID,                &    ! input: vector of ids in the routing layer
                runoff_data_in%hru_ix,  &    ! output: indices of hru ids in routing layer
                ierr, cmessage)              ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! check
   if(count(runoff_data_in%hru_ix/=integerMissing)/=size(basinID))then
    message=trim(message)//'unable to identify all polygons in the mapping file'
    ierr=20; return
   endif

   ! check that the basins match
   do iHRU = 1, size(runoff_data_in%hru_ix)
     jHRU = runoff_data_in%hru_ix(iHRU)
     if (jHRU == integerMissing) cycle
     if( runoff_data_in%hru_id(iHRU) /= basinID(jHRU) )then
      message=trim(message)//'mismatch in HRU ids for basins in the routing layer'
      ierr=20; return
     endif
   end do

 endif

 END SUBROUTINE init_runoff


 ! *****
 ! private subroutine: get indices of mapping points within runoff file...
 ! ***********************************************************************
 SUBROUTINE get_qix(qid,qidMaster,qix,ierr,message)

 USE nr_utility_module, ONLY: indexx  ! get rank of data value

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

 ! initialize error control
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

 write(iulog,*) 'nMissing = ', count(qix==integerMissing)

 ! check again
 do ix=1,size(qid)
  if(qix(ix) /= integerMissing)then
   if(qid(ix) /= qidMaster( qix(ix) ) )then
    message=trim(message)//'unable to find the match'
    ierr=20; return
   endif
  endif
 end do

 END SUBROUTINE get_qix


END MODULE model_setup
