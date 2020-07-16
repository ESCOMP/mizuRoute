MODULE model_setup

USE nrtype,    ONLY: i4b,dp,lgt          ! variable types, etc.
USE nrtype,    ONLY: strLen              ! length of characters

USE public_var, ONLY: iulog              ! i/o logical unit number
USE public_var, ONLY: debug
USE public_var, ONLY: integerMissing
USE public_var, ONLY: realMissing
USE public_var, ONLY: charMissing

USE init_model_data,  ONLY: init_ntopo_data
USE init_model_data,  ONLY: init_state_data
USE init_model_data,  ONLY: get_mpi_omp
USE init_model_data,  ONLY: update_time

USE nr_utility_module, ONLY : unique  ! get unique element array
USE nr_utility_module, ONLY : indexx  ! get rank of data value

implicit none

! privacy -- everything private unless declared explicitly
private
public :: init_mpi
public :: init_data
public :: infile_name

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

   ! runoff input files initialization
   call inFile_pop(ierr, message)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! time initialization
   call init_time(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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

   ! initialize error control
   ierr=0; message='init_runoff_data/'

   if (pid==0) then
     ! runoff and remap data initialization (TO DO: split runoff and remap initialization)
     call init_runoff(is_remap,        & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                      remap_data,      & ! output: data structure to remap data
                      runoff_data,     & ! output: data structure for runoff
                      ierr, cmessage)    ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if  ! if processor=0 (root)

 END SUBROUTINE init_runoff_data


 ! *********************************************************************
 ! public subroutine: read the name of the netcdf that is specified
 ! in a text file, populates the filed of inFiledata dataType
 ! *********************************************************************
 SUBROUTINE inFile_pop(ierr, message)  ! output

  ! data types
  USE dataTypes, ONLY: infileinfo               ! the data type for storing the infromation of the nc files and its attributes

  ! subroutines
  USE process_time_module, ONLY: process_time   ! process time information
  USE ascii_util_module,   ONLY: file_open      ! open file (performs a few checks as well)
  USE ascii_util_module,   ONLY: get_vlines     ! get a list of character strings from non-comment lines
  USE io_netcdf,           ONLY: get_nc         ! get the
  USE io_netcdf,           ONLY: get_var_attr   ! get the attributes interface
  USE io_netcdf,           ONLY: get_nc_dim_len ! get the nc dimension length

  ! Shared data
  USE public_var, ONLY: input_dir               ! directory containing input data
  USE public_var, ONLY: fname_qsim              ! simulated runoff txt file that includes the NetCDF file names
  USE public_var, ONLY: vname_time              ! variable name for time
  USE public_var, ONLY: dname_time              !
  USE globalData, ONLY: infileinfo_data         ! the information of the input files

  ! output: error control
  integer(i4b),         intent(out)    :: ierr             ! error code
  character(*),         intent(out)    :: message          ! error message

  ! local varibales
  integer(i4b)                         :: unit             ! file unit (free unit output from file_open)
  character(len=7)                     :: t_unit           ! time units. "<time_step> since yyyy-MM-dd hh:mm:ss"
  integer(i4b)                         :: iFile            ! counter for forcing files
  integer(i4b)                         :: nFile            ! number of nc files identified in the text file
  integer(i4b)                         :: nTime            ! hard coded for now
  real(dp)                             :: convTime2Days    ! conversion of the day to the local time
  character(len=strLen)                :: infilename       ! input filename
  character(len=strLen),allocatable    :: dataLines(:)     ! vector of lines of information (non-comment lines)
  character(len=strLen)                :: filenameData     ! name of forcing datafile
  character(len=strLen)                :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='inFile_pop/'

  ! build filename and its path containing list of NetCDF files
  infilename = trim(input_dir)//trim(fname_qsim)

  ! open the text file
  call file_open(trim(infilename),unit,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! get a list of character strings from non-commented lines
  call get_vlines(unit,dataLines,ierr,cmessage)
  if(ierr/=0)then; ierr=20; message=trim(message)//trim(cmessage); return; end if
  nFile = size(dataLines) ! get the name of the lines in the file

  ! allocate space for forcing information
  allocate(infileinfo_data(nFile), stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for forcFileInfo'; return; end if

  ! poputate the forcingInfo structure with filenames, and time variables/attributes
  do iFile=1,nFile

   ! split the line into "words" (expect one word: the file describing forcing data for that index)
   read(dataLines(iFile),*,iostat=ierr) filenameData
   if(ierr/=0)then; message=trim(message)//'problem reading a line of data from file ['//trim(infilename)//']'; return; end if

   ! set forcing file name
   infileinfo_data(iFile)%infilename = trim(filenameData)

   ! get the time units
   call get_var_attr(trim(input_dir)//trim(infileinfo_data(iFile)%infilename), &
                     trim(vname_time), 'units', infileinfo_data(iFile)%unit, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! get the calendar
   call get_var_attr(trim(input_dir)//trim(infileinfo_data(iFile)%infilename), &
                     trim(vname_time), 'calendar', infileinfo_data(iFile)%calendar, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! get the dimension of the time to populate nTime and pass it to the get_nc file
   call get_nc_dim_len(trim(input_dir)//trim(infileinfo_data(iFile)%infilename), &
                       trim(dname_time), infileinfo_data(iFile)%nTime, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   nTime = infileinfo_data(iFile)%nTime ! the length of time varibale for each nc file

   ! allocate space for time varibale of each file
   allocate(infileinfo_data(iFile)%timeVar(nTime))
   if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for infileinfo_data(:)%timeVar'; return; end if

   ! get the time varibale
   call get_nc(trim(input_dir)//trim(infileinfo_data(iFile)%infilename), &
               vname_time, infileinfo_data(iFile)%timeVar, 1, nTime, ierr, cmessage) ! does it needs timeVar(:)

   ! get the time multiplier needed to convert time to units of days for each nc file
   t_unit = trim( infileinfo_data(iFile)%unit(1:index(infileinfo_data(iFile)%unit,' ')) )
   select case( trim(t_unit) )
    case('seconds','second','sec','s'); convTime2Days=86400._dp
    case('minutes','minute','min','m'); convTime2Days=1440._dp
    case('hours'  ,'hour'  ,'hr' ,'h'); convTime2Days=24._dp
    case('days'   ,'day'   ,'d');       convTime2Days=1._dp
    case default
      ierr=20; message=trim(message)//'<time_units>= '//trim(t_unit)//': <time_units> must be seconds, minutes, hours or days.'; return
   end select
   infileinfo_data(iFile)%convTime2Days = convTime2Days

   ! get the reference julian day from the nc file
   call process_time(trim(infileinfo_data(iFile)%unit),infileinfo_data(iFile)%calendar,infileinfo_data(iFile)%ncrefjulday,ierr, cmessage)
   if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [ncrefjulday]'; return; endif

   ! populated the index of the iTimebound for each nc file
   if (iFile==1) then
    infileinfo_data(iFile)%iTimebound(1) = 1
    infileinfo_data(iFile)%iTimebound(2) = nTime
   else ! if multiple files specfied in the txt file
    infileinfo_data(iFile)%iTimebound(1) = infileinfo_data(iFile-1)%iTimebound(2) + 1 ! the last index from the perivous nc file + 1
    infileinfo_data(iFile)%iTimebound(2) = infileinfo_data(iFile-1)%iTimebound(2) + nTime ! the last index from the perivous nc file + 1
   endif

  enddo

  ! close ascii file
  close(unit=unit,iostat=ierr); if(ierr/=0)then;message=trim(message)//'problem closing forcing file list'; return; end if

  ! passing the first nc file as global file name to read
  fname_qsim = trim(infileinfo_data(1)%infilename)

 END SUBROUTINE inFile_pop


 ! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 SUBROUTINE init_time(ierr, message)  ! output

  ! subroutines:
  USE process_time_module, ONLY: process_time   ! process time information
  USE process_time_module, ONLY: process_calday ! compute data and time from julian day
  USE io_netcdf,           ONLY: get_nc         ! netcdf input
  ! derived datatype
  USE dataTypes,  ONLY: time                    ! time data type
  ! public data
  USE public_var, ONLY: time_units              ! time units (seconds, hours, or days)
  USE public_var, ONLY: simStart                ! date string defining the start of the simulation
  USE public_var, ONLY: simEnd                  ! date string defining the end of the simulation
  USE public_var, ONLY: calendar                ! calendar name
  USE public_var, ONLY: restart_write           ! restart write option
  USE public_var, ONLY: restart_date            ! restart date
  ! saved time variables
  USE globalData, ONLY: timeVar                 ! time variables (unit given by runoff data)
  USE globalData, ONLY: iTime                   ! time index at simulation time step
  USE globalData, ONLY: refJulday               ! julian day: reference
  USE globalData, ONLY: roJulday                ! julian day: runoff input time
  USE globalData, ONLY: startJulday             ! julian day: start of routing simulation
  USE globalData, ONLY: endJulday               ! julian day: end of routing simulation
  USE globalData, ONLY: modJulday               ! julian day: at model time step
  USE globalData, ONLY: restartJulday           ! julian day: at restart
  USE globalData, ONLY: modTime                 ! model time data (yyyy:mm:dd:hh:mm:ss)
  USE globalData, ONLY: infileinfo_data         ! the information of the input files

  implicit none

  ! output: error control
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer(i4b)                             :: ix
  integer(i4b)                             :: counter
  integer(i4b)                             :: nTime
  integer(i4b)                             :: nt
  integer(i4b)                             :: nFile ! number of nc files
  integer(i4b)                             :: iFile ! for loop over the nc files
  type(time)                               :: rofCal
  type(time)                               :: simCal
  character(len=strLen)                    :: cmessage         ! error message of downwind routine
  character(len=50)                        :: fmt1='(a,I4,a,I2.2,a,I2.2,x,I2.2,a,I2.2,a,F5.2)'

  ! initialize error control
  ierr=0; message='init_time/'

  ! Set time attributes for continuous time variables (saved in globalData to use for output)
  calendar   = infileinfo_data(1)%calendar
  time_units = infileinfo_data(1)%unit
  refJulday  = infileinfo_data(1)%ncrefjulday

  ! get the number of the total time length of all the nc files
  nFile = size(infileinfo_data)
  nTime = 0
  do iFile=1,nFile
   nTime = nTime + infileinfo_data(iFile)%nTime
  enddo

  ! Define time varialbes: timeVar and roJulday
  allocate(timeVar(nTime), roJulday(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! roJulday: Julian day of concatenated netCDF
  counter = 1;
  do iFile=1,nFile
    nt = infileinfo_data(iFile)%nTime
    roJulday(counter:counter+nt-1) = &
    infileinfo_data(iFile)%timeVar(1:nt)/infileinfo_data(iFile)%convTime2Days+infileinfo_data(iFile)%ncrefjulday
    counter = counter + infileinfo_data(iFile)%nTime
  end do

  ! timeVar: time variable in unit given by netCDF
  timeVar(1:nTime) = (roJulday(1:nTime) - refJulday)*infileinfo_data(1)%convTime2Days

  call process_time(trim(simStart),calendar, startJulday, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startJulday]'; return; endif
  call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endJulday]'; return; endif

  ! check that the dates are aligned
  if(endJulday<startJulday) then
    write(cmessage,'(7a)') 'simulation end is before simulation start:', new_line('a'), '<sim_start>= ', trim(simStart), new_line('a'), '<sim_end>= ', trim(simEnd)
    ierr=20; message=trim(message)//trim(cmessage); return
  endif

  ! check sim_start is before the last time step in runoff data
  if(startJulday>roJulday(nTime)) then
    call process_calday(roJulday(nTime), calendar, rofCal, ierr, cmessage)
    call process_calday(startJulday, calendar, simCal, ierr, cmessage)
    write(iulog,'(2a)') new_line('a'),'ERROR: <sim_start> is after the first time step in input runoff'
    write(iulog,fmt1)  ' runoff_end  : ', rofCal%iy,'-',rofCal%im,'-',rofCal%id, rofCal%ih,':', rofCal%imin,':',rofCal%dsec
    write(iulog,fmt1)  ' <sim_start> : ', simCal%iy,'-',simCal%im,'-',simCal%id, simCal%ih,':', simCal%imin,':',simCal%dsec
    ierr=20; message=trim(message)//'check <sim_start> against runoff input time'; return
  endif

  ! Compare sim_start vs. time at first time step in runoff data
  if (startJulday < roJulday(1)) then
    call process_calday(roJulday(1), calendar, rofCal, ierr, cmessage)
    call process_calday(startJulday, calendar, simCal, ierr, cmessage)
    write(iulog,'(2a)') new_line('a'),'WARNING: <sim_start> is before the first time step in input runoff'
    write(iulog,fmt1)  ' runoff_start: ', rofCal%iy,'-',rofCal%im,'-',rofCal%id, rofCal%ih,':', rofCal%imin,':',rofCal%dsec
    write(iulog,fmt1)  ' <sim_start> : ', simCal%iy,'-',simCal%im,'-',simCal%id, simCal%ih,':', simCal%imin,':',simCal%dsec
    write(iulog,'(a)') ' Reset <sim_start> to runoff_start'
    startJulday = roJulday(1)
  endif

  ! Compare sim_end vs. time at last time step in runoff data
  if (endJulday > roJulday(nTime)) then
    call process_calday(roJulday(nTime), calendar, rofCal, ierr, cmessage)
    call process_calday(endJulday,       calendar, simCal, ierr, cmessage)
    write(iulog,'(2a)')  new_line('a'),'WARNING: <sim_end> is after the last time step in input runoff'
    write(iulog,fmt1)   ' runoff_end: ', rofCal%iy,'-',rofCal%im,'-',rofCal%id, rofCal%ih,':', rofCal%imin,':',rofCal%dsec
    write(iulog,fmt1)   ' <sim_end> : ', simCal%iy,'-',simCal%im,'-',simCal%id, simCal%ih,':', simCal%imin,':',simCal%dsec
    write(iulog,'(a)')  ' Reset <sim_end> to runoff_end'
    endJulday = roJulday(nTime)
  endif

  ! fast forward time to time index at simStart and save iTime and modJulday
  do ix = 1, nTime
    modJulday = roJulday(ix)
    if( modJulday < startJulday ) cycle
    exit
  enddo
  iTime = ix

  ! initialize previous model time
  modTime(0) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

  ! restart drop off time
  select case(trim(restart_write))
    case('last','Last')
      call process_time(trim(simEnd), calendar, restartJulday, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [restartDate]'; return; endif
    case('never','Never')
      restartJulday = 0.0_dp
    case('specified','Specified')
      if (trim(restart_date) == charMissing) then
        ierr=20; message=trim(message)//'<restart_date> must be provided when <restart_write> option is "specified"'; return
      end if
      call process_time(trim(restart_date),calendar, restartJulday, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [restartDate]'; return; endif
    case default
      ierr=20; message=trim(message)//'Current accepted <restart_write> options: last[Last], never[Never], specified[Specified]'; return
  end select

 END SUBROUTINE init_time


 ! *********************************************************************
 ! public subroutine: get the name of input file based on iTime, will be called
 ! in get_hru_runoff to ajust for file name given iTime
 ! *********************************************************************
 SUBROUTINE infile_name(ierr, message)  ! output

  ! subroutines:
  USE process_time_module, ONLY: process_time  ! process time information
  USE io_netcdf,           ONLY: get_nc        ! netcdf input
  ! derived datatype
  USE dataTypes, ONLY: time           ! time data type
  ! Shared data
  USE public_var, ONLY: fname_qsim     ! simulated runoff netCDF name
  USE globalData, ONLY: iTime           ! time index at simulation time step
  USE globalData, ONLY: infileinfo_data ! the information of the input files
  USE globalData, ONLY: iTime_local     ! iTime index for the given netcdf file

  implicit none

  ! output:
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer(i4b)                             :: ix
  !character(len=strLen)                    :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='infile_name/'

  ! fast forward time to time index at simStart and save iTime and modJulday
  ixloop: do ix = 1, size(infileinfo_data) !loop over number of file
   if ((iTime >= infileinfo_data(ix)%iTimebound(1)).and.(iTime <= infileinfo_data(ix)%iTimebound(2))) then
    iTime_local = iTime - infileinfo_data(ix)%iTimebound(1) + 1
    fname_qsim = trim(infileinfo_data(ix)%infilename)
    exit ixloop
   endif
  enddo ixloop

 END SUBROUTINE infile_name

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
 integer(i4b), allocatable          :: unq_qhru_id(:)
 integer(i4b), allocatable          :: unq_idx(:)
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

   if (debug) then
     write(iulog,'(2a)') new_line('a'), 'DEBUG: Corresponding between River-Network(RN) hru in mapping data and RN hru in river network data'
     write(iulog,'(2x,a,I15)') '(1) number of RN hru in river-network = ', size(basinID)
     write(iulog,'(2x,a,I15)') '(2) number of RN hru in mapping       = ', size(remap_data_in%hru_id)
     write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two  = ', count(remap_data_in%hru_ix/=integerMissing)
     if(count(remap_data_in%hru_ix/=integerMissing)/=size(basinID))then
       message=trim(message)//'(1) not equal (2)'
       ierr=20; return
     endif
   end if

   if ( runoff_data_in%nSpace(2) == integerMissing ) then
     ! get indices of the "overlap HRUs" (the runoff input) in the runoff vector
     call get_qix(remap_data_in%qhru_id, &  ! input: vector of ids in mapping file
                  runoff_data_in%hru_id, &  ! input: vector of ids in runoff file
                  remap_data_in%qhru_ix, &  ! output: indices of mapping ids in runoff file
                  ierr, cmessage)           ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     if (debug) then
       call unique(remap_data_in%qhru_id, unq_qhru_id, unq_idx)
       write(iulog,'(2a)') new_line('a'),'DEBUG: corresponding between Hydro-Model (HM) hru in mapping data and HM hru in runoff data'
       write(iulog,'(2x,a,I15)') '(1) number of HM hru in hyrdo-model  = ', size(runoff_data_in%hru_id)
       write(iulog,'(2x,a,I15)') '(2) number of HM hru in mapping      = ', size(unq_qhru_id)
       write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two = ', count(remap_data_in%qhru_ix(unq_idx)/=integerMissing)
     end if
   end if

 else ! if runoff given in RN_HRU

   allocate(runoff_data_in%hru_ix(size(runoff_data_in%hru_id)), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data_in%hru_ix'; return; endif

   ! get indices of the HRU ids in the runoff file in the routing layer
   call get_qix(runoff_data_in%hru_id,  &    ! input: vector of ids in mapping file
                basinID,                &    ! input: vector of ids in the routing layer
                runoff_data_in%hru_ix,  &    ! output: indices of hru ids in routing layer
                ierr, cmessage)              ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (debug) then
     write(iulog,'(2a)') new_line('a'), 'DEBUG: corresponding between River-Network (RN) hru in runoff data and RN hru in river network data'
     write(iulog,'(2x,a,I15)') '(1) number of RN hru in river-network = ', size(basinID)
     write(iulog,'(2x,a,I15)') '(2) number of RN hru in hyrdo-model   = ', size(runoff_data_in%hru_id)
     write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two  = ', count(runoff_data_in%hru_ix/=integerMissing)
     if(count(runoff_data_in%hru_ix/=integerMissing)/=size(basinID))then
       message=trim(message)//'(1) not equal (2)'
       ierr=20; return
     endif
   end if

 endif

 END SUBROUTINE init_runoff


 ! *****
 ! private subroutine: get indices of mapping points within runoff file...
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
