MODULE model_setup

! data types
USE nrtype,    ONLY : i4b,i8b,dp,lgt          ! variable types, etc.
USE nrtype,    ONLY : strLen              ! length of characters
USE dataTypes, ONLY : var_ilength         ! integer type:          var(:)%dat
USE dataTypes, ONLY : var_clength         ! integer type:          var(:)%dat
USE dataTypes, ONLY : var_dlength,dlength ! double precision type: var(:)%dat, or dat

USE public_var, ONLY : iulog
USE public_var, ONLY : debug
USE public_var, ONLY : verySmall
USE public_var, ONLY : integerMissing
USE public_var, ONLY : realMissing
USE public_var, ONLY : charMissing

USE io_netcdf, ONLY : close_nc         ! close netcdf

USE nr_utility_module, ONLY : findIndex ! get array index of matching element
USE nr_utility_module, ONLY : unique  ! get unique element array
USE nr_utility_module, ONLY : indexx  ! get rank of data value

implicit none

private

public :: init_model
public :: init_data
public :: update_time

CONTAINS

 ! *********************************************************************
 ! public subroutine: model setup
 ! *********************************************************************
 SUBROUTINE init_model(cfile_name, ierr, message)

  ! used to read control files and namelist and broadcast to all processors

  ! shared data used
  USE public_var,          ONLY : ancil_dir
  USE public_var,          ONLY : param_nml
  USE globalData,          ONLY : nThreads         ! a number of threads
  ! subroutines: populate metadata
  USE popMetadat_module,   ONLY : popMetadat       ! populate metadata
  ! subroutines: model control
  USE read_control_module, ONLY : read_control     ! read the control file
  USE read_param_module,   ONLY : read_param       ! read the routing parameters

  implicit none

  character(*), intent(in)    :: cfile_name          ! name of the control file
  ! output: error control
  integer(i4b), intent(out)   :: ierr                ! error code
  character(*), intent(out)   :: message             ! error message
  ! local variables
  integer(i4b)                :: omp_get_num_threads ! number of threads used for openMP
  character(len=strLen)       :: cmessage            ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_model/'

  ! Get number of threads
  nThreads = 1
  !$OMP PARALLEL
  !$ nThreads = omp_get_num_threads()
  !$OMP END PARALLEL

  ! populate the metadata files
  call popMetadat(ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! read the control file
  call read_control(trim(cfile_name), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! read the routing parameter namelist
  call read_param(trim(ancil_dir)//trim(param_nml),ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE init_model


 ! *********************************************************************
 ! public subroutine: initialize river network, runoff, and runoff-mapping data
 ! *********************************************************************
 subroutine init_data(ierr, message)

  USE public_var,  ONLY : is_remap               ! logical whether or not runnoff needs to be mapped to river network HRU
  USE public_var,  ONLY : ntopAugmentMode        ! River network augmentation mode
  USE public_var,  ONLY : idSegOut               ! outlet segment ID (-9999 => no outlet segment specified)
  USE public_var,  ONLY : desireId               ! ID of reach to be checked by on-screen printing
  USE var_lookup,  ONLY : ixHRU2SEG              ! index of variables for data structure
  USE var_lookup,  ONLY : ixNTOPO                ! index of variables for data structure
  USE globalData,  ONLY : RCHFLX                 ! Reach flux data structures (entire river network)
  USE globalData,  ONLY : RCHSTA                 ! Reach state structures (entire river network)

  USE globalData,  ONLY : nHRU, nRch             ! number of HRUs and Reaches in the whole network
  USE globalData,  ONLY : nEns                   ! number of ensembles
  USE globalData,  ONLY : nRoutes                ! number of active routing methods
  USE globalData,  ONLY : basinID                ! HRU id vector
  USE globalData,  ONLY : reachID                ! reach ID vector
  USE globalData,  ONLY : ixPrint                ! reach index to be examined by on-screen printing
  USE globalData,  ONLY : runoff_data            ! runoff data structure
  USE globalData,  ONLY : remap_data             ! runoff mapping data structure

   implicit none
   ! input:
   ! None
   ! output: error control
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variable
   ! river network data structures for the entire domain
   type(var_dlength), allocatable           :: structHRU(:)     ! HRU properties
   type(var_dlength), allocatable           :: structSeg(:)     ! stream segment properties
   type(var_ilength), allocatable           :: structHRU2SEG(:) ! HRU-to-segment mapping
   type(var_ilength), allocatable           :: structNTOPO(:)   ! network topology
   type(var_clength), allocatable           :: structPFAF(:)    ! pfafstetter code
   ! others
   integer(i4b)                             :: iHRU, iRch       ! loop index
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ! initialize error control
   ierr=0; message='init_data/'

   ! populate various river network data strucutures for each proc
   ! read the river network data and compute additonal network attributes (inncludes spatial decomposition)
   call init_ntopo(nHRU, nRch,                                       &             ! output: number of HRU and Reaches
                   structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                   ierr, cmessage)                                                 ! output: error controls
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! check if network topology write option is on. If so, terminate the program
   if (ntopAugmentMode .or. idSegOut>0) stop

   ! allocate space for the entire river network
   allocate(RCHFLX(nEns,nRch), RCHSTA(nEns,nRch), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [RCHFLX, RCHSTA]'; return; endif

   do iRch = 1,nRch
     allocate(RCHFLX(nEns,iRch)%ROUTE(nRoutes))
   end do

   ! populate basiID and reachID vectors for output (in ONLY master processor)
   ! populate runoff data structure (only meta, no runoff values)
   ! populate remap data structure

   allocate(basinID(nHRU), reachID(nRch), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [basinID, reachID]'; return; endif

   do iHRU = 1,nHRU
     basinID(iHRU) = structHRU2SEG(iHRU)%var(ixHRU2SEG%hruId)%dat(1)
   enddo
   do iRch = 1,nRch
     reachID(iRch) = structNTOPO(iRch)%var(ixNTOPO%segId)%dat(1)
   end do

   ! get reach index to be examined by on-screen printing
   if (desireId/=integerMissing) ixPrint = findIndex(reachID, desireId, integerMissing)

   ! runoff and remap data initialization (TO DO: split runoff and remap initialization)
   call init_runoff(is_remap,        & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                    remap_data,      & ! output: data structure to remap data
                    runoff_data,     & ! output: data structure for runoff
                    ierr, cmessage)    ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! DateTime initialization
   call init_time(runoff_data%ntime, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! channel state initialization
   call init_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine init_data


 ! *********************************************************************
 ! public subroutine: update time to next time step
 ! *********************************************************************
 SUBROUTINE update_time(finished, ierr, message)

  USE public_var, ONLY : dt
  USE public_var, ONLY : calendar
  USE globalData, ONLY : TSEC          ! beginning/ending of simulation time step [sec]
  USE globalData, ONLY : iTime         ! time index at simulation time step
  USE globalData, ONLY : simout_nc     ! netCDF meta data
  USE globalData, ONLY : endCal        ! model ending datetime
  USE globalData, ONLY : modTime       ! model datetime

   implicit none
   ! output
   logical(lgt),              intent(out)   :: finished
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variables
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ! initialize error control
   ierr=0; message='update_time/'

   finished = .false.
   if (modTime(1)==endCal) then
     finished=.true.
     if (simout_nc%status == 2) then
       call close_nc(simout_nc%ncid, ierr, cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end if
     return
   endif

   ! update model time step bound
   TSEC(0) = TSEC(0) + dt
   TSEC(1) = TSEC(0) + dt

   ! update time index
   iTime=iTime+1

   ! increment model calendar
   modTime(0) = modTime(1)
   modTime(1) = modTime(1)%add_sec(dt, calendar, ierr, cmessage)

 END SUBROUTINE update_time


 ! *********************************************************************
 ! private subroutine: initialize channel state data
 ! *********************************************************************
 SUBROUTINE init_state(ierr, message)
  ! subroutines
  USE ascii_util_module, ONLY : lower             ! convert string to lower case
  USE read_restart,      ONLY : read_state_nc     ! read netcdf state output file
  ! global data
  USE public_var,    ONLY : dt                    ! simulation time step (seconds)
  USE public_var,    ONLY : impulseResponseFunc   ! IRF routing ID = 1
  USE public_var,    ONLY : kinematicWaveTracking ! KWT routing ID = 2
  USE public_var,    ONLY : kinematicWave         ! KW routing ID = 3
  USE public_var,    ONLY : muskingumCunge        ! MC routing ID = 4
  USE public_var,    ONLY : diffusiveWave         ! DW routing ID = 5
  USE public_var,    ONLY : fname_state_in        ! name of state input file
  USE public_var,    ONLY : restart_dir           ! directory containing output data
  USE globalData,    ONLY : nRoutes               !
  USE globalData,    ONLY : routeMethods          ! ID of active routing method
  USE globalData,    ONLY : RCHFLX                ! reach flux structure
  USE globalData,    ONLY : RCHSTA                ! reach restart state structure
  USE globalData,    ONLY : nMolecule             ! computational molecule
  USE globalData,    ONLY : TSEC                  ! begining/ending of simulation time step [sec]

  implicit none

  ! output: error control
  integer(i4b),        intent(out) :: ierr             ! error code
  character(*),        intent(out) :: message          ! error message
  ! local variable
  real(dp)                         :: T0,T1            ! begining/ending of simulation time step [sec]
  integer(i4b)                     :: iens             ! ensemble index (currently only 1)
  integer(i4b)                     :: ix,iRoute        ! loop index
  character(len=strLen)            :: cmessage         ! error message of downwind routine

  ierr=0; message='init_state/'

  iens = 1_i4b

  ! read restart file and initialize states
  if (trim(fname_state_in)==charMissing .or. lower(trim(fname_state_in))=='none' .or. lower(trim(fname_state_in))=='coldstart') then
    ! Cold start .......
    ! initialize flux structures
    RCHFLX(:,:)%BASIN_QI = 0._dp
    RCHFLX(:,:)%BASIN_QR(0) = 0._dp
    RCHFLX(:,:)%BASIN_QR(1) = 0._dp

    do iRoute = 1, nRoutes
      if (routeMethods(iRoute)==impulseResponseFunc) then
        do ix = 1, size(RCHSTA(1,:))
          RCHFLX(iens,ix)%ROUTE(iRoute)%REACH_VOL(0:1) = 0._dp
        end do
      else if (routeMethods(iRoute)==kinematicWaveTracking) then
        do ix = 1, size(RCHSTA(1,:))
          RCHFLX(iens,ix)%ROUTE(iRoute)%REACH_VOL(0:1) = 0._dp
        end do
      else if (routeMethods(iRoute)==kinematicWave) then
        nMolecule%KW_ROUTE = 2
        do ix = 1, size(RCHSTA(1,:))
          RCHFLX(iens,ix)%ROUTE(iRoute)%REACH_VOL(0:1) = 0._dp
          allocate(RCHSTA(iens,ix)%KW_ROUTE%molecule%Q(nMolecule%KW_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA]'; return; endif
          RCHSTA(iens,ix)%KW_ROUTE%molecule%Q(:) = 0._dp
        end do
      else if (routeMethods(iRoute)==muskingumCunge) then
        nMolecule%MC_ROUTE = 2
        do ix = 1, size(RCHSTA(1,:))
          RCHFLX(iens,ix)%ROUTE(iRoute)%REACH_VOL(0:1) = 0._dp
          allocate(RCHSTA(iens,ix)%MC_ROUTE%molecule%Q(nMolecule%MC_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA]'; return; endif
          RCHSTA(iens,ix)%MC_ROUTE%molecule%Q(:) = 0._dp
        end do
      else if (routeMethods(iRoute)==diffusiveWave) then
        nMolecule%DW_ROUTE = 5
        do ix = 1, size(RCHSTA(1,:))
          RCHFLX(iens,ix)%ROUTE(iRoute)%REACH_VOL(0:1) = 0._dp
          allocate(RCHSTA(iens,ix)%DW_ROUTE%molecule%Q(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA]'; return; endif
          RCHSTA(iens,ix)%DW_ROUTE%molecule%Q(:) = 0._dp
        end do
      end if
    end do

    ! initialize time
    TSEC(0)=0._dp; TSEC(1)=dt
  else
    call read_state_nc(trim(restart_dir)//trim(fname_state_in), T0, T1, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    TSEC(0)=T0; TSEC(1)=T1
  endif

 END SUBROUTINE init_state

! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 SUBROUTINE init_time(nTime,     &    ! input: number of time steps
                      ierr, message)  ! output

  ! subroutines:
  USE ascii_util_module,   ONLY : lower           ! convert string to lower case
  USE io_netcdf,           ONLY : open_nc         ! netcdf input
  USE io_netcdf,           ONLY : close_nc        ! netcdf input
  USE io_netcdf,           ONLY : get_nc          ! netcdf input
  ! derived datatype
  USE date_time,           ONLY : datetime        ! time data type
  ! public data
  USE public_var,          ONLY : input_dir     ! directory containing input data
  USE public_var,          ONLY : fname_qsim    ! simulated runoff netCDF name
  USE public_var,          ONLY : vname_time    ! variable name for time
  USE public_var,          ONLY : time_units    ! time units (seconds, hours, or days)
  USE public_var,          ONLY : simStart      ! date string defining the start of the simulation
  USE public_var,          ONLY : simEnd        ! date string defining the end of the simulation
  USE public_var,          ONLY : calendar      ! calendar name
  USE public_var,          ONLY : dt
  USE public_var,          ONLY : secprday
  USE public_var,          ONLY : restart_write ! restart write option
  USE public_var,          ONLY : restart_date  ! restart date
  USE public_var,          ONLY : restart_month !
  USE public_var,          ONLY : restart_day   !
  USE public_var,          ONLY : restart_hour  !
  ! saved time variables
  USE globalData,          ONLY : timeVar       ! time variables (unit given by runoff data)
  USE globalData,          ONLY : iTime         ! time index at runoff input time step
  USE globalData,          ONLY : modTime       ! model time data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,          ONLY : endCal        ! simulation end time data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,          ONLY : restCal       ! restart time data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,          ONLY : dropCal       ! restart dropoff calendar date/time

  implicit none

  ! input:
  integer(i4b),              intent(in)    :: nTime
  ! output: error control
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer(i4b)                             :: ncidRunoff
  integer(i4b)                             :: ix
  type(datetime)                           :: refCal
  type(datetime)                           :: roCal(nTime)
  type(datetime)                           :: startCal
  type(datetime)                           :: dummyCal
  integer(i4b)                             :: nDays          ! number of days in a month
  real(dp)                                 :: convTime2sec
  real(dp)                                 :: sec(nTime)
  character(len=7)                         :: t_unit
  character(len=strLen)                    :: cmessage         ! error message of downwind routine
  character(len=50)                        :: fmt1='(a,I4,a,I2.2,a,I2.2,x,I2.2,a,I2.2,a,F5.2)'

  ! initialize error control
  ierr=0; message='init_time/'

  ! time initialization
  allocate(timeVar(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time data
  call open_nc(trim(input_dir)//trim(fname_qsim), 'r', ncidRunoff, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  call get_nc(ncidRunoff, vname_time, timeVar, 1, nTime, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  call close_nc(ncidRunoff, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time multiplier needed to convert time to units of days
  t_unit = trim( time_units(1:index(time_units,' ')) )
  select case( trim(t_unit) )
   case('seconds','second','sec','s'); convTime2sec=1._dp
   case('minutes','minute','min');     convTime2sec=60._dp
   case('hours','hour','hr','h');      convTime2sec=3600._dp
   case('days','day','d');             convTime2sec=86400._dp
   case default
     ierr=20; message=trim(message)//'<time_units>= '//trim(t_unit)//': <time_units> must be seconds, minutes, hours or days.'; return
  end select

  ! extract datetime from the control information
  call refCal%str2datetime(time_units, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refCal]'; return; endif
  call startCal%str2datetime(simStart, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startCal]'; return; endif
  call endCal%str2datetime(simEnd, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endCal]'; return; endif

  ! calendar in runoff data
  sec(:) = timeVar(:)*convTime2sec
  do ix = 1, nTime
    roCal(ix) = refCal%add_sec(sec(ix), calendar, ierr, cmessage)
  end do

  ! check that the dates are aligned
  if(endCal<startCal) then
    write(cmessage,'(7a)') 'simulation end is before simulation start:', new_line('a'), '<sim_start>= ', trim(simStart), new_line('a'), '<sim_end>= ', trim(simEnd)
    ierr=20; message=trim(message)//trim(cmessage); return
  endif

  ! check sim_start is before the last time step in runoff data
  if(startCal > roCal(nTime)) then
    write(iulog,'(2a)') new_line('a'),'ERROR: <sim_start> is after the first time step in input runoff'
    write(iulog,fmt1) ' runoff_end  : ', roCal(nTime)%year(),'-',roCal(nTime)%month(),'-',roCal(nTime)%day(),roCal(nTime)%hour(),':',roCal(nTime)%minute(),':',roCal(nTime)%sec()
    write(iulog,fmt1) ' <sim_start> : ', startCal%year(),'-',startCal%month(),'-',startCal%day(), startCal%hour(),':', startCal%minute(),':',startCal%sec()
    ierr=20; message=trim(message)//'check <sim_start> against runoff input time'; return
  endif

  ! Compare sim_start vs. time at first time step in runoff data
  if (startCal < roCal(1)) then
    write(iulog,'(2a)') new_line('a'),'WARNING: <sim_start> is before the first time step in input runoff'
    write(iulog,fmt1) ' runoff_start : ', roCal(1)%year(),'-',roCal(1)%month(),'-',roCal(1)%day(), roCal(1)%hour(),':', roCal(1)%minute(),':',roCal(1)%sec()
    write(iulog,fmt1) ' <sim_start>  : ', startCal%year(),'-',startCal%month(),'-',startCal%day(), startCal%hour(),':', startCal%minute(),':',startCal%sec()
    write(iulog,'(a)') ' Reset <sim_start> to runoff_start'
    startCal = roCal(1)
  endif

  ! Compare sim_end vs. time at last time step in runoff data
  if (endCal > roCal(nTime)) then
    write(iulog,'(2a)')  new_line('a'),'WARNING: <sim_end> is after the last time step in input runoff'
    write(iulog,fmt1) ' runoff_end : ', roCal(nTime)%year(),'-',roCal(nTime)%month(),'-',roCal(nTime)%day(),roCal(nTime)%hour(),':',roCal(nTime)%minute(),':',roCal(nTime)%sec()
    write(iulog,fmt1) ' <sim_end>  : ', endCal%year(),'-',endCal%month(),'-',endCal%day(), endCal%hour(),':', endCal%minute(),':',endCal%sec()
    write(iulog,'(a)')  ' Reset <sim_end> to runoff_end'
    endCal = roCal(nTime)
  endif

  ! fast forward time to time index at simStart and save iTime and modTime(1)
  do ix = 1, nTime
    modTime(1) = roCal(ix)
    if( modTime(1) < startCal ) cycle
    exit
  enddo
  iTime = ix

 ! Set restart calendar date/time and dropoff calendar date/time and
 ! -- For periodic restart options  ---------------------------------------------------------------------
 ! Ensure that user-input restart month, day are valid.
 ! "Annual" option:  if user input day exceed number of days given user input month, set to last day
 ! "Monthly" option: use 2000-01 as template calendar yr/month
 ! "Daily" option:   use 2000-01-01 as template calendar yr/month/day
 select case(lower(trim(restart_write)))
   case('yearly')
     call dummyCal%set_datetime(2000, restart_month, 1, 0, 0, 0.0_dp)
     nDays = dummyCal%ndays_month(calendar, ierr, cmessage)
     if(ierr/=0) then; message=trim(message)//trim(cmessage); return; endif
     if (restart_day > nDays) restart_day=nDays
   case('monthly'); restart_month = 1
   case('daily');     restart_month = 1; restart_day = 1
 end select

  select case(lower(trim(restart_write)))
    case('last')
      dropCal = endCal
      restart_month = dropCal%month(); restart_day = dropCal%day(); restart_hour = dropCal%hour()
    case('specified')
      if (trim(restart_date) == charMissing) then
        ierr=20; message=trim(message)//'<restart_date> must be provided when <restart_write> option is "specified"'; return
      end if
      call restCal%str2datetime(restart_date, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [restart_date]'; return; endif
      dropCal = restCal%add_sec(-dt, calendar, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [restCal->dropCal]'; return; endif
      restart_month = dropCal%month(); restart_day = dropCal%day(); restart_hour = dropCal%hour()
    case('yearly','monthly','daily')
      call restCal%set_datetime(2000, restart_month, restart_day, restart_hour, 0, 0._dp)
      dropCal = restCal%add_sec(-dt, calendar, ierr, cmessage)
      if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [ dropCal for periodical restart]'; return; endif
      restart_month = dropCal%month(); restart_day = dropCal%day(); restart_hour = dropCal%hour()
    case('never')
      call dropCal%set_datetime(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
    case default
      ierr=20; message=trim(message)//'Current accepted <restart_write> options: last, never, specified, yearly, monthly, daily'; return
  end select

 END SUBROUTINE init_time


 ! *********************************************************************
 ! private subroutine: initialize river network data
 ! *********************************************************************
 SUBROUTINE init_ntopo(nHRU_out, nRch_out,                                           & ! output: number of HRU and Reaches
                       structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                       ierr, message)                                                  ! output: error controls
  ! public vars
  USE public_var,           ONLY : ancil_dir                ! name of the ancillary directory
  USE public_var,           ONLY : fname_ntopOld            ! name of the old network topology file
  USE public_var,           ONLY : fname_ntopNew            ! name of the new network topology file
  USE public_var,           ONLY : dname_nhru               ! dimension name for HRUs
  USE public_var,           ONLY : dname_sseg               ! dimension name for stream segments
  USE public_var,           ONLY : maxPfafLen               ! maximum digit of pfafstetter code (default 32)
  ! options
  USE public_var,           ONLY : ntopAugmentMode          ! River network augmentation mode
  USE public_var,           ONLY : idSegOut                 ! River network subset mode (idSegOut > 0)
  ! global data
  USE globalData,           ONLY : meta_PFAF                ! meta for pfafstetter code
  USE globalData,           ONLY : NETOPO, RPARAM           ! network and parameter data structure used in routing routine
  USE globalData,           ONLY : river_basin              ! OMP domain decompostion data strucuture
  ! variable index
  USE var_lookup,           ONLY : ixPFAF                   ! index of variables for the pfafstetter code
  ! external subroutines
  USE read_streamSeg,       ONLY : getData                  ! get the ancillary data
  USE write_streamSeg,      ONLY : writeData                ! write the ancillary data
  USE process_ntopo,        ONLY : check_river_properties   ! check if river network data is physically valid
  USE io_netcdf,            ONLY : get_var_dims
  USE process_ntopo,        ONLY : augment_ntopo            ! compute all the additional network topology (only compute option = on)
  USE process_ntopo,        ONLY : put_data_struct          ! populate NETOPO and RPARAM data structure
  USE domain_decomposition, ONLY : omp_domain_decomposition     ! domain decomposition for omp
!  USE domain_decomposition, ONLY : omp_domain_decomposition &    ! domain decomposition for omp
!                                => omp_domain_decomposition_stro
  implicit none
  ! input: None
  ! output (river network data structures for the entire domain)
  integer(i4b)                  , intent(out) :: nHRU_out                 ! number of HRUs
  integer(i4b)                  , intent(out) :: nRch_out                 ! number of reaches
  type(var_dlength), allocatable, intent(out) :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(out) :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(out) :: structHRU2SEG(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable, intent(out) :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable, intent(out) :: structPFAF(:)            ! pfafstetter code
  ! output: error control
  integer(i4b)      , intent(out)             :: ierr                     ! error code
  character(*)      , intent(out)             :: message                  ! error message
  ! local variable
  integer(i4b)                                :: tot_upstream             ! total number of all of the upstream stream segments for all stream segments
  integer(i4b)                                :: tot_upseg                ! total number of immediate upstream segments for all  stream segments
  integer(i4b)                                :: tot_hru                  ! total number of all the upstream hrus for all stream segments
  integer(i4b)                                :: tot_uh                   ! total number of unit hydrograph from all the stream segments
  integer(i4b),      allocatable              :: ixHRU_desired(:)         ! indices of desired hrus
  integer(i4b),      allocatable              :: ixSeg_desired(:)         ! indices of desired reaches
  integer(i4b)                                :: dummy(2)                 ! dummy variable to hold dimension length for 2D variables in netCDF
  integer(i4b)   , parameter                  :: maxUpstreamFile=10000000 ! 10 million: maximum number of upstream reaches to enable writing
  character(len=strLen)                       :: cmessage                 ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_ntopo/'

  if (meta_PFAF(ixPFAF%code)%varFile) then
    ! get the variable dimensions
    ! NOTE: need to update maxPfafLen to the exact character size for pfaf code in netCDF
    call get_var_dims(trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
                      trim(meta_PFAF(ixPFAF%code)%varName), & ! input: pfaf code variable name in netcdf
                      ierr, cmessage,                       & ! output: error control
                      dlen=dummy)                             ! output optional: dimension length
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    maxPfafLen = dummy(1)
  end if

  ! read river network data
  call getData(&
               ! input
               trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
               dname_nhru,   & ! input: dimension name of the HRUs
               dname_sseg,   & ! input: dimension name of the stream segments
               ! output: model control
               nHRU_out,      & ! output: number of HRUs
               nRch_out,      & ! output: number of stream segments
               ! output: populate data structures
               structHRU,    & ! ancillary data for HRUs
               structSeg,    & ! ancillary data for stream segments
               structHRU2seg,& ! ancillary data for mapping hru2basin
               structNTOPO,  & ! ancillary data for network topology
               structPFAF,   & ! ancillary data for pfafstetter code
               ! output: error control
               ierr,cmessage) ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call check_river_properties(structNTOPO, structHRU, structSEG, ierr, cmessage) ! input: data structure for physical river network data
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute additional network attributes
  call augment_ntopo(&
                     ! input: model control
                     nHRU_out,                         & ! number of HRUs
                     nRch_out,                         & ! number of stream segments
                     ! inout: populate data structures
                     structHRU,                        & ! ancillary data for HRUs
                     structSeg,                        & ! ancillary data for stream segments
                     structHRU2seg,                    & ! ancillary data for mapping hru2basin
                     structNTOPO,                      & ! ancillary data for network toopology
                     ! output:
                     ierr, cmessage,                   & ! error control
                     ! optional output
                     tot_hru       = tot_hru,          & ! total number of all the upstream hrus for all stream segments
                     tot_upseg     = tot_upseg,        & ! total number of all the immediate upstream segments for all stream segments
                     tot_upstream  = tot_upstream,     & ! total number of all the upstream segments for all stream segments
                     tot_uh        = tot_uh,           & ! total number of unit hydrograph for all stream segments
                     ixHRU_desired = ixHRU_desired,    & ! indices of desired hrus
                     ixSeg_desired = ixSeg_desired)      ! indices of desired reaches
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write network topology (if augment mode or subset mode)
  if(ntopAugmentMode .or. idSegOut>0)then

    ! disable the dimension containing all upstream reaches
    ! NOTE: For the CONUS this is 1,872,516,819 reaches !!
    !        --> it will always be quicker to recompute than read+write
    !        --> users can modify the hard-coded parameter "maxUpstreamFile" if desired
    if(tot_upstream > maxUpstreamFile) tot_upstream=0

    ! remove file if it exists
    call system('rm -f '//trim(ancil_dir)//trim(fname_ntopNew))

    ! write data
    call writeData(&
                   ! input
                   trim(ancil_dir)//trim(fname_ntopNew), & ! file name
                   ! input: total elements
                   tot_hru,                              & ! total number of all the upstream hrus for all stream segments
                   tot_upseg,                            & ! total number of immediate upstream segments for all  stream segments
                   tot_upstream,                         & ! total number of all of the upstream stream segments for all stream segments
                   tot_uh,                               & ! total number of unit hydrograph for all stream segments
                   ! input: reach masks
                   ixHRU_desired,                        & ! indices of desired hrus
                   ixSeg_desired,                        & ! indices of desired reaches
                   ! input: data structures
                   structHRU,                            & ! ancillary data for HRUs
                   structSeg,                            & ! ancillary data for stream segments
                   structHRU2seg,                        & ! ancillary data for mapping hru2basin
                   structNTOPO,                          & ! ancillary data for network topology
                   structPFAF,                           & ! ancillary data for pfafstetter code
                   ! output: error control
                   ierr,cmessage)                          ! error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    if (idSegOut>0)      write(iulog,'(2a)') new_line('a'), 'Running in river network subset mode'
    if (ntopAugmentMode) write(iulog,'(2a)') new_line('a'), 'Running in river network augmentation mode'
    write(iulog,'(x,a)') 'Created a new network topology file '//trim(fname_ntopNew)
    write(iulog,'(x,a)') '--> Run again using the new network topology file '
    write(iulog,'(x,a)') 'SUCCESSFUL EXECUTION '
    return
  endif

  ! copy data to the RPARAM and NETOPO structures
  call put_data_struct(nRch_out, structSEG, structNTOPO, & ! input
                       RPARAM, NETOPO,                   & ! output
                       ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! spatial domain decomposition for OMP parallelization
  call omp_domain_decomposition(nRch_out, structNTOPO, river_basin, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE init_ntopo

 ! *****
 ! public subroutine: get mapping data between runoff hru and river network hru
 ! *********************************************************************
 SUBROUTINE init_runoff(remap_flag,      & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                        remap_data_in,   & ! output: data structure to remap data
                        runoff_data_in,  & ! output: data structure for runoff
                        ierr, message)     ! output: error control

 USE public_var,  ONLY : ancil_dir              ! name of the ancillary directory
 USE public_var,  ONLY : input_dir              ! name of the runoff input directory
 USE public_var,  ONLY : fname_qsim             ! name of simulated runoff netCDF
 USE public_var,  ONLY : fname_remap            ! name of runoff mapping netCDF name
 USE public_var,  ONLY : calendar               ! name of calendar
 USE public_var,  ONLY : time_units             ! time units
 USE dataTypes,   ONLY : remap                  ! remapping data type
 USE dataTypes,   ONLY : runoff                 ! runoff data type
 USE read_runoff, ONLY : read_runoff_metadata   ! read meta data from runoff data
 USE read_remap,  ONLY : get_remap_data         ! read remap data
 USE globalData,  ONLY : basinID                ! basin ID

 implicit none
 ! data structures
 logical(lgt),      intent(in)      :: remap_flag       ! logical whether or not runnoff needs to be mapped to river network HRU
 type(remap)  , intent(out)         :: remap_data_in    ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
 type(runoff) , intent(out)         :: runoff_data_in   ! runoff for one time step for all HRUs
 ! error control
 integer(i4b), intent(out)          :: ierr             ! error code
 character(*), intent(out)          :: message          ! error message
 ! local variables
 integer(i8b), allocatable          :: unq_qhru_id(:)
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

 !print*, 'runoff_data_in%nSpace, nTime, trim(time_units) = ', runoff_data_in%nSpace(:), runoff_data_in%nTime, trim(time_units)

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
 integer(i8b), intent(in)  :: qid(:)                       ! ID of input vector
 integer(i8b), intent(in)  :: qidMaster(:)                 ! ID of master vector
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

 !print*, 'rankId = ', rankId(1:10)
 !print*, 'qId( rankId(1:10) ) = ', qId( rankId(1:10) )
 qix(1:size(qid)) = integerMissing
 nx=0
 jx=1
 ! loop through id vector
 do ix=1,size(qid)

  ! find match
  do ixMaster=jx,size(qidMaster) ! normally a very short loop

   ! keep track of trials
   nx=nx+1
   !print*, 'qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) ) = ', qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) )

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
   print*, trim(message)//'matching ids: ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) ) = ', &
                                         ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) )
  endif

 end do  ! looping through the vector

 ! check
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
