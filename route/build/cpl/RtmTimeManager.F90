MODULE RtmTimeManager

  USE ESMF
  USE shr_kind_mod, ONLY: r8 => shr_kind_r8
  USE shr_sys_mod , ONLY: shr_sys_abort, shr_sys_flush
  USE public_var  , ONLY: iulog
  USE public_var  , ONLY: debug
  USE public_var  , ONLY: integerMissing
  USE public_var  , ONLY: realMissing
  USE globalData  , ONLY: masterproc

  implicit none

  private
  public :: init_time                 ! setup startup values
  public :: shr_timeStr

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  ! Input from CESM driver
  integer,  save      :: nelapse               = integerMissing, & ! number of timesteps (or days if negative) to extend a run
                         start_ymd             = integerMissing, & ! starting date for run in yearmmdd format
                         start_tod             = 0,              & ! starting time of day for run in seconds
                         stop_ymd              = integerMissing, & ! stopping date for run in yearmmdd format
                         stop_tod              = 0,              & ! stopping time of day for run in seconds
                         ref_ymd               = integerMissing, & ! reference date for time coordinate in yearmmdd format
                         ref_tod               = 0                 ! reference time of day for time coordinate in seconds
  logical,  save      :: tm_first_restart_step = .false.           ! true for first step of a restart or branch run

CONTAINS

 ! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 SUBROUTINE init_time(ierr, message)

  ! Initialize mizuRoute time based on coupler imported clock

  USE public_var, ONLY: time_units             ! time units (seconds, hours, or days)
  USE public_var, ONLY: simStart               ! date string defining the start of the simulation
  USE public_var, ONLY: simEnd                 ! date string defining the end of the simulation
  USE public_var, ONLY: calendar               ! calendar name
  USE public_var, ONLY: dt                     ! calendar name
  USE globalData, ONLY: timeVar                ! time variables (unit given by runoff data)
  USE globalData, ONLY: sec2tunit              ! time unit conversion: seconds per t_units
  USE globalData, ONLY: iTime                  ! time index at simulation time step
  USE globalData, ONLY: begDatetime            ! datetime: start of routing simulation
  USE globalData, ONLY: endDatetime            ! datetime: end of routing simulation
  USE globalData, ONLY: simDatetime            ! current (0) and previous (1) simulation datetime (yyyy:mm:dd:hh:mm:ss)
  USE datetime_data, ONLY: datetime            ! datetime object definition

  implicit none
  ! Argument variables
  integer,                   intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! Local variables
  integer                                  :: nTime
  type(datetime)                           :: refDatetime      ! datetime: reference
  character(len=256)                       :: t_unit           ! unit of time
  character(len=256)                       :: cmessage         ! error message of downwind routine

  ierr=0; message='init_time/'

  ! get the time multiplier needed to convert time to units of days
  t_unit =  time_units(1:index(time_units,' '))
  select case( trim(t_unit)  )
    case('seconds','second','sec','s'); sec2tunit=1._r8
    case('minutes','minute','min');     sec2tunit=60._r8
    case('hours','hour','hr','h');      sec2tunit=3600._r8
    case('days','day','d');             sec2tunit=86400._r8
    case default
      ierr=20; message=trim(message)//'<time_units>= '//trim(time_units)//': <time_units> must be seconds, minutes, hours or days.'; return
  end select

  ! obtain reference, simulation start and end datetimes
  call refDatetime%str2datetime(time_units, calendar, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refDatetime]'; return; endif

  call begDatetime%str2datetime(simStart, calendar, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startDatetime]'; return; endif

  call endDatetime%str2datetime(simEnd, calendar, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endDatetime]'; return; endif

  ! number of time step from starting time to simulation end time
  nTime = int((endDatetime - begDatetime)/dt)

  ! Initialize timeVar : model time (time step endpoints) in model time unit (t_unit), used for time output
  timeVar(1) = begDatetime - refDatetime  ! second since reference datetime
  timeVar(2) = timeVar(1) + dt

  ! check that the dates are aligned
  if(endDatetime < begDatetime) then; ierr=20; message=trim(message)//'simulation end is before simulation start'; return; endif

  ! initialize model time at first time step (1) and previous time step (0)
  iTime = 1
  simDatetime(0) = datetime(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
  simDatetime(1) = begDatetime
  simDatetime(2) = simDatetime(1)%add_sec(dt, calendar, ierr, cmessage)

  if (masterproc .and. debug) then
    write(iulog,*) 'simStart datetime     = ', trim(simStart)
    write(iulog,*) 'simEnd   datetime     = ', trim(simEnd)
    write(iulog,*) 'reference datetime    = ', refDatetime%year(), refDatetime%month(), refDatetime%day(), refDatetime%hour(), refDatetime%minute(), refDatetime%sec()
    write(iulog,*) 'dt [sec]              = ', dt
    write(iulog,*) 'nTime                 = ', nTime
    write(iulog,*) 'iTime                 = ', iTime
    write(iulog,*) 'timeVar(1),timeVar(2) = ', timeVar(1), timeVar(2)
    call shr_sys_flush(iulog)
  end if

 END SUBROUTINE init_time


 ! Public subroutine:
 SUBROUTINE shr_timeStr(esmfTime, timeStr)

   implicit none

   ! Arguments
   type(ESMF_Time), intent(in)  :: esmfTime
   character(len=*),intent(out) :: timeStr
   ! Local variables
   integer                     :: yy,mm,dd              ! Temporaries for time query
   integer                     :: hr,mn,sec             ! Temporaries for time query
   integer                      :: rc

   call ESMF_TimeGet(esmfTime , yy=yy, mm=mm, dd=dd, h=hr, m=mn, s=sec, rc=rc )

   write(timeStr,'(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)'), yy,'-',mm,'-',dd,' ',hr,':',mn,':',sec

 END SUBROUTINE shr_timeStr


END MODULE RtmTimeManager
