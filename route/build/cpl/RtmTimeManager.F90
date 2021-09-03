MODULE RtmTimeManager

  USE ESMF
  USE shr_kind_mod, ONLY: r8 => shr_kind_r8
  USE shr_sys_mod , ONLY: shr_sys_abort, shr_sys_flush
  USE public_var  , ONLY: iulog
  USE public_var  , ONLY: integerMissing
  USE public_var  , ONLY: realMissing
  USE globalData  , ONLY: masterproc

  implicit none
  private ! except

  public :: init_time                 ! setup startup values
  public :: shr_timeStr

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  logical,  parameter :: debug_write = .true.

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

  ! Shared data
  USE public_var, ONLY: time_units             ! time units (seconds, hours, or days)
  USE public_var, ONLY: simStart               ! date string defining the start of the simulation
  USE public_var, ONLY: simEnd                 ! date string defining the end of the simulation
  USE public_var, ONLY: calendar               ! calendar name
  USE public_var, ONLY: dt                     ! calendar name
  USE public_var, ONLY: secprday               ! number of second per a day
  USE globalData, ONLY: timeVar                ! time variables (unit given by runoff data)
  USE globalData, ONLY: iTime                  ! time index at simulation time step
  USE globalData, ONLY: begDatetime            ! datetime: start of routing simulation
  USE globalData, ONLY: endDatetime            ! datetime: end of routing simulation
  USE globalData, ONLY: simDatetime            ! current (0) and previous (1) simulation datetime (yyyy:mm:dd:hh:mm:ss)
  USE datetime_data, ONLY: datetime            ! datetime object definition
  implicit none

  ! input:
  ! output: error control
  integer,                   intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer                                  :: nTime
  integer                                  :: ix
  real(r8)                                 :: dt_day           ! simulation time step in day
  real(r8)                                 :: timePerDay       ! number of time-unit per a day. time-unit is from t_unit
  real(r8)                                 :: secPerTime       ! number of sec per time-unit. time-unit is from t_unit
  real(r8)                                 :: refJulday        ! reference julian day
  real(r8)                                 :: begJulday        ! simulation start julian day
  real(r8)                                 :: endJulday        ! simulation end julian day
  type(datetime)                           :: refDatetime      ! datetime: reference
  character(len=256)                       :: t_unit           ! unit of time
  character(len=256)                       :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_time/'

  ! get the time multiplier needed to convert time to units of days
  t_unit =  time_units(1:index(time_units,' '))
  select case( trim(t_unit)  )
    case('seconds','second','sec','s'); secPerTime=1._r8;     timePerDay=86400._r8
    case('minutes','minute','min');     secPerTime=60._r8;    timePerDay=1440._r8
    case('hours','hour','hr','h');      secPerTime=3600._r8;  timePerDay=24._r8
    case('days','day','d');             secPerTime=86400._r8; timePerDay=1._r8
    case default
      ierr=20; message=trim(message)//'<time_units>= '//trim(time_units)//': <time_units> must be seconds, minutes, hours or days.'; return
  end select

  dt_day = dt/secprday  ! dt [sec] -> dt_day

  ! obtain reference, simulation start and end datetimes
  call refDatetime%str2datetime(time_units, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refDatetime]'; return; endif

  call begDatetime%str2datetime(simStart, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startDatetime]'; return; endif

  call endDatetime%str2datetime(simEnd, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endDatetime]'; return; endif

  ! obtain reference, simulation start and end julian days
  call refDatetime%julianday(calendar, refJulday, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refDatetime]'; return; endif

  call begDatetime%julianday(calendar, begJulday, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [begDatetime]'; return; endif

  call endDatetime%julianday(calendar, endJulday, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [begDatetime]'; return; endif

  ! number of time step from reference time to simulation end time
  nTime = int(endJulday - begJulday/dt_day) + 1

  ! Create timeVar array: starting with 0 and increment of model time step in model unit (t_unit)
  allocate(timeVar(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  timeVar(1) = (begJulday - refJulday)*timePerDay
  if (nTime>1) then
    do ix = 2, nTime
      timeVar(ix) = timeVar(ix-1) + dt/secPerTime
    end do
  end if

  ! check that the dates are aligned
  if(endDatetime < begDatetime) then; ierr=20; message=trim(message)//'simulation end is before simulation start'; return; endif

  ! initialize model time at first time step (1) and previous time step (0)
  iTime = 1
  simDatetime(0) = datetime(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
  simDatetime(1) = begDatetime

  if (masterproc .and. debug_write) then
    write(iulog,*) 'simStart datetime     = ', trim(simStart)
    write(iulog,*) 'simEnd   datetime     = ', trim(simEnd)
    write(iulog,*) 'reference datetime    = ', refDatetime%year(), refDatetime%month(), refDatetime%day(), refDatetime%hour(), refDatetime%minute(), refDatetime%sec()
    write(iulog,*) 'simDatetime           = ', simDatetime(1)%year(), simDatetime(1)%month(), simDatetime(1)%day(), simDatetime(1)%hour(), simDatetime(1)%minute(), simDatetime(1)%sec()
    write(iulog,*) 'dt [sec]              = ', dt
    write(iulog,*) 'nTime                 = ', nTime
    write(iulog,*) 'iTime, timeVar(iTime) = ', iTime, timeVar(iTime)
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
