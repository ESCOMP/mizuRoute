MODULE RtmTimeManager

  USE ESMF
  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE shr_sys_mod , only: shr_sys_abort, shr_sys_flush
  USE public_var  , only: iulog
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

  ! subroutines:
  USE process_time_module, ONLY: process_time  ! process time information
  ! Shared data
  USE dataTypes,  ONLY: time                   ! time data type
  USE public_var, ONLY: time_units             ! time units (seconds, hours, or days)
  USE public_var, ONLY: simStart               ! date string defining the start of the simulation
  USE public_var, ONLY: simEnd                 ! date string defining the end of the simulation
  USE public_var, ONLY: calendar               ! calendar name
  USE globalData, ONLY: timeVar                ! time variables (unit given by runoff data)
  USE globalData, ONLY: roJulday               !
  USE globalData, ONLY: iTime                  ! time index at simulation time step
  USE globalData, ONLY: refJulday              ! julian day: reference
  USE globalData, ONLY: roJulday               ! julian day: runoff input time
  USE globalData, ONLY: startJulday            ! julian day: start of routing simulation
  USE globalData, ONLY: endJulday              ! julian day: end of routing simulation
  USE globalData, ONLY: modJulday              ! julian day: at model time step
  USE globalData, ONLY: modTime                ! model time data (yyyy:mm:dd:hh:mm:ss)

  implicit none

  ! input:
  ! output: error control
  integer,                   intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer                                  :: nTime
  integer                                  :: ix
  real(r8)                                 :: convTime2Days
  character(len=256)                       :: t_unit           ! unit of time
  character(len=256)                       :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_time/'

  ! get the time multiplier needed to convert time to units of days
  t_unit =  time_units(1:index(time_units,' '))
  select case( trim(t_unit)  )
    case('seconds','second','sec','s'); convTime2Days=86400._r8
    case('minutes','minute','min');     convTime2Days=1440._r8
    case('hours','hour','hr','h');      convTime2Days=24._r8
    case('days','day','d');             convTime2Days=1._r8
    case default
      ierr=20; message=trim(message)//'<time_units>= '//trim(time_units)//': <time_units> must be seconds, minutes, hours or days.'; return
  end select

  ! extract time information from the control information
  call process_time(trim(time_units), calendar, refJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refJulday]'; return; endif
  call process_time(trim(simStart),calendar, startJulday, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startJulday]'; return; endif
  call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endJulday]'; return; endif

  nTime = int((endJulday - refJulday)*convTime2Days) + 1

  allocate(timeVar(nTime), roJulday(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! Create timeVar array: starting with 0 and increment of model time step in model unit
  timeVar(1) = (startJulday - refJulday)*convTime2Days
  if (nTime>1) then
    do ix = 2, nTime
      timeVar(ix) = timeVar(ix-1) + convTime2Days
    end do
  end if

  do ix = 1, nTime
    roJulday(ix) = refJulday + timeVar(ix)/convTime2Days
  end do

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

  if (masterproc .and. debug_write) then
    write(iulog,*) 'simStart (startJulDay) = ', trim(simStart),'(',startJulday,')'
    write(iulog,*) 'simEnd   (endJulDay)   = ', trim(simEnd),  '(',endJulday,  ')'
    write(iulog,*) 'refJulDay              = ', refJulday
    write(iulog,*) 'modJulDay              = ', modJulday
    write(iulog,*) 'nTime                  = ', nTime
    write(iulog,*) 'iTime, timeVar(iTime)  = ', iTime, timeVar(iTime)
    call shr_sys_flush(iulog)
  end if

  ! initialize previous model time
  modTime(0) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

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
