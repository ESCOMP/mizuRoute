MODULE datetime_data

! datetime class
!
! usage:
! initialize datetime
!   datetime1 = datetime(1969, 7,16, 7,32,  0.0)
!

USE nrtype
USE public_var, ONLY: realMissing, integerMissing
USE public_var, ONLY: secprday, secprhour, secprmin
USE public_var, ONLY: months_per_yr, days_per_yr
USE public_var, ONLY: hr_per_day, min_per_hour       ! hours per day and minutes per hour

implicit none

type, public :: datetime

  private
  integer(i4b)           :: iy       = integerMissing
  integer(i4b)           :: im       = integerMissing
  integer(i4b)           :: id       = integerMissing
  integer(i4b)           :: ih       = integerMissing
  integer(i4b)           :: imin     = integerMissing
  real(dp)               :: dsec     = realMissing
  character(20)          :: calendar = 'empty'

CONTAINS

  procedure, public :: jul2datetime  => sub_jul2datetime
  procedure, public :: str2datetime  => sub_str2datetime
  procedure, public :: year          => fn_get_year
  procedure, public :: month         => fn_get_month
  procedure, public :: day           => fn_get_day
  procedure, public :: hour          => fn_get_hour
  procedure, public :: minute        => fn_get_min
  procedure, public :: sec           => fn_get_sec
  procedure, public :: cal           => fn_get_calendar
  procedure, public :: is_leap_year  => fn_is_leap_year
  procedure, public :: ndays_month   => fn_ndays_month
  procedure, public :: julianday     => sub_julian_day
  procedure, public :: dayofyear     => fn_dayofyear
  procedure, public :: add_mon       => fn_add_months
  procedure, public :: add_day       => fn_add_days
  procedure, public :: add_hr        => fn_add_hours
  procedure, public :: add_sec       => fn_add_sec
  procedure, public :: is_equal_mon  => fn_is_equal_month
  procedure, public :: is_equal_day  => fn_is_equal_day
  procedure, public :: is_equal_time => fn_is_equal_time

  procedure, private :: fn_is_equal
  procedure, private :: fn_is_gt
  procedure, private :: fn_is_ge
  procedure, private :: fn_is_lt
  procedure, private :: fn_is_le
  procedure, private :: sub_assign
  procedure, private :: fn_delta
  procedure, private :: fn_delta_fast
  generic            :: operator(==)  => fn_is_equal
  generic            :: operator(>)   => fn_is_gt
  generic            :: operator(>=)  => fn_is_ge
  generic            :: operator(<)   => fn_is_lt
  generic            :: operator(<=)  => fn_is_le
  generic            :: assignment(=) => sub_assign
  generic            :: operator(-)   => fn_delta
  generic            :: operator(.minus.) => fn_delta_fast ! this works only standard calendar (but faster than fn_delta)

end type datetime

private :: sub_jul2datetime, sub_str2datetime
private :: fn_get_year, fn_get_month, fn_get_day, fn_get_hour, fn_get_min, fn_get_sec
private :: fn_get_calendar
private :: fn_dayofyear
private :: fn_is_leap_year, sub_julian_day, fn_ndays_month
private :: fn_add_months, fn_add_days, fn_add_hours, fn_add_sec
private :: fn_is_equal_month, fn_is_equal_day, fn_is_equal_time

INTERFACE datetime
  module procedure constructor
END INTERFACE datetime

CONTAINS

  FUNCTION constructor(iyr, imo, ida, ihr, imin, dsec, calendar) RESULT(instDatetime)

    implicit none
    type(datetime)                     :: instDatetime
    integer(i4b),           intent(in) :: iyr
    integer(i4b),           intent(in) :: imo
    integer(i4b),           intent(in) :: ida
    integer(i4b),           intent(in) :: ihr
    integer(i4b),           intent(in) :: imin
    real(dp),               intent(in) :: dsec
    character(*), optional, intent(in) :: calendar

    instDatetime%iy       = iyr
    instDatetime%im       = imo
    instDatetime%id       = ida
    instDatetime%ih       = ihr
    instDatetime%imin     = imin
    instDatetime%dsec     = dsec
    instDatetime%calendar = 'standard'
    if (present(calendar)) instDatetime%calendar = calendar

  END FUNCTION constructor

  SUBROUTINE sub_jul2datetime(this, julday, calendar, ierr, message)

    implicit none
    class(datetime)                    :: this
    real(dp),             intent(in)   :: julday
    character(*),         intent(in)   :: calendar
    integer(i4b),         intent(out)  :: ierr
    character(len=strLen),intent(out)  :: message
    ! local variable
    character(len=strLen)              :: cmessage

    ierr=0; message='sub_jul2datetime/'

    select case(trim(calendar))
     case ('noleap','365_day')
      call compCalday_noleap(julday,this%iy,this%im,this%id,this%ih,this%imin,this%dsec, ierr, cmessage)
     case ('standard','gregorian','proleptic_gregorian')
      call compCalday(julday,this%iy,this%im,this%id,this%ih,this%imin,this%dsec, ierr, cmessage)
     case default; ierr=20; message=trim(message)//'calendar name: '//trim(calendar)//' invalid'; return
    end select
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    this%calendar = calendar

  END SUBROUTINE sub_jul2datetime

  SUBROUTINE sub_str2datetime(this, str, calendar, ierr, message)

    implicit none
    class(datetime)                    :: this
    character(*),         intent(in)   :: str
    character(*),         intent(in)   :: calendar
    integer(i4b),         intent(out)  :: ierr
    character(len=strLen),intent(out)  :: message
    ! local variable
    character(len=strLen)              :: cmessage

    ierr=0; message='sub_str2datetime/'

    this%calendar = calendar
    call extractTime(str,this%iy,this%im,this%id,this%ih,this%imin,this%dsec, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  END SUBROUTINE sub_str2datetime

  SUBROUTINE sub_assign(this, that)
    implicit none
    class(datetime),intent(out)  :: this
    type(datetime), intent(in)   :: that
    this%iy   = that%iy
    this%im   = that%im
    this%id   = that%id
    this%ih   = that%ih
    this%imin = that%imin
    this%dsec = that%dsec
    this%calendar = that%calendar
  END SUBROUTINE sub_assign


  pure elemental integer(i4b) FUNCTION fn_get_year(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_year = this%iy
  END FUNCTION fn_get_year

  pure elemental integer(i4b) FUNCTION fn_get_month(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_month = this%im
  END FUNCTION fn_get_month

  pure elemental integer(i4b) FUNCTION fn_get_day(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_day = this%id
  END FUNCTION fn_get_day

  pure elemental integer(i4b) FUNCTION fn_get_hour(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_hour = this%ih
  END FUNCTION fn_get_hour

  pure elemental integer(i4b) FUNCTION fn_get_min(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_min = this%imin
  END FUNCTION fn_get_min

  pure elemental real(dp) FUNCTION fn_get_sec(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_sec = this%dsec
  END FUNCTION fn_get_sec

  pure elemental character(20) FUNCTION fn_get_calendar(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_calendar = this%calendar
  END FUNCTION fn_get_calendar


  integer(i4b) FUNCTION fn_dayofyear(this)
    ! Returns the integer day of the year (ordinal date).
    implicit none
    class(datetime), intent(in) :: this
    type(datetime)              :: jan1datetime

    jan1datetime = datetime(this%iy, 1, 1, 0, 0, 0.0, calendar=this%calendar)
    fn_dayofyear = int((this - jan1datetime)/secprday) + 1

  END FUNCTION fn_dayofyear


  pure elemental logical(lgt) FUNCTION fn_is_leap_year(this)
    implicit none
    class(datetime), intent(in)  :: this
    if (mod(this%iy, 4) == 0) then
      if (mod(this%iy, 100) == 0) then
        if (mod(this%iy, 400) == 0) then
          fn_is_leap_year = .True.
        else
          fn_is_leap_year = .False.
        end if
      else
        fn_is_leap_year = .True.
      end if
    else
      fn_is_leap_year = .False.
    end if
  END FUNCTION fn_is_leap_year

  pure elemental integer(i4b) FUNCTION fn_ndays_month(this)
    ! return total numbers of days in current month
    implicit none
    class(datetime),      intent(in)  :: this
    ! local variables
    integer(i4b)                      :: nmonths(12)

    select case(trim(this%calendar))
      case ('standard','gregorian','proleptic_gregorian')
        if ( fn_is_leap_year(this)) then
          nmonths = [31,29,31,30,31,30,31,31,30,31,30,31]
        else
          nmonths = [31,28,31,30,31,30,31,31,30,31,30,31]
        end if
      case('noleap')
        nmonths = [31,28,31,30,31,30,31,31,30,31,30,31]
      case default; nmonths=integerMissing
    end select
    fn_ndays_month = nmonths(this%im)
  END FUNCTION fn_ndays_month


  SUBROUTINE sub_julian_day(this, julianDay, ierr, message)
    implicit none
    class(datetime),      intent(in)  :: this
    real(dp),             intent(out) :: julianDay
    integer(i4b),         intent(out) :: ierr
    character(len=strLen),intent(out) :: message
    ! local variables
    character(len=strLen)             :: cmessage

    ierr=0; message='sub_julian_day/'
    select case(trim(this%calendar))
      case ('standard','gregorian','proleptic_gregorian')
        call compJulday(this%iy, this%im, this%id, this%ih, this%imin, this%dsec, julianDay, ierr, cmessage)
      case('noleap')
        call compJulday_noleap(this%iy, this%im, this%id, this%ih, this%imin, this%dsec, julianDay, ierr, cmessage)
      case default; ierr=20; message=trim(message)//'calendar name: '//trim(this%calendar)//' invalid'; return
    end select
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  END SUBROUTINE sub_julian_day

  ! -----
  ! increment(+)/decrement(+) datetime
  ! -----

  type(datetime) FUNCTION fn_add_months(this, nmonths)
    implicit none
    class(datetime),     intent(in)    :: this
    integer(i4b),        intent(in)    :: nmonths
    fn_add_months%iy   = this%iy + (nmonths/12_i4b)
    fn_add_months%im   = this%im + mod(nmonths,12)
    fn_add_months%id   = this%id
    fn_add_months%ih   = this%ih
    fn_add_months%imin = this%imin
    fn_add_months%dsec = this%dsec
  END FUNCTION fn_add_months

  type(datetime) FUNCTION fn_add_days(this, days, ierr, message)
    implicit none
    class(datetime),      intent(in)    :: this
    integer(i4b),         intent(in)    :: days
    integer(i4b),         intent(out)   :: ierr
    character(len=strLen),intent(out)   :: message
    ! local variables
    real(dp)                            :: julday
    character(len=strLen)               :: cmessage

    ierr=0; message='fn_add_days/'

    call sub_julian_day(this, julday, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

    julday = julday + real(days,dp)
    call sub_jul2datetime(fn_add_days, julday, this%calendar, ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

  END FUNCTION fn_add_days

  type(datetime) FUNCTION fn_add_hours(this, hrs, ierr, message)
    implicit none
    class(datetime),      intent(in)  :: this
    integer(i4b),         intent(in)  :: hrs
    integer(i4b),         intent(out) :: ierr
    character(len=strLen),intent(out) :: message
    ! local variables
    real(dp)                          :: julday
    character(len=strLen)             :: cmessage

    ierr=0; message='fn_add_hours/'
    call sub_julian_day(this, julday, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

    julday = julday + real(hrs,dp)/hr_per_day
    call sub_jul2datetime(fn_add_hours, julday, this%calendar, ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

  END FUNCTION fn_add_hours

  type(datetime) FUNCTION fn_add_sec(this, sec, ierr, message)
    implicit none
    class(datetime),      intent(in)   :: this
    real(dp),             intent(in)   :: sec
    integer(i4b),         intent(out)  :: ierr
    character(len=strLen),intent(out)  :: message
    ! local variables
    real(dp)                          :: julday
    character(len=strLen)             :: cmessage

    ierr=0; message='fn_add_sec/'
    call sub_julian_day(this, julday, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

    julday = julday + sec/secprday
    call sub_jul2datetime(fn_add_sec, julday, this%calendar, ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

  END FUNCTION fn_add_sec

  ! -----
  ! check relational logic between two datetime
  ! -----

  logical(lgt) FUNCTION fn_is_equal_month(this, that)
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    fn_is_equal_month = (this%im==that%im)
  END FUNCTION fn_is_equal_month

  logical(lgt) FUNCTION fn_is_equal_day(this, that)
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    fn_is_equal_day = (this%id==that%id)
  END FUNCTION fn_is_equal_day

  logical(lgt) FUNCTION fn_is_equal_time(this, that)
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    fn_is_equal_time = (this%ih==that%ih .and. this%imin==that%imin .and. abs(this%dsec-that%dsec)<=epsilon(this%dsec))
  END FUNCTION fn_is_equal_time

  logical(lgt) FUNCTION fn_is_equal(this, that)
    ! if this == that T, otherwise F
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    fn_is_equal = (this%iy==that%iy .and. this%im==that%im .and. this%id==that%id .and. &
                   this%ih==that%ih .and. this%imin==that%imin .and. abs(this%dsec-that%dsec)<=epsilon(this%dsec))
  END FUNCTION fn_is_equal

  logical(lgt) FUNCTION fn_is_gt(this, that)
    ! if this > that T, otherwise F
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    if (this%iy > that%iy) then
      fn_is_gt = .true.
    else if (this%iy < that%iy) then
      fn_is_gt = .false.
    else
      if (this%im > that%im) then
        fn_is_gt = .true.
      else if (this%im < that%im) then
        fn_is_gt = .false.
      else
        if (this%id > that%id) then
          fn_is_gt = .true.
        else if (this%id < that%id) then
          fn_is_gt = .false.
        else
          if (this%ih > that%ih) then
            fn_is_gt = .true.
          else if (this%ih < that%ih) then
            fn_is_gt = .false.
          else
            if (this%imin > that%imin) then
              fn_is_gt = .true.
            else if (this%imin < that%imin) then
              fn_is_gt = .false.
            else
              if (this%dsec > that%dsec) then
                fn_is_gt = .true.
              else
                fn_is_gt = .false.
              endif ! dsec
            endif ! min
          endif ! hr
        endif ! day
      endif ! mon
    endif ! yr
  END FUNCTION fn_is_gt

  logical(lgt) FUNCTION fn_is_ge(this, that)
    ! if this >= that T, otherwise F
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    fn_is_ge = .false.
    if ( fn_is_gt(this, that) .or. fn_is_equal(this, that) ) then
      fn_is_ge = .true.
    end if
  END FUNCTION fn_is_ge

  logical(lgt) FUNCTION fn_is_lt(this, that)
    ! if this < that T, otherwise F
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    fn_is_lt = .false.
    if (.not.( fn_is_ge(this, that))) then
      fn_is_lt = .true.
    end if
  END FUNCTION fn_is_lt

  logical(lgt) FUNCTION fn_is_le(this, that)
    ! if this <= that T, otherwise F
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    fn_is_le = .false.
    if (.not.(fn_is_gt(this, that))) then
      fn_is_le = .true.
    end if
  END FUNCTION fn_is_le

  real(dp) FUNCTION fn_delta(this, that) result(ans)
    ! this - that give in time difference in second
    implicit none
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    ! local variables
    real(dp)                         :: sec_in_day1
    real(dp)                         :: sec_in_day2
    real(dp)                         :: diffDay
    integer(i4b)                     :: ierr
    type(datetime)                   :: date1
    type(datetime)                   :: date2
    real(dp)                         :: julday1
    real(dp)                         :: julday2
    character(len=strLen)            :: cmessage

    sec_in_day1 = real(this%ih*secprhour +this%imin*secprmin, kind=dp) +this%dsec
    sec_in_day2 = real(that%ih*secprhour +that%imin*secprmin, kind=dp) +that%dsec

    date1 = datetime(this%iy, this%im, this%id, 0, 0, 0.0, calendar=this%calendar)
    date2 = datetime(that%iy, that%im, that%id, 0, 0, 0.0, calendar=that%calendar)

    call date1%julianday(julday1, ierr, cmessage)
    call date2%julianday(julday2, ierr, cmessage)

    diffDay=julday1-julday2

    ans = diffDay*secprday + sec_in_day1 - sec_in_day2
  END FUNCTION fn_delta

  real(dp) FUNCTION fn_delta_fast(this, that) result(elapsedSec)
    implicit none
    ! Argument variables
    class(datetime),     intent(in)  :: this
    class(datetime),     intent(in)  :: that
    ! local variables
    integer(i4b)                   :: elapsedDay                         ! elapsed full days
    integer(i4b)                   :: yy                                 ! index of year
    ! number of days of each month
    integer(i4b)                   :: days1(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer(i4b)                   :: days2(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    ! calculate the elapsed time smaller than a day
    elapsedSec = (this%dsec-that%dsec) + (this%imin-that%imin)*secprmin + (this%ih-that%ih)*secprhour

    ! check if the run is within the same day otherwise calculate how many days
    if (this%iy>that%iy .or. this%im>that%im .or. this%id>that%id) then
      elapsedDay = 0
      ! diffenece in year
      do yy = that%iy, this%iy-1
        elapsedDay = elapsedDay + 365
        if ((mod(yy,4)==0 .and. .not. mod(yy,100)==0) .or. (mod(yy,400)==0)) elapsedDay = elapsedDay + 1
      end do
      if ((mod(that%iy,4)==0 .and. .not. mod(that%iy,100)==0) .or. (mod(that%iy,400)==0)) days1(2) = 29
      if ((mod(this%iy,4)==0 .and. .not. mod(this%iy,100)==0) .or. (mod(this%iy,400)==0)) days2(2) = 29
      ! difference in month
      if (that%im > 1) elapsedDay = elapsedDay - sum(days1(1:(that%im-1)))
      elapsedDay = elapsedDay - that%id
      ! difference in day
      if (this%im > 1) elapsedDay = elapsedDay + sum(days2(1:(this%im-1)))
      elapsedDay = elapsedDay + this%id
      ! convert to seconds
      elapsedSec = elapsedSec + elapsedDay*secprday
    end if
  END FUNCTION

  ! -------------------------------
  ! no type-bound procedure
  ! -------------------------------
  ! ***************************************************************************************
  ! public subroutine compJulday: convert datetime to julian day (units of days)
  ! ***************************************************************************************
  SUBROUTINE compJulday(iyyy,mm,id,ih,imin,dsec,&  ! input
                        juldayss,err,message)      ! output
    implicit none
    ! Argument variables
    integer(i4b), intent(in)   :: iyyy,mm,id   ! year, month, day
    integer(i4b), intent(in)   :: ih,imin      ! hour, minute
    real(dp),     intent(in)   :: dsec         ! seconds
    real(dp),     intent(out)  :: juldayss
    integer(i4b), intent(out)  :: err          ! error code
    character(*), intent(out)  :: message      ! error message
    ! local variables
    integer(i4b)               :: julday       ! julian day
    integer(i4b),parameter     :: igreg=15+31*(10+12*1582)  !IGREG = 588829
    integer(i4b)               :: ja,jm,jy
    real(dp)                   :: jfrac        ! fraction of julian day

    err=0; message="juldayss"

    ! compute julian day
    jy=iyyy
    if (jy.eq.0) then; err=10; message=trim(message)//"noYearZero/"; return; end if
    if (jy.lt.0) jy=jy+1
    if (mm.gt.2) then
      jm=mm+1
    else
      jy=jy-1
      jm=mm+13
    end if
    julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
    if (id+31*(mm+12*iyyy).ge.IGREG) then
      ja=int(0.01*jy)
      julday=julday+2-ja+int(0.25*ja)
    end if

    ! compute fraction of the day
    jfrac = (real(ih,kind(dp))*secprhour + real(imin,kind(dp))*secprmin + dsec) / secprday

    ! and return the julian day, expressed in fraction of a day
    juldayss = real(julday,kind(dp)) + jfrac

  END SUBROUTINE compJulday

  ! ***************************************************************************************
  ! public subroutine compJulday: convert datetime to julian day (units of days) for noleap calendar
  ! reference: https://github.com/nmizukami/VIC/blob/VIC.5.0.0/vic/drivers/shared_all/src/vic_time.c
  ! ***************************************************************************************
  SUBROUTINE compJulday_noleap(iyyy,mm,id,ih,imin,dsec,&   ! input
                               juldayss,err,message)       ! output
    implicit none
    ! Argument variables
    integer(i4b), intent(in)   :: iyyy,mm,id   ! year, month, day
    integer(i4b), intent(in)   :: ih,imin      ! hour, minute
    real(dp),     intent(in)   :: dsec         ! seconds
    real(dp),     intent(out)  :: juldayss
    integer(i4b), intent(out)  :: err          ! error code
    character(*), intent(out)  :: message      ! error message
    ! local variables
    integer(i4b)               :: mm_tmp       ! adjusted month
    integer(i4b)               :: iyyy_tmp     ! adjusted year
    real(dp)                   :: dfrac        ! fraction of day

    err=0; message="compJulday_noleap"

    ! compute fraction of the day
    dfrac = real(id,kind(dp))+(real(ih,kind(dp))*secprhour + real(imin,kind(dp))*secprmin + dsec) / secprday

    ! compute julian day
    ! Start Meeus algorithm (variables are in his notation)
    mm_tmp = mm; iyyy_tmp = iyyy
    if (mm < 3) then
      mm_tmp = mm + months_per_yr
      iyyy_tmp = iyyy - 1
    end if

    juldayss = real(days_per_yr*(iyyy_tmp + 4716)) + &
               real(floor(30.6001 * real(mm_tmp+1))) + dfrac - 1524.5

  END SUBROUTINE compJulday_noleap

  ! ***************************************************************************************
  ! public subroutine compgregcal: convert julian day (units of days) to calendar date
  ! source: https://en.wikipedia.org/wiki/Julian_day#Julian_or_Gregorian_calendar_from_Julian_day_number
  ! ***************************************************************************************
  subroutine compCalday(julday,                              & !input
                        iyyy,mm,id,ih,imin,dsec,err,message)   !output
    implicit none
    ! Argument variables
    real(dp), intent(in)          :: julday       ! julian day
    integer(i4b), intent(out)     :: iyyy         ! year
    integer(i4b), intent(out)     :: mm           ! month
    integer(i4b), intent(out)     :: id           ! day
    integer(i4b), intent(out)     :: ih           ! hour
    integer(i4b), intent(out)     :: imin         ! minute
    real(dp),     intent(out)     :: dsec         ! seconds
    integer(i4b), intent(out)     :: err          ! error code
    character(*), intent(out)     :: message      ! error message
    ! local parameters
    integer(i4b),parameter       :: y = 4716
    integer(i4b),parameter       :: j = 1401
    integer(i4b),parameter       :: m = 2
    integer(i4b),parameter       :: n = 12
    integer(i4b),parameter       :: r = 4
    integer(i4b),parameter       :: p = 1461
    integer(i4b),parameter       :: v = 3
    integer(i4b),parameter       :: u = 5
    integer(i4b),parameter       :: s = 153
    integer(i4b),parameter       :: w = 2
    integer(i4b),parameter       :: b = 274277
    integer(i4b),parameter       :: c = -38
    ! local variables
    integer(i4b)                 :: f,e,g,h                          ! various step variables from wikipedia
    integer(i4b)                 :: step_1a,step_1b,step_1c,step_1d  ! temporary variables for calendar calculations
    real(dp)                     :: frac_day                         ! fractional day
    real(dp)                     :: remainder                        ! remainder of modulus operation

    err=0; message="compCalday"

    if(julday<=0)then;err=10;message=trim(message)//"no negative julian days/"; return; end if

    ! step 1
    step_1a = 4*int(julday)+b
    step_1b = step_1a/146097
    step_1c = step_1b*3
    step_1d = step_1c/4

    f = int(julday)+j+step_1d+c

    ! step 2
    e = r * f + v

    ! step 3
    g = mod(e,p)/r

    ! step 4
    h = u * g + w

    ! find day
    id = (mod(h,s))/u + 1

    ! find month
    mm = mod(h/s+m,n)+1

    ! find year
    iyyy = (e/p)-y + (n+m-mm)/n

    ! now find hour,min,second

    frac_day = julday - floor(julday)
    ih = floor((frac_day+1e-9)*hr_per_day)

    remainder = (frac_day+1e-9)*hr_per_day - ih
    imin = floor(remainder*min_per_hour)

    remainder = remainder*min_per_hour - imin
    dsec = nint(remainder*secprmin)

  END SUBROUTINE compCalday

  ! ***************************************************************************************
  ! public subroutine compgregcal_noleap: compute yy,mm,dd,hr,min,hr from a noleap julian day
  ! source: https://github.com/nmizukami/VIC/blob/VIC.5.0.0/vic/drivers/shared_all/src/vic_time.c
  ! ***************************************************************************************
  SUBROUTINE compCalday_noleap(julday,                              & !input
                               iyyy,mm,id,ih,imin,dsec,err,message)   !output
    implicit none
    ! Argument variables
    real(dp), intent(in)         :: julday       ! julian day
    integer(i4b), intent(out)    :: iyyy         ! year
    integer(i4b), intent(out)    :: mm           ! month
    integer(i4b), intent(out)    :: id           ! day
    integer(i4b), intent(out)    :: ih           ! hour
    integer(i4b), intent(out)    :: imin         ! minute
    real(dp),     intent(out)    :: dsec         ! seconds
    integer(i4b), intent(out)    :: err          ! error code
    character(*), intent(out)    :: message      ! error message
    ! local variables
    integer(i4b)                 :: A,B,C,D,E    ! various step variable
    integer(i4b)                 :: nday         ! various step variable
    real(dp)                     :: F            ! various step variable
    integer(i4b)                 :: dayofyr      ! day of year
    real(dp)                     :: frac_day     ! fractional day
    real(dp)                     :: days         ! day with a fraction days
    real(dp)                     :: remainder    ! remainder of modulus operation

    err=0; message="compCalday_noleap"

    if(julday<=0)then;err=10;message=trim(message)//"no negative julian days/"; return; end if

    A = floor(julday+0.5_dp)
    F = julday + 0.5_dp - real(A,kind(dp))
    B = A + 1524_i4b
    C = int((real(B,kind(dp)) - 122.1_dp)/real(days_per_yr,kind(dp)))
    D = days_per_yr*C
    E = int(real(B-D,kind(dp))/30.6001_dp)

    ! compute day
    days = real(B-D - int(30.6001_dp*real(E,kind(dp))),kind(dp)) + F
    id = floor(days)

    ! compute day in a year
    nday = B - D - 123_i4b
    if (nday <= 305_i4b) then
       dayofyr = nday + 60_i4b
    else
       dayofyr = nday - 305_i4b
    endif

    ! compute month
    if (E < 14_i4b) then
      mm = E - 1_i4b
    else
      mm = E - 13_i4b
    endif

    ! compute year
    if (mm > 2_i4b) then
     iyyy = C - 4716_i4b
    else
     iyyy = C - 4715_i4b
    endif

    ! Convert fractions of a day to time
    ! now find hour,min,second
    frac_day = days - real(id, kind(dp))
    ih = floor((frac_day+1e-9)*hr_per_day)

    remainder = (frac_day+1e-9)*hr_per_day - ih
    imin = floor(remainder*min_per_hour)

    remainder = remainder*min_per_hour - imin
    dsec = nint(remainder*secprmin)

 END SUBROUTINE compCalday_noleap

 ! ******************************************************************************************
 ! public subroutine extractTime: extract year/month/day/hour/minute/second from units string
 ! ******************************************************************************************
 SUBROUTINE extractTime(refdate,iyyy,im,id,ih,imin,dsec,err,message)
   implicit none
   ! Argument variables
   character(*), intent(in)    :: refdate             ! units string (time since...)
   integer(i4b), intent(out)   :: iyyy,im,id,ih,imin  ! time (year/month/day/hour/minute)
   real(dp),     intent(out)   :: dsec                ! seconds
   integer(i4b), intent(out)   :: err                 ! error code
   character(*), intent(out)   :: message             ! error message
   ! local variables
   integer(i4b)               :: n                   ! length of the string
   integer(i4b)               :: istart,iend         ! position in string

   err=0; message="extractTime/"

   ! get the length of the string
   n      = len_trim(refdate)
   ! move to a position in string past the time units (seconds since , days since , hours since )
   istart = index(refdate,'since')  ! get the index at the beginning of the word "since"
   if (istart>0) then ! if the word "since" exists
     iend   = index(refdate(istart:n)," ")
     istart = istart+iend
   else
     istart=1
   end if

   ! get the year
   call extract(refdate(istart:n),"-",iend,iyyy,err,message); if (err/=0) return
   if(iyyy <    0)then; err=20; message=trim(message)//'year <    0'; return; end if
   if(iyyy > 3000)then; err=20; message=trim(message)//'year > 3000'; return; end if
   ! get the month
   istart=istart+iend
   call extract(refdate(istart:n),"-",iend,im,err,message);   if (err/=0) return
   if(im <  1)then; err=20; message=trim(message)//'month < 1'; return; end if
   if(im > 12)then; err=20; message=trim(message)//'month > 12'; return; end if
   ! get the day
   istart=istart+iend
   call extract(refdate(istart:n)," ",iend,id,err,message);   if (err/=0) return
   if(id <  1)then; err=20; message=trim(message)//'day < 1'; return; end if
   if(id > 31)then; err=20; message=trim(message)//'day > 31'; return; end if
   ! check if we are at the end of the string
   if (istart+(iend-2)==n) then
     ih=0; imin=0; dsec=0._dp; return
   end if

   ! get the hour (":" at end of hour)
   istart = istart+iend
   if(istart > len_trim(refdate))then; err=20; message=trim(message)//'string does not include hours'; return; end if
   call extract(refdate(istart:n),":",iend,ih,err,message);   if (err/=0) return
   if(ih <  0)then; err=20; message=trim(message)//'hour < 0'; return; end if
   if(ih > 24)then; err=20; message=trim(message)//'hour > 24'; return; end if
   ! get the minute (":" at end of minute)
   istart = istart+iend
   if(istart > len_trim(refdate))then; err=20; message=trim(message)//'string does not include minutes'; return; end if
   call extract(refdate(istart:n),":",iend,imin,err,message); if (err/=0) return
   if(imin <  0)then; err=20; message=trim(message)//'minute < 0'; return; end if
   if(imin > 60)then; err=20; message=trim(message)//'minute > 60'; return; end if

   ! get the second
   istart = istart+iend
   ! if second is missing (e.g., yyyy-mm-dd hh:mm)...
   if(istart > len_trim(refdate)) then
     dsec=0._dp; return
   end if
   iend   = index(refdate(istart:n)," ")
   read(refdate(istart:n),*) dsec
   !write(*,'(a,i4,1x,4(i2,1x))') 'refdate: iyyy, im, id, ih, imin = ', iyyy, im, id, ih, imin

 CONTAINS
   ! internal subroutine extract: extract substring
   SUBROUTINE extract(substring,cdelim,iend,itemp,err,message)
   implicit none
   ! Argument variables
   character(*),    intent(in)     :: substring  ! sub-string to process
   character(len=1),intent(in)     :: cdelim     ! string delimiter
   integer(i4b),    intent(out)    :: iend       ! index at the end of desired string
   integer(i4b),    intent(out)    :: itemp      ! output date
   integer(i4b),    intent(out)    :: err        ! error code
   character(*),    intent(out)    :: message    ! error message

   err=0; message="extract/"
   ! identify end-point of string
   iend = index(substring,cdelim)
   ! if sub-string does not exist, assume end is at end of string
   if (iend==0) iend=len_trim(substring)+1
   ! convert string to integer
   read(substring(1:iend-1),*,iostat=err) itemp
   ! read error
   if (err/=0) then
     err=20; message=trim(message)//"unexpected characters [string='"//trim(substring)//"']"; return
   end if
   END SUBROUTINE extract
 END SUBROUTINE extractTime

END MODULE datetime_data
