MODULE datetime_data

! datetime class
!
! usage:
! initialize datetime
!   datetime1 = datetime(1969, 7,16, 7,32,  0.0)
!

USE nrtype
USE public_var,        ONLY: realMissing, integerMissing
USE public_var,        ONLY: secprday, hr_per_day
USE time_utils_module, ONLY: extractTime                 !
USE time_utils_module, ONLY: compJulday,&                ! compute julian day
                             compJulday_noleap           ! compute julian day for noleap calendar
USE time_utils_module, ONLY: compCalday,&                ! compute calendar date and time
                             compCalday_noleap           ! compute calendar date and time for noleap calendar

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
  generic            :: operator(==)  => fn_is_equal
  generic            :: operator(>)   => fn_is_gt
  generic            :: operator(>=)  => fn_is_ge
  generic            :: operator(<)   => fn_is_lt
  generic            :: operator(<=)  => fn_is_le
  generic            :: assignment(=) => sub_assign
  generic            :: operator(-)   => fn_delta

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
    fn_dayofyear = int((this - jan1datetime)/86400._dp) + 1

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
    class(datetime),     intent(inout)  :: this
    integer(i4b),         intent(in)    :: nmonths
    fn_add_months%iy   = this%iy + (nmonths/12_i4b)
    fn_add_months%im   = this%im + mod(nmonths,12)
    fn_add_months%id   = this%id
    fn_add_months%ih   = this%ih
    fn_add_months%imin = this%imin
    fn_add_months%dsec = this%dsec
  END FUNCTION fn_add_months

  type(datetime) FUNCTION fn_add_days(this, days, ierr, message)
    implicit none
    class(datetime),      intent(inout) :: this
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

    sec_in_day1 = real(this%ih*3600+this%imin*60, kind=dp) +this%dsec
    sec_in_day2 = real(that%ih*3600+that%imin*60, kind=dp) +that%dsec

    date1 = datetime(this%iy, this%im, this%id, 0, 0, 0.0, calendar=this%calendar)
    date2 = datetime(that%iy, that%im, that%id, 0, 0, 0.0, calendar=that%calendar)

    call date1%julianday(julday1, ierr, cmessage)
    call date2%julianday(julday2, ierr, cmessage)

    diffDay=julday1-julday2

    ans = diffDay*86400._dp+sec_in_day1-sec_in_day2
  END FUNCTION fn_delta


END MODULE datetime_data
