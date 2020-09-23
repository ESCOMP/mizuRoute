MODULE date_time

! datetime class

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

CONTAINS

  procedure, public :: set_datetime  => sub_set_datetime
  procedure, public :: jul2datetime  => sub_jul2datetime
  procedure, public :: str2datetime  => sub_str2datetime
  procedure, public :: year          => fn_get_year
  procedure, public :: month         => fn_get_month
  procedure, public :: day           => fn_get_day
  procedure, public :: hour          => fn_get_hour
  procedure, public :: minute        => fn_get_min
  procedure, public :: sec           => fn_get_sec
  procedure, public :: is_leap_year  => fn_is_leap_year
  procedure, public :: ndays_month   => fn_ndays_month
  procedure, public :: julianday     => sub_julian_day
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
  generic            :: operator(==)  => fn_is_equal
  generic            :: operator(>)   => fn_is_gt
  generic            :: operator(>=)  => fn_is_ge
  generic            :: operator(<)   => fn_is_lt
  generic            :: operator(<=)  => fn_is_le
  generic            :: assignment(=) => sub_assign

end type datetime

private :: sub_set_datetime, sub_jul2datetime, sub_str2datetime
private :: fn_get_year, fn_get_month, fn_get_day, fn_get_hour, fn_get_min, fn_get_sec
private :: fn_is_leap_year, sub_julian_day, fn_ndays_month
private :: fn_add_months, fn_add_days, fn_add_hours, fn_add_sec
private :: fn_is_equal_month, fn_is_equal_day, fn_is_equal_time

CONTAINS

  SUBROUTINE sub_set_datetime(this, iy, im, id, ih, imin, dsec)

    implicit none
    class(datetime)              :: this
    integer(i4b),     intent(in) :: iy
    integer(i4b),     intent(in) :: im
    integer(i4b),     intent(in) :: id
    integer(i4b),     intent(in) :: ih
    integer(i4b),     intent(in) :: imin
    real(dp),         intent(in) :: dsec

    this%iy       = iy
    this%im       = im
    this%id       = id
    this%ih       = ih
    this%imin     = imin
    this%dsec     = dsec

  END SUBROUTINE sub_set_datetime

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

  END SUBROUTINE sub_jul2datetime

  SUBROUTINE sub_str2datetime(this, str, ierr, message)

    implicit none
    class(datetime)                    :: this
    character(*),         intent(in)   :: str
    integer(i4b),         intent(out)  :: ierr
    character(len=strLen),intent(out)  :: message
    ! local variable
    character(len=strLen)              :: cmessage

    ierr=0; message='sub_str2datetime/'

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
  END SUBROUTINE sub_assign


  integer(i4b) FUNCTION fn_get_year(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_year = this%iy
  END FUNCTION fn_get_year

  integer(i4b) FUNCTION fn_get_month(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_month = this%im
  END FUNCTION fn_get_month

  integer(i4b) FUNCTION fn_get_day(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_day = this%id
  END FUNCTION fn_get_day

  integer(i4b) FUNCTION fn_get_hour(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_hour = this%ih
  END FUNCTION fn_get_hour

  integer(i4b) FUNCTION fn_get_min(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_min = this%imin
  END FUNCTION fn_get_min

  real(dp) FUNCTION fn_get_sec(this)
    implicit none
    class(datetime), intent(in) :: this
    fn_get_sec = this%dsec
  END FUNCTION fn_get_sec


  logical(lgt) FUNCTION fn_is_leap_year(this)
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

  integer(i4b) FUNCTION fn_ndays_month(this, calendar, ierr, message)
    implicit none
    class(datetime),      intent(in)  :: this
    character(*),         intent(in)  :: calendar
    integer(i4b),         intent(out) :: ierr
    character(len=strLen),intent(out) :: message
    ! local variables
    integer(i4b)                      :: nmonths(12)

    ierr=0; message="ndays_month/"

    fn_ndays_month = integerMissing
    select case(trim(calendar))
      case ('standard','gregorian','proleptic_gregorian')
        if ( fn_is_leap_year(this)) then
          nmonths = [31,29,31,30,31,30,31,31,30,31,30,31]
        else
          nmonths = [31,28,31,30,31,30,31,31,30,31,30,31]
        end if
      case('noleap')
        nmonths = [31,28,31,30,31,30,31,31,30,31,30,31]
      case default; ierr=20; message=trim(message)//'calendar name: '//trim(calendar)//' invalid'; return
    end select
    fn_ndays_month = nmonths(this%im)

  END FUNCTION fn_ndays_month


  SUBROUTINE sub_julian_day(this, calendar, julianDay, ierr, message)
    implicit none
    class(datetime),      intent(in)  :: this
    character(*),         intent(in)  :: calendar
    real(dp),             intent(out) :: julianDay
    integer(i4b),         intent(out) :: ierr
    character(len=strLen),intent(out) :: message
    ! local variables
    character(len=strLen)             :: cmessage

    ierr=0; message='sub_julian_day/'
    select case(trim(calendar))
      case ('standard','gregorian','proleptic_gregorian')
        call compJulday(this%iy, this%im, this%id, this%ih, this%imin, this%dsec, julianDay, ierr, cmessage)
      case('noleap')
        call compJulday_noleap(this%iy, this%im, this%id, this%ih, this%imin, this%dsec, julianDay, ierr, cmessage)
      case default; ierr=20; message=trim(message)//'calendar name: '//trim(calendar)//' invalid'; return
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

  type(datetime) FUNCTION fn_add_days(this, days, calendar, ierr, message)
    implicit none
    class(datetime),      intent(inout) :: this
    integer(i4b),         intent(in)    :: days
    character(*),         intent(in)    :: calendar
    integer(i4b),         intent(out)   :: ierr
    character(len=strLen),intent(out)   :: message
    ! local variables
    real(dp)                            :: julday
    character(len=strLen)               :: cmessage

    ierr=0; message='fn_add_days/'

    call sub_julian_day(this, calendar, julday, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

    julday = julday + real(days,dp)
    call sub_jul2datetime(fn_add_days, julday, calendar, ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

  END FUNCTION fn_add_days

  type(datetime) FUNCTION fn_add_hours(this, hrs, calendar, ierr, message)
    implicit none
    class(datetime),      intent(in)  :: this
    integer(i4b),         intent(in)  :: hrs
    character(*),         intent(in)  :: calendar
    integer(i4b),         intent(out) :: ierr
    character(len=strLen),intent(out) :: message
    ! local variables
    real(dp)                          :: julday
    character(len=strLen)             :: cmessage

    ierr=0; message='fn_add_hours/'
    call sub_julian_day(this, calendar, julday, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

    julday = julday + real(hrs,dp)/hr_per_day
    call sub_jul2datetime(fn_add_hours, julday, calendar, ierr, message)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

  END FUNCTION fn_add_hours

  type(datetime) FUNCTION fn_add_sec(this, sec, calendar, ierr, message)
    implicit none
    class(datetime),      intent(in)   :: this
    real(dp),             intent(in)   :: sec
    character(*),         intent(in)   :: calendar
    integer(i4b),         intent(out)  :: ierr
    character(len=strLen),intent(out)  :: message
    ! local variables
    real(dp)                          :: julday
    character(len=strLen)             :: cmessage

    ierr=0; message='fn_add_sec/'
    call sub_julian_day(this, calendar, julday, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

    julday = julday + sec/secprday
    call sub_jul2datetime(fn_add_sec, julday, calendar, ierr, message)
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
    ! if this >= that T, otherwise F
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

END MODULE date_time
