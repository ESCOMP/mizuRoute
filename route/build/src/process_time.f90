MODULE process_time_module

! data types
USE nrtype,    only : i4b,dp              ! variable types, etc.
USE nrtype,    only : strLen              ! length of characters
USE dataTypes, only : time                ! time data type

! subroutines: model time info
USE time_utils_module, only : extractTime       ! get time from units string
USE time_utils_module, only : compJulday,&      ! compute julian day
                              compJulday_noleap ! compute julian day for noleap calendar
USE time_utils_module, only : compcalday,&      ! compute calendar date and time
                              compcalday_noleap ! compute calendar date and time for noleap calendar

implicit none

! privacy -- everything private unless declared explicitly
private
public::process_time
public::conv_cal2julian
public::conv_julian2cal

CONTAINS

 ! *********************************************************************
 ! public subroutine: extract time information from the control information
 ! *********************************************************************
 SUBROUTINE process_time(&
                         ! input
                         timeUnits,        & ! time units string
                         calendar,         & ! calendar
                         ! output
                         julianDate,       & ! julian date
                         ierr, message)
 implicit none
 ! input
 character(*)      , intent(in)               :: timeUnits        ! time units string
 character(*)      , intent(in)               :: calendar         ! calendar string
 ! output
 real(dp)          , intent(out)              :: julianDate       ! julian date
 integer(i4b)      , intent(out)              :: ierr             ! error code
 character(*)      , intent(out)              :: message          ! error message
 ! --------------------------------------------------------------------------------------------------------------
 ! local variables
 type(time)                      :: timeStruct       ! time data structure
 character(len=strLen)           :: cmessage         ! error message of downwind routine

 ierr=0; message='process_time/'

 ! extract time from the units string
 call extractTime(timeUnits,timeStruct%iy,timeStruct%im,timeStruct%id,timeStruct%ih,timeStruct%imin,timeStruct%dsec,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call conv_cal2julian(timeStruct, calendar, julianDate, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE process_time

 ! *********************************************************************
 ! public subroutine: convert julian date to calendar date/time
 ! *********************************************************************
 SUBROUTINE conv_julian2cal(julianDate,       & ! input: julian date
                            calendar,         & ! input: calendar
                            datetime,         & ! output: calendar date/time
                            ierr, message)
 implicit none
 ! input
 real(dp)          , intent(in)               :: julianDate       ! julian date
 character(*)      , intent(in)               :: calendar         ! calendar string
 ! output
 type(time)        , intent(out)              :: datetime         ! time data structure
 integer(i4b)      , intent(out)              :: ierr             ! error code
 character(*)      , intent(out)              :: message          ! error message
 ! --------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=strLen)           :: cmessage         ! error message of downwind routine

 ierr=0; message='conv_julian2cal/'

 select case(trim(calendar))
  case ('noleap','365_day')
   call compcalday_noleap(julianDate, datetime%iy,datetime%im,datetime%id,datetime%ih,datetime%imin,datetime%dsec, ierr, cmessage)
  case ('standard','gregorian','proleptic_gregorian')
   call compcalday(julianDate, datetime%iy,datetime%im,datetime%id,datetime%ih,datetime%imin,datetime%dsec, ierr, cmessage)
  case default; ierr=20; message=trim(message)//trim(calendar)//': calendar invalid; accept either noleap, 365_day, standard, gregorian, or proleptic_gregorian'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE conv_julian2cal

 ! *********************************************************************
 ! public subroutine: convert calendar date/time to julian days
 ! *********************************************************************
 SUBROUTINE conv_cal2julian(datetime,         & ! input: calendar date
                            calendar,         & ! input: calendar
                            julianDate,       & ! output: julian date
                            ierr, message)
 implicit none
 ! input
 type(time)        , intent(in)               :: datetime         ! time data structure
 character(*)      , intent(in)               :: calendar         ! calendar string
 ! output
 real(dp)          , intent(out)              :: julianDate       ! julian date
 integer(i4b)      , intent(out)              :: ierr             ! error code
 character(*)      , intent(out)              :: message          ! error message
 ! --------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=strLen)           :: cmessage         ! error message of downwind routine

 ierr=0; message='conv_cal2julian/'

 select case(trim(calendar))
  case ('noleap','365_day')
   call compjulday_noleap(datetime%iy, datetime%im, datetime%id, datetime%ih, datetime%imin, datetime%dsec, julianDate, ierr, cmessage)
  case ('standard','gregorian','proleptic_gregorian')
   call compjulday(datetime%iy, datetime%im, datetime%id, datetime%ih, datetime%imin, datetime%dsec, julianDate, ierr, cmessage)
  case default; ierr=20; message=trim(message)//trim(calendar)//': calendar invalid; accept either noleap, 365_day, standard, gregorian, or proleptic_gregorian'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE conv_cal2julian


END MODULE process_time_module
