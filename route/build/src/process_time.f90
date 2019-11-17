module process_time_module

! data types
USE nrtype,    ONLY: i4b,dp              ! variable types, etc.
USE nrtype,    ONLY: strLen              ! length of characters
USE dataTypes, ONLY: time                ! time data type

! subroutines: model time info
USE time_utils_module, ONLY: extractTime       ! get time from units string
USE time_utils_module, ONLY: compJulday,&      ! compute julian day
                              compJulday_noleap ! compute julian day for noleap calendar
implicit none

! privacy -- everything private unless declared explicitly
private
public::process_time

contains

 ! *********************************************************************
 ! public subroutine: extract time information from the control information
 ! *********************************************************************
 subroutine process_time(&
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
 ! initialize error control
 ierr=0; message='process_time/'

 ! extract time from the units string
 call extractTime(timeUnits,timeStruct%iy,timeStruct%im,timeStruct%id,timeStruct%ih,timeStruct%imin,timeStruct%dsec,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! calculate the julian day for the start time
 select case(trim(calendar))
  case ('noleap')
   call compjulday_noleap(timeStruct%iy,timeStruct%im,timeStruct%id,timeStruct%ih,timeStruct%imin,timeStruct%dsec,julianDate,ierr,cmessage)
  case ('standard','gregorian','proleptic_gregorian')
   call compjulday(timeStruct%iy,timeStruct%im,timeStruct%id,timeStruct%ih,timeStruct%imin,timeStruct%dsec,julianDate,ierr,cmessage)
  case default; ierr=20; message=trim(message)//trim(calendar)//': calendar name invalid'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine process_time

end module process_time_module
