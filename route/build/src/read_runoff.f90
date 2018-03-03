module read_runoff
use nrtype
USE netcdf
use public_var
use read_netcdf, only:get_nc
use read_netcdf, only:get_nc_dim_len
use globalData,  only:runoff_data

implicit none

private
public::get_runoff

contains

 ! *********************************************************************
 ! new subroutine: read runoff data
 ! *********************************************************************
 subroutine get_runoff(fname,          &  ! input: filename
                       iTime,          &  ! input: time index
                       nHRU,           &  ! input:number of HRUs
                       dTime,          &  ! output: time
                       qSim,           &  ! output: runoff data
                       ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname         ! filename
 integer(i4b), intent(in)        :: iTime         ! index of time element
 integer(i4b), intent(in)        :: nHRU          ! number of HRUs
 ! output variables
 real(dp)    , intent(out)       :: dTime         ! time
 real(dp)    , intent(out)       :: qSim(:)       ! runoff for one time step for all HRUs
 integer(i4b), intent(out)       :: ierr          ! error code
 character(*), intent(out)       :: message       ! error message
 ! local variables
 real(dp)                        :: dummy(nHRU,1) ! data read
 character(len=strLen)           :: cmessage      ! error message from subroutine

 ! initialize error control
 ierr=0; message='get_runoff/'

 ! get the time data
 call get_nc(trim(fname), vname_time, dtime, iTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the simulated runoff data
 call get_nc(trim(fname),vname_qsim, dummy, (/1,iTime/), (/nHRU,1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! reshape
 qSim(1:nHRU)=dummy(1:nHRU,1)

 end subroutine

end module
