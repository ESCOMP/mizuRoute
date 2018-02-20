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
                       dTime,          &  ! output: time
                       ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 integer(i4b), intent(in)        :: iTime        ! index of time element
 ! output variables
 real(dp), intent(out)           :: dTime        ! time
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 real(dp),allocatable            :: dummy(:,:)   ! dummy 2d array
 integer(i4b)                    :: nHRU         ! number of hrus
 character(len=strLen)           :: cmessage     ! error message from subroutine

 ! initialize error control
 ierr=0; message='get_runoff/'

 call get_nc(trim(fname),vname_time, dTime, iTime, ierr, cmessage);
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call get_nc_dim_len(fname, dname_hruid, nHRU, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 allocate(dummy(nHRU,1),stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating dummy'; return; endif

 if (.not.(allocated(runoff_data%qsim))) then
   allocate(runoff_data%qsim(nHRU),stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data%qsim'; return; endif
 endif

 call get_nc(trim(fname),vname_qsim, dummy, (/1,iTime/), (/nHRU,1/), ierr, cmessage);
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 runoff_data%qsim = reshape(dummy,(/nHRU/))

 end subroutine

end module
