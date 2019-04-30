module read_runoff
USE nrtype
USE netcdf
USE public_var
USE read_netcdf, only:get_nc, &
                      get_var_attr_real
USE globalData,  only:runoff_data
USE dataTypes,   only:runoff                 ! runoff data type

implicit none

private
public::get_runoff

contains

 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************

 subroutine get_runoff(fname,          &  ! input: filename
                       iTime,          &  ! input: time index
                       nSpace,         &  ! input: size of HRUs
                       runoff_data,    &  ! inout: runoff data structure
                       ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)      :: fname              ! filename
 integer(i4b), intent(in)      :: iTime              ! index of time element
 integer(i4b), intent(in)      :: nSpace(1:2)        ! size of spatial dimensions
 ! input/output variables
 type(runoff), intent(inout)   :: runoff_data        ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)     :: ierr               ! error code
 character(*), intent(out)     :: message            ! error message
 ! local variables
 character(len=strLen)         :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='get_runoff/'

 if (nSpace(2) == integerMissing) then
  call get_1D_runoff(fname, iTime, nSpace(1), runoff_data, ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 else
  call get_2D_runoff(fname, iTime, nSpace, runoff_data, ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif

 end subroutine get_runoff

 ! *********************************************************************
 ! private subroutine: read 2D runoff data
 ! *********************************************************************
 subroutine get_1D_runoff(fname,          &  ! input: filename
                          iTime,          &  ! input: time index
                          nSpace,         &  ! input: size of HRUs
                          runoff_data,    &  ! inout: runoff data structure
                          ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)      :: fname              ! filename
 integer(i4b), intent(in)      :: iTime              ! index of time element
 integer(i4b), intent(in)      :: nSpace             ! size of spatial dimensions
 ! input/output variables
 type(runoff), intent(inout)   :: runoff_data        ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)     :: ierr               ! error code
 character(*), intent(out)     :: message            ! error message
 ! local variables
 real(dp)                      :: fill_value         ! fill_value
 real(dp)                      :: dummy(nSpace,1)    ! data read
 character(len=strLen)         :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='get_1D_runoff/'

 ! get the time data
 call get_nc(trim(fname), vname_time, runoff_data%time, iTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the simulated runoff data
 call get_nc(trim(fname),vname_qsim, dummy, (/1,iTime/), (/nSpace,1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the _fill_values for runoff variable
 call get_var_attr_real(trim(fname), vname_qsim, '_FillValue', fill_value, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! replace _fill_value with -999 for dummy
 where ( abs(dummy - fill_value) < verySmall ) dummy = realMissing

 ! reshape
 runoff_data%qsim(1:nSpace) = dummy(1:nSpace,1)

 end subroutine get_1D_runoff

 ! *********************************************************************
 ! private subroutine: read 2D runoff data
 ! *********************************************************************
 subroutine get_2D_runoff(fname,          &  ! input: filename
                          iTime,          &  ! input: time index
                          nSpace,         &  ! input: size of grid dimensions
                          runoff_data,    &  ! output: runoff data structure
                          ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)    :: fname            ! filename
 integer(i4b), intent(in)    :: iTime            ! index of time element
 integer(i4b), intent(in)    :: nSpace(1:2)      ! size of spatial dimensions
 ! input/output variables
 type(runoff), intent(inout) :: runoff_data      ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)   :: ierr             ! error code
 character(*), intent(out)   :: message          ! error message
 ! local variables
 real(dp)                   :: fill_value                   ! fill_value
 real(dp)                   :: dummy(nSpace(2),nSpace(1),1) ! data read
 character(len=strLen)      :: cmessage                     ! error message from subroutine

 ! initialize error control
 ierr=0; message='get_2D_runoff/'

 ! get the time data
 call get_nc(trim(fname), vname_time, runoff_data%time, iTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the simulated runoff data
 call get_nc(trim(fname), vname_qsim, dummy, (/1,1,iTime/), (/nSpace(2), nSpace(1), 1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the _fill_values for runoff variable
 call get_var_attr_real(trim(fname), vname_qsim, '_FillValue', fill_value, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! replace _fill_value with -999 for dummy
 where ( abs(dummy - fill_value) < verySmall ) dummy = realMissing

 ! reshape
 runoff_data%qsim2d(1:nSpace(2),1:nSpace(1)) = dummy(1:nSpace(2),1:nSpace(1),1)

 end subroutine get_2D_runoff

end module
