module read_runoff

USE netcdf
USE nrtype
USE public_var
USE io_netcdf, only:open_nc
USE io_netcdf, only:get_nc
USE io_netcdf, only:get_var_attr
USE io_netcdf, only:get_nc_dim_len
USE dataTypes, only:runoff                 ! runoff data type


implicit none

private
public::read_runoff_metadata
public::read_runoff_data

contains

 ! *****
 ! public subroutine: get runoff  metadata...
 ! ******************************************
 subroutine read_runoff_metadata(&
                                ! input
                                fname          , &  ! filename including directory
                                var_name       , &  ! varibale name
                                var_time_name  , &  ! name of varibale time
                                dim_time_name  , &  ! name of dimension time
                                var_hru_name   , &  ! name of varibale HRUs
                                dim_hru_name   , &  ! name of dimension HRUs
                                dim_ylat_name  , &  ! name of dimension lat in case of a 2D input varibale
                                dim_xlon_name  , &  ! name of dimension lon in case of a 2D input varibale
                                ! output
                                runoff_data_in , &  ! runoff data structure
                                timeUnits      , &  ! time units
                                calendar       , &  ! calendar
                                ! error control
                                ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname           ! filename
 character(*), intent(in)        :: var_name        ! name of the varibale for simulated runoff or abstraction/injection
 character(*), intent(in)        :: var_time_name   ! name of the varibale time
 character(*), intent(in)        :: dim_time_name   ! name of dimension for time
 character(*), intent(in)        :: var_hru_name    ! name of the varibale hru
 character(*), intent(in)        :: dim_hru_name    ! name of dimension for hydrological HRUs
 character(*), intent(in)        :: dim_ylat_name   ! name of dimension along lat
 character(*), intent(in)        :: dim_xlon_name   ! name of dimension along lon
 ! output variables
 type(runoff), intent(out)       :: runoff_data_in  ! runoff for one time step for all HRUs
 character(*), intent(out)       :: timeUnits       ! time units
 character(*), intent(out)       :: calendar        ! calendar
 ! error control
 integer(i4b), intent(out)       :: ierr            ! error code
 character(*), intent(out)       :: message         ! error message
 ! local variables
 integer(i4b)                    :: ncid            ! netcdf id
 integer(i4b)                    :: ivarID          ! variable id
 integer(i4b)                    :: nDims           ! number of dimension in runoff file
 character(len=strLen)           :: cmessage        ! error message from subroutine
 ! initialize error control
 ierr=0; message='read_runoff_metadata/'

 ! open NetCDF file
 call open_nc(trim(fname), 'r', ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the ID of runoff variable
 ierr = nf90_inq_varid(ncid, trim(var_name), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the number of dimensions - must be 2D(hru, time) or 3D(y, x, time)
 ierr= nf90_inquire_variable(ncid, ivarID, ndims = nDims)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get runoff metadata
 select case( nDims )
  case(2); call read_1D_runoff_metadata(fname,            &
                                        var_time_name,    &
                                        dim_time_name,    &
                                        var_hru_name,     &
                                        dim_hru_name,     &
                                        runoff_data_in,   &
                                        timeUnits,        &
                                        calendar,         &
                                        ierr, cmessage)
  case(3); call read_2D_runoff_metadata(fname,            &
                                        var_time_name,    &
                                        dim_time_name,    &
                                        dim_ylat_name,    &
                                        dim_xlon_name,    &
                                        runoff_data_in,   &
                                        timeUnits,        &
                                        calendar,         &
                                        ierr, cmessage)
  case default; ierr=20; message=trim(message)//'runoff input must be 2-dimension (e.g, [time, hru]) or 3-dimension (e.g., [time, lat, lon]'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine read_runoff_metadata

 ! *****
 ! private subroutine: get 2D runoff (hru, time) metadata...
 ! ******************************************
 subroutine read_1D_runoff_metadata(&
                                   ! input
                                   fname          , & ! filename
                                   var_time_name  , & ! name of varibale time
                                   dim_time_name  , & ! name of dimension time
                                   var_hru_name   , & ! name of varibale HRUs
                                   dim_hru_name   , & ! name of dimension HRUs
                                   ! output
                                   runoff_data_in , & ! runoff data structure
                                   timeUnits      , & ! time units
                                   calendar       , & ! calendar
                                   ! error control
                                   ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 character(*), intent(in)                :: var_time_name   ! name of the varibale time
 character(*), intent(in)                :: dim_time_name   ! name of dimension for time
 character(*), intent(in)                :: var_hru_name    ! name of the varibale hru
 character(*), intent(in)                :: dim_hru_name    ! name of dimension for hydrological HRUs
 ! output variables
 type(runoff), intent(out)               :: runoff_data_in  ! runoff for one time step for all HRUs
 character(*), intent(out)               :: timeUnits       ! time units
 character(*), intent(out)               :: calendar        ! calendar
 ! error control
 integer(i4b), intent(out)               :: ierr            ! error code
 character(*), intent(out)               :: message         ! error message
 ! local variables
 character(len=strLen)                   :: cmessage        ! error message from subroutine
 ! initialize error control
 ierr=0; message='read_1D_runoff_metadata/'

 runoff_data_in%nSpace(2) = integerMissing

 ! get the number of HRUs
 call get_nc_dim_len(fname, trim(dim_hru_name), runoff_data_in%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get number of time steps from the runoff file
 call get_nc_dim_len(fname, trim(dim_time_name), runoff_data_in%nTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time units
 if (trim(timeUnits) == charMissing) then
   call get_var_attr(fname, trim(var_time_name), 'units', timeUnits, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! get the calendar
 if (trim(calendar) == charMissing) then
   call get_var_attr(fname, trim(var_time_name), 'calendar', calendar, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! allocate space for hru_id
 allocate(runoff_data_in%hru_id(runoff_data_in%nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data_in%hruId'; return; endif

 ! allocate space for simulated runoff
 allocate(runoff_data_in%sim(runoff_data_in%nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data_in%sim'; return; endif

 ! get HRU ids from the runoff file
 call get_nc(fname, var_hru_name, runoff_data_in%hru_id, 1, runoff_data_in%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine read_1D_runoff_metadata

 ! *****
 ! private subroutine: get 3D runoff (lat, lon, time) metadata...
 ! ******************************************
 subroutine read_2D_runoff_metadata(&
                                    ! input
                                    fname          , & ! filename
                                    var_time_name  , & ! name of varibale time
                                    dim_time_name  , & ! name of dimension time
                                    dim_ylat_name  , & ! name of varibale HRUs
                                    dim_xlon_name  , & ! name of dimension HRUs
                                    ! output
                                    runoff_data_in , & ! runoff data structure
                                    timeUnits      , & ! time units
                                    calendar       , & ! calendar
                                    ! error control
                                    ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 character(*), intent(in)                :: var_time_name   ! name of the varibale time
 character(*), intent(in)                :: dim_time_name   ! name of dimension for time
 character(*), intent(in)                :: dim_ylat_name   ! name of dimension along lat
 character(*), intent(in)                :: dim_xlon_name   ! name of dimension along lon
 ! output variables
 type(runoff), intent(out)               :: runoff_data_in  ! runoff for one time step for all HRUs
 character(*), intent(out)               :: timeUnits       ! time units
 character(*), intent(out)               :: calendar        ! calendar
 ! error control
 integer(i4b), intent(out)               :: ierr            ! error code
 character(*), intent(out)               :: message         ! error message
 ! local variables
 character(len=strLen)                   :: cmessage        ! error message from subroutine
 ! initialize error control
 ierr=0; message='read_2D_runoff_metadata/'

 ! get number of time steps from the runoff file
 call get_nc_dim_len(fname, trim(dim_time_name), runoff_data_in%nTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time units
 if (trim(timeUnits) == charMissing) then
   call get_var_attr(fname, trim(var_time_name), 'units', timeUnits, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! get the calendar
 if (trim(calendar) == charMissing) then
   call get_var_attr(fname, trim(var_time_name), 'calendar', calendar, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! get size of ylat dimension
 call get_nc_dim_len(fname, trim(dim_ylat_name), runoff_data_in%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get size of xlon dimension
 call get_nc_dim_len(fname, trim(dim_xlon_name), runoff_data_in%nSpace(2), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for simulated runoff. sim2d = runoff(lon, lat)
 allocate(runoff_data_in%sim2d(runoff_data_in%nSpace(2),runoff_data_in%nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating sim'; return; endif

 end subroutine read_2D_runoff_metadata


 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 subroutine read_runoff_data(fname,          &  ! input: filename
                             var_name,       &  ! input: varibale name
                             iTime,          &  ! input: time index
                             runoff_data_in, &  ! inout: runoff data structure
                             ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)           :: fname              ! filename
 character(*), intent(in)           :: var_name           ! variable name
 integer(i4b), intent(in)           :: iTime              ! index of time element
 ! input/output variables
 type(runoff), intent(inout)        :: runoff_data_in     ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)          :: ierr               ! error code
 character(*), intent(out)          :: message            ! error message
 ! local variables
 character(len=strLen)              :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_runoff_data/'

 if (runoff_data_in%nSpace(2) == integerMissing) then
  call read_1D_runoff(fname, var_name, iTime, runoff_data_in%nSpace(1), runoff_data_in, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 else
  call read_2D_runoff(fname, var_name, iTime, runoff_data_in%nSpace, runoff_data_in, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif

 end subroutine read_runoff_data

 ! *********************************************************************
 ! private subroutine: read 2D runoff data
 ! *********************************************************************
 subroutine read_1D_runoff(fname,          &  ! input: filename
                           var_name,       &  ! input: variable name
                           iTime,          &  ! input: time index
                           nSpace,         &  ! input: size of HRUs
                           runoff_data_in, &  ! inout: runoff data structure
                           ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)           :: fname              ! filename
 character(*), intent(in)           :: var_name           ! variable name
 integer(i4b), intent(in)           :: iTime              ! index of time element
 integer(i4b), intent(in)           :: nSpace             ! size of spatial dimensions
 ! input/output variables
 type(runoff), intent(inout)        :: runoff_data_in     ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)          :: ierr               ! error code
 character(*), intent(out)          :: message            ! error message
 ! local variables
 real(dp)                           :: fill_value         ! fill_value
 real(dp)                           :: dummy(nSpace,1)    ! data read
 character(len=strLen)              :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_1D_runoff/'

 ! get the simulated runoff data
 call get_nc(trim(fname), trim(var_name), dummy, (/1,iTime/), (/nSpace,1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the _fill_values for runoff variable
 call get_var_attr(trim(fname), trim(var_name), '_FillValue', fill_value, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! replace _fill_value with -999 for dummy
 where ( abs(dummy - fill_value) < verySmall ) dummy = realMissing

 ! reshape
 runoff_data_in%sim(1:nSpace) = dummy(1:nSpace,1)

 end subroutine read_1D_runoff

 ! *********************************************************************
 ! private subroutine: read 2D runoff data
 ! *********************************************************************
 subroutine read_2D_runoff(fname,          &  ! input: filename
                           var_name,       &  ! input: variable name
                           iTime,          &  ! input: time index
                           nSpace,         &  ! input: size of HRUs
                           runoff_data_in, &  ! output: runoff data structure
                           ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)    :: fname            ! filename
 character(*), intent(in)    :: var_name         ! variable name
 integer(i4b), intent(in)    :: iTime            ! index of time element
 integer(i4b), intent(in)    :: nSpace(1:2)      ! size of spatial dimensions
 ! input/output variables
 type(runoff), intent(inout) :: runoff_data_in   ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)   :: ierr             ! error code
 character(*), intent(out)   :: message          ! error message
 ! local variables
 real(dp)                    :: fill_value                   ! fill_value
 real(dp)                    :: dummy(nSpace(2),nSpace(1),1) ! data read
 character(len=strLen)       :: cmessage                     ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_2D_runoff/'

 ! get the simulated runoff data
 call get_nc(trim(fname), trim(var_name), dummy, (/1,1,iTime/), (/nSpace(2), nSpace(1), 1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the _fill_values for runoff variable
 call get_var_attr(trim(fname), trim(var_name), '_FillValue', fill_value, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! replace _fill_value with -999 for dummy
 where ( abs(dummy - fill_value) < verySmall ) dummy = realMissing

 ! reshape
 runoff_data_in%sim2d(1:nSpace(2),1:nSpace(1)) = dummy(1:nSpace(2),1:nSpace(1),1)

 end subroutine read_2D_runoff

end module read_runoff
