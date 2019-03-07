module read_runoff

USE nrtype
USE netcdf
USE public_var
USE read_netcdf, only:get_nc
USE read_netcdf, only:get_var_attr_real
USE read_netcdf, only:get_nc
USE read_netcdf, only:get_nc_dim_len
USE read_netcdf, only:get_var_attr_char
USE dataTypes,   only:runoff                 ! runoff data type

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
                                fname        , & ! filename
                                ! output
                                _runoff_data , & ! runoff data structure
                                timeUnits    , & ! time units
                                calendar     , & ! calendar
                                ! error control
                                ierr, message)  ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname           ! filename
 ! output variables
 type(runoff), intent(out)       :: _runoff_data    ! runoff for one time step for all HRUs
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
 ierr = nf90_open(trim(fname), nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//'['//trim(nf90_strerror(ierr))//'; file='//trim(fname)//']'; return; endif

 ! get the ID of runoff variable
 ierr = nf90_inq_varid(ncid, trim(vname_qsim), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the number of dimensions - must be 2D(hru, time) or 3D(y, x, time)
 ierr= nf90_inquire_variable(ncid, ivarID, ndims = nDims)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get runoff metadata
 select case( nDims )
  case(2); call read_1D_runoff_metadata(fname, _runoff_data, timeUnits, calendar, ierr, cmessage)
  case(3); call read_2D_runoff_metadata(fname, _runoff_data, timeUnits, calendar, ierr, cmessage)
  case default; ierr=20; message=trim(message)//'runoff array nDimensions must be 2 or 3'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine read_runoff_metadata

 ! *****
 ! private subroutine: get 2D runoff (hru, time) metadata...
 ! ******************************************
 subroutine read_1D_runoff_metadata(&
                                   ! input
                                   fname        , & ! filename
                                   ! output
                                   _runoff_data , & ! runoff data structure
                                   timeUnits    , & ! time units
                                   calendar     , & ! calendar
                                   ! error control
                                   ierr, message)  ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 ! output variables
 type(runoff), intent(out)               :: _runoff_data     ! runoff for one time step for all HRUs
 character(*), intent(out)               :: timeUnits       ! time units
 character(*), intent(out)               :: calendar        ! calendar
 ! error control
 integer(i4b), intent(out)               :: ierr            ! error code
 character(*), intent(out)               :: message         ! error message
 ! local variables
 character(len=strLen)                   :: cmessage        ! error message from subroutine
 ! initialize error control
 ierr=0; message='read_1D_runoff_metadata/'

 _runoff_data%nSpace(2) = integerMissing

 ! get the number of HRUs
 call get_nc_dim_len(fname, trim(dname_hruid), _runoff_data%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get number of time steps from the runoff file
 call get_nc_dim_len(fname, trim(dname_time), _runoff_data%nTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time units
 call get_var_attr_char(fname, trim(vname_time), 'units', timeUnits, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the calendar
 call get_var_attr_char(fname, trim(vname_time), 'calendar', calendar, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for hru_id
 allocate(_runoff_data%hru_id(_runoff_data%nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating _runoff_data%hruId'; return; endif

 ! allocate space for simulated runoff
 allocate(_runoff_data%qSim(_runoff_data%nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating _runoff_data%qsim'; return; endif

 ! get HRU ids from the runoff file
 call get_nc(fname, vname_hruid, _runoff_data%hru_id, 1, _runoff_data%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine read_1D_runoff_metadata

 ! *****
 ! private subroutine: get 3D runoff (lat, lon, time) metadata...
 ! ******************************************
 subroutine read_2D_runoff_metadata(&
                                    ! input
                                    fname       , & ! filename
                                    ! output
                                    _runoff_data , & ! runoff data structure
                                    timeUnits   , & ! time units
                                    calendar    , & ! calendar
                                    ! error control
                                    ierr, message)  ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 ! output variables
 type(runoff), intent(out)               :: _runoff_data     ! runoff for one time step for all HRUs
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
 call get_nc_dim_len(fname, trim(dname_time), _runoff_data%nTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time units
 call get_var_attr_char(fname, trim(vname_time), 'units', timeUnits, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the calendar
 call get_var_attr_char(fname, trim(vname_time), 'calendar', calendar, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get size of ylat dimension
 call get_nc_dim_len(fname, trim(dname_ylat), _runoff_data%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get size of xlon dimension
 call get_nc_dim_len(fname, trim(dname_xlon), _runoff_data%nSpace(2), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for simulated runoff. qSim2d = runoff(lon, lat)
 allocate(_runoff_data%qSim2d(_runoff_data%nSpace(2),_runoff_data%nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating qsim'; return; endif

 end subroutine read_2D_runoff_metadata


 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 subroutine read_runoff_data(fname,          &  ! input: filename
                             iTime,          &  ! input: time index
                             _runoff_data,   &  ! inout: runoff data structure
                             ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)      :: fname              ! filename
 integer(i4b), intent(in)      :: iTime              ! index of time element
 ! input/output variables
 type(runoff), intent(inout)   :: _runoff_data       ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)     :: ierr               ! error code
 character(*), intent(out)     :: message            ! error message
 ! local variables
 character(len=strLen)         :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='get_runoff/'

 if (_runoff_data%nSpace(2) == integerMissing) then
  call get_1D_runoff(fname, iTime, _runoff_data%nSpace(1), _runoff_data, ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 else
  call get_2D_runoff(fname, iTime, _runoff_data%nSpace, _runoff_data, ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif

 end subroutine read_runoff_data

 ! *********************************************************************
 ! private subroutine: read 2D runoff data
 ! *********************************************************************
 subroutine read_1D_runoff(fname,          &  ! input: filename
                           iTime,          &  ! input: time index
                           nSpace,         &  ! input: size of HRUs
                           _runoff_data,   &  ! inout: runoff data structure
                           ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)      :: fname              ! filename
 integer(i4b), intent(in)      :: iTime              ! index of time element
 integer(i4b), intent(in)      :: nSpace             ! size of spatial dimensions
 ! input/output variables
 type(runoff), intent(inout)   :: _runoff_data       ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)     :: ierr               ! error code
 character(*), intent(out)     :: message            ! error message
 ! local variables
 real(dp)                      :: fill_value         ! fill_value
 real(dp)                      :: dummy(nSpace,1)    ! data read
 character(len=strLen)         :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_1D_runoff/'

 ! get the time data
 call get_nc(trim(fname), vname_time, _runoff_data%time, iTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the simulated runoff data
 call get_nc(trim(fname),vname_qsim, dummy, (/1,iTime/), (/nSpace,1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the _fill_values for runoff variable
 call get_var_attr_real(trim(fname), vname_qsim, '_FillValue', fill_value, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! replace _fill_value with -999 for dummy
 where ( abs(dummy - fill_value) < verySmall ) dummy =realMissing

 ! reshape
 _runoff_data%qsim(1:nSpace) = dummy(1:nSpace,1)

 end subroutine read_1D_runoff

 ! *********************************************************************
 ! private subroutine: read 2D runoff data
 ! *********************************************************************
 subroutine read_2D_runoff(fname,          &  ! input: filename
                           iTime,          &  ! input: time index
                           nSpace,         &  ! input: size of HRUs
                           _runoff_data,   &  ! output: runoff data structure
                           ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)    :: fname            ! filename
 integer(i4b), intent(in)    :: iTime            ! index of time element
 integer(i4b), intent(in)    :: nSpace(1:2)      ! size of spatial dimensions
 ! input/output variables
 type(runoff), intent(inout) :: _runoff_data     ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)   :: ierr             ! error code
 character(*), intent(out)   :: message          ! error message
 ! local variables
 real(dp)                   :: fill_value                   ! fill_value
 real(dp)                   :: dummy(nSpace(2),nSpace(1),1) ! data read
 character(len=strLen)      :: cmessage                     ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_2D_runoff/'

 ! get the time data
 call get_nc(trim(fname), vname_time, _runoff_data%time, iTime, ierr, cmessage)
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
 _runoff_data%qsim2d(1:nSpace(2),1:nSpace(1)) = dummy(1:nSpace(2),1:nSpace(1),1)

 end subroutine read_2D_runoff

end module read_runoff
