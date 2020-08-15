module read_runoff

USE netcdf
USE nrtype
USE public_var
USE io_netcdf, only:open_nc
USE io_netcdf, only:get_nc
USE io_netcdf, only:get_var_attr
USE io_netcdf, only:check_attr
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
                                fname          , & ! filename
                                ! output
                                runoff_data_in , & ! runoff data structure
                                timeUnits      , & ! time units
                                calendar       , & ! calendar
                                ! error control
                                ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname           ! filename
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
 ierr = nf90_inq_varid(ncid, trim(vname_qsim), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the number of dimensions - must be 2D(hru, time) or 3D(y, x, time)
 ierr= nf90_inquire_variable(ncid, ivarID, ndims = nDims)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get runoff metadata
 select case( nDims )
  case(2); call read_1D_runoff_metadata(fname, runoff_data_in, timeUnits, calendar, ierr, cmessage)
  case(3); call read_2D_runoff_metadata(fname, runoff_data_in, timeUnits, calendar, ierr, cmessage)
  case default; ierr=20; message=trim(message)//'runoff input must be 2-dimension (e.g, [time, hru]) or 3-dimension (e.g., [time, lat, lon]'; return
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
                                   runoff_data_in , & ! runoff data structure
                                   timeUnits    , & ! time units
                                   calendar     , & ! calendar
                                   ! error control
                                   ierr, message)  ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 ! output variables
 type(runoff), intent(out)               :: runoff_data_in     ! runoff for one time step for all HRUs
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
 call get_nc_dim_len(fname, trim(dname_hruid), runoff_data_in%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get number of time steps from the runoff file
 call get_nc_dim_len(fname, trim(dname_time), runoff_data_in%nTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time units
 if (trim(timeUnits) == charMissing) then
   call get_var_attr(fname, trim(vname_time), 'units', timeUnits, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! get the calendar
 if (trim(calendar) == charMissing) then
   call get_var_attr(fname, trim(vname_time), 'calendar', calendar, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! allocate space for hru_id
 allocate(runoff_data_in%hru_id(runoff_data_in%nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data_in%hruId'; return; endif

 ! allocate space for simulated runoff
 allocate(runoff_data_in%sim(runoff_data_in%nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data_in%sim'; return; endif

 ! get HRU ids from the runoff file
 call get_nc(fname, vname_hruid, runoff_data_in%hru_id, 1, runoff_data_in%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine read_1D_runoff_metadata

 ! *****
 ! private subroutine: get 3D runoff (lat, lon, time) metadata...
 ! ******************************************
 subroutine read_2D_runoff_metadata(&
                                    ! input
                                    fname       , & ! filename
                                    ! output
                                    runoff_data_in , & ! runoff data structure
                                    timeUnits   , & ! time units
                                    calendar    , & ! calendar
                                    ! error control
                                    ierr, message)  ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 ! output variables
 type(runoff), intent(out)               :: runoff_data_in     ! runoff for one time step for all HRUs
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
 call get_nc_dim_len(fname, trim(dname_time), runoff_data_in%nTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time units
 if (trim(timeUnits) == charMissing) then
   call get_var_attr(fname, trim(vname_time), 'units', timeUnits, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! get the calendar
 if (trim(calendar) == charMissing) then
   call get_var_attr(fname, trim(vname_time), 'calendar', calendar, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! get size of ylat dimension
 call get_nc_dim_len(fname, trim(dname_ylat), runoff_data_in%nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get size of xlon dimension
 call get_nc_dim_len(fname, trim(dname_xlon), runoff_data_in%nSpace(2), ierr, cmessage)
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
 integer(i4b)                       :: iStart(2)
 integer(i4b)                       :: iCount(2)
 logical(lgt)                       :: existFillVal
 real(dp)                           :: dummy(nSpace,1)    ! data read
 character(len=strLen)              :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_1D_runoff/'

 ! get the simulated runoff data
 iStart = [1,iTime]
 iCount = [nSpace,1]
 call get_nc(fname, var_name, dummy, iStart, iCount, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the _fill_values for runoff variable if exist
 existFillVal = check_attr(fname, var_name, '_FillValue')
 if (existFillval) then
   call get_var_attr(fname, var_name, '_FillValue', input_fillvalue, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! replace _fill_value with -999 for dummy
 where ( abs(dummy - input_fillvalue) < verySmall ) dummy = realMissing

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
 logical(lgt)                :: existFillVal
 integer(i4b)                :: iStart(3)
 integer(i4b)                :: iCount(3)
 real(dp)                    :: dummy(nSpace(2),nSpace(1),1) ! data read
 character(len=strLen)       :: cmessage                     ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_2D_runoff/'

 ! get the simulated runoff data
 iStart = [1,1,iTime]
 iCount = [nSpace(2),nSpace(1),1]
 call get_nc(fname, var_name, dummy, iStart, iCount, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the _fill_values for runoff variable
 existFillVal = check_attr(fname, var_name, '_FillValue')
 if (existFillval) then
   call get_var_attr(fname, var_name, '_FillValue', input_fillvalue, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! replace _fill_value with -999 for dummy
 where ( abs(dummy - input_fillvalue) < verySmall ) dummy = realMissing

 ! reshape
 runoff_data_in%sim2d(1:nSpace(2),1:nSpace(1)) = dummy(1:nSpace(2),1:nSpace(1),1)

 end subroutine read_2D_runoff

end module read_runoff
