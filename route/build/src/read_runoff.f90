module read_runoff
USE nrtype
USE netcdf
USE public_var
USE read_netcdf, only:get_nc
USE read_netcdf, only:get_var_attr_real
USE read_netcdf, only:get_nc
USE read_netcdf, only:get_nc_dim_len
USE read_netcdf, only:get_var_attr_char

USE globalData,  only:runoff_data
USE dataTypes,   only:runoff                 ! runoff data type

implicit none

private
public::get_runoff_metadata
public::get_runoff
public::get_qDims
public::get_qMeta
public::getVarQsim

contains

 ! *****
 ! public subroutine: get runoff  metadata...
 ! ******************************************
 subroutine get_runoff_metadata(&
                                ! input
                                fname       , & ! filename
                                ! output
                                runoff_data , & ! runoff data structure
                                nSpatial    , & ! number of spatial elements
                                nTime       , & ! number of time steps
                                timeUnits   , & ! time units
                                calendar    , & ! calendar
                                ! error control
                                ierr, message)  ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname           ! filename
 ! output variables
 type(runoff), intent(out)       :: runoff_data     ! runoff for one time step for all HRUs
 integer(i4b), intent(out)       :: nSpatial(1:2)   ! number of spatial elements
 integer(i4b), intent(out)       :: nTime           ! number of time steps
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
 ierr=0; message='get_runoff_metadata/'

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
  case(2); call get_1D_runoff_metadata(fname, runoff_data, nSpatial, nTime, timeUnits, calendar, ierr, cmessage)
  case(3); call get_2D_runoff_metadata(fname, runoff_data, nSpatial, nTime, timeUnits, calendar, ierr, cmessage)
  case default; ierr=20; message=trim(message)//'runoff array nDimensions must be 2 or 3'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine get_runoff_metadata

 ! *****
 ! private subroutine: get 2D runoff (hru, time) metadata...
 ! ******************************************
 subroutine get_1D_runoff_metadata(&
                                   ! input
                                   fname       , & ! filename
                                   ! output
                                   runoff_data , & ! runoff data structure
                                   nSpatial    , & ! number of spatial elements
                                   nTime       , & ! number of time steps
                                   timeUnits   , & ! time units
                                   calendar    , & ! calendar
                                   ! error control
                                   ierr, message)  ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 ! output variables
 type(runoff), intent(out)               :: runoff_data     ! runoff for one time step for all HRUs
 integer(i4b), intent(out)               :: nSpatial(1:2)   ! number of spatial elements
 integer(i4b), intent(out)               :: nTime           ! number of time steps
 character(*), intent(out)               :: timeUnits       ! time units
 character(*), intent(out)               :: calendar        ! calendar
 ! error control
 integer(i4b), intent(out)               :: ierr            ! error code
 character(*), intent(out)               :: message         ! error message
 ! local variables
 character(len=strLen)                   :: cmessage        ! error message from subroutine
 ! initialize error control
 ierr=0; message='get_1D_runoff_metadata/'

 nSpatial(2) = integerMissing

 ! get the number of HRUs
 call get_nc_dim_len(fname, trim(dname_hruid), nSpatial(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get number of time steps from the runoff file
 call get_nc_dim_len(fname, trim(dname_time), nTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time units
 call get_var_attr_char(fname, trim(vname_time), 'units', timeUnits, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the calendar
 call get_var_attr_char(fname, trim(vname_time), 'calendar', calendar, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for hru_id
 allocate(runoff_data%hru_id(nSpatial(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data%hruId'; return; endif

 ! allocate space for simulated runoff
 allocate(runoff_data%qSim(nSpatial(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data%qsim'; return; endif

 ! get HRU ids from the runoff file
 call get_nc(fname, vname_hruid, runoff_data%hru_id, 1, nSpatial(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine get_1D_runoff_metadata

 ! *****
 ! private subroutine: get 3D runoff (lat, lon, time) metadata...
 ! ******************************************
 subroutine get_2D_runoff_metadata(&
                                ! input
                                fname       , & ! filename
                                ! output
                                runoff_data , & ! runoff data structure
                                nSpatial    , & ! number of ylat dimensions
                                nTime       , & ! number of time steps
                                timeUnits   , & ! time units
                                calendar    , & ! calendar
                                ! error control
                                ierr, message)  ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 ! output variables
 type(runoff), intent(out)               :: runoff_data     ! runoff for one time step for all HRUs
 integer(i4b), intent(out)               :: nSpatial(1:2)   ! number of spatial elements (lat, lon) (y,x),(i,j)
 integer(i4b), intent(out)               :: nTime           ! number of time steps
 character(*), intent(out)               :: timeUnits       ! time units
 character(*), intent(out)               :: calendar        ! calendar
 ! error control
 integer(i4b), intent(out)               :: ierr            ! error code
 character(*), intent(out)               :: message         ! error message
 ! local variables
 character(len=strLen)                   :: cmessage        ! error message from subroutine
 ! initialize error control
 ierr=0; message='get_2D_runoff_metadata/'

 ! get number of time steps from the runoff file
 call get_nc_dim_len(fname, trim(dname_time), nTime, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time units
 call get_var_attr_char(fname, trim(vname_time), 'units', timeUnits, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the calendar
 call get_var_attr_char(fname, trim(vname_time), 'calendar', calendar, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get size of ylat dimension
 call get_nc_dim_len(fname, trim(dname_ylat), nSpatial(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get size of xlon dimension
 call get_nc_dim_len(fname, trim(dname_xlon), nSpatial(2), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for simulated runoff. qSim2d = runoff(lon, lat)
 allocate(runoff_data%qSim2d(nSpatial(2),nSpatial(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating qsim'; return; endif

 end subroutine get_2D_runoff_metadata


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
 where ( abs(dummy - fill_value) < verySmall ) dummy =realMissing

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

 ! *********************************************************************
 ! new subroutine: read dimensions from runoff file
 ! *********************************************************************
 subroutine get_qDims(fname,           &  ! input: filename
                      vname_hruid,     &  ! input: name of coordinate dimension HRUid
                      vname_time,      &  ! input: name of coordinate dimension time
                      units_time,      &  ! output: time units
                      nTime,           &  ! output: number of time elements
                      nHRU,            &  ! output: number of HRUs in the runoff data file
                      ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname_hruid  ! coordinate variable for HRUs
 character(*), intent(in)        :: vname_time   ! coordinate variable for time
 ! output variables
 character(*), intent(out)       :: units_time   ! time units
 integer(i4b), intent(out)       :: nTime        ! number of time elements
 integer(i4b), intent(out)       :: nHRU         ! number of HRUs
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 integer(i4b)                    :: idimID_time  ! dimension ID for time
 integer(i4b)                    :: idimID_hru   ! dimension ID for hru
 ! initialize error control
 ierr=0; message='get_qDims/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the ID of the time dimension
 ierr = nf90_inq_dimid(ncid, vname_time, idimID_time)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(vname_time); return; endif

 ! get the length of the time dimension
 ierr = nf90_inquire_dimension(ncid, idimID_time, len=nTime)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the ID of the time variable
 ierr = nf90_inq_varid(ncid, trim(vname_time), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the time units
 ierr = nf90_get_att(ncid, ivarID, 'units', units_time)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the ID of the HRU dimension
 ierr = nf90_inq_dimid(ncid, vname_hruID, idimID_hru)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(vname_hruID); return; endif

 ! get the length of the HRU dimension
 ierr = nf90_inquire_dimension(ncid, idimID_hru, len=nHRU)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_qDims

 ! *********************************************************************
 ! new subroutine: read metadata from runoff file
 ! *********************************************************************
 subroutine get_qMeta(fname,           &  ! input: filename
                      vname_hruid,     &  ! input: name of coordinate dimension HRUid
                      qsimHRUid,       &  ! output: HRUid in the simulations
                      ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname_hruid  ! coordinate variable for HRUs
 ! output variables
 integer(i4b), intent(out)       :: qsimHRUid(:) ! HRUid in the simulations
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 ! initialize error control
 ierr=0; message='get_qMeta/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the variable ID for the HRUid
 !print*, 'trim(vname_hruid) = ', trim(vname_hruid)
 ierr = nf90_inq_varid(ncid, trim(vname_hruid), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the data for the HRUid
 !print*, 'size(qsimHRUid) = ', size(qsimHRUid)
 ierr = nf90_get_var(ncid, ivarID, qsimHRUid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_qMeta


 ! *********************************************************************
 ! new subroutine: read runoff data
 ! *********************************************************************
 subroutine getVarQsim(fname,          &  ! input: filename
                       vname_time,     &  ! input: name of coordinate dimension time
                       vname_qsim,     &  ! input: name of runoff variable
                       iTime,          &  ! input: time index
                       dTime,          &  ! output: time
                       qsim_hru,       &  ! output: simulated runoff
                       ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname_time   ! name of coordinate dimension time
 character(*), intent(in)        :: vname_qsim   ! name of runoff variable
 integer(i4b), intent(in)        :: iTime        ! index of time element
 ! output variables
 real(dp), intent(out)           :: dTime        ! time
 real(dp), intent(out)           :: qsim_HRU(:)  ! simulated runoff for each HRU
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 real(dp),dimension(1)           :: tempTime     ! temporarary time vector (length=1)
 ! initialize error control
 ierr=0; message='getVarQsim/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the variable ID for time
 ierr = nf90_inq_varid(ncid, trim(vname_time), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the time
 ierr = nf90_get_var(ncid, ivarID, tempTime, start=(/iTime/), count=(/1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 dTime = tempTime(1)

 ! get the variable ID for the simulated runoff
 ierr = nf90_inq_varid(ncid, trim(vname_qsim), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the simulated runoff data
 ierr = nf90_get_var(ncid, ivarID, qsim_HRU, start=(/1,iTime/), count=(/size(qsim_HRU),1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine getVarQsim


end module
