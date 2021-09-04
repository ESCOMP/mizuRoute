MODULE read_runoff

USE netcdf
USE nrtype
USE public_var
USE ncio_utils, ONLY: open_nc
USE ncio_utils, ONLY: get_nc
USE ncio_utils, ONLY: get_var_attr
USE ncio_utils, ONLY: check_attr
USE ncio_utils, ONLY: get_nc_dim_len


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
                                 fname            , & ! filename including directory
                                 var_name         , & ! varibale name
                                 var_time_name    , & ! name of varibale time
                                 dim_time_name    , & ! name of dimension time
                                 var_hru_name     , & ! name of varibale HRUs
                                 dim_hru_name     , & ! name of dimension HRUs
                                 dim_ylat_name    , & ! name of dimension lat in case of a 2D input varibale
                                 dim_xlon_name    , & ! name of dimension lon in case of a 2D input varibale
                                 ! output
                                 nSpace           , & ! nSpace of the input in runoff or wm strcuture
                                 nTime            , & ! nTime of the input in runoff or wm strcuture
                                 sim              , & ! 1D simulation
                                 sim2D            , & ! 2D simulation
                                 ID_array         , & ! ID of seg or hru in data
                                 timeUnits        , & ! time units
                                 calendar         , & ! calendar
                                 ! error control
                                 ierr, message)     ! output: error control
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
 integer(i4b),                intent(out)   :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
 integer(i4b),                intent(out)   :: nTime           ! nTime of the input in runoff or wm strcuture
 real(dp),      allocatable,  intent(out)   :: sim(:)          ! 1D simulation
 real(dp),      allocatable,  intent(out)   :: sim2D(:,:)      ! 2D simulation
 integer(i4b),  allocatable,  intent(out)   :: ID_array(:)     ! ID of seg or hru in data
 character(*),                intent(out)   :: timeUnits       ! time units
 character(*),                intent(out)   :: calendar        ! calendar
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
  case(2); call read_1D_runoff_metadata(fname,            & ! name of dile including its directory
                                        var_time_name,    & ! name of varibale time
                                        dim_time_name,    & ! dimension of varibale time
                                        var_hru_name,     & ! name of varibale hrus
                                        dim_hru_name,     & ! dimension of varibales hrus
                                        nSpace        , & ! nSpace of the input in runoff or wm strcuture
                                        nTime           , & ! nTime of the input in runoff or wm strcuture
                                        sim             , & ! 1D simulation
                                        sim2D           , & ! 2D simulation
                                        ID_array        , & ! ID of seg or hru in data
                                        timeUnits,        &
                                        calendar,         &
                                        ierr, cmessage)
  case(3); call read_2D_runoff_metadata(fname,            & ! name of file including its directory
                                        var_time_name,    & ! name of varibale time
                                        dim_time_name,    & ! dimension of varibale time
                                        dim_ylat_name,    & ! dimesnion of lat
                                        dim_xlon_name,    & ! dimension of lon
                                        nSpace        , & ! nSpace of the input in runoff or wm strcuture
                                        nTime           , & ! nTime of the input in runoff or wm strcuture
                                        sim             , & ! 1D simulation
                                        sim2D           , & ! 2D simulation
                                        ID_array        , & ! ID of seg or hru in data
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
                                   fname            , & ! filename
                                   var_time_name    , & ! name of varibale time
                                   dim_time_name    , & ! name of dimension time
                                   var_hru_name     , & ! name of varibale HRUs
                                   dim_hru_name     , & ! name of dimension HRUs
                                   ! output
                                   nSpace           , & ! nSpace of the input in runoff or wm strcuture
                                   nTime            , & ! nTime of the input in runoff or wm strcuture
                                   sim              , & ! 1D simulation
                                   sim2D            , & ! 2D simulation
                                   ID_array         , & ! ID of seg or hru in data
                                   timeUnits        , & ! time units
                                   calendar         , & ! calendar
                                   ! error control
                                   ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*),               intent(in)          :: fname           ! filename
 character(*),               intent(in)          :: var_time_name   ! name of the varibale time
 character(*),               intent(in)          :: dim_time_name   ! name of dimension for time
 character(*),               intent(in)          :: var_hru_name    ! name of the varibale hru
 character(*),               intent(in)          :: dim_hru_name    ! name of dimension for hydrological HRUs
 ! output variables
 integer(i4b),               intent(out)         :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
 integer(i4b),               intent(out)         :: nTime           ! nTime of the input in runoff or wm strcuture
 real(dp),     allocatable,  intent(out)         :: sim(:)          ! 1D simulation
 real(dp),     allocatable,  intent(out)         :: sim2D(:,:)      ! 2D simulation
 integer(i4b), allocatable,  intent(out)         :: ID_array(:)     ! ID of seg or hru in data
 character(*),               intent(out)         :: timeUnits       ! time units
 character(*),               intent(out)         :: calendar        ! calendar
 ! error control
 integer(i4b),               intent(out)         :: ierr            ! error code
 character(*),               intent(out)         :: message         ! error message
 ! local variables
 character(len=strLen)                           :: cmessage        ! error message from subroutine


 ! initialize error control
 ierr=0; message='read_1D_runoff_metadata/'

 nSpace(2) = integerMissing

 ! get the number of HRUs
 call get_nc_dim_len(fname, trim(dim_hru_name), nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get number of time steps from the runoff file
 call get_nc_dim_len(fname, trim(dim_time_name), nTime, ierr, cmessage)
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
 allocate(ID_array(nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating ID'; return; endif

 ! allocate space for simulated runoff
 allocate(sim(nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating read_struct_data%sim'; return; endif

 ! get HRU ids from the runoff file
 call get_nc(fname, var_hru_name, ID_array, 1, nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine read_1D_runoff_metadata

 ! *****
 ! private subroutine: get 3D runoff (lat, lon, time) metadata...
 ! ******************************************
 subroutine read_2D_runoff_metadata(&
                                    ! input
                                    fname            , & ! filename
                                    var_time_name    , & ! name of varibale time
                                    dim_time_name    , & ! name of dimension time
                                    dim_ylat_name    , & ! name of varibale HRUs
                                    dim_xlon_name    , & ! name of dimension HRUs
                                    ! output
                                    nSpace           , & ! nSpace of the input in runoff or wm strcuture
                                    nTime            , & ! nTime of the input in runoff or wm strcuture
                                    sim              , & ! 1D simulation
                                    sim2D            , & ! 2D simulation
                                    ID_array         , & ! ID of seg or hru in data
                                    timeUnits        , & ! time units
                                    calendar         , & ! calendar
                                    ! error control
                                    ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*),                  intent(in)   :: fname           ! filename
 character(*),                  intent(in)   :: var_time_name   ! name of the varibale time
 character(*),                  intent(in)   :: dim_time_name   ! name of dimension for time
 character(*),                  intent(in)   :: dim_ylat_name   ! name of dimension along lat
 character(*),                  intent(in)   :: dim_xlon_name   ! name of dimension along lon
 ! output variables
 integer(i4b),                  intent(out)  :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
 integer(i4b),                  intent(out)  :: nTime           ! nTime of the input in runoff or wm strcuture
 real(dp),     allocatable,     intent(out)  :: sim(:)          ! 1D simulation
 real(dp),     allocatable,     intent(out)  :: sim2D(:,:)      ! 2D simulation
 integer(i4b), allocatable,     intent(out)  :: ID_array(:)     ! ID of seg or hru in data
 character(*),                  intent(out)  :: timeUnits       ! time units
 character(*),                  intent(out)  :: calendar        ! calendar
 ! error control
 integer(i4b),                  intent(out)  :: ierr            ! error code
 character(*),                  intent(out)  :: message         ! error message
 ! local variables
 character(len=strLen)                       :: cmessage        ! error message from subroutine


 ! initialize error control
 ierr=0; message='read_2D_runoff_metadata/'

 ! get number of time steps from the runoff file
 call get_nc_dim_len(fname, trim(dim_time_name), nTime, ierr, cmessage)
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
 call get_nc_dim_len(fname, trim(dim_ylat_name), nSpace(1), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get size of xlon dimension
 call get_nc_dim_len(fname, trim(dim_xlon_name), nSpace(2), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for simulated runoff. sim2d = runoff(lon, lat)
 allocate(sim2d(nSpace(2),nSpace(1)), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating sim'; return; endif

 end subroutine read_2D_runoff_metadata


 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 subroutine read_runoff_data(fname,            &  ! input: filename
                             var_name,         &  ! input: varibale name
                             time_index,       &  ! input: time index
                             nSpace,           &  ! input: dimension of data to be read
                             sim,              &  ! input/output: read data 1D sim
                             sim2D,            &  ! input/output: read data 2D sim
                             ierr, message)       ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname              ! filename
 character(*), intent(in)                :: var_name           ! variable name
 integer(i4b), intent(in)                :: time_index         ! index of time element
 integer(i4b), intent(in)                :: nSpace(1:2)        ! dimension of data for one time step
 ! input/output variables
 real(dp), allocatable,  intent(inout)   :: sim(:)             ! runoff for one time step for all spatial dimension
 real(dp), allocatable,  intent(inout)   :: sim2D(:,:)         ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)               :: ierr               ! error code
 character(*), intent(out)               :: message            ! error message
 ! local variables
 character(len=strLen)                   :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_runoff_data/'

 if (nSpace(2) == integerMissing) then
  call read_1D_runoff(fname, var_name, time_index, nSpace(1), sim, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 else
  call read_2D_runoff(fname, var_name, time_index, nSpace, sim2D, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif

 end subroutine read_runoff_data

 ! *********************************************************************
 ! private subroutine: read 2D runoff data
 ! *********************************************************************
 subroutine read_1D_runoff(fname,          &  ! input: filename
                           var_name,       &  ! input: variable name
                           time_index,     &  ! input: time index
                           nSpace,         &  ! input: size of HRUs
                           sim,            &  ! input/output: runoff data structure
                           ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname              ! filename
 character(*), intent(in)                :: var_name           ! variable name
 integer(i4b), intent(in)                :: time_index         ! index of time element
 integer(i4b), intent(in)                :: nSpace             ! size of spatial dimensions
 ! input/output variables
 real(dp), allocatable,  intent(inout)   :: sim(:)             ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)               :: ierr               ! error code
 character(*), intent(out)               :: message            ! error message
 ! local variables
 integer(i4b)                            :: iStart(2)
 integer(i4b)                            :: iCount(2)
 logical(lgt)                            :: existFillVal
 real(dp)                                :: dummy(nSpace,1)    ! data read
 character(len=strLen)                   :: cmessage           ! error message from subroutine

 ! initialize error control
 ierr=0; message='read_1D_runoff/'

 ! get the simulated runoff data
 iStart = [1,time_index]
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
 sim(1:nSpace) = dummy(1:nSpace,1)

 end subroutine read_1D_runoff

 ! *********************************************************************
 ! private subroutine: read 2D runoff data
 ! *********************************************************************
 subroutine read_2D_runoff(fname,            &  ! input: filename
                           var_name,         &  ! input: variable name
                           time_index,       &  ! input: time index
                           nSpace,           &  ! input: size of HRUs
                           sim2D,            &  ! input/output: runoff data structure
                           ierr, message)       ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                :: fname            ! filename
 character(*), intent(in)                :: var_name         ! variable name
 integer(i4b), intent(in)                :: time_index       ! index of time element
 integer(i4b), intent(in)                :: nSpace(1:2)      ! size of spatial dimensions
 ! input/output variables
 real(dp), allocatable,  intent(inout)   :: sim2D(:,:)       ! runoff for one time step for all spatial dimension
 ! output variables
 integer(i4b), intent(out)               :: ierr             ! error code
 character(*), intent(out)               :: message          ! error message
 ! local variables
 logical(lgt)                            :: existFillVal
 integer(i4b)                            :: iStart(3)
 integer(i4b)                            :: iCount(3)
 real(dp)                                :: dummy(nSpace(2),nSpace(1),1) ! data read
 character(len=strLen)                   :: cmessage                     ! error message from subroutine


 ! initialize error control
 ierr=0; message='read_2D_runoff/'

 ! get the simulated runoff data
 iStart = [1,1,time_index]
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
 sim2d(1:nSpace(2),1:nSpace(1)) = dummy(1:nSpace(2),1:nSpace(1),1)

 end subroutine read_2D_runoff

END MODULE read_runoff
