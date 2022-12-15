MODULE read_runoff

USE netcdf
USE nrtype
USE public_var
USE ncio_utils, ONLY: open_nc
USE ncio_utils, ONLY: get_nc
USE ncio_utils, ONLY: get_var_attr
USE ncio_utils, ONLY: check_attr
USE ncio_utils, ONLY: get_nc_dim_len
USE dataTypes,  ONLY: inputData           ! input data class (runoff and wm  inheritated)
USE dataTypes,  ONLY: map_time            ! time step mapping data type
USE dataTypes,  ONLY: inFileInfo          ! input file metadata


implicit none

private
public::read_runoff_metadata
public::read_forcing_data

CONTAINS

 ! *****
 ! public subroutine: get runoff  metadata...
 ! ******************************************
 SUBROUTINE read_runoff_metadata(fname            , & ! input: filename including directory
                                 var_name         , & ! input: varibale name
                                 var_time_name    , & ! input: name of varibale time
                                 dim_time_name    , & ! input: name of dimension time
                                 var_hru_name     , & ! input: name of varibale HRUs
                                 dim_hru_name     , & ! input: name of dimension HRUs
                                 dim_ylat_name    , & ! input: name of dimension lat in case of a 2D input varibale
                                 dim_xlon_name    , & ! input: name of dimension lon in case of a 2D input varibale
                                 nSpace           , & ! output: nSpace of the input in runoff or wm strcuture
                                 nTime            , & ! output: nTime of the input in runoff or wm strcuture
                                 sim              , & ! output: 1D simulation
                                 sim2D            , & ! output: 2D simulation
                                 ID_array         , & ! output: ID of seg or hru in data
                                 timeUnits        , & ! output: time units
                                 calendar         , & ! output: calendar
                                 ierr, message)       ! output: error control
 implicit none
 ! argument variables
 character(*),              intent(in)    :: fname           ! filename
 character(*),              intent(in)    :: var_name        ! name of the varibale for simulated runoff or abstraction/injection
 character(*),              intent(in)    :: var_time_name   ! name of the varibale time
 character(*),              intent(in)    :: dim_time_name   ! name of dimension for time
 character(*),              intent(in)    :: var_hru_name    ! name of the varibale hru
 character(*),              intent(in)    :: dim_hru_name    ! name of dimension for hydrological HRUs
 character(*),              intent(in)    :: dim_ylat_name   ! name of dimension along lat
 character(*),              intent(in)    :: dim_xlon_name   ! name of dimension along lon
 integer(i4b),              intent(out)   :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
 integer(i4b),              intent(out)   :: nTime           ! nTime of the input in runoff or wm strcuture
 real(dp),     allocatable, intent(out)   :: sim(:)          ! 1D simulation
 real(dp),     allocatable, intent(out)   :: sim2D(:,:)      ! 2D simulation
 integer(i4b), allocatable, intent(out)   :: ID_array(:)     ! ID of seg or hru in data
 character(*),              intent(out)   :: timeUnits       ! time units
 character(*),              intent(out)   :: calendar        ! calendar
 integer(i4b),              intent(out)   :: ierr            ! error code
 character(*),              intent(out)   :: message         ! error message
 ! local variables
 integer(i4b)                             :: ncid            ! netcdf id
 integer(i4b)                             :: ivarID          ! variable id
 integer(i4b)                             :: nDims           ! number of dimension in runoff file
 character(len=strLen)                    :: cmessage        ! error message from subroutine

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
                                        nSpace,           & ! nSpace of the input in runoff or wm strcuture
                                        nTime,            & ! nTime of the input in runoff or wm strcuture
                                        sim,              & ! 1D simulation
                                        ID_array,         & ! ID of seg or hru in data
                                        timeUnits,        &
                                        calendar,         &
                                        ierr, cmessage)
  case(3); call read_2D_runoff_metadata(fname,            & ! name of file including its directory
                                        var_time_name,    & ! name of varibale time
                                        dim_time_name,    & ! dimension of varibale time
                                        dim_ylat_name,    & ! dimesnion of lat
                                        dim_xlon_name,    & ! dimension of lon
                                        nSpace,           & ! nSpace of the input in runoff or wm strcuture
                                        nTime,            & ! nTime of the input in runoff or wm strcuture
                                        sim2D,            & ! 2D simulation
                                        ID_array,         & ! ID of seg or hru in data
                                        timeUnits,        &
                                        calendar,         &
                                        ierr, cmessage)
  case default; ierr=20; message=trim(message)//'runoff input must be 2-dimension (e.g, [time, hru]) or 3-dimension (e.g., [time, lat, lon]'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE read_runoff_metadata

 ! *****
 ! private subroutine: get 2D runoff (hru, time) metadata...
 ! ******************************************
 SUBROUTINE read_1D_runoff_metadata(fname           , & ! input: filename
                                   var_time_name    , & ! input: name of varibale time
                                   dim_time_name    , & ! input: name of dimension time
                                   var_hru_name     , & ! input: name of varibale HRUs
                                   dim_hru_name     , & ! input: name of dimension HRUs
                                   nSpace           , & ! output: nSpace of the input in runoff or wm strcuture
                                   nTime            , & ! output: nTime of the input in runoff or wm strcuture
                                   sim              , & ! output: 1D simulation
                                   ID_array         , & ! output: ID of seg or hru in data
                                   timeUnits        , & ! output: time units
                                   calendar         , & ! output: calendar
                                   ierr, message)       ! output: error control
 implicit none
 ! argument variables
 character(*),               intent(in)          :: fname           ! filename
 character(*),               intent(in)          :: var_time_name   ! name of the varibale time
 character(*),               intent(in)          :: dim_time_name   ! name of dimension for time
 character(*),               intent(in)          :: var_hru_name    ! name of the varibale hru
 character(*),               intent(in)          :: dim_hru_name    ! name of dimension for hydrological HRUs
 integer(i4b),               intent(out)         :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
 integer(i4b),               intent(out)         :: nTime           ! nTime of the input in runoff or wm strcuture
 real(dp),     allocatable,  intent(out)         :: sim(:)          ! 1D simulation
 integer(i4b), allocatable,  intent(out)         :: ID_array(:)     ! ID of seg or hru in data
 character(*),               intent(out)         :: timeUnits       ! time units
 character(*),               intent(out)         :: calendar        ! calendar
 integer(i4b),               intent(out)         :: ierr            ! error code
 character(*),               intent(out)         :: message         ! error message
 ! local variables
 character(len=strLen)                           :: cmessage        ! error message from subroutine

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

 END SUBROUTINE read_1D_runoff_metadata

 ! *****
 ! private subroutine: get 3D runoff (lat, lon, time) metadata...
 ! ******************************************
 SUBROUTINE read_2D_runoff_metadata(fname            , & ! input: filename
                                    var_time_name    , & ! input: name of varibale time
                                    dim_time_name    , & ! input: name of dimension time
                                    dim_ylat_name    , & ! input: name of varibale HRUs
                                    dim_xlon_name    , & ! input: name of dimension HRUs
                                    nSpace           , & ! output: nSpace of the input in runoff or wm strcuture
                                    nTime            , & ! output: nTime of the input in runoff or wm strcuture
                                    sim2D            , & ! output: 2D simulation
                                    ID_array         , & ! output: ID of seg or hru in data
                                    timeUnits        , & ! output: time units
                                    calendar         , & ! output: calendar
                                    ierr, message)       ! output: error control
 implicit none
 ! argument variables
 character(*),                  intent(in)   :: fname           ! filename
 character(*),                  intent(in)   :: var_time_name   ! name of the varibale time
 character(*),                  intent(in)   :: dim_time_name   ! name of dimension for time
 character(*),                  intent(in)   :: dim_ylat_name   ! name of dimension along lat
 character(*),                  intent(in)   :: dim_xlon_name   ! name of dimension along lon
 integer(i4b),                  intent(out)  :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
 integer(i4b),                  intent(out)  :: nTime           ! nTime of the input in runoff or wm strcuture
 real(dp),     allocatable,     intent(out)  :: sim2D(:,:)      ! 2D simulation
 integer(i4b), allocatable,     intent(out)  :: ID_array(:)     ! ID of seg or hru in data
 character(*),                  intent(out)  :: timeUnits       ! time units
 character(*),                  intent(out)  :: calendar        ! calendar
 integer(i4b),                  intent(out)  :: ierr            ! error code
 character(*),                  intent(out)  :: message         ! error message
 ! local variables
 character(len=strLen)                       :: cmessage        ! error message from subroutine

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

 END SUBROUTINE read_2D_runoff_metadata

 ! *********************************************************************
 ! public subroutine: main interface for forcing data reading
 ! *********************************************************************
 SUBROUTINE read_forcing_data(indir,            &  ! input: forcing input directory
                              inFileInfo_in,    &  ! input: forcing input file metadata
                              var_name,         &  ! input: varibale name
                              tmap_sim_forc_in, &  ! input: time-step mapping between model and forcing
                              forcing_data_in,  &  ! inout: forcing data structure
                              ierr, message)       ! output: error control
 implicit none
 ! argument variables
 character(*),     intent(in)      :: indir              ! forcing input directory
 type(inFileInfo), intent(in)      :: inFileInfo_in(:)   ! input file (forcing or water-management) metadata
 character(*),     intent(in)      :: var_name           ! variable name
 type(map_time),   intent(in)      :: tmap_sim_forc_in   ! time-step mapping between model and forcing
 class(inputData), intent(inout)   :: forcing_data_in    ! forcing data structure
 integer(i4b),     intent(out)     :: ierr               ! error code
 character(*),     intent(out)     :: message            ! error message
 ! local variables
 character(len=strLen)             :: cmessage           ! error message from subroutine

 ierr=0; message='read_forcing_data/'

 if (forcing_data_in%nSpace(2) == integerMissing) then
  call read_1D_forcing(indir, inFileInfo_in, var_name, tmap_sim_forc_in, forcing_data_in%nSpace(1), forcing_data_in, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 else
  call read_2D_forcing(indir, inFileInfo_in, var_name, tmap_sim_forc_in, forcing_data_in%nSpace, forcing_data_in, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif

 END SUBROUTINE read_forcing_data

 ! *********************************************************************
 ! private subroutine: read 1D forcing data
 ! *********************************************************************
 SUBROUTINE read_1D_forcing(indir,            &  ! input: forcing input directory
                            inFileInfo_in,    &  ! input: forcing input file metadata
                            var_name,         &  ! input: variable name
                            tmap_sim_forc_in, &  ! input: time-step mapping between model and forcing
                            nSpace,           &  ! input: size of spatial elements (e.g., HRU or reach)
                            forc_data_in,     &  ! inout: forcing data structure
                            ierr, message)       ! output: error control

 ! handle time step aggregation if forcing time step is finer than simulation time step

 implicit none
 ! Argument variables
 character(*),       intent(in)    :: indir             ! input directory
 type(inFileInfo),   intent(in)    :: inFileInfo_in(:)  ! input file (forcing or water-management) metadata
 character(*),       intent(in)    :: var_name          ! variable name
 type(map_time),     intent(in)    :: tmap_sim_forc_in  ! time-step mapping between model and forcing
 integer(i4b),       intent(in)    :: nSpace            ! size of spatial dimensions
 class(inputData),   intent(inout) :: forc_data_in      ! forcing data structure
 integer(i4b),       intent(out)   :: ierr              ! error code
 character(*),       intent(out)   :: message           ! error message
 ! Local variables
 character(len=strLen)             :: fname             ! filename
 integer(i4b)                      :: ix,it             ! loop index
 integer(i4b)                      :: nTime             ! number of forcing time step within a simulation time-step
 integer(i4b)                      :: iStart(2)         ! first indices in the variable to be read
 integer(i4b)                      :: iCount(2)         ! numbers of elements to be read
 logical(lgt)                      :: existFillVal      ! logical to indicate whether fillvalue exist in the variable
 real(dp)                          :: sumWeights        ! sum of time weight
 real(dp), allocatable             :: dummy(:,:)        ! array storing the read variable
 character(len=strLen)             :: cmessage          ! error message from subroutine

 ierr=0; message='read_1D_forcing/'

 nTime = size(tmap_sim_forc_in%iTime)

 allocate(dummy(nSpace, nTime), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the forcing data
 forc_data_in%sim(1:nSpace) = 0._dp
 do it = 1, nTime
   fname = trim(indir)//trim(inFileInfo_in(tmap_sim_forc_in%iFile(it))%infilename)

   iStart = [1, tmap_sim_forc_in%iTime(it)]
   iCount = [nSpace,1]
   call get_nc(fname, var_name, dummy(1:nSpace,it:it), iStart, iCount, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end do

 ! get the _fill_values for forcing variable if exist
 existFillVal = check_attr(fname, var_name, '_FillValue')
 if (existFillval) then
   call get_var_attr(fname, var_name, '_FillValue', input_fillvalue, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! finalize
 if (nTime>1) then
   do ix = 1, nSpace
     sumWeights = 0._dp
     do it = 1, nTime
       if (abs(dummy(ix,it)-input_fillvalue) < verySmall) cycle
       sumWeights = sumWeights + tmap_sim_forc_in%frac(it)
       forc_data_in%sim(ix) = forc_data_in%sim(ix) + dummy(ix,it)*tmap_sim_forc_in%frac(it)
     end do
    if(abs(0._dp - sumWeights)<verySmall) forc_data_in%sim(ix) = realMissing
    if(sumWeights > 0._dp .and. sumWeights < 1.0_dp) forc_data_in%sim(ix) = forc_data_in%sim(ix) / sumWeights
   end do
 else
   where ( abs(dummy - input_fillvalue) < verySmall ) dummy = realMissing
   forc_data_in%sim(1:nSpace) = dummy(1:nSpace,1)
 end if

 END SUBROUTINE read_1D_forcing

 ! *********************************************************************
 ! private subroutine: read 2D forcing data
 ! *********************************************************************
 SUBROUTINE read_2D_forcing(indir,            &  ! input: input directory
                            inFileInfo_in,    &  ! input: meta for input file
                            var_name,         &  ! input: variable name
                            tmap_sim_forc_in, &  ! input: time-step mapping between model and forcing
                            nSpace,           &  ! input: size of HRUs
                            forc_data_in,     &  ! inout: forcing data structure
                            ierr, message)       ! output: error control

 ! handle time step aggregation if forcing time step is finer than simulation time step

 implicit none
 ! Argument variables
 character(*),       intent(in)      :: indir             ! input directory
 type(inFileInfo),   intent(in)      :: inFileInfo_in(:)  ! input file (forcing or water-management) meta data
 character(*),       intent(in)      :: var_name          ! variable name
 type(map_time),     intent(in)      :: tmap_sim_forc_in  ! time-step mapping between model and forcing
 integer(i4b),       intent(in)      :: nSpace(1:2)       ! size of spatial dimensions
 class(inputData),   intent(inout)   :: forc_data_in      ! forcing data structure
 integer(i4b),       intent(out)     :: ierr              ! error code
 character(*),       intent(out)     :: message           ! error message
 ! local variables
 character(len=strLen)               :: fname             ! filename
 integer(i4b)                        :: ix1,ix2,it        ! loop index
 integer(i4b)                        :: nTime             ! number of forcing time step within a simulation time-step
 integer(i4b)                        :: iStart(3)         ! first indices in the variable to be read
 integer(i4b)                        :: iCount(3)         ! numbers of elements to be read
 logical(lgt)                        :: existFillVal      ! logical to indicate whether fillvalue exist in the variable
 real(dp)                            :: sumWeights        ! sum of time weight
 real(dp), allocatable               :: dummy(:,:,:)      ! array storing the read variable
 character(len=strLen)               :: cmessage          ! error message from subroutine

 ierr=0; message='read_2D_forcing/'

 nTime = size(tmap_sim_forc_in%iTime)

 allocate(dummy(nSpace(2), nSpace(1), nTime), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the forcing data
 forc_data_in%sim2d(1:nSpace(2),1:nSpace(1)) = 0._dp
 do it = 1, nTime
   fname = trim(indir)//trim(inFileInfo_in(tmap_sim_forc_in%iFile(it))%infilename)

   iStart = [1,1,tmap_sim_forc_in%iTime(it)]
   iCount = [nSpace(2),nSpace(1),1]
   call get_nc(fname, var_name, dummy(1:nSpace(2),1:nSpace(1),it:it), iStart, iCount, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do

 ! get the _fill_values for forcing variable
 existFillVal = check_attr(fname, var_name, '_FillValue')
 if (existFillval) then
   call get_var_attr(fname, var_name, '_FillValue', input_fillvalue, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! finalize
 if (nTime>1) then ! simulation time step includes multiple forcing time steps
   do ix2 = 1, nSpace(2)
     do ix1 = 1, nSpace(1)
       sumWeights = 0._dp
       do it = 1, nTime
         if (abs(dummy(ix2,ix1,it)-input_fillvalue) < verySmall) cycle
         sumWeights = sumWeights + tmap_sim_forc_in%frac(it)
         forc_data_in%sim2d(ix2,ix1) = forc_data_in%sim2d(ix2,ix1) + dummy(ix2,ix1,it)*tmap_sim_forc_in%frac(it)
       end do
       if(abs(0._dp - sumWeights)<verySmall) forc_data_in%sim2d(ix2,ix1) = realMissing
       if(sumWeights > 0._dp .and. sumWeights < 1.0_dp) forc_data_in%sim2d(ix2,ix1) = forc_data_in%sim2d(ix2,ix1) / sumWeights
     end do
   end do
 else ! if simulation time step include one forcing time step
   where ( abs(dummy - input_fillvalue) < verySmall ) dummy = realMissing
   forc_data_in%sim2d(1:nSpace(2),1:nSpace(1)) = dummy(1:nSpace(2),1:nSpace(1),1)
 end if

 END SUBROUTINE read_2D_forcing

END MODULE read_runoff
