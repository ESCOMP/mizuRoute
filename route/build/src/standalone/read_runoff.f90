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

INTERFACE read_forcing_metadata
  module procedure read_1D_forcing_i4b_metadata
  module procedure read_1D_forcing_i8b_metadata
  module procedure read_2D_forcing_metadata
END INTERFACE

private
public::read_forcing_metadata
public::read_forcing_data

CONTAINS

! *****
! private subroutine: get 2D runoff (hru, time) with long integer runoff hru_id  metadata...
! ******************************************
  SUBROUTINE read_1D_forcing_i4b_metadata(fname            , & ! input: filename
                                          var_name         , & ! input: varibale name
                                          var_hru_name     , & ! input: name of varibale HRUs
                                          dim_hru_name     , & ! input: name of dimension HRUs
                                          dim_ylat_name    , & ! input: name of dimension lat in case of a 2D input varibale
                                          dim_xlon_name    , & ! input: name of dimension lon in case of a 2D input varibale
                                          nSpace           , & ! output: nSpace of the input in runoff or wm strcuture
                                          sim              , & ! output: 1D simulation
                                          ID_array         , & ! output: ID of seg or hru in data
                                          fillvalue        , & ! output: fillvalue for data
                                          ierr, message)       ! output: error control
    implicit none
    ! argument variables
    character(*),              intent(in)    :: fname           ! filename
    character(*),              intent(in)    :: var_name        ! name of the varibale for simulated runoff or abstraction/injection
    character(*),              intent(in)    :: var_hru_name    ! name of the varibale hru
    character(*),              intent(in)    :: dim_hru_name    ! name of dimension for hydrological HRUs
    character(*),              intent(in)    :: dim_ylat_name   ! name of dimension along lat
    character(*),              intent(in)    :: dim_xlon_name   ! name of dimension along lon
    integer(i4b),              intent(out)   :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
    real(dp),     allocatable, intent(out)   :: sim(:)          ! input forcing data - 1D
    integer(i4b), allocatable, intent(out)   :: ID_array(:)     ! ID of seg or hru in data
    real(dp),                  intent(out)   :: fillvalue       ! fillvalue for data
    integer(i4b),              intent(out)   :: ierr            ! error code
    character(*),              intent(out)   :: message         ! error message
    ! local variables
    character(len=strLen)                    :: cmessage        ! error message from subroutine

    ierr=0; message='read_1D_forcing_metadata/'

    call read_forcing_attrs(fname, var_name, fillvalue, ierr, cmessage)

    ! get the number of HRUs
    call get_nc_dim_len(fname, trim(dim_hru_name), nSpace(1), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! allocate space for hru_id
    allocate(ID_array(nSpace(1)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating ID'; return; endif

    ! allocate space for simulated runoff
    allocate(sim(nSpace(1)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating read_struct_data%sim'; return; endif

    ! get HRU ids from the runoff file
    call get_nc(fname, var_hru_name, ID_array, 1, nSpace(1), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  END SUBROUTINE read_1D_forcing_i4b_metadata

! *****
! private subroutine: get 2D runoff (hru, time) with integer runoff hru_id  metadata...
! ******************************************
  SUBROUTINE read_1D_forcing_i8b_metadata(fname            , & ! input: filename
                                          var_name         , & ! input: varibale name
                                          var_hru_name     , & ! input: name of varibale HRUs
                                          dim_hru_name     , & ! input: name of dimension HRUs
                                          dim_ylat_name    , & ! input: name of dimension lat in case of a 2D input varibale
                                          dim_xlon_name    , & ! input: name of dimension lon in case of a 2D input varibale
                                          nSpace           , & ! output: nSpace of the input in runoff or wm strcuture
                                          sim              , & ! output: 1D simulation
                                          ID_array         , & ! output: ID of seg or hru in data
                                          fillvalue        , & ! output: fillvalue for data
                                          ierr, message)       ! output: error control
    implicit none
    ! argument variables
    character(*),              intent(in)    :: fname           ! filename
    character(*),              intent(in)    :: var_name        ! name of the varibale for simulated runoff or abstraction/injection
    character(*),              intent(in)    :: var_hru_name    ! name of the varibale hru
    character(*),              intent(in)    :: dim_hru_name    ! name of dimension for hydrological HRUs
    character(*),              intent(in)    :: dim_ylat_name   ! name of dimension along lat
    character(*),              intent(in)    :: dim_xlon_name   ! name of dimension along lon
    integer(i4b),              intent(out)   :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
    real(dp),     allocatable, intent(out)   :: sim(:)          ! 1D simulation
    integer(i8b), allocatable, intent(out)   :: ID_array(:)     ! ID of seg or hru in data
    real(dp),                  intent(out)   :: fillvalue       ! fillvalue for data
    integer(i4b),              intent(out)   :: ierr            ! error code
    character(*),              intent(out)   :: message         ! error message
    ! local variables
    character(len=strLen)                    :: cmessage        ! error message from subroutine

    ierr=0; message='read_1D_forcing_i8b_metadata/'

    call read_forcing_attrs(fname, var_name, fillvalue, ierr, cmessage)

    ! get the number of HRUs
    call get_nc_dim_len(fname, trim(dim_hru_name), nSpace(1), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! allocate space for hru_id
    allocate(ID_array(nSpace(1)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating ID'; return; endif

    ! allocate space for simulated runoff
    allocate(sim(nSpace(1)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating read_struct_data%sim'; return; endif

    ! get HRU ids from the runoff file
    call get_nc(fname, var_hru_name, ID_array, 1, nSpace(1), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  END SUBROUTINE read_1D_forcing_i8b_metadata

! *****
! private subroutine: get 3D runoff (lat, lon, time) metadata...
! ******************************************
  SUBROUTINE read_2D_forcing_metadata(fname            , & ! input: filename
                                      var_name         , & ! input: varibale name
                                      var_hru_name     , & ! input: name of varibale HRUs
                                      dim_hru_name     , & ! input: name of dimension HRUs
                                      dim_ylat_name    , & ! input: name of dimension lat in case of a 2D input varibale
                                      dim_xlon_name    , & ! input: name of dimension lon in case of a 2D input varibale
                                      nSpace           , & ! output: nSpace of the input in runoff or wm strcuture
                                      sim              , & ! output: 1D simulation
                                      ID_array         , & ! output: ID of seg or hru in data
                                      fillvalue        , & ! output: fillvalue for data
                                      ierr, message)       ! output: error control
    implicit none
    ! argument variables
    character(*),              intent(in)    :: fname           ! filename
    character(*),              intent(in)    :: var_name        ! name of the varibale for simulated runoff or abstraction/injection
    character(*),              intent(in)    :: var_hru_name    ! name of the varibale hru
    character(*),              intent(in)    :: dim_hru_name    ! name of dimension for hydrological HRUs
    character(*),              intent(in)    :: dim_ylat_name   ! name of dimension along lat
    character(*),              intent(in)    :: dim_xlon_name   ! name of dimension along lon
    integer(i4b),              intent(out)   :: nSpace(1:2)     ! nSpace of the input in runoff or wm strcuture
    real(dp),     allocatable, intent(out)   :: sim(:,:)        ! 2D simulation
    integer(i8b), allocatable, intent(out)   :: ID_array(:)     ! ID of seg or hru in data
    real(dp),                  intent(out)   :: fillvalue       ! fillvalue for data
    integer(i4b),              intent(out)   :: ierr            ! error code
    character(*),              intent(out)   :: message         ! error message
    ! local variables
    character(len=strLen)                    :: cmessage        ! error message from subroutine

    ierr=0; message='read_2D_forcing_metadata/'

    call read_forcing_attrs(fname, var_name, fillvalue, ierr, cmessage)

    ! get size of ylat dimension
    call get_nc_dim_len(fname, trim(dim_ylat_name), nSpace(1), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! get size of xlon dimension
    call get_nc_dim_len(fname, trim(dim_xlon_name), nSpace(2), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! allocate space for simulated runoff. sim2d = runoff(lon, lat)
    allocate(sim(nSpace(2),nSpace(1)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating sim'; return; endif

  END SUBROUTINE read_2D_forcing_metadata

! *****
! private subroutine: get forcing attributes ...
! ******************************************
  SUBROUTINE read_forcing_attrs(fname            , & ! input: filename including directory
                                var_name         , & ! input: varibale name
                                fillvalue        , & ! output: fillvalue for data
                                ierr, message)       ! output: error control
    implicit none
    ! argument variables
    character(*),              intent(in)    :: fname           ! filename
    character(*),              intent(in)    :: var_name        ! name of the varibale for simulated runoff or abstraction/injection
    real(dp),                  intent(out)   :: fillvalue       ! fillvalue for data
    integer(i4b),              intent(out)   :: ierr            ! error code
    character(*),              intent(out)   :: message         ! error message
    ! local variables
    integer(i4b)                             :: ncid            ! netcdf id
    integer(i4b)                             :: ivarID          ! variable id
    logical(lgt)                             :: existAttr       ! logical to indicate whether attribute exists in the variable
    character(len=strLen)                    :: cmessage        ! error message from subroutine

    ierr=0; message='read_forcing_attrs/'

    ! open NetCDF file
    call open_nc(trim(fname), 'r', ncid, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! get the ID of runoff variable
    ierr = nf90_inq_varid(ncid, trim(var_name), ivarID)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

    ! get the _fill_values for forcing variable if exist
    existAttr = check_attr(fname, var_name, '_FillValue')
    if (existAttr) then
      call get_var_attr(fname, var_name, '_FillValue', fillvalue, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    else
      if (input_fillvalue==realMissing) then
        write(iulog, '(3A)')       'WARNING:', trim(var_name), '. No _FillValue exist in attributes nor provided by user in input_fillvalue'
        write(iulog, '(A,x,G15.4)') '         _FillValue is set to default:', realMissing
      end if
      fillvalue = input_fillvalue
    end if

  END SUBROUTINE read_forcing_attrs

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

 ! finalize
 if (nTime>1) then
   do ix = 1, nSpace
     sumWeights = 0._dp
     do it = 1, nTime
       if (abs(dummy(ix,it)-forc_data_in%fillvalue) < verySmall) cycle
       sumWeights = sumWeights + tmap_sim_forc_in%frac(it)
       forc_data_in%sim(ix) = forc_data_in%sim(ix) + dummy(ix,it)*tmap_sim_forc_in%frac(it)
     end do
    if(abs(0._dp - sumWeights)<verySmall) forc_data_in%sim(ix) = realMissing
    if(sumWeights > 0._dp .and. sumWeights < 1.0_dp) forc_data_in%sim(ix) = forc_data_in%sim(ix) / sumWeights
   end do
 else
   where ( abs(dummy - forc_data_in%fillvalue) < verySmall ) dummy = realMissing
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
 real(dp)                            :: sumWeights        ! sum of time weight
 real(dp), allocatable               :: dummy(:,:,:)      ! array storing the read variable
 character(len=strLen)               :: cmessage          ! error message from subroutine

 ierr=0; message='read_2D_forcing/'

 nTime = size(tmap_sim_forc_in%iTime)

 allocate(dummy(nSpace(2), nSpace(1), nTime), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the forcing data
 do it = 1, nTime
   fname = trim(indir)//trim(inFileInfo_in(tmap_sim_forc_in%iFile(it))%infilename)

   iStart = [1,1,tmap_sim_forc_in%iTime(it)]
   iCount = [nSpace(2),nSpace(1),1]
   call get_nc(fname, var_name, dummy(1:nSpace(2),1:nSpace(1),it:it), iStart, iCount, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do

 forc_data_in%sim2d(1:nSpace(2),1:nSpace(1)) = 0._dp
 ! finalize
 if (nTime>1) then ! simulation time step includes multiple forcing time steps
   do ix2 = 1, nSpace(2)
     do ix1 = 1, nSpace(1)
       sumWeights = 0._dp
       do it = 1, nTime
         if (abs(dummy(ix2,ix1,it)-forc_data_in%fillvalue) < verySmall) cycle
         sumWeights = sumWeights + tmap_sim_forc_in%frac(it)
         forc_data_in%sim2d(ix2,ix1) = forc_data_in%sim2d(ix2,ix1) + dummy(ix2,ix1,it)*tmap_sim_forc_in%frac(it)
       end do
       if(abs(0._dp - sumWeights)<verySmall) forc_data_in%sim2d(ix2,ix1) = realMissing
       if(sumWeights > 0._dp .and. sumWeights < 1.0_dp) forc_data_in%sim2d(ix2,ix1) = forc_data_in%sim2d(ix2,ix1) / sumWeights
     end do
   end do
 else ! if simulation time step include one forcing time step
   where ( abs(dummy - forc_data_in%fillvalue) < verySmall ) dummy = realMissing
   forc_data_in%sim2d(1:nSpace(2),1:nSpace(1)) = dummy(1:nSpace(2),1:nSpace(1),1)
 end if

 END SUBROUTINE read_2D_forcing

END MODULE read_runoff
