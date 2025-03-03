MODULE ncio_utils

USE nrtype
USE netcdf

implicit none

private

public::get_nc
public::get_var_dims
public::get_nc_dim_len
public::get_var_attr
public::check_attr
public::check_variable
public::put_global_attr
public::def_nc
public::def_dim
public::def_var
public::end_def
public::write_nc
public::open_nc
public::close_nc

INTERFACE get_nc
  module procedure get_scalar
  module procedure get_vec
  module procedure get_charvec
  module procedure get_2d_array
  module procedure get_3d_array
  module procedure get_4d_array
END INTERFACE

INTERFACE write_nc
  module procedure write_vec
  module procedure write_2d_array
  module procedure write_3d_array
END INTERFACE

! public netCDF parameter
integer(i4b),parameter,public :: ncd_short     = nf90_short
integer(i4b),parameter,public :: ncd_int       = nf90_int
integer(i4b),parameter,public :: ncd_float     = nf90_float
integer(i4b),parameter,public :: ncd_double    = nf90_double
integer(i4b),parameter,public :: ncd_char      = nf90_char
integer(i4b),parameter,public :: ncd_unlimited = nf90_unlimited
integer(i4b),parameter,public :: ncd_global    = nf90_global

CONTAINS

 ! *********************************************************************
 ! subroutine: Define netcdf file
 ! *********************************************************************
 SUBROUTINE def_nc(fname,          &  ! input: file name
                   ncid,           &  ! error control
                   ierr,message,   &  ! error control
                   nctype)            ! netCDF type
 implicit none
 ! input variables
 character(*),           intent(in)  :: fname          ! filename
 character(*), optional, intent(in)  :: nctype         ! netCDF type
 ! output variables
 integer(i4b),           intent(out) :: ncid           ! netcdf id
 integer(i4b),           intent(out) :: ierr           ! error code
 character(*),           intent(out) :: message        ! error message
 ! local variables
 character(len=strLen)               :: nctype_local   ! local string for netCDF type name
 integer(i4b)                        :: nctypeID       ! netCDF type ID

 ! initialize error control
 ierr=0; message='def_nc/'

 if (present(nctype)) then
   nctype_local = nctype
 else
   nctype_local = '64bit_offset'
 end if

 select case(trim(nctype_local))
   case('64bit_offset'); nctypeID = nf90_64bit_offset
   case('netcdf4');      nctypeID = nf90_netcdf4
   case('classic');      nctypeID = nf90_classic_model
   case default; ierr=20; message=trim(message)//'unable to identify netCDF type'; return
 end select

 ierr = nf90_create(trim(fname), nctypeID, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 END SUBROUTINE def_nc


 ! *********************************************************************
 ! subroutine: get dimension name and length for given 2D variable
 ! *********************************************************************
 subroutine get_var_dims(fname,          &  ! input: file name
                         vname,          &  ! input: variable name
                         ierr,message,   &  ! error control
                         dlen,           &  ! dimension length
                         dname)             ! dimension name
 implicit none
 ! input variables
 character(*), intent(in)                :: fname           ! filename
 character(*), intent(in)                :: vname           ! dimension name
 ! output variables
 integer(i4b), optional, intent(out)     :: dlen(:)         ! dimension length
 character(*), optional, intent(out)     :: dname(:)        ! dimension name
 integer(i4b), intent(out)               :: ierr            ! error code
 character(*), intent(out)               :: message         ! error message
 ! local variables
 integer(i4b)                            :: ncid            ! NetCDF file ID
 integer(i4b)                            :: ivarID          ! variable ID
 integer(i4b), allocatable               :: ncDimIDs(:)     ! dimension IDs for a given variable
 integer(i4b)                            :: nDims           ! number of dimensions in a variable
 integer(i4b)                            :: ii              ! loop indix

 ! initialize error control
 ierr=0; message='get_var_dims/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//'['//trim(nf90_strerror(ierr))//'; file='//trim(fname)//']'; return; endif

 ! get the ID of the variable
 ierr = nf90_inq_varid(ncid, trim(vname), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the number of dimesion in variable
 ierr = nf90_inquire_variable(ncid, ivarID, ndims=nDims)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! check size of output arguments and number of dimensions in the variable
 if (present(dlen)) then
   if(nDims/=size(dlen)) then; ierr=10;  message=trim(message)//': '//trim(vname)//' number of dimensions mismatch with dlen size'; return; endif
 end if
 if (present(dname)) then
   if(nDims/=size(dname)) then; ierr=10; message=trim(message)//': '//trim(vname)//' number of dimensions mismatch with dname size'; return; endif
 end if

 allocate(ncDimIDs(nDims), stat=ierr)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating spatial dimension for data structures'; return; endif

 ! get the dimension IDs
 ierr = nf90_inquire_variable(ncid, ivarID, dimids=ncDimIDs(:nDims))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the dimension name and length
 do ii = 1,nDims
  if ( present(dlen) .and. present(dname) )then
    ierr = nf90_inquire_dimension(ncid, ncDimIDs(ii), len=dlen(ii), name=dname(ii))
  else if ( present(dlen) .and. .not. present(dname) )then
    ierr = nf90_inquire_dimension(ncid, ncDimIDs(ii), len=dlen(ii))
  else if ( .not. present(dlen) .and. present(dname) )then
    ierr = nf90_inquire_dimension(ncid, ncDimIDs(ii), name=dname(ii))
  endif
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; varname='//trim(vname); return; endif
 end do

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_var_dims

 ! *********************************************************************
 ! subroutine: get vector dimension from netCDF
 ! *********************************************************************
 subroutine get_nc_dim_len(fname,           &  ! input: filename
                           dname,           &  ! input: variable name
                           nDim,            &  ! output: Size of dimension
                           ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: dname        ! dimension name
  ! output variables
  integer(i4b), intent(out)       :: nDim         ! size of dimension
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iDimID       ! NetCDF dimension ID
  ! initialize error control
  ierr=0; message='get_nc_dim_len/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//'['//trim(nf90_strerror(ierr))//'; file='//trim(fname)//']'; return; endif

  ! get the ID of the dimension
  ierr = nf90_inq_dimid(ncid, dname, iDimID)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname); return; endif

  ! get the length of the dimension
  ierr = nf90_inquire_dimension(ncid, iDimID, len=nDim)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: check if variable exists
 ! *********************************************************************
 function check_variable(ncid, vname)
   implicit none
   ! input
   integer(i4b), intent(in)        :: ncid         ! NetCDF file ID
   character(*), intent(in)        :: vname        ! variable name
   logical(lgt)                    :: check_variable
   ! local
   integer(i4b)                    :: ierr         ! error code
   integer(i4b)                    :: iVarID       ! variable ID

   ! get the ID of the variable
   ierr = nf90_inq_varid(ncid, trim(vname), iVarID)
   check_variable = (ierr == nf90_noerr)

 end function check_variable

 ! *********************************************************************
 ! subroutine: check if attribute exists
 ! *********************************************************************
 FUNCTION check_attr(fname, vname, attr_name)
   implicit none
   ! input
   character(*), intent(in)        :: fname        ! filename
   character(*), intent(in)        :: vname        ! variable name
   character(*), intent(in)        :: attr_name    ! attribute name
   logical(lgt)                    :: check_attr
   ! local
   integer(i4b)                    :: ierr         ! error code
   integer(i4b)                    :: ncid         ! NetCDF file ID
   integer(i4b)                    :: iVarID       ! variable ID

   ! open file for reading
   ierr = nf90_open(fname, nf90_nowrite, ncid)

   ! get the ID of the variable
   ierr = nf90_inq_varid(ncid, trim(vname), iVarID)

   ierr = nf90_inquire_attribute(ncid, iVarID, attr_name)
   check_attr = (ierr == nf90_noerr)

  ! close output file
  ierr = nf90_close(ncid)

 END FUNCTION check_attr


 ! *********************************************************************
 ! subroutine: get attribute values for a variable
 ! *********************************************************************
 subroutine get_var_attr(fname,           &  ! input: filename
                         vname,           &  ! input: variable name
                         attr_name,       &  ! inpu: attribute name
                         attr_value,      &  ! output: attribute value
                         ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname        ! variable name
 character(*), intent(in)        :: attr_name    ! attribute name
 ! output variables
 class(*),     intent(out)       :: attr_value   ! attribute value
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID

 ierr=0; message='get_var_attr/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the ID of the variable
 ierr = nf90_inq_varid(ncid, trim(vname), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the attribute value
 select type (attr_value)
   type is (integer(i4b))
     ierr = nf90_get_att(ncid, ivarID, attr_name, attr_value)
   type is (integer(i8b))
     ierr = nf90_get_att(ncid, ivarID, attr_name, attr_value)
   type is (real(sp))
     ierr = nf90_get_att(ncid, ivarID, attr_name, attr_value)
   type is (real(dp))
     ierr = nf90_get_att(ncid, ivarID, attr_name, attr_value)
   type is (character(len=*))
     ierr = nf90_get_att(ncid, ivarID, attr_name, attr_value)
   class default
     ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_var_attr


 ! *********************************************************************
 ! subroutine: double precision scalar value from netCDF
 ! *********************************************************************
 subroutine get_scalar(fname,           &  ! input:  filename
                       vname,           &  ! input:  variable name
                       array,           &  ! output: variable data
                       iStart,          &  ! input:  start index
                       ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname     ! filename
  character(*), intent(in)        :: vname     ! variable name
  integer(i4b), intent(in)        :: iStart    ! start index
  ! output variables
  class(*),     intent(out)       :: array     ! output variable data
  integer(i4b), intent(out)       :: ierr      ! error code
  character(*), intent(out)       :: message   ! error message
  ! local variables
  character(len=strLen)           :: cmessage   ! error message of downwind routine
  integer(i4b)                    :: vec_14b(1) ! 32bit int variable data
  integer(i8b)                    :: vec_18b(1) ! 64bit int variable data
  real(sp)                        :: vec_sp(1)  ! double precision float variable data
  real(dp)                        :: vec_dp(1)  ! single precision float variable data

 ierr=0; message='get_scalar/'

 select type (array)
   type is (integer(i4b))
     call get_vec(fname, vname, vec_14b, iStart, 1, ierr, cmessage)
     array = vec_14b(1)
   type is (integer(i8b))
     call get_vec(fname, vname, vec_18b, iStart, 1, ierr, cmessage)
     array = vec_18b(1)
   type is (real(sp))
     call get_vec(fname, vname, vec_sp, iStart, 1, ierr, cmessage)
     array = vec_sp(1)
   type is (real(dp))
     call get_vec(fname, vname, vec_dp, iStart, 1, ierr, cmessage)
     array = vec_dp(1)
   class default
     ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
 end select

 end subroutine get_scalar

 ! *********************************************************************
 ! subroutine: get vector (1D aray) value from netCDF
 ! *********************************************************************
 subroutine get_vec(fname,           &  ! input:  filename
                    vname,           &  ! input:  variable name
                    array,           &  ! output: variable data
                    iStart,          &  ! input:  start index
                    iCount,          &  ! input:  length of vector
                    ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname     ! filename
  character(*), intent(in)        :: vname     ! variable name
  integer(i4b), intent(in)        :: iStart    ! start index
  integer(i4b), intent(in)        :: iCount    ! length of vector to be read in
  ! output variables
  class(*),     intent(out)       :: array(:)  ! output variable data
  integer(i4b), intent(out)       :: ierr      ! error code
  character(*), intent(out)       :: message   ! error message
  ! local variables
  integer(i4b)                    :: ncid      ! NetCDF file ID
  integer(i4b)                    :: iVarID    ! NetCDF variable ID

 ierr=0; message='get_vec/'

 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the data
 select type (array)
   type is (integer(i4b))
    ierr = nf90_get_var(ncid, iVarID, array, start=(/iStart/), count=(/iCount/))
   type is (integer(i8b))
    ierr = nf90_get_var(ncid, iVarID, array, start=(/iStart/), count=(/iCount/))
   type is (real(sp))
     ierr = nf90_get_var(ncid, iVarID, array, start=(/iStart/), count=(/iCount/))
   type is (real(dp))
     ierr = nf90_get_var(ncid, iVarID, array, start=(/iStart/), count=(/iCount/))
   class default
     ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_vec

 ! *********************************************************************
 ! subroutine: read a integer 2D array
 ! *********************************************************************
 subroutine get_2d_array(fname,           &  ! input: filename
                         vname,           &  ! input: variable name
                         array,           &  ! output: variable data
                         iStart,          &  ! input: start index
                         iCount,          &  ! input: length of vector
                         ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  integer(i4b), intent(in)        :: iStart(1:2)  ! start indices
  integer(i4b), intent(in)        :: iCount(1:2)  ! length of vector
  ! output variables
  class(*),     intent(out)       :: array(:,:)   ! variable data
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID

  ! initialize error control
  ierr=0; message='get_2d_array/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  select type (array)
    type is (integer(i4b))
      ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (integer(i8b))
      ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(sp))
      ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(dp))
      ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    class default
      ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_2d_array


 ! *********************************************************************
 ! subroutine: read a 3D array
 ! *********************************************************************
 subroutine get_3d_array(fname,          &  ! input: filename
                         vname,          &  ! input: variable name
                         array,          &  ! output: variable data
                         iStart,         &  ! input: start index
                         iCount,         &  ! input: length of vector
                         ierr, message)     ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  integer(i4b), intent(in)        :: iStart(1:3)  ! start indices
  integer(i4b), intent(in)        :: iCount(1:3)  ! length of vector
  ! output variables
  class(*),     intent(out)       :: array(:,:,:) ! variable data
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID

  ierr=0; message='get_3d_array/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  select type (array)
    type is (integer(i4b))
      ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (integer(i8b))
      ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(sp))
      ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(dp))
      ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    class default
      ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_3d_array

 ! *********************************************************************
 ! subroutine: read a 4D array
 ! *********************************************************************
 subroutine get_4d_array(fname,          &  ! input: filename
                         vname,          &  ! input: variable name
                         array,          &  ! output: variable data
                         iStart,         &  ! input: start index
                         iCount,         &  ! input: length of vector
                         ierr, message)     ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  integer(i4b), intent(in)        :: iStart(1:4)  ! start indices
  integer(i4b), intent(in)        :: iCount(1:4)  ! length of vector
  ! output variables
  class(*),     intent(out)       :: array(:,:,:,:) ! variable data
  integer(i4b), intent(out)       :: ierr           ! error code
  character(*), intent(out)       :: message        ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID

  ierr=0; message='get_4d_array/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  select type (array)
    type is (integer(i4b))
     ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (integer(i8b))
     ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(sp))
     ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(dp))
     ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
    class default
      ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_4d_array

  ! *********************************************************************
  ! subroutine: read a character vector
  ! *********************************************************************
  subroutine get_charvec(fname,           &  ! input: netcdf id
                         vname,           &  ! input: variable name
                         array,           &  ! output: variable data
                         iStart,          &  ! input: start index
                         iCount,          &  ! input: length of vector
                         ierr, message)      ! output: error control
    implicit none
    ! input variables
    character(*), intent(in)        :: fname        ! filename
    character(*), intent(in)        :: vname       ! variable name
    integer(i4b), intent(in)        :: iStart(1:2) ! start index
    integer(i4b), intent(in)        :: iCount(1:2) ! length of vector to be read in
    ! output variables
    character(*), intent(out)       :: array(:)    ! output variable data
    integer(i4b), intent(out)       :: ierr        ! error code
    character(*), intent(out)       :: message     ! error message
    ! local variables
    integer(i4b)                    :: ncid        ! NetCDF file ID
    integer(i4b)                    :: iVarID      ! NetCDF variable ID

    ierr=0; message='get_charvec/'

    ! open NetCDF file
    ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

    ! get variable ID
    ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

    ! get the data
    ierr = nf90_get_var(ncid, iVarID, array, start=iStart, count=iCount)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

    array = adjustl(array)

    ! close output file
    ierr = nf90_close(ncid)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

  ! *********************************************************************
  ! subroutine: write global attribute
  ! *********************************************************************
  subroutine put_global_attr(ncid,          & ! input: netCDF ID
                             attName,       & ! input: global attribute name
                             attValue,      & ! input: global attribute values
                             ierr, message)   ! output: error control
  implicit none
  ! input variables
  integer(i4b), intent(in)        :: ncid      ! Input: netcdf fine ID
  character(*), intent(in)        :: attName   ! attribute name
  character(*), intent(in)        :: attValue  ! attribute values
  ! output variables
  integer(i4b), intent(out)       :: ierr          ! error code
  character(*), intent(out)       :: message       ! error message

  ! initialize error control
  ierr=0; message='put_global_attr/'

  ierr = nf90_put_att(ncid, ncd_global, trim(attName), trim(attValue))
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

  ! *********************************************************************
  ! subroutine: write an vector (1D array)
  ! *********************************************************************
  SUBROUTINE write_vec(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        array,           &  ! input: variable data
                        iStart,          &  ! input: start index
                        iCount,          &  ! input: length of vector
                        ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  class(*),     intent(in)        :: array(:)     ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start index
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID

  ierr=0; message='write_vec/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  select type (array)
    type is (integer(i4b))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (integer(i8b))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(sp))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(dp))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (character(len=*))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    class default
      ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  END SUBROUTINE write_vec

  ! *********************************************************************
  ! subroutine: write a 2D array
  ! *********************************************************************
  SUBROUTINE write_2d_array(fname,           &  ! input: filename
                            vname,           &  ! input: variable name
                            array,           &  ! input: variable data
                            iStart,          &  ! input: start index
                            iCount,          &  ! input: length of vector
                            ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  class(*),     intent(in)        :: array(:,:)   ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start indices
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID

  ierr=0; message='write_2d_array/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  select type (array)
    type is (integer(i4b))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (integer(i8b))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(sp))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(dp))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (character(len=*))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    class default
      ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  END SUBROUTINE write_2d_array

  ! *********************************************************************
  ! subroutine: write a 3D array
  ! *********************************************************************
  SUBROUTINE write_3d_array(fname,          &  ! input: filename
                            vname,          &  ! input: variable name
                            array,          &  ! input: variable data
                            iStart,         &  ! input: start index
                            iCount,         &  ! input: length of vector
                            ierr, message)     ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  class(*),     intent(in)        :: array(:,:,:) ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start indices
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID

  ierr=0; message='write_3d_array/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  select type (array)
    type is (integer(i4b))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (integer(i8b))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(sp))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (real(dp))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    type is (character(len=*))
      ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
    class default
      ierr = 1; message=trim(message)//'ERROR: invalid array type'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  END SUBROUTINE write_3d_array


 ! *********************************************************************
 ! Public subroutine: define variable attributes NetCDF file
 ! *********************************************************************
 SUBROUTINE def_var(ncid, vname, dimNames, ivtype, ierr, message, vdesc, vunit, vcal)

   implicit none
   ! input
   integer(i4b), intent(in)             :: ncid                   ! Input: netcdf ID
   character(*), intent(in)             :: vname                  ! Input: variable name
   character(*), intent(in)             :: dimNames(:)            ! Input: variable dimension names
   integer(i4b), intent(in)             :: ivtype                 ! Input: variable type
   character(*), intent(in), optional   :: vdesc                  ! Input: variable description
   character(*), intent(in), optional   :: vunit                  ! Input: variable units
   character(*), intent(in), optional   :: vcal                   ! Input: calendar (if time variable)
   ! output
   integer(i4b), intent(out)            :: ierr                   ! error code
   character(*), intent(out)            :: message                ! error message
   ! local
   integer(i4b)                         :: id                     ! loop through dimensions
   integer(i4b)                         :: dimIDs(size(dimNames)) ! vector of dimension IDs
   integer(i4b)                         :: iVarId                 ! variable ID

   ! initialize error control
   ierr=0; message='def_var/'

   ! define dimension IDs
   do id=1,size(dimNames)
    ierr=nf90_inq_dimid(ncid,trim(dimNames(id)),dimIDs(id))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
   end do

   ! define variable
   ierr = nf90_def_var(ncid,trim(vname),ivtype,dimIds,iVarId)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

   if (present(vdesc)) then ! add long_name
     ierr = nf90_put_att(ncid,iVarId,'long_name',trim(vdesc))
     if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
   end if

   if (present(vunit)) then ! add variable unit
     ierr = nf90_put_att(ncid,iVarId,'units',trim(vunit))
     if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
   end if

   if (present(vcal)) then ! add time calendar
     ierr = nf90_put_att(ncid,iVarId,'calendar',trim(vcal))
     if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
   end if

 END SUBROUTINE def_var


 ! *********************************************************************
 ! Public subroutine: define dimension NetCDF file
 ! *********************************************************************
 SUBROUTINE def_dim(ncid, dimName, dimLen, dimID, ierr, message)

   implicit none
   ! input
   integer(i4b), intent(in)  :: ncid     ! Input: netcdf fine ID
   character(*), intent(in)  :: dimName  ! Input: variable dimension name
   integer(i4b), intent(in)  :: dimLen   ! Input: dimension length
   ! output
   integer(i4b), intent(out) :: dimID    ! dimension id
   integer(i4b), intent(out) :: ierr     ! error code
   character(*), intent(out) :: message  ! error message

   ! initialize error control
   ierr=0; message='def_dim/'

   ierr = nf90_def_dim(ncid, trim(dimName), dimLen, dimId)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 END SUBROUTINE def_dim


 ! *********************************************************************
 ! Public subroutine: End defining netCDF file
 ! *********************************************************************
 SUBROUTINE end_def(ncid, ierr, message)

   implicit none
   ! input
   integer(i4b), intent(in)  :: ncid     ! Input: netcdf fine ID
   ! output
   integer(i4b), intent(out) :: ierr     ! error code
   character(*), intent(out) :: message  ! error message

   ! initialize error control
   ierr=0; message='end_def/'

   ierr = nf90_enddef(ncid)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 END SUBROUTINE end_def


 ! *********************************************************************
 ! Public subroutine: close netcdf
 ! *********************************************************************
 SUBROUTINE open_nc(fname, mode, ncid, ierr, message)

   implicit none
   ! input variables
   character(*), intent(in)     :: fname     ! filename
   character(1), intent(in)     :: mode      ! Input: mode ID
   ! output
   integer(i4b), intent(out)    :: ncid      ! output: netcdf ID
   integer(i4b), intent(out)    :: ierr      ! error code
   character(*), intent(out)    :: message   ! error message
   ! local
   integer(i4b)                 :: modeId    ! mode flag

   ! initialize error control
   ierr=0; message='open_nc/'

   select case(trim(mode))
     case('r','R'); modeId = NF90_NOWRITE
     case('w','W'); modeId = NF90_WRITE
     case('s','S'); modeId = NF90_SHARE
     case default; ierr=20; message=trim(message)//'netCDF open mode not recongnized; chose "r", "w", "s"'; return
   end select

   ierr = nf90_open(trim(fname), modeId, ncid)
   if(ierr/=0)then; message=trim(message)//'['//trim(nf90_strerror(ierr))//'; file='//trim(fname)//']'; return; endif

 END SUBROUTINE open_nc


 ! *********************************************************************
 ! Public subroutine: close netcdf
 ! *********************************************************************
 SUBROUTINE close_nc(ncid, ierr, message)

   implicit none
   ! input
   integer(i4b), intent(in)     :: ncid      ! Input: netcdf fine ID
   ! output
   integer(i4b), intent(out)    :: ierr      ! error code
   character(*), intent(out)    :: message   ! error message

   ! initialize error control
   ierr=0; message='close_nc/'

   ierr = nf90_close(ncid)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); endif

 END SUBROUTINE close_nc

END MODULE ncio_utils
