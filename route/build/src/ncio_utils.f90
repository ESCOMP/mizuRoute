module io_netcdf

USE nrtype
USE netcdf

implicit none

private

public::get_nc
public::get_var_dims
public::get_nc_dim_len
public::get_var_attr_char
public::get_var_attr_real
public::defvar
public::write_nc

interface get_nc
  module procedure get_iscalar
  module procedure get_dscalar
  module procedure get_ivec
  module procedure get_dvec
  module procedure get_2d_iarray
  module procedure get_2d_darray
  module procedure get_3d_iarray
  module procedure get_3d_darray
  module procedure get_4d_iarray
  module procedure get_4d_darray
end interface

interface write_nc
  module procedure write_ivec
  module procedure write_dvec
  module procedure write_charvec
  module procedure write_2d_iarray
  module procedure write_2d_darray
  module procedure write_3d_iarray
  module procedure write_3d_darray
end interface

contains

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
 ! subroutine: get attribute values for a variable
 ! *********************************************************************
 subroutine get_var_attr_char(fname,           &  ! input: filename
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
 character(*), intent(out)       :: attr_value   ! attribute value
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: var_type     ! attribute variable type
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 ! initialize error control
 ierr=0; message='get_var_attr_char/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the ID of the variable
 ierr = nf90_inq_varid(ncid, trim(vname), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! Inquire attribute type, NF90_CHAR(=2)
 ierr = nf90_inquire_attribute(ncid, ivarID, attr_name, xtype=var_type)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 if (var_type /= 2)then; ierr=20; message=trim(message)//'attribute type is not character'; return; endif

 ! get the attribute value
 ierr = nf90_get_att(ncid, ivarID, attr_name, attr_value)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_var_attr_char

 ! *********************************************************************
 ! subroutine: get attribute values for a real variable
 ! *********************************************************************
 subroutine get_var_attr_real(fname,           &  ! input: filename
                              vname,           &  ! input: variable name
                              attr_name,       &  ! inpu: attribute name
                              attr_value,      &  ! output: attribute value in real
                              ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname        ! variable name
 character(*), intent(in)        :: attr_name    ! attribute name
 ! output variables
 real(dp), intent(out)           :: attr_value   ! attribute value in real
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: var_type     ! attribute variable type
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 ! initialize error control
 ierr=0; message='get_var_attr_real/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the ID of the variable
 ierr = nf90_inq_varid(ncid, trim(vname), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! Inquire attribute type, NF90_CHAR(=2)
 ierr = nf90_inquire_attribute(ncid, ivarID, attr_name, xtype=var_type)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 if (var_type /= 5 .and. var_type /= 6)then; ierr=20; message=trim(message)//'attribute type is not real'; return; endif

 ! get the attribute value
 ierr = nf90_get_att(ncid, ivarID, attr_name, attr_value)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_var_attr_real


 ! *********************************************************************
 ! subroutine: get integer scalar value from netCDF
 ! *********************************************************************
 subroutine get_iscalar(fname,           &  ! input:  filename
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
  integer(i4b), intent(out)       :: array     ! output variable data
  integer(i4b), intent(out)       :: ierr      ! error code
  character(*), intent(out)       :: message   ! error message
  ! local variables
  character(len=strLen)           :: cmessage     ! error message of downwind routine
  integer(i4b)                    :: array_vec(1) ! output variable data

 ! initialize error control
 ierr=0; message='get_iscalar/'

 call get_ivec(fname, vname, array_vec, iStart, 1, ierr, cmessage)
 array = array_vec(1)

 end subroutine

 ! *********************************************************************
 ! subroutine: double precision scalar value from netCDF
 ! *********************************************************************
 subroutine get_dscalar(fname,           &  ! input:  filename
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
  real(dp), intent(out)           :: array     ! output variable data
  integer(i4b), intent(out)       :: ierr      ! error code
  character(*), intent(out)       :: message   ! error message
  ! local variables
  character(len=strLen)           :: cmessage     ! error message of downwind routine
  real(dp)                        :: array_vec(1) ! output variable data

 ! initialize error control
 ierr=0; message='get_dscalar/'

 call get_dvec(fname, vname, array_vec, iStart, 1, ierr, cmessage)
 array = array_vec(1)

 end subroutine

 ! *********************************************************************
 ! subroutine: get integer vector value from netCDF
 ! *********************************************************************
 subroutine get_ivec(fname,           &  ! input:  filename
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
  integer(i4b), intent(out)       :: array(:)  ! output variable data
  integer(i4b), intent(out)       :: ierr      ! error code
  character(*), intent(out)       :: message   ! error message
  ! local variables
  integer(i4b)                    :: ncid      ! NetCDF file ID
  integer(i4b)                    :: iVarID    ! NetCDF variable ID

 ! initialize error control
 ierr=0; message='get_ivec/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the data
 ierr = nf90_get_var(ncid, iVarID, array, start=(/iStart/), count=(/iCount/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: read a double precision vector
 ! *********************************************************************
 subroutine get_dvec(fname,           &  ! input: filename
                     vname,           &  ! input: variable name
                     array,           &  ! output: variable data
                     iStart,          &  ! input: start index
                     iCount,          &  ! input: length of vector
                     ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname     ! filename
  character(*), intent(in)        :: vname     ! variable name
  integer(i4b), intent(in)        :: iStart    ! start index
  integer(i4b), intent(in)        :: iCount    ! length of vector to be read in
  ! output variables
  real(dp), intent(out)           :: array(:)  ! output variable data
  integer(i4b), intent(out)       :: ierr      ! error code
  character(*), intent(out)       :: message   ! error message
  ! local variables
  integer(i4b)                    :: ncid      ! NetCDF file ID
  integer(i4b)                    :: iVarID    ! NetCDF variable ID

  ! initialize error control
  ierr=0; message='get_dvec/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get the data
  ierr = nf90_get_var(ncid, iVarID, array, start=(/iStart/), count=(/iCount/))
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: read a integer 2D array
 ! *********************************************************************
 subroutine get_2d_iarray(fname,           &  ! input: filename
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
  integer(i4b), intent(out)       :: array(:,:)   ! variable data
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID

  ! initialize error control
  ierr=0; message='get_2d_iarray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: read a integer 3D array
 ! *********************************************************************
 subroutine get_3d_iarray(fname,          &  ! input: filename
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
  integer(i4b), intent(out)       :: array(:,:,:) ! variable data
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='get_3d_iarray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: read a integer 4D array
 ! *********************************************************************
 subroutine get_4d_iarray(fname,          &  ! input: filename
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
  integer(i4b), intent(out)       :: array(:,:,:,:) ! variable data
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='get_4d_iarray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine


 ! *********************************************************************
 ! subroutine: read a double precision 2D array
 ! *********************************************************************
 subroutine get_2d_darray(fname,          &  ! input: filename
                          vname,          &  ! input: variable name
                          array,          &  ! output: variable data
                          iStart,         &  ! input: start index
                          iCount,         &  ! input: length of vector
                          ierr, message)     ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  integer(i4b), intent(in)        :: iStart(1:2)  ! start indices
  integer(i4b), intent(in)        :: iCount(1:2)  ! length of vector
  ! output variables
  real(dp), intent(out)           :: array(:,:)   ! variable data
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='get_2d_darray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: read a double precision 3D array
 ! *********************************************************************
 subroutine get_3d_darray(fname,          &  ! input: filename
                          vname,          &  ! input: variable name
                          array,          &  ! output: variable data
                          iStart,         &  ! input: start index
                          iCount,         &  ! input: length of vector
                          ierr, message)     ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname         ! filename
  character(*), intent(in)        :: vname         ! variable name
  integer(i4b), intent(in)        :: iStart(1:3)   ! start indices
  integer(i4b), intent(in)        :: iCount(1:3)   ! length of vector
  ! output variables
  real(dp), intent(out)           :: array(:,:,:)  ! variable data
  integer(i4b), intent(out)       :: ierr          ! error code
  character(*), intent(out)       :: message       ! error message
  ! local variables
  integer(i4b)                    :: ncid          ! NetCDF file ID
  integer(i4b)                    :: iVarId        ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='get_3d_darray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: read a double precision 4D array
 ! *********************************************************************
 subroutine get_4d_darray(fname,          &  ! input: filename
                          vname,          &  ! input: variable name
                          array,          &  ! output: variable data
                          iStart,         &  ! input: start index
                          iCount,         &  ! input: length of vector
                          ierr, message)     ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname         ! filename
  character(*), intent(in)        :: vname         ! variable name
  integer(i4b), intent(in)        :: iStart(1:4)   ! start indices
  integer(i4b), intent(in)        :: iCount(1:4)   ! length of vector
  ! output variables
  real(dp), intent(out)           :: array(:,:,:,:)  ! variable data
  integer(i4b), intent(out)       :: ierr          ! error code
  character(*), intent(out)       :: message       ! error message
  ! local variables
  integer(i4b)                    :: ncid          ! NetCDF file ID
  integer(i4b)                    :: iVarId        ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='get_4d_darray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! read data
  ierr = nf90_get_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

  ! *********************************************************************
  ! subroutine: write an integer vector
  ! *********************************************************************
  subroutine write_ivec(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        array,           &  ! input: variable data
                        iStart,          &  ! input: start index
                        iCount,          &  ! input: length of vector
                        ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  integer(i4b), intent(in)        :: array(:)     ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start index
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='write_ivec/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

  ! *********************************************************************
  ! subroutine: write a double precision vector
  ! *********************************************************************
  subroutine write_dvec(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        array,           &  ! input: variable data
                        iStart,          &  ! input: start index
                        iCount,          &  ! input: length of vector
                        ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  real(dp), intent(in)            :: array(:)     ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start indices
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='write_dvec/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

  ! *********************************************************************
  ! subroutine: write a character vector
  ! *********************************************************************
  subroutine write_charvec(fname,           &  ! input: filename
                           vname,           &  ! input: variable name
                           array,           &  ! input: variable data
                           iStart,          &  ! input: start index
                           iCount,          &  ! input: length of vector
                           ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  character(*), intent(in)        :: array(:)     ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start indices
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='write_charvec/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

  ! *********************************************************************
  ! subroutine: write a double precision 2D array
  ! *********************************************************************
  subroutine write_2d_iarray(fname,           &  ! input: filename
                             vname,           &  ! input: variable name
                             array,           &  ! input: variable data
                             iStart,          &  ! input: start index
                             iCount,          &  ! input: length of vector
                             ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  integer(i4b), intent(in)        :: array(:,:)   ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start indices
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='write_2d_iarray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

  ! *********************************************************************
  ! subroutine: write a double precision 3D array
  ! *********************************************************************
  subroutine write_3d_iarray(fname,          &  ! input: filename
                             vname,          &  ! input: variable name
                             array,          &  ! input: variable data
                             iStart,         &  ! input: start index
                             iCount,         &  ! input: length of vector
                             ierr, message)     ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  integer(i4b), intent(in)        :: array(:,:,:) ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start indices
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='write_3d_iarray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

  ! *********************************************************************
  ! subroutine: write a double precision 2D array
  ! *********************************************************************
  subroutine write_2d_darray(fname,          &  ! input: filename
                             vname,          &  ! input: variable name
                             array,          &  ! input: variable data
                             iStart,         &  ! input: start index
                             iCount,         &  ! input: length of vector
                             ierr, message)     ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: vname        ! variable name
  real(dp), intent(in)            :: array(:,:)   ! variable data
  integer(i4b), intent(in)        :: iStart(:)    ! start indices
  integer(i4b), intent(in)        :: iCount(:)    ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: iVarId       ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='write_2d_darray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

  ! *********************************************************************
  ! subroutine: write a double precision 3D array
  ! *********************************************************************
  subroutine write_3d_darray(fname,          &  ! input: filename
                             vname,          &  ! input: variable name
                             array,          &  ! input: variable data
                             iStart,         &  ! input: start index
                             iCount,         &  ! input: length of vector
                             ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname         ! filename
  character(*), intent(in)        :: vname         ! variable name
  real(dp), intent(in)            :: array(:,:,:)  ! variable data
  integer(i4b), intent(in)        :: iStart(:)     ! start indices
  integer(i4b), intent(in)        :: iCount(:)     ! length of vector
  ! output variables
  integer(i4b), intent(out)       :: ierr          ! error code
  character(*), intent(out)       :: message       ! error message
  ! local variables
  integer(i4b)                    :: ncid          ! NetCDF file ID
  integer(i4b)                    :: iVarId        ! NetCDF variable ID
  ! initialize error control
  ierr=0; message='write_3d_darray/'

  ! open NetCDF file
  ierr = nf90_open(trim(fname),nf90_write,ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get variable ID
  ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! write data
  ierr = nf90_put_var(ncid,iVarId,array,start=iStart,count=iCount)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

 ! *********************************************************************
 ! private subroutine: define variable attributes NetCDF file
 ! *********************************************************************
 subroutine defvar(ncid, vname, dimNames, ivtype, ierr, message, vdesc, vunit, vcal)
  ! input
  integer(i4b), intent(in)             :: ncid                   ! Input: netcdf fine ID
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

 end subroutine defvar

end module io_netcdf
