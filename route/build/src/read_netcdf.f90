module read_netcdf

USE nrtype
USE netcdf
use public_var

implicit none

private

public::get_nc
public::get_nc_dim_len

interface get_nc
  module procedure get_iscalar
  module procedure get_dscalar
  module procedure get_ivec
  module procedure get_dvec
  module procedure get_2d_iarray
  module procedure get_2d_darray
  module procedure get_3d_iarray
  module procedure get_3d_darray
end interface

contains

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
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get the ID of the time dimension
  ierr = nf90_inq_dimid(ncid, dname, iDimID)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname); return; endif

  ! get the length of the time dimension
  ierr = nf90_inquire_dimension(ncid, iDimID, len=nDim)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close output file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

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
 ! subroutine: read a double precision 2D array
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
 ! subroutine: read a double precision 3D array
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

end module read_netcdf
