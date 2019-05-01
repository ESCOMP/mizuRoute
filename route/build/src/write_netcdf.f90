module write_netcdf

USE nrtype
USE netcdf

implicit none

private

public::write_nc

interface write_nc
  module procedure write_ivec
  module procedure write_dvec
  module procedure write_2d_iarray
  module procedure write_2d_darray
  module procedure write_3d_iarray
  module procedure write_3d_darray
end interface

contains

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

end module write_netcdf
