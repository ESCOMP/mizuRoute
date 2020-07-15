MODULE example_netcdf_module

! test program for Lieke
USE netcdf                                  ! use netCDF libraries

! parameter
integer,parameter  :: str_len=64

! dimension lengths
integer  :: nx = 3  ! length of the x dimension
integer  :: ny = 5  ! length of the y dimension
integer  :: nt = 10 ! length of time dimension

character(len=str_len) :: cmode_name = 'netcdf4'

CONTAINS

  SUBROUTINE example_netcdf()
    implicit none
    ! error code
    integer  :: ierr

    ! NetCDF file id
    integer  :: ncid

    ! NetCDF dimension ids
    integer  :: ix_dim,iy_dim,it_dim

    ! NetCDF variable IDs
    integer  :: varid

    ! NetCDF creation mode name and flag
    integer  :: cmode

    ! loop counter
    integer  :: ix

    ! some array
    real     :: stuff(nx,ny)

    ! character strings
    character(len=str_len) :: filename  ! filename
    character(len=str_len) :: var_name  ! variable name
    character(len=str_len) :: var_desc  ! variable description

    ! -------------- end of declaration ----------------------------------

    ! define some stuff
    filename   = 'stuff.nc'
    var_name   = 'stuff'
    var_desc   = 'just some random crap'

    select case(trim(cmode_name))
      case('64bit_offset'); cmode = nf90_64bit_offset
      case('netcdf4');      cmode = nf90_netcdf4
      case default;         cmode = nf90_clobber
    end select

    ! create file
    ierr = nf90_create(trim(filename), cmode, ncid)
     call handle_err(ierr,trim(nf90_strerror(ierr)))

    ! define dimensions
    ierr = nf90_def_dim(ncid,'x',nx,ix_dim)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ierr = nf90_def_dim(ncid,'y',ny,iy_dim)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ierr = nf90_def_dim(ncid,'t',nf90_unlimited,it_dim)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ! define a variable
    ierr = nf90_def_var(ncid,'time',nf90_double,it_dim, varid)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ierr = nf90_put_att(ncid,varid,'long_name','time')
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ierr = nf90_def_var(ncid,'stuff',nf90_double,(/ix_dim,iy_dim,it_dim/),varid)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ierr = nf90_put_att(ncid,varid,'long_name',trim(var_desc))
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ! end definitions and close file
    ierr = nf90_enddef(ncid)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ierr = nf90_close(ncid)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ! -------------- end of netcdf definition  ----------------------------------

    ! write stuff
    ierr = nf90_open(trim(filename),nf90_write, ncid)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    ! inquire variable
    ierr = nf90_inq_varid(ncid,'stuff',varid)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

    do ix = 1,nt
      stuff = ix
      ierr = nf90_put_var(ncid, varid, stuff, start=[1,1,ix], count=[nx,ny,1])
      call handle_err(ierr,trim(nf90_strerror(ierr)))
    end do

    ! close
    ierr = nf90_close(ncid)
    call handle_err(ierr,trim(nf90_strerror(ierr)))

  END SUBROUTINE example_netcdf


  SUBROUTINE handle_err(err,message)
    implicit none
    integer,     intent(in)::err             ! error code
    character(*),intent(in)::message         ! error message
    if(err/=0)then
      print*,'FATAL ERROR: '//trim(message)
      call flush(6)
      stop
    endif
  END SUBROUTINE handle_err

END MODULE example_netcdf_module


PROGRAM test_netcdf
  USE example_netcdf_module, only:example_netcdf
  ! call subroutine
  call example_netcdf()
  stop
END
