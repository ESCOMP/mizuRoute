module example_netcdf_module
! test program for Lieke
USE netcdf                                  ! use netCDF libraries
contains


subroutine example_netcdf()

! error code
integer  :: ierr

! NetCDF file id
integer  :: ncid

! NetCDF dimension ids
integer  :: ix_dim,iy_dim,it_dim

! NetCDF variable IDs
integer  :: ivar_id

! dimension lengths
integer  :: nx,ny

! character strings
integer,parameter  :: str_len=64
character(len=str_len) :: filename  ! filename
character(len=str_len) :: var_name  ! variable name
character(len=str_len) :: var_desc  ! variable description

! -------------- end of definitions ----------------------------------

! define some stuff

nx = 2  ! length of the x dimension
ny = 3  ! length of the y dimension

filename = 'stuff.nc'
var_name = 'stuff'
var_desc = 'just some random crap'

! create a file

! create file
ierr = nf90_create(trim(filename),nf90_clobber,ncid); call handle_err(ierr)

 ! define dimensions
 ierr = nf90_def_dim(ncid,'x',nx,ix_dim); call handle_err(ierr)
 ierr = nf90_def_dim(ncid,'y',ny,iy_dim); call handle_err(ierr)
 ierr = nf90_def_dim(ncid,'t',nf90_unlimited,it_dim); call handle_err(ierr)

 ! define a variable
 ierr = nf90_def_var(ncid,'stuff',nf90_real,(/ix_dim,iy_dim,it_dim/),ivar_id); call handle_err(ierr)
 ierr = nf90_put_att(ncid,ivar_id,'long_name',trim(var_desc)); call handle_err(ierr)

! end definitions and close file
ierr = nf90_enddef(ncid); call handle_err(ierr)
ierr = nf90_close(ncid); call handle_err(ierr)


end subroutine example_netcdf

subroutine handle_err(err)
! used to handle errors for NetCDF calls
implicit none
! declare dummies
integer, intent(in)        :: err
! start procedure here
if (err/=nf90_noerr) then
 print*, '['//trim(nf90_strerror(err))//']'
 stop '[fatal error in netcdf]'
endif
end subroutine handle_err



end module example_netcdf_module
