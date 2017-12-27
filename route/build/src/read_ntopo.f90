module read_ntopo

USE nrtype
USE netcdf
implicit none
private
public::get_vec_dim
public::get_vec_ivar
public::get_scl_ivar
public::get_vec_dvar
public::get_scl_dvar

contains

! *********************************************************************
! subroutine: get vector dimension from netCDF 
! *********************************************************************
subroutine get_vec_dim(fname,           &  ! input: filename
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
 ierr=0; message='get_vec_dim/'

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

end subroutine get_vec_dim

! *********************************************************************
! new subroutine: get vector value from netCDF
! *********************************************************************
subroutine get_vec_ivar(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        iVec,            &  ! output: outputvariable data
                        iStart,          &  ! input: start index
                        iCount,          &  ! input: length of vector
                        ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                           :: fname     ! filename
 character(*), intent(in)                           :: vname     ! variable name
 integer(i4b), intent(in)                           :: iStart(:) ! start index
 integer(i4b), intent(in)                           :: iCount(:) ! length of vector to be read in
 ! output variables
 !integer(i4b), intent(out),dimension(:),allocatable :: iVec      ! output variable data
 integer(i4b), intent(out),dimension(:)             :: iVec      ! output variable data
 integer(i4b), intent(out)                          :: ierr      ! error code
 character(*), intent(out)                          :: message   ! error message
 ! local variables
 integer(i4b)                                       :: ncid      ! NetCDF file ID
 integer(i4b)                                       :: iVarId    ! NetCDF variable ID
 ! initialize error control
 ierr=0; message='get_vec_ivar/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! allocate space for the output 
 !allocate(iVec(iCount),stat=ierr)
 !if(ierr/=0)then; message=trim(message)//'problem allocating space for iVec'; return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the data
 ierr = nf90_get_var(ncid, ivarID, iVec, start=iStart, count=iCount)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

end subroutine get_vec_ivar


! *********************************************************************
! new subroutine: get vector value from netCDF
! *********************************************************************
subroutine get_scl_ivar(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        iScl,            &  ! output: outputvariable data
                        iStart,          &  ! input: start index
                        ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                           :: fname     ! filename
 character(*), intent(in)                           :: vname     ! variable name
 integer(i4b), intent(in)                           :: iStart    ! start index
 ! output variables
 integer(i4b), intent(out)                          :: iScl      ! output variable data
 integer(i4b), intent(out)                          :: ierr      ! error code
 character(*), intent(out)                          :: message   ! error message
 ! local variables
 integer(i4b),dimension(1)                          :: iDummy    ! temporary vector of length 1
 integer(i4b)                                       :: ncid      ! NetCDF file ID
 integer(i4b)                                       :: iVarId    ! NetCDF variable ID
 ! initialize error control
 ierr=0; message='get_scl_ivar/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the data
 ierr = nf90_get_var(ncid, ivarID, iDummy, start=(/iStart/), count=(/1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! save the output
 iScl = iDummy(1)

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

end subroutine get_scl_ivar

! *********************************************************************
! new subroutine: get vector value from netCDF
! *********************************************************************
subroutine get_scl_dvar(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        dScl,            &  ! output: outputvariable data
                        iStart,          &  ! input: start index
                        ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)               :: fname     ! filename
 character(*), intent(in)               :: vname     ! variable name
 integer(i4b), intent(in)               :: iStart    ! start index
 ! output variables
 real(dp), intent(out)                  :: dScl      ! output variable data
 integer(i4b), intent(out)              :: ierr      ! error code
 character(*), intent(out)              :: message   ! error message
 ! local variables
 real(dp),dimension(1)                  :: dDummy    ! temporary vector of length 1
 integer(i4b)                           :: ncid      ! NetCDF file ID
 integer(i4b)                           :: iVarId    ! NetCDF variable ID
 ! initialize error control
 ierr=0; message='get_scl_dvar/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the data
 ierr = nf90_get_var(ncid, ivarID, dDummy, start=(/iStart/), count=(/1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! save the output
 dScl = dDummy(1)

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

end subroutine get_scl_dvar


 ! *********************************************************************
 ! new subroutine: read a double precision vector
 ! *********************************************************************
 subroutine get_vec_dvar(fname,           &  ! input: filename
                         vname,           &  ! input: variable name
                         dVec,            &  ! input: variable data
                         iStart,          &  ! input: start index
                         iCount,          &  ! input: length of vector
                         ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                           :: fname     ! filename
 character(*), intent(in)                           :: vname     ! variable name
 integer(i4b), intent(in)                           :: iStart(:) ! start index
 integer(i4b), intent(in)                           :: iCount(:) ! length of vector to be read in
 ! output variables
 !real(dp), intent(out), dimension(:), allocatable   :: dVec      ! output variable data
 real(dp), intent(out), dimension(:)                :: dVec      ! output variable data
 integer(i4b), intent(out)                          :: ierr      ! error code
 character(*), intent(out)                          :: message   ! error message
 ! local variables
 integer(i4b)                                       :: ncid      ! NetCDF file ID
 integer(i4b)                                       :: iVarId    ! NetCDF variable ID

 ! initialize error control
 ierr=0; message='get_vec_dVec/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! allocate space for the output 
! allocate(dVec(iCount),stat=ierr)
! if(ierr/=0)then; message=trim(message)//'problem allocating space for dVec'; return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the data
 ierr = nf90_get_var(ncid, ivarID, dVec, start=iStart, count=iCount)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_vec_dvar

end module read_ntopo
