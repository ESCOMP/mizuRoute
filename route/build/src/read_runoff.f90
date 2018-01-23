module read_runoff
use nrtype
USE netcdf
use public_var
use read_netcdf, only:get_nc
use data_runoff, only:runoff_data

implicit none

private
public::get_runoff_meta
public::get_runoff_hru
public::get_runoff

contains

  ! *********************************************************************
  ! subroutine: get metadata runoff file
  ! *********************************************************************
  subroutine get_runoff_meta(fname,           &  ! input: filename
                             ierr, message,   &  ! output: error control
                             n_hru,           &  ! output(optional): number of hru in runoff data
                             n_time,          &  ! output(optional): number of time in runoff data
                             t_unit)             ! output(optional): time units
    implicit none
    ! input variables
    character(*), intent(in)              :: fname        ! filename
    ! output variables
    integer(i4b), optional, intent(out)   :: n_hru        ! number of runoff hrus
    integer(i4b), optional, intent(out)   :: n_time       ! number of time steps
    character(*), optional, intent(out)   :: t_unit       ! time units
    integer(i4b), intent(out)             :: ierr         ! error code
    character(*), intent(out)             :: message      ! error message
    ! local variables
    character(len=strLen)                 :: cmessage     ! error message from subroutine

    ! initialize error control
    ierr=0; message='get_runoff_meta/'

    if (present(n_hru)) then
      call get_qDims(fname, dname_hruid, n_hru, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end if
    if (present(n_time)) then
      call get_qDims(fname, dname_time, n_time, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end if
    if (present(t_unit)) then
      call get_time_unit(fname, dname_time, t_unit, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

    return
  end subroutine

  ! *********************************************************************
  ! subroutine: populate runoff hru
  ! *********************************************************************
  subroutine get_runoff_hru(fname,           &  ! input: filename
                            ierr, message)      ! output: error control
    implicit none
    ! input variables
    character(*), intent(in)              :: fname        ! filename
    ! output variables
    integer(i4b), intent(out)             :: ierr         ! error code
    character(*), intent(out)             :: message      ! error message
    ! local variables
    character(len=strLen)                 :: cmessage     ! error message from subroutine
    integer(i4b)                          :: nHRU        ! number of runoff hrus

    ! initialize error control
    ierr=0; message='get_runoff_hru'

    call get_qDims(fname, dname_hruid, nHRU, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    allocate(runoff_data%hru_id(nHRU), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data%hru_id'; return; endif

    call get_nc(trim(fname),vname_hruid, runoff_data%hru_id, 1, nHRU, ierr, cmessage);

    return
  end subroutine

 ! *********************************************************************
 ! subroutine: get unit from runoff file
 ! *********************************************************************
 subroutine get_time_unit(fname,           &  ! input: filename
                          vname_time,      &  ! input: name of coordinate dimension time
                          units_time,      &  ! output: time units
                          ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname_time   ! coordinate variable for time
 ! output variables
 character(*), intent(out)       :: units_time   ! time units
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 ! initialize error control
 ierr=0; message='get_time_unit/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the ID of the time variable
 ierr = nf90_inq_varid(ncid, trim(vname_time), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the time units
 ierr = nf90_get_att(ncid, ivarID, 'units', units_time)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: read dimensions from runoff file
 ! *********************************************************************
 subroutine get_qDims(fname,           &  ! input: filename
                      dname,           &  ! input: dimension name
                      nDim,            &  ! output: dimension size
                      ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: dname        ! dimension name
 ! output variables
 integer(i4b), intent(out)       :: nDim         ! diminesion size
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: idimID       ! dimension ID

 ! initialize error control
 ierr=0; message='get_qDims/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the ID of the dimension
 ierr = nf90_inq_dimid(ncid, dname, idimID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname); return; endif

 ! get the length of the dimension
 ierr = nf90_inquire_dimension(ncid, idimID, len=nDim)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine


 ! *********************************************************************
 ! new subroutine: read runoff data
 ! *********************************************************************
 subroutine get_runoff(fname,          &  ! input: filename
                       iTime,          &  ! input: time index
                       dTime,          &  ! output: time
                       ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 integer(i4b), intent(in)        :: iTime        ! index of time element
 ! output variables
 real(dp), intent(out)           :: dTime        ! time
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 real(dp),allocatable            :: dummy(:,:)   ! dummy 2d array
 integer(i4b)                    :: nHRU         ! number of hrus
 character(len=strLen)           :: cmessage     ! error message from subroutine

 ! initialize error control
 ierr=0; message='get_runoff/'

 call get_nc(trim(fname),vname_time, dTime, iTime, ierr, cmessage);
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call get_qDims(fname, dname_hruid, nHRU, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 allocate(dummy(nHRU,1),stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating dummy'; return; endif

 if (.not.(allocated(runoff_data%qsim))) then
   allocate(runoff_data%qsim(nHRU),stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data%qsim'; return; endif
 endif

 call get_nc(trim(fname),vname_qsim, dummy, (/1,iTime/), (/nHRU,1/), ierr, cmessage);
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 runoff_data%qsim = reshape(dummy,(/nHRU/))

 end subroutine

end module
