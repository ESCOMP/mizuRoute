module read_remap

use nrtype
use netcdf
use public_var
use read_netcdf, only:get_nc
use data_remap, only:remap_data

implicit none

private

public::get_remap_data

contains

  ! *********************************************************************
  ! public subroutine: get mapping data between runoff hru and river network hru
  ! *********************************************************************
  subroutine get_remap_data(fname,         &   ! input: file name
                            ierr, message)     ! output: error control
    implicit none
    ! input variables
    character(*), intent(in)           :: fname           ! filename
    ! input-output
    ! output variables
    integer(i4b), intent(out)          :: ierr            ! error code
    character(*), intent(out)          :: message         ! error message
    ! local variables
    integer(i4b)                       :: nHRU            ! number of HRU in mapping files (this should match up with river network hru)
    integer(i4b)                       :: nData           ! number of data (weight, runoff hru id) in mapping files
    character(len=strLen)              :: cmessage        ! error message from subroutine

    ! initialize error control
    ierr=0; message='get_remap_data/'

    call get_map_dims(fname, dname_hru_remap, dname_data_remap, nData, nHRU, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    allocate(remap_data%hru_id(nHRU), stat=ierr)
    allocate(remap_data%num_qhru(nHRU), stat=ierr)
    allocate(remap_data%qhru_id(nData), stat=ierr)
    allocate(remap_data%weight(nData), stat=ierr)

    call get_nc(trim(fname),vname_hruid_in_remap, remap_data%hru_id, 1, nHRU, ierr, cmessage);
    call get_nc(trim(fname),vname_weight, remap_data%weight, 1, nData, ierr, cmessage);
    call get_nc(trim(fname),vname_qhruid, remap_data%qhru_id, 1, nData, ierr, cmessage);
    call get_nc(trim(fname),vname_num_qhru, remap_data%num_qhru, 1, nHRU, ierr, cmessage);

    return
  end subroutine

  ! *********************************************************************
  ! subroutine: read dimensions from runoff file
  ! *********************************************************************
  subroutine get_map_dims(fname,           &  ! input: filename
                          dname_hru,       &  ! input: name of coordinate dimension HRUid
                          dname_data,      &  ! input: name of coordinate dimension time
                          nData,           &  ! output: number of time elements
                          nHRU,            &  ! output: number of HRUs in the runoff data file
                          ierr, message)      ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)        :: fname        ! filename
  character(*), intent(in)        :: dname_hru    ! coordinate variable for HRUs
  character(*), intent(in)        :: dname_data   ! coordinate variable for time
  ! output variables
  integer(i4b), intent(out)       :: nData        ! number of time elements
  integer(i4b), intent(out)       :: nHRU         ! number of HRUs
  integer(i4b), intent(out)       :: ierr         ! error code
  character(*), intent(out)       :: message      ! error message
  ! local variables
  integer(i4b)                    :: ncid         ! NetCDF file ID
  integer(i4b)                    :: dimID_data  ! dimension ID for time
  integer(i4b)                    :: dimID_hru   ! dimension ID for hru
  ! initialize error control
  ierr=0; message='get_mapDims/'

  ! open file for reading
  ierr = nf90_open(fname, nf90_nowrite, ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

  ! get the ID of the data dimension
  ierr = nf90_inq_dimid(ncid, dname_data, dimID_Data)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(vname_time); return; endif

  ! get the length of the data dimension
  ierr = nf90_inquire_dimension(ncid, dimID_data, len=nData)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get the ID of the HRU dimension
  ierr = nf90_inq_dimid(ncid, dname_hru, dimID_hru)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(vname_hruID); return; endif

  ! get the length of the HRU dimension
  ierr = nf90_inquire_dimension(ncid, dimID_hru, len=nHRU)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! close the NetCDF file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine

end module
