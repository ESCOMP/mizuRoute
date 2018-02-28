module read_remap

use nrtype
use netcdf
use public_var
use read_netcdf, only:get_nc
use read_netcdf, only:get_nc_dim_len
use read_netcdf, only:get_units

implicit none

private

public::get_remap_data
public::get_runoff_metadata

contains

  ! *****
  ! public subroutine: get runoff metadata...
  ! *****************************************
  subroutine get_runoff_metadata(&
                                 ! input
                                 fname       , & ! filename
                                 ! output
                                 nSpatial    , & ! number of spatial elements
                                 hruId       , & ! vector of HRU Ids in the output file
                                 nTime       , & ! number of time steps
                                 timeUnits   , & ! time units
                                 ! error control
                                 ierr, message)  ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)                :: fname           ! filename
  ! output variables
  integer(i4b), intent(out)               :: nSpatial        ! number of spatial elements
  integer(i4b), intent(out) , allocatable :: hruId(:)        ! vector of HRU Ids in the output file
  integer(i4b), intent(out)               :: nTime           ! number of time steps
  character(*), intent(out)               :: timeUnits       ! time units
  ! error control
  integer(i4b), intent(out)               :: ierr            ! error code
  character(*), intent(out)               :: message         ! error message
  ! local variables
  character(len=strLen)                   :: cmessage        ! error message from subroutine
  ! initialize error control
  ierr=0; message='get_runoff_metadata/'

  ! get the number of HRUs
  call get_nc_dim_len(fname, trim(dname_hruid), nSpatial, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get number of time steps from the runoff file
  call get_nc_dim_len(fname, trim(dname_time), nTime, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time units
  call get_units(fname, trim(vname_time), timeUnits, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate space
  allocate(hruId(nSpatial), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating hruId'; return; endif

  ! get HRU ids from the runoff file
  call get_nc(fname, vname_hruid, hruId, 1, nSpatial, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine get_runoff_metadata

  ! *****
  ! public subroutine: get mapping data between runoff hru and river network hru...
  ! ********************************************************************************
  subroutine get_remap_data(fname,         &   ! input: file name
                            remap_data,    &   ! output: data structure to remap data from a polygon
                            ierr, message)     ! output: error control
  USE dataTypes,  only : remap                 ! remapping data type
  implicit none
  ! input variables
  character(*), intent(in)           :: fname           ! filename
  ! output variables
  type(remap),  intent(out)          :: remap_data      ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
  integer(i4b), intent(out)          :: ierr            ! error code
  character(*), intent(out)          :: message         ! error message
  ! local variables
  integer(i4b)                       :: iVar            ! index of variables
  integer(i4b)                       :: nHRU            ! number of HRU in mapping files (this should match up with river network hru)
  integer(i4b)                       :: nData           ! number of data (weight, runoff hru id) in mapping files
  character(len=strLen)              :: cmessage        ! error message from subroutine

  ! initialize error control
  ierr=0; message='get_remap_data/'

  ! get the number of HRUs
  call get_nc_dim_len(fname, dname_hru_remap, nHRU, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the number of spatial elements in the runoff file
  call get_nc_dim_len(fname, dname_data_remap, nData, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate space for structure components
  allocate(remap_data%hru_id(nHRU), remap_data%num_qhru(nHRU), remap_data%qhru_id(nData), remap_data%weight(nData), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space'; return; endif

  ! get data
  do iVar=1,4
   select case(iVar)
    case(1); call get_nc(fname, vname_hruid_in_remap, remap_data%hru_id,   1, nHRU,  ierr, cmessage)
    case(2); call get_nc(fname, vname_num_qhru,       remap_data%num_qhru, 1, nHRU,  ierr, cmessage)
    case(3); call get_nc(fname, vname_weight,         remap_data%weight,   1, nData, ierr, cmessage)
    case(4); call get_nc(fname, vname_qhruid,         remap_data%qhru_id,  1, nData, ierr, cmessage)
    case default; ierr=20; message=trim(message)//'unable to find variable'; return
   end select
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end do  ! looping through variables

  end subroutine

end module
