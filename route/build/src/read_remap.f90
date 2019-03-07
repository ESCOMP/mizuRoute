module read_remap

! data types
use nrtype

! global data
USE public_var

! Netcdf
use netcdf
use read_netcdf, only:get_nc
use read_netcdf, only:get_nc_dim_len
use read_netcdf, only:get_var_attr_char

implicit none

! privacy
private
public::get_remap_data

contains

 ! *****
 ! public subroutine: get mapping data between runoff hru and river network hru...
 ! ********************************************************************************
 subroutine get_remap_data(fname,         &   ! input: file name
                           nSpatial,      &   ! input: number of spatial elements
                           _remap_data,   &   ! output: data structure to remap data from a polygon
                           ierr, message)     ! output: error control
 USE dataTypes,  only : remap                 ! remapping data type
 implicit none
 ! input variables
 character(*), intent(in)           :: fname           ! filename
 integer(i4b), intent(in)           :: nSpatial(1:2)   ! number of spatial elements
 ! output variables
 type(remap),  intent(out)          :: _remap_data     ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
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

 ! allocate space for info in mapping file
 allocate(_remap_data%hru_id(nHRU), _remap_data%num_qhru(nHRU), &
          _remap_data%hru_ix(nHRU), _remap_data%weight(nData), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for common mapping info '; return; endif

 ! if runoff input is hru vector...
 if (nSpatial(2) == integerMissing) then
  allocate(_remap_data%qhru_id(nData), _remap_data%qhru_ix(nData),  stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for mapping info for vector runoff'; return; endif
 ! if runoff input is grid...
 else
  allocate(_remap_data%i_index(nData),_remap_data%j_index(nData), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for mapping info'; return; endif
 endif

 ! get data
 do iVar=1,6
  if (nSpatial(2) == integerMissing .and. iVar >= 5) cycle
  if (nSpatial(2) /= integerMissing .and. iVar == 4) cycle
  select case(iVar)
   case(1); call get_nc(fname, vname_hruid_in_remap, _remap_data%hru_id,   1, nHRU,  ierr, cmessage)
   case(2); call get_nc(fname, vname_num_qhru,       _remap_data%num_qhru, 1, nHRU,  ierr, cmessage)
   case(3); call get_nc(fname, vname_weight,         _remap_data%weight,   1, nData, ierr, cmessage)
   case(4); call get_nc(fname, vname_qhruid,         _remap_data%qhru_id,  1, nData, ierr, cmessage)
   case(5); call get_nc(fname, vname_i_index,        _remap_data%i_index,  1, nData, ierr, cmessage)
   case(6); call get_nc(fname, vname_j_index,        _remap_data%j_index,  1, nData, ierr, cmessage)
   case default; ierr=20; message=trim(message)//'unable to find variable'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through variables

 end subroutine

end module read_remap
