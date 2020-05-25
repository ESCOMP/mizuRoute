module read_remap
! data types
use nrtype

! global data
USE public_var

! Netcdf
use io_netcdf, only:get_nc
use io_netcdf, only:get_nc_dim_len

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
                           remap_data_in, &   ! output: data structure to remap data from a polygon
                           ierr, message)     ! output: error control
 USE dataTypes,  only : remap                 ! remapping data type
 implicit none
 ! input variables
 character(*), intent(in)           :: fname           ! filename
 integer(i4b), intent(in)           :: nSpatial(1:2)   ! number of spatial elements
 ! output variables
 type(remap),  intent(out)          :: remap_data_in   ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
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
 allocate(remap_data_in%hru_id(nHRU), remap_data_in%num_qhru(nHRU), &
          remap_data_in%hru_ix(nHRU), remap_data_in%weight(nData), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for common mapping info '; return; endif

 ! if runoff input is hru vector...
 if (nSpatial(2) == integerMissing) then
  allocate(remap_data_in%qhru_id(nData), remap_data_in%qhru_ix(nData),  stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for mapping info for vector runoff'; return; endif
 ! if runoff input is grid...
 else
  allocate(remap_data_in%i_index(nData),remap_data_in%j_index(nData), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for mapping info'; return; endif
 endif

 ! get data
 do iVar=1,6
  if (nSpatial(2) == integerMissing .and. iVar >= 5) cycle
  if (nSpatial(2) /= integerMissing .and. iVar == 4) cycle
  select case(iVar)
   case(1); call get_nc(fname, vname_hruid_in_remap, remap_data_in%hru_id,   1, nHRU,  ierr, cmessage)
   case(2); call get_nc(fname, vname_num_qhru,       remap_data_in%num_qhru, 1, nHRU,  ierr, cmessage)
   case(3); call get_nc(fname, vname_weight,         remap_data_in%weight,   1, nData, ierr, cmessage)
   case(4); call get_nc(fname, vname_qhruid,         remap_data_in%qhru_id,  1, nData, ierr, cmessage)
   case(5); call get_nc(fname, vname_i_index,        remap_data_in%i_index,  1, nData, ierr, cmessage)
   case(6); call get_nc(fname, vname_j_index,        remap_data_in%j_index,  1, nData, ierr, cmessage)
   case default; ierr=20; message=trim(message)//'unable to find variable'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through variables

 end subroutine

end module read_remap
