module read_remap
! data types
use nrtype

! global data
USE public_var

! Netcdf
USE io_netcdf, ONLY: open_nc
USE io_netcdf, ONLY: close_nc
USE io_netcdf, ONLY: get_nc
USE io_netcdf, ONLY: get_nc_dim_len

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
 integer(i4b)                       :: ncidMapping     ! mapping netcdf id
 integer(i4b)                       :: iVar            ! index of variables
 integer(i4b)                       :: nHRU            ! number of HRU in mapping files (this should match up with river network hru)
 integer(i4b)                       :: nData           ! number of data (weight, runoff hru id) in mapping files
 character(len=strLen)              :: cmessage        ! error message from subroutine

 ierr=0; message='get_remap_data/'

 call open_nc(fname, 'r', ncidMapping, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the number of HRUs
 call get_nc_dim_len(ncidMapping, dname_hru_remap, nHRU, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the number of spatial elements in the runoff file
 call get_nc_dim_len(ncidMapping, dname_data_remap, nData, ierr, cmessage)
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
   case(1); call get_nc(ncidMapping, vname_hruid_in_remap, remap_data_in%hru_id,   1, nHRU,  ierr, cmessage)
   case(2); call get_nc(ncidMapping, vname_num_qhru,       remap_data_in%num_qhru, 1, nHRU,  ierr, cmessage)
   case(3); call get_nc(ncidMapping, vname_weight,         remap_data_in%weight,   1, nData, ierr, cmessage)
   case(4); call get_nc(ncidMapping, vname_qhruid,         remap_data_in%qhru_id,  1, nData, ierr, cmessage)
   case(5); call get_nc(ncidMapping, vname_i_index,        remap_data_in%i_index,  1, nData, ierr, cmessage)
   case(6); call get_nc(ncidMapping, vname_j_index,        remap_data_in%j_index,  1, nData, ierr, cmessage)
   case default; ierr=20; message=trim(message)//'unable to find variable'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through variables

 ! check data
 call check_remap_data(remap_data_in, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call close_nc(ncidMapping, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine

 ! *****
 ! public subroutine: check spatial weight files
 ! ********************************************************************************
 subroutine check_remap_data(remap_data_in, &   ! inout: data structure to remap data
                             ierr, message)     ! output: error control
 USE dataTypes,          ONLY : remap           ! remapping data type
 USE nr_utility_module,  ONLY : arth
 implicit none
 ! input/output variables
 type(remap),  intent(inout)        :: remap_data_in     ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
 ! output variables
 integer(i4b), intent(out)          :: ierr              ! error code
 character(*), intent(out)          :: message           ! error message
 ! local variables
 integer(i4b)                       :: ix,jx             ! index of variables
 integer(i4b), allocatable          :: idxZero(:)        ! indices where array variables are zero
 integer(i4b)                       :: nHRU_map          ! number of HRU in mapping files (this should match up with river network hru)
 integer(i4b)                       :: nData             ! number of data (weight, runoff hru id) in mapping files
 integer(i4b)                       :: total_intersects  ! sum of intersected hydrologic model hrus over the entire routing HRUs
 integer(i4b)                       :: nZero             ! number of array variables with zero
 logical(lgt), allocatable          :: logical_array(:)  !
 real(dp),     allocatable          :: real_array(:)     !
 integer(i8b), allocatable          :: int_array(:)      !

 ! initialize error control
 ierr=0; message='check_remap_data/'

 total_intersects = sum(remap_data_in%num_qhru)
 nData            = size(remap_data_in%weight)
 nHRU_map         = size(remap_data_in%num_qhru)

 write(iulog,'(2a)') new_line('a'), 'Check total intersected hydrologic model hrus'
 write(iulog,'(3a,I7)') ' 1) sum of ', trim(vname_num_qhru),    '= ', total_intersects
 write(iulog,'(3a,I7)') ' 2) size of ', trim(dname_data_remap), '= ', nData

 if (total_intersects == nData) then
   write(iulog,'(a)') ' 1) and 2) identical. Should be good to go'
 else
   write(iulog,'(a)')    ' 1) and 2) differ'
   write(iulog,'(3a)')   ' Correct ',trim(dname_data_remap),' dimension variables at routing HRUs without any interesected hydrologic model hrus'

   nZero=count(remap_data_in%num_qhru==0)

   allocate(logical_array(nData))
   logical_array = .true.

   if (nZero>0) then
     allocate(idxZero(nZero))
     idxZero = pack(arth(1,1,nHRU_map), remap_data_in%num_qhru==0)
     do ix = 1,nZero
       jx = sum(remap_data_in%num_qhru(1:idxZero(ix)))+1
       logical_array(jx) = .false.
     end do

     allocate(real_array(total_intersects), int_array(total_intersects))
     if (allocated(remap_data_in%qhru_id)) then
       real_array = pack(remap_data_in%weight, logical_array)
       deallocate(remap_data_in%weight)
       allocate(remap_data_in%weight(total_intersects))
       remap_data_in%weight = real_array
     end if
     if (allocated(remap_data_in%qhru_id)) then
       int_array = pack(remap_data_in%qhru_id, logical_array)
       deallocate(remap_data_in%qhru_id)
       allocate(remap_data_in%qhru_id(total_intersects))
       remap_data_in%qhru_id = int_array
     end if
     if (allocated(remap_data_in%i_index)) then
       int_array = pack(remap_data_in%i_index, logical_array)
       deallocate(remap_data_in%i_index)
       allocate(remap_data_in%i_index(total_intersects))
       remap_data_in%i_index = int_array
     end if
     if (allocated(remap_data_in%j_index)) then
       int_array = pack(remap_data_in%j_index, logical_array)
       deallocate(remap_data_in%j_index)
       allocate(remap_data_in%j_index(total_intersects))
       remap_data_in%j_index = int_array
     end if
   end if
 end if

 end subroutine

end module read_remap
