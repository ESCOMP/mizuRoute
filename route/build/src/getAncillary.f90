module getAncillary_module

! data types
use nrtype
USE dataTypes, only : remap                  ! remapping data type
USE dataTypes, only : runoff                 ! runoff data type
USE dataTypes, only : var_ilength            ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength            ! double precision type: var(:)%dat

! lookup variables
USE var_lookup,only:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup,only:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup,only:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup,only:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology

! global data
USE public_var

! Netcdf
use netcdf
use read_netcdf, only:get_nc
use read_netcdf, only:get_nc_dim_len
use read_netcdf, only:get_units

implicit none

! privacy
private
public::getAncillary

contains

 ! *****
 ! public subroutine: get mapping data between runoff hru and river network hru
 ! *********************************************************************
 subroutine getAncillary(&
                         ! data structures
                         nHRU,            & ! input:  number of HRUs in the routing layer
                         structHRU2seg,   & ! input:  ancillary data for mapping hru2basin
                         remap_data,      & ! output: data structure to remap data
                         runoff_data,     & ! output: data structure for runoff
                         ! dimensions
                         nSpatial,        & ! output: number of spatial elements in runoff data
                         nTime,           & ! output: number of time steps
                         time_units,      & ! output: time units
                         ! error control
                         ierr, message)     ! output: error control
 implicit none
 ! data structures
 integer(i4b)     , intent(in)      :: nHRU             ! number of HRUs
 type(var_ilength), intent(in)      :: structHRU2seg(:) ! HRU-to-segment mapping
 type(remap)  , intent(out)         :: remap_data       ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
 type(runoff) , intent(out)         :: runoff_data      ! runoff for one time step for all HRUs
 ! ancillary data
 integer(i4b) , intent(out)         :: nSpatial         ! number of spatial elements
 integer(i4b) , intent(out)         :: nTime            ! number of time steps
 character(*) , intent(out)         :: time_units       ! time units
 ! error control
 integer(i4b), intent(out)          :: ierr             ! error code
 character(*), intent(out)          :: message          ! error message
 ! -----------------------------------------------------------------------------------------------------------------

 ! local variables
 integer(i4b)                       :: iHRU             ! HRU array index
 integer(i4b)                       :: route_id(nHRU)   ! ID of routing layer
 character(len=strLen)              :: cmessage         ! error message from subroutine

 ! initialize error control
 ierr=0; message='getAncillary/'

 ! get runoff metadata
 call get_runoff_metadata(trim(input_dir)//trim(fname_qsim), & ! input: filename
                          nSpatial, runoff_data%hru_id,      & ! output: number spatial elements and HRU ids
                          nTime,    time_units,              & ! output: number of time steps and time units
                          ierr,     cmessage)                  ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for simulated runoff
 allocate(runoff_data%qSim(nSpatial), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating qsim'; return; endif

 !print*, 'nSpatial, nTime, trim(time_units) = ', nSpatial, nTime, trim(time_units)
 !print*, trim(message)//'PAUSE : '; read(*,*)

 ! get runoff mapping file
 call get_remap_data(trim(ancil_dir)//trim(fname_remap),     & ! input: file name
                     remap_data,                             & ! output: data structure to remap data from a polygon
                     ierr, cmessage)                           ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get vector of HRU ids in the routing layer
 forall(iHRU=1:nHRU) route_id(iHRU) = structHRU2seg(iHRU)%var(ixHRU2seg%hruId)%dat(1)

 ! get indices of the HRU ids in the mapping file in the routing layer
 call get_qix(remap_data%hru_id,  &    ! input: vector of ids in mapping file
              route_id,           &    ! input: vector of ids in the routing layer
              remap_data%hru_ix,  &    ! output: indices of hru ids in routing layer
              ierr, cmessage)          ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get indices of the "overlap HRUs" (the runoff input) in the runoff vector
 call get_qix(remap_data%qhru_id, &    ! input: vector of ids in mapping file
              runoff_data%hru_id, &    ! input: vector of ids in runoff file
              remap_data%qhru_ix, &    ! output: indices of mapping ids in runoff file
              ierr, cmessage)          ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! check
 if(count(remap_data%hru_ix/=integerMissing)/=nHRU)then
  message=trim(message)//'unable to identify all polygons in the mapping file'
  ierr=20; return
 endif

 !print*, trim(message)//'PAUSE : '; read(*,*)

 end subroutine getAncillary

 ! ============================================================================================
 ! ============================================================================================
 ! ============================================================================================
 ! ============================================================================================
 ! ============================================================================================
 ! ============================================================================================

 ! *****
 ! private subroutine: get runoff metadata...
 ! ******************************************
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
 integer(i4b), intent(out) , allocatable :: hruId(:)        ! vector of HRU Ids
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
 ! private subroutine: get mapping data between runoff hru and river network hru...
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

 ! allocate space for info in mapping file
 allocate(remap_data%hru_id(nHRU), remap_data%num_qhru(nHRU), &
          remap_data%qhru_id(nData), remap_data%weight(nData), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for mapping info'; return; endif

 ! allocate space for index vectors
 allocate(remap_data%hru_ix(nHRU), remap_data%qhru_ix(nData), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for index vectors'; return; endif

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

 ! *****
 ! private subroutine: get indices of mapping points within runoff file...
 ! ***********************************************************************
 subroutine get_qix(qid,qidMaster,qix,ierr,message)
 USE nr_utility_module, ONLY: indexx  ! get rank of data value
 implicit none
 ! input
 integer(i4b), intent(in)  :: qid(:)                       ! ID of input vector
 integer(i4b), intent(in)  :: qidMaster(:)                 ! ID of master vector
 ! output
 integer(i4b), intent(out) :: qix(:)                       ! index within master vector
 integer(i4b), intent(out) :: ierr                         ! error code
 character(*), intent(out) :: message                      ! error message
 ! local
 integer(i4b)             :: rankID( size(qid) )           ! rank of input vector
 integer(i4b)             :: rankMaster( size(qidMaster) ) ! rank of master vector
 integer(i4b)             :: ix,jx,ixMaster                ! array indices
 integer(i4b)             :: nx                            ! counter
 ! initialize error control
 ierr=0; message='get_qix/'

 ! sort the data vector from smallest to largest
 call indexx(qid,       rankID)
 call indexx(qidMaster, rankMaster)

 !print*, 'rankId = ', rankId(1:10)
 !print*, 'qId( rankId(1:10) ) = ', qId( rankId(1:10) )
 !print*, 'PAUSE: '; read(*,*)
 qix(1:size(qid)) = integerMissing
 nx=0
 jx=1
 ! loop through id vector
 do ix=1,size(qid)

  ! find match
  do ixMaster=jx,size(qidMaster) ! normally a very short loop

   ! keep track of trials
   nx=nx+1
   !print*, 'qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) ) = ', qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) )

   ! find match
   if( qid( rankId(ix) ) == qidMaster( rankMaster(ixMaster) ) )then
    qix( rankId(ix) ) = rankMaster(ixMaster)
    jx = ixMaster
    exit
   endif

   ! unable to find match
   if( qidMaster( rankMaster(ixMaster) ) > qid( rankId(ix) ) )then
    qix( rankId(ix) ) = integerMissing
    jx = ixMaster
    exit
   endif

  end do  ! ixMaster

  ! print progress
  if(qix( rankId(ix) )/=integerMissing .and. mod(ix,1000000)==0)then
   print*, trim(message)//'matching ids: ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) ) = ', &
                                         ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) )
   !print*, 'PAUSE : '; read(*,*)
  endif

 end do  ! looping through the vector

 print*, 'nMissing = ', count(qix==integerMissing)

 ! check again
 do ix=1,size(qid)
  if(qix(ix) /= integerMissing)then
   if(qid(ix) /= qidMaster( qix(ix) ) )then
    message=trim(message)//'unable to find the match'
    ierr=20; return
   endif
  endif
 end do

 end subroutine get_qix

end module getAncillary_module
