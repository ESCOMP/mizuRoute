module data_initialization

! data types
USE nrtype,    only : i4b,dp,lgt          ! variable types, etc.
USE nrtype,    only : strLen              ! length of characters
USE dataTypes, only : var_ilength         ! integer type:          var(:)%dat
USE dataTypes, only : var_clength         ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength,dlength ! double precision type: var(:)%dat, or dat

USE public_var, only : integerMissing

implicit none

! privacy -- everything private unless declared explicitly
private
public::init_data

contains

 ! *********************************************************************
 ! private subroutine: initialize river network data
 ! *********************************************************************
 subroutine init_data(nNodes, pid, nEns, ixSubHRU, ixSubSEG, remap_data, runoff_data, ierr, message)

  USE public_var,  only : is_remap               ! logical whether or not runnoff needs to be mapped to river network HRU
  USE var_lookup,  only : ixHRU2SEG              ! index of variables for data structure
  USE dataTypes,   only : remap                  ! remapping data type
  USE dataTypes,   only : runoff                 ! runoff data type

   implicit none
   ! input:
   integer(i4b),              intent(in)    :: nNodes           ! number of procs
   integer(i4b),              intent(in)    :: pid              ! proc id
   integer(i4b),              intent(in)    :: nEns             ! number of ensembles
   ! out:
   integer(i4b), allocatable, intent(out)   :: ixSubHRU(:)      ! global HRU index in the order of domains
   integer(i4b), allocatable, intent(out)   :: ixSubSEG(:)      ! global reach index in the order of domains
   type(remap),               intent(out)   :: remap_data       ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
   type(runoff),              intent(out)   :: runoff_data      ! runoff for one time step for all HRUs
   ! output: error control
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variable
   ! river network data structures for the entire domain
   type(var_dlength), allocatable           :: structHRU(:)     ! HRU properties
   type(var_dlength), allocatable           :: structSeg(:)     ! stream segment properties
   type(var_ilength), allocatable           :: structHRU2SEG(:) ! HRU-to-segment mapping
   type(var_ilength), allocatable           :: structNTOPO(:)   ! network topology
   type(var_clength), allocatable           :: structPFAF(:)    ! pfafstetter code
   ! others
   character(len=strLen)                    :: cmessage         ! error message of downwind routine
   integer(i4b)                             :: iHRU             ! loop index
   integer(i4b)     , allocatable           :: hruId_network(:) ! hruID of river network layer
   integer(i4b)                             :: nHRU_network     ! number of HRU in river network
   integer(i4b)                             :: nSpatial(1:2)    ! number of spatial elements in runoff data
   integer(i4b)                             :: nTime            ! number of time steps
   character(len=strLen)                    :: time_units       ! time units
   character(len=strLen)                    :: calendar         ! calendar


   ! initialize error control
   ierr=0; message='init_data/'

   call init_ntopo(nNodes, pid, nEns, structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, ixSubHRU, ixSubSEG, &
                   ierr, message)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (pid==0) then

     ! get vector of HRU ids in the river network layer
     nHRU_network = size(structHRU2SEG)
     allocate(hruId_network(nHRU_network), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [hruId_network]'; return; endif
     forall(iHRU=1:nHRU_network) hruId_network(iHRU) = structHRU2SEG(iHRU)%var(ixHRU2SEG%hruId)%dat(1)

     call init_runoff(&
                      ! data structures
                      hruId_network,   & ! input:  ancillary data for mapping hru2basin
                      is_remap,        & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                      remap_data,      & ! output: data structure to remap data
                      runoff_data,     & ! output: data structure for runoff
                      ! dimensions
                      nSpatial,        & ! output: number of spatial elements in runoff data
                      nTime,           & ! output: number of time steps
                      time_units,      & ! output: time units
                      calendar,        & ! output: calendar
                      ! error control
                      ierr, message)     ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end if

 end subroutine init_data

 ! *********************************************************************
 ! private subroutine: initialize river network data
 ! *********************************************************************
 subroutine init_ntopo(nNodes, pid, nEns, structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, ixSubHRU, ixSubSEG, ierr, message)
  ! public vars
  USE public_var,           only : ancil_dir                ! name of the ancillary directory
  USE public_var,           only : fname_ntopOld            ! name of the old network topology file
  USE public_var,           only : fname_ntopNew            ! name of the new network topology file
  USE public_var,           only : dname_nhru               ! dimension name for HRUs
  USE public_var,           only : dname_sseg               ! dimension name for stream segments
  USE public_var,           only : idSegOut                 ! ID for stream segment at the bottom of the subset
  USE public_var,           only : maxPfafLen               ! maximum digit of pfafstetter code (default 32)
  USE public_var,           only : nThresh                  ! maximum number of segments for each sub-basin
  ! options
  USE public_var,           only : ntopWriteOption          ! option to write updated network topology
  ! common variables
  USE public_var,           only : realMissing              ! missing value for real
  USE public_var,           only : integerMissing           ! missing value for integers
  ! global data
  USE globalData,           only : meta_PFAF                ! meta for pfafstetter code
  ! variable index
  USE var_lookup,           only : ixPFAF                   ! index of variables for the pfafstetter code
  ! external subroutines
  USE read_streamSeg,       only : getData                  ! get the ancillary data
  USE write_streamSeg,      only : writeData                ! write the ancillary data
  USE read_netcdf,          only : get_var_dims
  USE process_ntopo,        only : augment_ntopo            ! compute all the additional network topology (only compute option = on)
  USE domain_decomposition, only : classify_river_basin_mpi
  USE mpi_routine,          only : comm_ntopo_data
  implicit none
  ! input:
  integer(i4b)      , intent(in)              :: nNodes                   ! number of procs
  integer(i4b)      , intent(in)              :: pid                      ! proc id
  integer(i4b)      , intent(in)              :: nEns                     ! number of ensembles
  ! output (river network data structures for the entire domain)
  type(var_dlength), allocatable, intent(out) :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(out) :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(out) :: structHRU2SEG(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable, intent(out) :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable, intent(out) :: structPFAF(:)            ! pfafstetter code
  integer(i4b),      allocatable, intent(out) :: ixSubHRU(:)              ! global HRU index in the order of domains
  integer(i4b),      allocatable, intent(out) :: ixSubSEG(:)              ! global reach index in the order of domains
  ! output: error control
  integer(i4b)      , intent(out)             :: ierr                     ! error code
  character(*)      , intent(out)             :: message                  ! error message
  ! local variable
  integer(i4b)                                :: nHRU                     ! number of HRUs in network data
  integer(i4b)                                :: nSeg                     ! number of stream segments in network data
  integer(i4b)                                :: tot_upstream             ! total number of all of the upstream stream segments for all stream segments
  integer(i4b)                                :: tot_upseg                ! total number of immediate upstream segments for all  stream segments
  integer(i4b)                                :: tot_hru                  ! total number of all the upstream hrus for all stream segments
  integer(i4b)                                :: tot_uh                   ! total number of unit hydrograph from all the stream segments
  integer(i4b),      allocatable              :: ixHRU_desired(:)         ! indices of desired hrus
  integer(i4b),      allocatable              :: ixSeg_desired(:)         ! indices of desired reaches
  integer(i4b)                                :: dummy(2)                 ! dummy variable to hold dimension length for 2D variables in netCDF
  integer(i4b)   , parameter                  :: maxUpstreamFile=10000000 ! 10 million: maximum number of upstream reaches to enable writing
  character(len=strLen)                       :: cmessage                 ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_ntopo/'

  if (pid==0) then
  ! need to update maxPfafLen to the exact character size for pfaf code in netCDF
  call get_var_dims(trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
                    trim(meta_PFAF(ixPFAF%code)%varName), & ! input: pfaf code variable name in netcdf
                    ierr, cmessage,                       & ! output: error control
                    dlen=dummy)                             ! output optional: dimension length
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  maxPfafLen = dummy(1)

  ! read river network data
  call getData(&
               ! input
               trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
               dname_nhru,   & ! input: dimension name of the HRUs
               dname_sseg,   & ! input: dimension name of the stream segments
               ! output: model control
               nHRU,         & ! output: number of HRUs
               nSeg,         & ! output: number of stream segments
               ! output: populate data structures
               structHRU,    & ! ancillary data for HRUs
               structSeg,    & ! ancillary data for stream segments
               structHRU2seg,& ! ancillary data for mapping hru2basin
               structNTOPO,  & ! ancillary data for network topology
               structPFAF,   & ! ancillary data for pfafstetter code
               ! output: error control
               ierr,cmessage) ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call augment_ntopo(&
                  ! input: model control
                  nHRU,             & ! number of HRUs
                  nSeg,             & ! number of stream segments
                  ! inout: populate data structures
                  structHRU,        & ! ancillary data for HRUs
                  structSeg,        & ! ancillary data for stream segments
                  structHRU2seg,    & ! ancillary data for mapping hru2basin
                  structNTOPO,      & ! ancillary data for network toopology
                  ! output
                  tot_hru,          & ! total number of all the upstream hrus for all stream segments
                  tot_upseg,        & ! total number of all the immediate upstream segments for all stream segments
                  tot_upstream,     & ! total number of all the upstream segments for all stream segments
                  tot_uh,           & ! total number of unit hydrograph for all stream segments
                  ixHRU_desired,    & ! indices of desired hrus
                  ixSeg_desired,    & ! indices of desired reaches
                  ! output: error control
                  ierr, message)

  ! Write out augmented network topology if desired
  if(idSegOut>0) ntopWriteOption=.true.   ! ensure that localized network topology is written if a particular outlet is specified

  if(ntopWriteOption)then

    ! disable the dimension containing all upstream reaches
    ! NOTE: For the CONUS this is 1,872,516,819 reaches !!
    !        --> it will always be quicker to recompute than read+write
    !        --> users can modify the hard-coded parameter "maxUpstreamFile" if desired
    if(tot_upstream > maxUpstreamFile) tot_upstream=0

    ! remove file if it exists
    call system('rm -f '//trim(ancil_dir)//trim(fname_ntopNew))

    ! write data
    call writeData(&
                   ! input
                   trim(ancil_dir)//trim(fname_ntopNew), & ! input: file name
                   ! input: model control
                   tot_hru,       & ! input: total number of all the upstream hrus for all stream segments
                   tot_upseg,     & ! input: total number of immediate upstream segments for all  stream segments
                   tot_upstream,  & ! input: total number of all of the upstream stream segments for all stream segments
                   tot_uh,        & ! input: total number of unit hydrograph for all stream segments
                   ! input: reach masks
                   ixHRU_desired, & ! input: indices of desired hrus
                   ixSeg_desired, & ! input: indices of desired reaches
                   ! input: data structures
                   structHRU,     & ! input: ancillary data for HRUs
                   structSeg,     & ! input: ancillary data for stream segments
                   structHRU2seg, & ! input: ancillary data for mapping hru2basin
                   structNTOPO,   & ! input: ancillary data for network topology
                   structPFAF,    & ! input: ancillary data for pfafstetter code
                   ! output: error control
                   ierr,cmessage) ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif
  endif

  ! created a subset = sucessful execution: Need to run again with the subset
  if(idSegOut>0)then
   write(*,'(a)') 'Running in subsetting mode'
   write(*,'(a)') 'Created a subset network topology file '//trim(fname_ntopNew)
   write(*,'(a)') ' --> Run again using the new network topology file '
   write(*,'(a)') ' SUCCESSFUL EXECUTION '
   stop
  endif

  if (pid==0) then
    call classify_river_basin_mpi(nNodes, nSeg, structPFAF, structNTOPO, nThresh, nHRU, ierr, cmessage)     !Warning: nHRU may change
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  call comm_ntopo_data(pid, nNodes, nEns, nSeg, nHRU, structHRU, structSEG, structHRU2SEG, structNTOPO, &
                       ixSubHRU, ixSubSEG, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine init_ntopo

 ! *****
 ! public subroutine: get mapping data between runoff hru and river network hru
 ! *********************************************************************
 subroutine init_runoff(&
                         ! data structures
                         hruID_network,   & ! input:  ancillary data for mapping hru2basin
                         remap_flag,      & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                         remap_data,      & ! output: data structure to remap data
                         runoff_data,     & ! output: data structure for runoff
                         ! dimensions
                         nSpatial,        & ! output: number of spatial elements in runoff data
                         nTime,           & ! output: number of time steps
                         time_units,      & ! output: time units
                         calendar,        & ! output: calendar
                         ! error control
                         ierr, message)     ! output: error control

 USE public_var,  only : ancil_dir              ! name of the ancillary directory
 USE public_var,  only : input_dir              ! name of the runoff input directory
 USE public_var,  only : fname_qsim             ! name of simulated runoff netCDF
 USE public_var,  only : fname_remap            ! name of runoff mapping netCDF name
 USE dataTypes,   only : remap                  ! remapping data type
 USE dataTypes,   only : runoff                 ! runoff data type

 USE read_runoff, only : get_runoff_metadata  ! read meta data from runoff data
 USE read_remap,  only : get_remap_data       ! read remap data

 implicit none
 ! data structures
 integer(i4b)     , intent(in)      :: hruID_network(:) ! ID of routing layer
 logical(lgt),      intent(in)      :: remap_flag       ! logical whether or not runnoff needs to be mapped to river network HRU
 type(remap)  , intent(out)         :: remap_data       ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
 type(runoff) , intent(out)         :: runoff_data      ! runoff for one time step for all HRUs
 ! ancillary data
 integer(i4b) , intent(out)         :: nSpatial(1:2)    ! number of spatial elements
 integer(i4b) , intent(out)         :: nTime            ! number of time steps
 character(*) , intent(out)         :: time_units       ! time units
 character(*) , intent(out)         :: calendar         ! calendar
 ! error control
 integer(i4b), intent(out)          :: ierr             ! error code
 character(*), intent(out)          :: message          ! error message
 ! local variables
 character(len=strLen)              :: cmessage         ! error message from subroutine

 ! initialize error control
 ierr=0; message='init_runoff/'

 ! get runoff metadata
 call get_runoff_metadata(trim(input_dir)//trim(fname_qsim), & ! input: filename
                          runoff_data,                       & ! output: runoff data structure
                          nSpatial,                          & ! output: number spatial elements
                          nTime, time_units, calendar,       & ! output: number of time steps, time units, calendar
                          ierr, cmessage)                      ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 !print*, 'nSpatial, nTime, trim(time_units) = ', nSpatial(:), nTime, trim(time_units)
 !print*, trim(message)//'PAUSE : '; read(*,*)

 ! need to remap runoff to HRUs
 if (remap_flag) then

   ! get runoff mapping file
   call get_remap_data(trim(ancil_dir)//trim(fname_remap),     & ! input: file name
                       nSpatial,                               & ! input: number of spatial elements
                       remap_data,                             & ! output: data structure to remap data from a polygon
                       ierr, cmessage)                           ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! get indices of the HRU ids in the mapping file in the routing layer
   call get_qix(remap_data%hru_id,  &    ! input: vector of ids in mapping file
                hruID_network,      &    ! input: vector of ids in the routing layer
                remap_data%hru_ix,  &    ! output: indices of hru ids in routing layer
                ierr, cmessage)          ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if ( nSpatial(2) == integerMissing ) then
     ! get indices of the "overlap HRUs" (the runoff input) in the runoff vector
     call get_qix(remap_data%qhru_id, &    ! input: vector of ids in mapping file
                  runoff_data%hru_id, &    ! input: vector of ids in runoff file
                  remap_data%qhru_ix, &    ! output: indices of mapping ids in runoff file
                  ierr, cmessage)          ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   ! check
   if(count(remap_data%hru_ix/=integerMissing)/=size(hruID_network))then
    message=trim(message)//'unable to identify all polygons in the mapping file'
    ierr=20; return
   endif
 endif
 !print*, trim(message)//'PAUSE : '; read(*,*)

 end subroutine init_runoff

 ! ============================================================================================
 ! ============================================================================================
 ! ============================================================================================
 ! ============================================================================================
 ! ============================================================================================
 ! ============================================================================================

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

end module data_initialization
