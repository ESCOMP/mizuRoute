module data_initialization

! data types
USE nrtype,    only : i4b,dp,lgt          ! variable types, etc.
USE nrtype,    only : strLen              ! length of characters
USE dataTypes, only : var_ilength         ! integer type:          var(:)%dat
USE dataTypes, only : var_clength         ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength,dlength ! double precision type: var(:)%dat, or dat

implicit none

! privacy -- everything private unless declared explicitly
private
public::initial_ntopo

contains

 ! *********************************************************************
 ! public subroutine: initialize river network data
 ! *********************************************************************
 subroutine initial_ntopo(nNodes, pid, ierr, message)

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
  integer(i4b)      , intent(in)  :: nNodes                   ! number of procs
  integer(i4b)      , intent(in)  :: pid                      ! proc id
  ! output: error control
  integer(i4b)      , intent(out) :: ierr                     ! error code
  character(*)      , intent(out) :: message                  ! error message
  ! local variable
  integer(i4b)                    :: nHRU                     ! number of HRUs in network data
  integer(i4b)                    :: nSeg                     ! number of stream segments in network data
  type(var_dlength), allocatable  :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable  :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable  :: structHRU2seg(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable  :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable  :: structPFAF(:)            ! pfafstetter code
  integer(i4b)                    :: tot_upstream             ! total number of all of the upstream stream segments for all stream segments
  integer(i4b)                    :: tot_upseg                ! total number of immediate upstream segments for all  stream segments
  integer(i4b)                    :: tot_hru                  ! total number of all the upstream hrus for all stream segments
  integer(i4b)                    :: tot_uh                   ! total number of unit hydrograph from all the stream segments
  integer(i4b),      allocatable  :: ixHRU_desired(:)         ! indices of desired hrus
  integer(i4b),      allocatable  :: ixSeg_desired(:)         ! indices of desired reaches
  integer(i4b)                    :: dummy(2)                 ! dummy variable to hold dimension length for 2D variables in netCDF
  integer(i4b)   , parameter      :: maxUpstreamFile=10000000 ! 10 million: maximum number of upstream reaches to enable writing
  character(len=strLen)           :: cmessage                 ! error message of downwind routine

  ! initialize error control
  ierr=0; message='initial_ntopo/'

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
                  structPFAF,       & ! ancillary data for pfafstetter code
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

  call comm_ntopo_data(pid, nNodes, nSeg, nHRU, structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine initial_ntopo


end module data_initialization
