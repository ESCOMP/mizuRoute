module model_setup

! data types
USE nrtype,    only : i4b,dp,lgt          ! variable types, etc.
USE nrtype,    only : strLen              ! length of characters
USE dataTypes, only : var_ilength         ! integer type:          var(:)%dat
USE dataTypes, only : var_clength         ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength,dlength ! double precision type: var(:)%dat, or dat

USE public_var, only : integerMissing
USE public_var, only : realMissing

implicit none

! privacy -- everything private unless declared explicitly
private
public :: init_model
public :: init_data

contains

 ! *********************************************************************
 ! public subroutine: model setupt
 ! *********************************************************************
 subroutine init_model(pid, cfile_name, ierr, message)

  ! shared data used
  USE public_var,          only : ancil_dir
  USE public_var,          only : param_nml
  ! subroutines: populate metadata
  USE popMetadat_module,   only : popMetadat       ! populate metadata
  ! subroutines: model control
  USE read_control_module, only : read_control     ! read the control file
  USE read_param_module,   only : read_param       ! read the routing parameters

  implicit none

  integer(i4b), intent(in)    :: pid               ! process id
  character(*), intent(in)    :: cfile_name        ! name of the control file
  ! output: error control
  integer(i4b), intent(out)   :: ierr              ! error code
  character(*), intent(out)   :: message           ! error message
  ! local variables
  character(len=strLen)       :: cmessage          ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_model/'

  ! if the master node
  if (pid==0) then

   ! populate the metadata files
   call popMetadat(ierr,cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! read the control file
   call read_control(trim(cfile_name), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! read the routing parameter namelist
   call read_param(trim(ancil_dir)//trim(param_nml),ierr,cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif  ! if the master node
 end subroutine init_model


 ! *********************************************************************
 ! public subroutine: initialize river network, runoff, and runoff-mapping data
 ! *********************************************************************
 subroutine init_data(nNodes,        & ! input: number of procs
                      pid,           & ! input: proc id
                      ixSubHRU,      & ! output:
                      ixSubSEG,      & ! output:
                      hru_per_proc,  & ! output:
                      seg_per_proc,  & ! output:
                      ierr, message)

  USE public_var,  only : is_remap               ! logical whether or not runnoff needs to be mapped to river network HRU
  USE var_lookup,  only : ixHRU2SEG              ! index of variables for data structure
  USE var_lookup,  only : ixNTOPO                ! index of variables for data structure
  USE globalData,  only : nHRU, nRch             ! number of HRUs and Reaches in the whole network
  USE globalData,  only : basinID                ! HRU id vector
  USE globalData,  only : reachID                ! reach ID vector
  USE globalData,  only : runoff_data            ! runoff data structure
  USE globalData,  only : remap_data             ! runoff mapping data structure

   implicit none
   ! input:
   integer(i4b),              intent(in)    :: nNodes           ! number of procs
   integer(i4b),              intent(in)    :: pid              ! proc id
   ! out:
   integer(i4b), allocatable, intent(out)   :: ixSubHRU(:)      ! global HRU index in the order of domains
   integer(i4b), allocatable, intent(out)   :: ixSubSEG(:)      ! global reach index in the order of domains
   integer(i4b), allocatable, intent(out)   :: hru_per_proc(:)  ! number of hrus assigned to each proc (i.e., node)
   integer(i4b), allocatable, intent(out)   :: seg_per_proc(:)  ! number of reaches assigned to each proc (i.e., node)
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
   integer(i4b)                             :: iHRU, iRch       ! loop index
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ! initialize error control
   ierr=0; message='init_data/'

   call init_ntopo(nNodes, pid,                                                  & ! input:  number of Procs and proc id
                   nHRU, nRch,                                                   & ! output: number of HRU and Reaches
                   structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                   ixSubHRU, ixSubSEG, hru_per_proc, seg_per_proc,               & ! output: MPI domain decomposition data
                   ierr, cmessage)                                                 ! output: error controls
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (pid==0) then

     allocate(basinID(nHRU), reachID(nRch), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [basinID, reachID]'; return; endif

     forall(iHRU=1:nHRU) basinID(iHRU) = structHRU2SEG(iHRU)%var(ixHRU2SEG%hruId)%dat(1)
     forall(iRch=1:nRch) reachID(iRch) = structNTOPO(iRch)%var(ixNTOPO%segId)%dat(1)

     ! runoff and remap data initialization (TO DO: split runoff and remap initialization)
     call init_runoff(&
                      is_remap,        & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                      remap_data,      & ! output: data structure to remap data
                      runoff_data,     & ! output: data structure for runoff
                      ierr, message)     ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! DateTime initialization
     call init_time(runoff_data%ntime, ierr, message)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! channel state initialization
     call init_state(ierr, message)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end if

 end subroutine init_data

 ! *********************************************************************
 ! private subroutine: initialize channel state data
 ! *********************************************************************
 subroutine init_state(ierr, message)
  ! subroutines
  USE read_restart,  only : read_state_nc     ! read netcdf state output file
  USE write_restart, only : define_state_nc   ! define netcdf state output file
  ! global data
  USE public_var,    only : dt                ! simulation time step (seconds)
  USE public_var,    only : isRestart         ! restart option: True-> model run with restart, F -> model run with empty channels
  USE public_var,    only : routOpt           ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error
  USE public_var,    only : fname_state_in    ! name of state input file
  USE public_var,    only : fname_state_out   ! name of state output file
  USE public_var,    only : output_dir        ! directory containing output data
  USE public_var,    only : time_units        ! time units (seconds, hours, or days)
  USE globalData,    only : RCHFLX            ! reach flux structure
  USE globalData,    only : TSEC              ! begining/ending of simulation time step [sec]

  implicit none

  ! output: error control
  integer(i4b),        intent(out) :: ierr             ! error code
  character(*),        intent(out) :: message          ! error message
  ! local variable
  real(dp)                         :: T0,T1            ! begining/ending of simulation time step [sec]
  integer(i4b)                     :: ix               ! index for the stream segment
  character(len=strLen)            :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_state/'

  ! read restart file and initialize states
  if (isRestart) then

   call read_state_nc(trim(output_dir)//trim(fname_state_in), routOpt, T0, T1, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   TSEC(0)=T0; TSEC(1)=T1

  else

   ! Cold start .......
   ! initialize flux structures
   RCHFLX(:,:)%BASIN_QI = 0._dp
   forall(ix=0:1) RCHFLX(:,:)%BASIN_QR(ix) = 0._dp

   ! initialize time
   TSEC(0)=0._dp; TSEC(1)=dt

  endif

  ! Define output state netCDF
  call define_state_nc(trim(output_dir)//trim(fname_state_out), time_units, routOpt, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine init_state

! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 subroutine init_time(nTime,     &    ! input: number of time steps
                      ierr, message)  ! output

  ! subroutines:
  USE process_time_module, only : process_time  ! process time information
  USE read_netcdf,         only : get_nc        ! netcdf input
  ! derived datatype
  USE dataTypes,           only : time          ! time data type
  ! public data
  USE public_var,          only : input_dir     ! directory containing input data
  USE public_var,          only : fname_qsim    ! simulated runoff netCDF name
  USE public_var,          only : vname_time    ! variable name for time
  USE public_var,          only : time_units    ! time units (seconds, hours, or days)
  USE public_var,          only : simStart      ! date string defining the start of the simulation
  USE public_var,          only : simEnd        ! date string defining the end of the simulation
  USE public_var,          only : calendar      ! calendar name
  USE public_var,          only : startJulday   ! julian day: start
  USE public_var,          only : endJulday     ! julian day: end
  USE public_var,          only : refJulday     ! julian day: reference
  USE globalData,          only : timeVar       ! time variables (unit given by runoff data)
  USE globalData,          only : modTime       ! model time data (yyyy:mm:dd:hh:mm:ss)

  implicit none

  ! input:
  integer(i4b),              intent(in)    :: nTime
  ! output: error control
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  character(len=strLen)                    :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_time/'

  ! time initialization
  allocate(timeVar(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time data
  call get_nc(trim(input_dir)//trim(fname_qsim), vname_time, timeVar, 1, nTime, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! extract time information from the control information
  call process_time(time_units,    calendar, refJulday,   ierr, cmessage); if(ierr/=0) message=trim(message)//trim(cmessage)//' [refJulday]'; return
  call process_time(trim(simStart),calendar, startJulday, ierr, cmessage); if(ierr/=0) message=trim(message)//trim(cmessage)//' [startJulday]'; return
  call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage); if(ierr/=0) message=trim(message)//trim(cmessage)//' [endJulday]'; return

  ! check that the dates are aligned
  if(endJulday<startJulday) ierr=20; message=trim(message)//'simulation end is before simulation start'; return

  ! initialize previous model time
  modTime(0:1) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

 end subroutine init_time


 ! *********************************************************************
 ! private subroutine: initialize river network data
 ! *********************************************************************
 subroutine init_ntopo(nNodes, pid,                                                  & ! input:  number of Procs and proc id
                       nHRU_out, nRch_out,                                             & ! output: number of HRU and Reaches
                       structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                       ixSubHRU, ixSubSEG, hru_per_proc, seg_per_proc,               & ! output: MPI domain decomposition data
                       ierr, message)                                                  ! output: error controls
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
  USE domain_decomposition, only : classify_river_basin_mpi ! domain decomposition for mpi
  USE mpi_routine,          only : comm_ntopo_data
  implicit none
  ! input:
  integer(i4b)      , intent(in)              :: nNodes                   ! number of procs
  integer(i4b)      , intent(in)              :: pid                      ! proc id
  ! output (river network data structures for the entire domain)
  integer(i4b)                  , intent(out) :: nHRU_out                 ! number of HRUs
  integer(i4b)                  , intent(out) :: nRch_out                 ! number of reaches
  type(var_dlength), allocatable, intent(out) :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(out) :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(out) :: structHRU2SEG(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable, intent(out) :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable, intent(out) :: structPFAF(:)            ! pfafstetter code
  integer(i4b),      allocatable, intent(out) :: ixSubHRU(:)              ! global HRU index in the order of domains
  integer(i4b),      allocatable, intent(out) :: ixSubSEG(:)              ! global reach index in the order of domains
  integer(i4b),      allocatable, intent(out) :: hru_per_proc(:)          ! number of hrus assigned to each proc (i.e., node)
  integer(i4b),      allocatable, intent(out) :: seg_per_proc(:)          ! number of reaches assigned to each proc (i.e., node)
  ! output: error control
  integer(i4b)      , intent(out)             :: ierr                     ! error code
  character(*)      , intent(out)             :: message                  ! error message
  ! local variable
  integer(i4b)                                :: nContribHRU              ! number of HRUs that are connected to any reaches
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
               nHRU_out,      & ! output: number of HRUs
               nRch_out,      & ! output: number of stream segments
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
                  nHRU_out,         & ! number of HRUs
                  nRch_out,         & ! number of stream segments
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

  ! created a subset = sucessful execution: Need to run again with the subset
  if(idSegOut>0)then
   write(*,'(a)') 'Running in subsetting mode'
   write(*,'(a)') 'Created a subset network topology file '//trim(fname_ntopNew)
   write(*,'(a)') ' --> Run again using the new network topology file '
   write(*,'(a)') ' SUCCESSFUL EXECUTION '
   stop
  endif

  endif

  if (pid==0) then
    call classify_river_basin_mpi(nNodes, nRch_out, structPFAF, structNTOPO, nThresh, nContribHRU, ierr, cmessage)       !Warning: nHRU /= nContribHRU
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  call comm_ntopo_data(pid, nNodes, nRch_out, nContribHRU, structHRU, structSEG, structHRU2SEG, structNTOPO, & ! input
                       ixSubHRU, ixSubSEG, hru_per_proc, seg_per_proc, ierr, cmessage)                         ! output
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine init_ntopo

 ! *****
 ! public subroutine: get mapping data between runoff hru and river network hru
 ! *********************************************************************
 subroutine init_runoff(&
                         ! data structures
                         remap_flag,      & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                         remap_data_in,   & ! output: data structure to remap data
                         runoff_data_in,  & ! output: data structure for runoff
                         ! error control
                         ierr, message)     ! output: error control

 USE public_var,  only : ancil_dir              ! name of the ancillary directory
 USE public_var,  only : input_dir              ! name of the runoff input directory
 USE public_var,  only : fname_qsim             ! name of simulated runoff netCDF
 USE public_var,  only : fname_remap            ! name of runoff mapping netCDF name
 USE public_var,  only : calendar               ! name of calendar
 USE public_var,  only : time_units             ! time units
 USE dataTypes,   only : remap                  ! remapping data type
 USE dataTypes,   only : runoff                 ! runoff data type
 USE read_runoff, only : read_runoff_metadata   ! read meta data from runoff data
 USE read_remap,  only : get_remap_data         ! read remap data
 USE globalData,  only : basinID                ! basin ID

 implicit none
 ! data structures
 logical(lgt),      intent(in)      :: remap_flag       ! logical whether or not runnoff needs to be mapped to river network HRU
 type(remap)  , intent(out)         :: remap_data_in    ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
 type(runoff) , intent(out)         :: runoff_data_in   ! runoff for one time step for all HRUs
 ! error control
 integer(i4b), intent(out)          :: ierr             ! error code
 character(*), intent(out)          :: message          ! error message
 ! local variables
 integer(i4b)                       :: iHRU, jHRU       ! hru loop indices
 character(len=strLen)              :: cmessage         ! error message from subroutine

 ! initialize error control
 ierr=0; message='init_runoff/'

 ! get runoff metadata
 call read_runoff_metadata(trim(input_dir)//trim(fname_qsim), & ! input: filename
                          runoff_data_in,                     & ! output: runoff data structure
                          time_units, calendar,               & ! output: number of time steps, time units, calendar
                          ierr, cmessage)                      ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 !print*, 'runoff_data_in%nSpace, nTime, trim(time_units) = ', runoff_data_in%nSpace(:), runoff_data_in%nTime, trim(time_units)
 !print*, trim(message)//'PAUSE : '; read(*,*)

 ! need to remap runoff to HRUs
 if (remap_flag) then

   ! get runoff mapping file
   call get_remap_data(trim(ancil_dir)//trim(fname_remap),     & ! input: file name
                       runoff_data_in%nSpace,                  & ! input: number of spatial elements
                       remap_data_in,                          & ! output: data structure to remap data from a polygon
                       ierr, cmessage)                           ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! get indices of the HRU ids in the mapping file in the routing layer
   call get_qix(remap_data_in%hru_id, &  ! input: vector of ids in mapping file
                basinID,              &  ! input: vector of ids in the routing layer
                remap_data_in%hru_ix, &  ! output: indices of hru ids in routing layer
                ierr, cmessage)          ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if ( runoff_data_in%nSpace(2) == integerMissing ) then
     ! get indices of the "overlap HRUs" (the runoff input) in the runoff vector
     call get_qix(remap_data_in%qhru_id, &  ! input: vector of ids in mapping file
                  runoff_data_in%hru_id, &  ! input: vector of ids in runoff file
                  remap_data_in%qhru_ix, &  ! output: indices of mapping ids in runoff file
                  ierr, cmessage)           ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   ! check
   if(count(remap_data_in%hru_ix/=integerMissing)/=size(basinID))then
    message=trim(message)//'unable to identify all polygons in the mapping file'
    ierr=20; return
   endif

   ! check that the basins match
   do iHRU = 1, size(remap_data_in%hru_ix)
     jHRU = remap_data_in%hru_ix(iHRU)
     if (jHRU == integerMissing) cycle
     if( remap_data_in%hru_id(iHRU) /= basinID(jHRU) )then
      message=trim(message)//'mismatch in HRU ids for basins in the routing layer'
      ierr=20; return
     endif
   end do

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

end module model_setup
