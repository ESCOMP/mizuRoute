MODULE init_model_data

! data types
USE nrtype,    ONLY: i4b,dp,lgt          ! variable types, etc.
USE nrtype,    ONLY: strLen              ! length of characters
USE dataTypes, ONLY: var_ilength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_clength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_dlength         ! double precision type: var(:)%dat, or dat

USE var_lookup, ONLY: ixNTOPO            ! index of variables for the network topology
USE var_lookup, ONLY: ixHRU2SEG          ! index of variables for data structure
USE var_lookup, ONLY: ixPFAF             ! index of variables for the pfafstetter code

! Shared data
USE public_var, ONLY: iulog              ! i/o logical unit number
USE public_var, ONLY: verySmall
USE public_var, ONLY: integerMissing
USE public_var, ONLY: realMissing
USE public_var, ONLY: charMissing

implicit none

integer(i4b),parameter  :: upstream_size=1
integer(i4b),parameter  :: stream_order=2

! privacy -- everything private unless declared explicitly
private
public :: get_mpi_omp
public :: init_model
public :: init_ntopo_data
public :: init_state_data
public :: update_time

CONTAINS

 ! *********************************************************************
 ! public subroutine: get mpi and omp info
 ! *********************************************************************
 SUBROUTINE get_mpi_omp(comm)

  ! Obtain mpi rank/ntasks and omp thread number

  ! shared data used
  USE public_var, ONLY: root         ! root proce id
  USE globalData, ONLY: nNodes       ! number of tasks
  USE globalData, ONLY: masterproc   ! root proc logical
  USE globalData, ONLY: multiProcs   ! mpi multi-procs logical (.true. -> use more than 1 processors)
  USE globalData, ONLY: pid          ! procs id (rank)
  USE globalData, ONLY: nThreads     ! number of OMP threads
  ! subroutines: populate metadata
  USE mpi_mod, ONLY: shr_mpi_commsize
  USE mpi_mod, ONLY: shr_mpi_commrank

  implicit none

  ! input:  None
  integer(i4b),  intent(in)  :: comm      ! communicator
  ! output: None
  ! local variables
  character(len=strLen)      :: message             ! error message
  integer(i4b)               :: omp_get_num_threads ! number of threads used for openMP

  ! initialize error control
  message='get_mpi_omp/'

  ! Get the number of processes
  call shr_mpi_commsize(comm, nNodes, message)

  ! Get the individual process ID
  call shr_mpi_commrank(comm, pid, message)

  if (pid == root) then
     masterproc = .true.
  else
     masterproc = .false.
  end if

  if (nNodes > 1) then
     multiProcs = .true.
  else
     multiProcs = .false.
  endif

  !  Get number of threads
  nThreads = 1
  !$OMP PARALLEL
  !$ nThreads = omp_get_num_threads()
  !$OMP END PARALLEL

 END SUBROUTINE get_mpi_omp


 ! *********************************************************************
 ! public subroutine: model setup
 ! *********************************************************************
 SUBROUTINE init_model(cfile_name, ierr, message)

  ! used to read control files and namelist and broadcast to all processors

  ! shared data used
  USE public_var, ONLY: ancil_dir
  USE public_var, ONLY: input_dir
  USE public_var, ONLY: param_nml
  ! subroutines: populate metadata
  USE popMetadat_module, ONLY: popMetadat       ! populate metadata
  ! subroutines: model control
  USE read_control_module, ONLY: read_control     ! read the control file
  USE read_param_module,   ONLY: read_param       ! read the routing parameters

  implicit none

  character(*), intent(in)    :: cfile_name        ! name of the control file
  ! output: error control
  integer(i4b), intent(out)   :: ierr              ! error code
  character(*), intent(out)   :: message           ! error message
  ! local variables
  character(len=strLen)       :: cmessage          ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_model/'

  ! populate the metadata files
  call popMetadat(ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! read the control file
  call read_control(trim(cfile_name), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! read the routing parameter namelist
  call read_param(trim(input_dir)//trim(param_nml),ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE init_model


 ! *********************************************************************
 ! public subroutine: initialize river network and state data
 ! *********************************************************************
 SUBROUTINE init_ntopo_data(pid,           & ! input: proc id
                            nNodes,        & ! input: number of procs
                            comm,          & ! input: communicator
                            ierr, message)   ! output: error control

  ! Shared data
  USE public_var,  ONLY: ntopAugmentMode        ! River network augmentation mode
  USE public_var,  ONLY: idSegOut               ! outlet segment ID (-9999 => no outlet segment specified)
  USE globalData,  ONLY: masterproc             ! root proc logical
  USE globalData,  ONLY: multiProcs             ! mpi multi-procs logical (.true. -> use more than 1 processors)
  USE globalData,  ONLY: nHRU, nRch             ! number of HRUs and Reaches in the whole network
  USE globalData,  ONLY: nRch_mainstem          ! number of mainstem reaches
  USE globalData,  ONLY: nHRU_mainstem          ! number of mainstem HRUs
  USE globalData,  ONLY: RCHFLX_main            ! Reach flux data structures (master proc, mainstem)
  USE globalData,  ONLY: RCHSTA_main            ! Reach state data structures (master proc, mainstem)
  USE globalData,  ONLY: NETOPO_main            !
  USE globalData,  ONLY: RPARAM_main            !
  USE globalData,  ONLY: nContribHRU            ! number of HRUs that are connected to any reaches
  USE globalData,  ONLY: nEns                   ! number of ensembles
  USE globalData,  ONLY: basinID                ! HRU id vector
  USE globalData,  ONLY: reachID                ! reach ID vector
  ! external subroutines
  USE model_utils, ONLY: model_finalize
  USE mpi_routine, ONLY: comm_ntopo_data        ! mpi routine: initialize river network data in slave procs (incl. river data transfer from root proc)
  USE process_ntopo, ONLY: put_data_struct      ! populate NETOPO and RPARAM data structure
  USE mpi_mod      , ONLY: shr_mpi_initialized  ! If MPI is being used
  USE domain_decomposition, ONLY: mpi_domain_decomposition ! domain decomposition for mpi

   implicit none
   ! input:
   integer(i4b),              intent(in)    :: pid              ! proc id
   integer(i4b),              intent(in)    :: nNodes           ! number of procs
   integer(i4b),              intent(in)    :: comm             ! communicator
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
   ierr=0; message='init_ntopo_data/'

   ! populate various river network data strucutures for each proc
   if (masterproc) then

     ! read the river network data and compute additonal network attributes (inncludes spatial decomposition)
     call init_ntopo(nHRU, nRch,                                                   & ! output: number of HRU and Reaches
                     structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                     ierr, cmessage)                                                 ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end if

   ! check if network topology write option is on. If so, terminate the program
   if (ntopAugmentMode .or. idSegOut>0) then
    call model_finalize(comm)
   end if

   if (masterproc) then

     ! populate basiID and reachID vectors for output (in only master processor)

     allocate(basinID(nHRU), reachID(nRch), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [basinID, reachID]'; return; endif

     do iHRU = 1,nHRU
      basinID(iHRU) = structHRU2SEG(iHRU)%var(ixHRU2SEG%hruId)%dat(1)
     enddo
     do iRch = 1,nRch
      reachID(iRch) = structNTOPO(iRch)%var(ixNTOPO%segId)%dat(1)
     enddo

   end if  ! if processor=0 (root)

   if (multiProcs) then

     ! spatial domain decomposition for MPI parallelization
     if (masterproc) then
       call mpi_domain_decomposition(nNodes, nRch,               & ! input:
                                     structNTOPO, structHRU2SEG, & ! input:  input data structures
                                     nContribHRU,                & ! output: number of HRUs that are connected to any reaches
                                     ierr, cmessage)               ! output: error controls
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end if

     ! distribute network topology data and network parameters to the different processors
     call comm_ntopo_data(pid, nNodes, comm,                                    & ! input: proc id, # of procs and commnicator
                          nRch, nHRU,                                           & ! input: number of reach and HRUs that contribut to any reaches
                          structHRU, structSEG, structHRU2SEG, structNTOPO,     & ! input: river network data structures for the entire network
                          ierr, cmessage)                                         ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   else

     ! Count a number of contributory HRUs
     nContribHRU = 0
     do iRch = 1, nRch
       nContribHRU = nContribHRU + structNTOPO(iRch)%var(ixNTOPO%nHRU)%dat(1)
     enddo

     nRch_mainstem = nRch
     nHRU_mainstem = nHRU

     allocate(RCHFLX_main(nEns, nRch_mainstem), RCHSTA_main(nEns, nRch_mainstem), stat=ierr, errmsg=cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     call put_data_struct(nRch_mainstem, structSEG, structNTOPO, & ! input
                          RPARAM_main, NETOPO_main,              & ! output
                          ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! Setup openmp for whole domain
     call all_domain_omp_decomp(structNTOPO, structHRU2SEG, & ! input: river network data structures for the entire network
                                ierr, cmessage)               ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   endif

 END SUBROUTINE init_ntopo_data


 ! *********************************************************************
 ! public subroutine: update time to next time
 ! *********************************************************************
 SUBROUTINE update_time(finished, ierr, message)

  USE public_var, ONLY: dt            ! time step [sec]
  USE globalData, ONLY: TSEC          ! beginning/ending of simulation time step [sec]
  USE globalData, ONLY: iTime         ! time index at simulation time step
  USE globalData, ONLY: roJulday      ! julian day: runoff input time
  USE globalData, ONLY: modJulday     ! julian day: at model time step
  USE globalData, ONLY: endJulday     ! julian day: at end of simulation
  ! external routine
  USE write_simoutput_pio, ONLY: close_output_nc

   implicit none
   ! output: error control
   logical(lgt),              intent(out)   :: finished
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message

   ! initialize error control
   ierr=0; message='update_time/'

   if (abs(modJulday-endJulday)<verySmall) then
     call close_output_nc()
     finished=.true.;return
   endif

   ! update model time step bound
   TSEC(0) = TSEC(0) + dt
   TSEC(1) = TSEC(0) + dt

   ! update time index
   iTime=iTime+1

   ! update the julian day of the model simulation
   modJulday = roJulday(iTime)

 END SUBROUTINE update_time


 ! *********************************************************************
 ! public subroutine: initialize channel state data
 ! *********************************************************************
 SUBROUTINE init_state_data(pid, nNodes, comm, ierr, message)

  ! external routines
  USE read_restart,      ONLY: read_state_nc     ! read netcdf state output file
  USE mpi_routine,       ONLY: mpi_restart
  ! shared data
  USE public_var, ONLY: dt                ! simulation time step (seconds)
  USE public_var, ONLY: routOpt           ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error
  USE public_var, ONLY: output_dir        ! directory containing output data
  USE public_var, ONLY: fname_state_in    ! name of state input file
  USE public_var, ONLY: kinematicWaveEuler!
  USE globalData, ONLY: masterproc        ! root proc logical
  USE globalData, ONLY: nRch_mainstem     ! number of mainstem reaches
  USE globalData, ONLY: rch_per_proc      ! number of tributary reaches
  USE globalData, ONLY: RCHFLX_main       ! reach flux structure
  USE globalData, ONLY: RCHFLX_trib       ! reach flux structure
  USE globalData, ONLY: RCHSTA_main       ! reach flux structure
  USE globalData, ONLY: RCHSTA_trib       ! reach flux structure
  USE globalData, ONLY: TSEC              ! begining/ending of simulation time step [sec]

  implicit none
  ! input:
  integer(i4b),        intent(in)  :: pid              ! proc id
  integer(i4b),        intent(in)  :: nNodes           ! number of procs
  integer(i4b),        intent(in)  :: comm             ! communicator
  ! output: error control
  integer(i4b),        intent(out) :: ierr             ! error code
  character(*),        intent(out) :: message          ! error message
  ! local variable
  real(dp)                         :: T0,T1            ! begining/ending of simulation time step [sec]
  integer(i4b)                     :: iens             ! ensemble index (currently only 1)
  integer(i4b)                     :: ix               ! loop index
  character(len=strLen)            :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_state_data/'

  iens = 1_i4b

  ! read restart file and initialize states
  if (trim(fname_state_in)/=charMissing) then

   if (masterproc) then
    call read_state_nc(trim(output_dir)//trim(fname_state_in), routOpt, T0, T1, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! time bound [sec] is at previous time step, so need to add dt for curent time step
    TSEC(0)=T0+dt; TSEC(1)=T1+dt

   end if

   if (nNodes>0) then
     call mpi_restart(pid, nNodes, comm, iens, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

  else
   ! Cold start .......
    ! initialize flux structures

   if (masterproc) then

     if (nRch_mainstem > 0) then
       RCHFLX_main(:,:)%BASIN_QI = 0._dp
       RCHFLX_main(:,:)%BASIN_QR(0) = 0._dp
       RCHFLX_main(:,:)%BASIN_QR(1) = 0._dp

       if (routOpt==kinematicWaveEuler) then
         do ix = 1,4
           RCHSTA_main(:,:)%EKW_ROUTE%A(ix) = 0._dp
           RCHSTA_main(:,:)%EKW_ROUTE%Q(ix) = 0._dp
         end do
       end if

     end if

   end if

   if (rch_per_proc(pid) > 0) then

     RCHFLX_trib(:,:)%BASIN_QI = 0._dp
     RCHFLX_trib(:,:)%BASIN_QR(0) = 0._dp
     RCHFLX_trib(:,:)%BASIN_QR(1) = 0._dp

     if (routOpt==kinematicWaveEuler) then
       do ix = 1,4
         RCHSTA_trib(:,:)%EKW_ROUTE%A(ix) = 0._dp
         RCHSTA_trib(:,:)%EKW_ROUTE%Q(ix) = 0._dp
       end do
     end if

   end if

   ! initialize time
   TSEC(0)=0._dp; TSEC(1)=dt

  endif

 END SUBROUTINE init_state_data


 ! *********************************************************************
 ! private subroutine: initialize river network data
 ! *********************************************************************
 SUBROUTINE init_ntopo(nHRU_out, nRch_out,                                           & ! output: number of HRU and Reaches
                       structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                       ierr, message)                                                  ! output: error controls
  ! Shared data
  USE public_var, ONLY: ancil_dir                ! name of the ancillary directory
  USE public_var, ONLY: fname_ntopOld            ! name of the old network topology file
  USE public_var, ONLY: fname_ntopNew            ! name of the new network topology file
  USE public_var, ONLY: dname_nhru               ! dimension name for HRUs
  USE public_var, ONLY: dname_sseg               ! dimension name for stream segments
  USE public_var, ONLY: maxPfafLen               ! maximum digit of pfafstetter code (default 32)
  ! options
  USE public_var, ONLY: ntopAugmentMode          ! River network augmentation mode
  USE public_var, ONLY: idSegOut                 ! River network subset mode (idSegOut > 0)
  ! common variables
  USE public_var, ONLY: realMissing              ! missing value for real
  USE public_var, ONLY: integerMissing           ! missing value for integers
  ! global data
  USE globalData, ONLY: meta_PFAF                ! meta for pfafstetter code
  ! external subroutines
  USE read_streamSeg,       ONLY: getData                  ! get the ancillary data
  USE write_streamSeg,      ONLY: writeData                ! write the ancillary data
  USE process_ntopo,        ONLY: check_river_properties   ! check if river network data is physically valid
  USE io_netcdf,            ONLY: get_var_dims
  USE process_ntopo,        ONLY: augment_ntopo            ! compute all the additional network topology (only compute option = on)

  implicit none
  ! input: None
  ! output (river network data structures for the entire domain)
  integer(i4b)                  , intent(out) :: nHRU_out                 ! number of HRUs
  integer(i4b)                  , intent(out) :: nRch_out                 ! number of reaches
  type(var_dlength), allocatable, intent(out) :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(out) :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(out) :: structHRU2SEG(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable, intent(out) :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable, intent(out) :: structPFAF(:)            ! pfafstetter code
  ! output: error control
  integer(i4b)      , intent(out)             :: ierr                     ! error code
  character(*)      , intent(out)             :: message                  ! error message
  ! local variable
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

  ! get the variable dimensions
  ! NOTE: need to update maxPfafLen to the exact character size for pfaf code in netCDF
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

  call check_river_properties(structNTOPO, structHRU, structSEG, ierr, cmessage) ! input: data structure for physical river network data
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute additional network attributes
  call augment_ntopo(&
                     ! input: model control
                     nHRU_out,                         & ! number of HRUs
                     nRch_out,                         & ! number of stream segments
                     ! inout: populate data structures
                     structHRU,                        & ! ancillary data for HRUs
                     structSeg,                        & ! ancillary data for stream segments
                     structHRU2seg,                    & ! ancillary data for mapping hru2basin
                     structNTOPO,                      & ! ancillary data for network toopology
                     ! output:
                     ierr, cmessage,                   & ! error control
                     ! optional output
                     tot_hru       = tot_hru,          & ! total number of all the upstream hrus for all stream segments
                     tot_upseg     = tot_upseg,        & ! total number of all the immediate upstream segments for all stream segments
                     tot_upstream  = tot_upstream,     & ! total number of all the upstream segments for all stream segments
                     tot_uh        = tot_uh,           & ! total number of unit hydrograph for all stream segments
                     ixHRU_desired = ixHRU_desired,    & ! indices of desired hrus
                     ixSeg_desired = ixSeg_desired)      ! indices of desired reaches
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write network topology (if augment mode or subset mode)
  if(ntopAugmentMode .or. idSegOut>0)then

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

    if (idSegOut>0) write(iulog,'(a)') 'Running in river network subset mode'
    if (ntopAugmentMode) write(iulog,'(a)') 'Running in river network augmentation mode'
    write(iulog,'(a)') 'Created a new network topology file '//trim(fname_ntopNew)
    write(iulog,'(a)') ' --> Run again using the new network topology file '
    return
  endif

 END SUBROUTINE init_ntopo


 ! *********************************************************************
 ! private subroutine: openmp domain decomposition for entire network
 ! *********************************************************************
 SUBROUTINE all_domain_omp_decomp(structNTOPO,structHRU2SEG, & ! input: data structure for network toopology
                                  ierr,message)                ! output: error control

  ! Perform the entire network openmp domain decomposition
  ! This sub-routine is called only with single mpi proc use

  USE public_var,          ONLY: desireId
  USE globalData,          ONLY: ixPrint
  USE globalData,          ONLY: river_basin_main         ! OMP domain data structure for mainstem
  USE globalData,          ONLY: rch_per_proc             ! number of reaches assigned to each proc (size = num of procs+1)
  USE globalData,          ONLY: hru_per_proc             ! number of hrus assigned to each proc (size = num of procs+1)
  USE globalData,          ONLY: ixRch_order              ! global reach index in the order of proc assignment (size = total number of reaches in the entire network)
  USE globalData,          ONLY: ixHRU_order              ! global HRU index in the order of proc assignment (size = total number of HRUs contributing to any reaches, nContribHRU)
  USE domain_decomposition,ONLY: omp_domain_decomposition ! domain decomposition for omp
  USE nr_utility_module,   ONLY: findIndex                ! find index within a vector

  implicit none
  ! Input variables
  type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)           ! network topology
  type(var_ilength), allocatable, intent(in)  :: structHRU2SEG(:)         ! HRU-to-segment mapping
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                  ! error message
  ! Local variables
  integer(i4b)                                :: nRch_tmp                 ! number of total reaches
  integer(i4b)                                :: nHRU_tmp                 ! number of total HRUs
  integer(i4b), allocatable                   :: segId(:)                 ! reach id
  integer(i4b)                                :: iSeg, iHRU               ! reach and hru loop indices
  character(len=strLen)                       :: cmessage                 ! error message from subroutine

  ierr=0; message='all_domain_omp_decomp/'

  nRch_tmp = size(structNTOPO)
  nHRU_tmp = size(structHRU2SEG)

  ! allocate local and global indices
  allocate(rch_per_proc(-1:0), hru_per_proc(-1:0), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  allocate(ixRch_order(nRch_tmp), segId(nRch_tmp), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  allocate(ixHRU_order(nHRU_tmp), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! Count the number of reaches and hrus in each node (order of hru and reach should be the same as input data
  rch_per_proc(-1) = nRch_tmp
  rch_per_proc(0)  = 0
  do iSeg =1,nRch_tmp
    segId(iSeg)       = structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1)
    ixRch_order(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
  end do

  hru_per_proc(-1) = nHRU_tmp
  hru_per_proc(0)  = 0
  do iHRU =1,nHRU_tmp
    ixHRU_order(iHRU) = structHRU2SEG(iHRU)%var(ixHRU2SEG%HRUindex)%dat(1)
  end do

  ! OMP domain decomposition
  call omp_domain_decomposition(stream_order, rch_per_proc(-1), structNTOPO, river_basin_main, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (desireId/=integerMissing) then
    ixPrint(1) = findIndex(segId, desireId, integerMissing);
    print*, 'desireId=',desireId, 'ixPrint(1)=', ixPrint(1)
  end if

  ! -------------------
  ! print*,'segid,branch,order'
  ! do ix = 1, size(river_basin_main)
  !   do ixx = 1, size(river_basin_main(ix)%branch)
  !     do iSeg = 1, river_basin_main(ix)%branch(ixx)%nRch
  !       associate (idx_tmp => river_basin_main(ix)%branch(ixx)%segIndex(iSeg))
  !       write(*,"(I15,A,I9,A,I9)") structNTOPO(idx_tmp)%var(ixNTOPO%segId)%dat(1),',',ixx,',',ix
  !       end associate
  !     end do
  !   end do
  ! enddo
  ! -------------------

 END SUBROUTINE all_domain_omp_decomp


END MODULE init_model_data
