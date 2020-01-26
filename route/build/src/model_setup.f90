MODULE model_setup

! data types
USE nrtype,    ONLY: i4b,dp,lgt          ! variable types, etc.
USE nrtype,    ONLY: strLen              ! length of characters
USE dataTypes, ONLY: var_ilength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_clength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_dlength         ! double precision type: var(:)%dat, or dat

! Shared data
USE public_var, ONLY: iulog              ! i/o logical unit number
USE public_var, ONLY: verySmall
USE public_var, ONLY: integerMissing
USE public_var, ONLY: realMissing

implicit none

integer(i4b),parameter  :: upstream_size=1
integer(i4b),parameter  :: stream_order=2

! privacy -- everything private unless declared explicitly
private
public :: init_mpi           ! For stand-alone only
public :: get_mpi_omp
public :: init_model
public :: init_data          ! For stand-alone only
public :: init_ntopo_data
public :: init_state_data
public :: update_time
public :: model_finalize

CONTAINS

 ! *********************************************************************
 ! public subroutine: initialize MPI for stand-alone program
 ! *********************************************************************
 SUBROUTINE init_mpi()

  ! Initialize MPI and get OMP thread

  ! shared data used
  USE globalData, ONLY: mpicom_route ! communicator id
  ! subroutines: populate metadata
  USE mpi_mod, ONLY: shr_mpi_init

  implicit none

  ! input:  None
  ! output: None
  ! local variables
  character(len=strLen)       :: message             ! error message

  ! initialize error control
  message='init_mpi/'

  call shr_mpi_init(mpicom_route, message)

  call get_mpi_omp(mpicom_route)

 END SUBROUTINE init_mpi


 ! *********************************************************************
 ! public subroutine: get mpi and omp info
 ! *********************************************************************
 SUBROUTINE get_mpi_omp(comm)

  ! Obtain mpi rank/ntasks and omp thread number

  ! shared data used
  USE public_var, ONLY: root         ! root proce id
  USE globalData, ONLY: nNodes       ! number of tasks
  USE globalData, ONLY: masterproc   ! root proc logical
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
  call read_param(trim(ancil_dir)//trim(param_nml),ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE init_model


 ! *********************************************************************
 ! public subroutine: initialize all the data - For stand-alone
 ! *********************************************************************
 SUBROUTINE init_data(pid,           & ! input: proc id
                      nNodes,        & ! input: number of procs
                      comm,          & ! input: communicator
                      ierr, message)   ! output: error control

   USE globalData,  ONLY: isFileOpen             ! file open/close status
   ! external subroutines
   USE mpi_routine, ONLY: pass_global_data       ! mpi globaldata copy to slave proc

   IMPLICIT NONE
   integer(i4b),              intent(in)    :: pid              ! proc id
   integer(i4b),              intent(in)    :: nNodes           ! number of procs
   integer(i4b),              intent(in)    :: comm             ! communicator
   ! output: error control
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variable
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ! initialize error control
   ierr=0; message='init_data/'

   ! network topology data initialization
   call init_ntopo_data(pid, nNodes, comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! runoff and mapping data initialization
   call init_runoff_data(pid, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! broadcast public and some global variables
   call pass_global_data(comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! channel state initialization
   call init_state_data(pid, nNodes, comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   isFileOpen=.false.

 END SUBROUTINE init_data


 ! *********************************************************************
 ! public subroutine: initialize river network and state data
 ! *********************************************************************
 SUBROUTINE init_ntopo_data(pid,           & ! input: proc id
                            nNodes,        & ! input: number of procs
                            comm,          & ! input: communicator
                            ierr, message)   ! output: error control

  USE var_lookup,  ONLY: ixHRU2SEG              ! index of variables for data structure
  USE var_lookup,  ONLY: ixNTOPO                ! index of variables for data structure
  ! Shared data
  USE public_var,  ONLY: ntopAugmentMode        ! River network augmentation mode
  USE public_var,  ONLY: idSegOut               ! outlet segment ID (-9999 => no outlet segment specified)
  USE globalData,  ONLY: RCHFLX                 ! Reach flux data structures (entire river network)
  USE globalData,  ONLY: RCHSTA                 ! Reach state data structures (entire river network)
  USE globalData,  ONLY: nHRU, nRch             ! number of HRUs and Reaches in the whole network
  USE globalData,  ONLY: nContribHRU            ! number of HRUs that are connected to any reaches
  USE globalData,  ONLY: nEns                   ! number of ensembles
  USE globalData,  ONLY: basinID                ! HRU id vector
  USE globalData,  ONLY: reachID                ! reach ID vector
  ! external subroutines
  USE mpi_routine, ONLY: comm_ntopo_data        ! mpi routine: initialize river network data in slave procs (incl. river data transfer from root proc)

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
   if (pid==0) then

     ! read the river network data and compute additonal network attributes (inncludes spatial decomposition)
     call init_ntopo(nNodes,                                                       & ! input:  number of nodes
                     nHRU, nRch,                                                   & ! output: number of HRU and Reaches
                     structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                     nContribHRU,                                                  & ! output: MPI domain decomposition data
                     ierr, cmessage)                                                 ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end if

   ! check if network topology write option is on. If so, terminate the program
   if (ntopAugmentMode .or. idSegOut>0) then
    call model_finalize(comm)
   end if

   if (pid==0) then

     ! allocate space for the entire river network
     allocate(RCHFLX(nEns,nRch), RCHSTA(nEns,nRch), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [RCHFLX, RCHSTA]'; return; endif

     ! populate basiID and reachID vectors for output (in only master processor)
     ! populate runoff data structure (only meta, no runoff values)
     ! populate remap data structure

     allocate(basinID(nHRU), reachID(nRch), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [basinID, reachID]'; return; endif

     do iHRU = 1,nHRU
      basinID(iHRU) = structHRU2SEG(iHRU)%var(ixHRU2SEG%hruId)%dat(1)
     enddo
     do iRch = 1,nRch
      reachID(iRch) = structNTOPO(iRch)%var(ixNTOPO%segId)%dat(1)
     enddo

   end if  ! if processor=0 (root)

   if (nNodes>1) then
     ! distribute network topology data and network parameters to the different processors
     call comm_ntopo_data(pid, nNodes, comm,                                    & ! input: proc id, # of procs and commnicator
                          nRch, nContribHRU,                                    & ! input: number of reach and HRUs that contribut to any reaches
                          structHRU, structSEG, structHRU2SEG, structNTOPO,     & ! input: river network data structures for the entire network
                          ierr, cmessage)                                         ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else
     ! Setup openmp for whole domain
     call all_domain_omp_decomp(nRch,           & ! input: number of reach and HRUs that contribut to any reaches
                                structNTOPO,    & ! input: river network data structures for the entire network
                                ierr, cmessage)   ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   endif

 END SUBROUTINE init_ntopo_data


 ! ********************************************************************************
 ! public subroutine: initialize runoff, and runoff-mapping data - For stand-alone
 ! ********************************************************************************
 SUBROUTINE init_runoff_data(pid,           & ! input: proc id
                             ierr, message)   ! output: error control

   USE public_var, ONLY: is_remap             ! logical whether or not runnoff needs to be mapped to river network HRU
   USE globalData, ONLY: remap_data           ! runoff mapping data structure
   USE globalData, ONLY: runoff_data          ! runoff data structure

   implicit none
   ! input:
   integer(i4b),              intent(in)    :: pid              ! proc id
   ! output: error control
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local:
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   if (pid==0) then
     ! runoff and remap data initialization (TO DO: split runoff and remap initialization)
     call init_runoff(is_remap,        & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                      remap_data,      & ! output: data structure to remap data
                      runoff_data,     & ! output: data structure for runoff
                      ierr, cmessage)    ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! DateTime initialization
     call init_time(runoff_data%ntime, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end if  ! if processor=0 (root)

 END SUBROUTINE init_runoff_data


 ! *********************************************************************
 ! public subroutine: update time to next time
 ! *********************************************************************
 SUBROUTINE update_time(finished, ierr, message)

  ! Shared data
  USE public_var, ONLY: dt
  USE globalData, ONLY: TSEC          ! beginning/ending of simulation time step [sec]
  USE globalData, ONLY: timeVar       ! time variables (unit given by runoff data)
  USE globalData, ONLY: iTime         ! time index at simulation time step
  USE globalData, ONLY: convTime2Days ! conversion multipliers for time unit of runoff input to day
  USE globalData, ONLY: refJulday     ! julian day: reference
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

   ! update model time step bound
   TSEC(0) = TSEC(0) + dt
   TSEC(1) = TSEC(0) + dt

   if (abs(modJulday-endJulday)<verySmall) then
     call close_output_nc()
     finished=.true.;return
   endif

   ! update time index
   iTime=iTime+1

   ! update the julian day of the model simulation
   modJulday = refJulday + timeVar(iTime)/convTime2Days

 END SUBROUTINE update_time


 ! *********************************************************************
 ! public subroutine: initialize channel state data
 ! *********************************************************************
 SUBROUTINE init_state_data(pid, nNodes, comm, ierr, message)

  ! external routines
  USE read_restart,      ONLY: read_state_nc     ! read netcdf state output file
  USE mpi_routine,       ONLY: mpi_restart
  ! global data
  USE public_var, ONLY: dt                ! simulation time step (seconds)
  USE public_var, ONLY: isRestart         ! restart option: True-> model run with restart, F -> model run with empty channels
  USE public_var, ONLY: routOpt           ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error
  USE public_var, ONLY: fname_state_in    ! name of state input file
  USE public_var, ONLY: output_dir        ! directory containing output data
  USE globalData, ONLY: nRch_mainstem     ! number of mainstem reaches
  USE globalData, ONLY: RCHFLX_main       ! reach flux structure
  USE globalData, ONLY: RCHFLX_trib       ! reach flux structure
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
  character(len=strLen)            :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_state_data/'

  iens = 1_i4b

  ! read restart file and initialize states
  if (isRestart) then

   if (pid==0) then
    call read_state_nc(trim(output_dir)//trim(fname_state_in), routOpt, T0, T1, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    TSEC(0)=T0; TSEC(1)=T1

   end if

   if (nNodes>0) then
     call mpi_restart(pid, nNodes, comm, iens, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

  else
   ! Cold start .......
    ! initialize flux structures

   if (pid==0) then
     if (nRch_mainstem > 0) then
       RCHFLX_main(:,:)%BASIN_QI = 0._dp
       RCHFLX_main(:,:)%BASIN_QR(0) = 0._dp
       RCHFLX_main(:,:)%BASIN_QR(1) = 0._dp
     end if
   end if

   RCHFLX_trib(:,:)%BASIN_QI = 0._dp
   RCHFLX_trib(:,:)%BASIN_QR(0) = 0._dp
   RCHFLX_trib(:,:)%BASIN_QR(1) = 0._dp

   ! initialize time
   TSEC(0)=0._dp; TSEC(1)=dt

  endif

 END SUBROUTINE init_state_data


 ! *********************************************************************
 ! public subroutine: finalize model
 ! *********************************************************************
 SUBROUTINE model_finalize(comm)
  USE globalData, ONLY: masterproc       ! root proc logical
  USE mpi_mod,    ONLY: shr_mpi_finalize ! mpi utilities: shut down mpi
  implicit none
  integer(i4b), intent(in) :: comm   ! communicator

  if (masterproc) then
    write(iulog,*) '----------------------'
    write(iulog,*) ' SUCCESSFUL EXECUTION '
    write(iulog,*) '----------------------'
  end if

  call shr_mpi_finalize(comm)

  stop

 END SUBROUTINE model_finalize


 ! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 SUBROUTINE init_time(nTime,     &    ! input: number of time steps
                      ierr, message)  ! output

  ! subroutines:
  USE process_time_module, ONLY: process_time  ! process time information
  USE io_netcdf,           ONLY: get_nc        ! netcdf input
  ! derived datatype
  USE dataTypes, ONLY: time           ! time data type
  ! Shared data
  USE public_var, ONLY: input_dir     ! directory containing input data
  USE public_var, ONLY: fname_qsim    ! simulated runoff netCDF name
  USE public_var, ONLY: vname_time    ! variable name for time
  USE public_var, ONLY: time_units    ! time units (seconds, hours, or days)
  USE public_var, ONLY: simStart      ! date string defining the start of the simulation
  USE public_var, ONLY: simEnd        ! date string defining the end of the simulation
  USE public_var, ONLY: calendar      ! calendar name
  USE globalData, ONLY: timeVar       ! time variables (unit given by runoff data)
  USE globalData, ONLY: iTime         ! time index at simulation time step
  USE globalData, ONLY: convTime2Days ! conversion multipliers for time unit of runoff input to day
  USE globalData, ONLY: refJulday     ! julian day: reference
  USE globalData, ONLY: startJulday   ! julian day: start of routing simulation
  USE globalData, ONLY: endJulday     ! julian day: end of routing simulation
  USE globalData, ONLY: modJulday     ! julian day: at model time step
  USE globalData, ONLY: modTime       ! model time data (yyyy:mm:dd:hh:mm:ss)

  implicit none

  ! input:
  integer(i4b),              intent(in)    :: nTime
  ! output: error control
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer(i4b)                             :: ix
  character(len=strLen)                    :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_time/'

  ! time initialization
  allocate(timeVar(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time data
  call get_nc(trim(input_dir)//trim(fname_qsim), vname_time, timeVar, 1, nTime, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time multiplier needed to convert time to units of days
  select case( trim( time_units(1:index(time_units,' ')) ) )
   case('seconds'); convTime2Days=86400._dp
   case('hours');   convTime2Days=24._dp
   case('days');    convTime2Days=1._dp
   case default;    ierr=20; message=trim(message)//'unable to identify time units'; return
  end select

  ! extract time information from the control information
  call process_time(time_units,    calendar, refJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refJulday]'; return; endif
  call process_time(trim(simStart),calendar, startJulday, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startJulday]'; return; endif
  call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endJulday]'; return; endif

  ! check that the dates are aligned
  if(endJulday<startJulday) then; ierr=20; message=trim(message)//'simulation end is before simulation start'; return; endif

  ! fast forward time to time index at simStart and save iTime and modJulday
  ! need to convert time unit in timeVar to day
  do ix = 1, nTime
    modJulday = refJulday + timeVar(ix)/convTime2Days
    if( modJulday < startJulday ) cycle
    exit
  enddo
  iTime = ix

  ! initialize previous model time
  !modTime(0:1) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
  modTime(0) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

 END SUBROUTINE init_time


 ! *********************************************************************
 ! private subroutine: initialize river network data
 ! *********************************************************************
 SUBROUTINE init_ntopo(nNodes,                                                       & ! input:  number of nodes
                       nHRU_out, nRch_out,                                           & ! output: number of HRU and Reaches
                       structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                       nContribHRU,                                                  & ! output: number of HRUs that are connected to any reaches
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
  USE globalData, ONLY: NETOPO, RPARAM           !
  ! variable index
  USE var_lookup, ONLY: ixPFAF                   ! index of variables for the pfafstetter code
  ! external subroutines
  USE read_streamSeg,       ONLY: getData                  ! get the ancillary data
  USE write_streamSeg,      ONLY: writeData                ! write the ancillary data
  USE process_ntopo,        ONLY: check_river_properties   ! check if river network data is physically valid
  USE io_netcdf,            ONLY: get_var_dims
  USE process_ntopo,        ONLY: augment_ntopo            ! compute all the additional network topology (only compute option = on)
  USE process_ntopo,        ONLY: put_data_struct          ! populate NETOPO and RPARAM data structure
  USE domain_decomposition, ONLY: mpi_domain_decomposition ! domain decomposition for mpi

  implicit none
  ! input: None
  integer(i4b),                   intent(in)  :: nNodes                   ! number of procs
  ! output (river network data structures for the entire domain)
  integer(i4b)                  , intent(out) :: nHRU_out                 ! number of HRUs
  integer(i4b)                  , intent(out) :: nRch_out                 ! number of reaches
  type(var_dlength), allocatable, intent(out) :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(out) :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(out) :: structHRU2SEG(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable, intent(out) :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable, intent(out) :: structPFAF(:)            ! pfafstetter code
  integer(i4b),                   intent(out) :: nContribHRU              ! number of HRUs that are connected to any reaches
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

  ! copy data to the RPARAM and NETOPO structures
  call put_data_struct(nRch_out, structSEG, structNTOPO, & ! input
                       RPARAM, NETOPO,                   & ! output
                       ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! spatial domain decomposition for MPI parallelization
  call mpi_domain_decomposition(nNodes, nRch_out, structNTOPO, nContribHRU, ierr, cmessage)     !Warning: nHRU /= nContribHRU
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE init_ntopo

 ! *****
 ! private subroutine: get mapping data between runoff hru and river network hru
 ! *********************************************************************
 SUBROUTINE init_runoff(remap_flag,      & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                        remap_data_in,   & ! output: data structure to remap data
                        runoff_data_in,  & ! output: data structure for runoff
                        ierr, message)     ! output: error control

 USE public_var,  ONLY: ancil_dir              ! name of the ancillary directory
 USE public_var,  ONLY: input_dir              ! name of the runoff input directory
 USE public_var,  ONLY: fname_qsim             ! name of simulated runoff netCDF
 USE public_var,  ONLY: fname_remap            ! name of runoff mapping netCDF name
 USE public_var,  ONLY: calendar               ! name of calendar
 USE public_var,  ONLY: time_units             ! time units
 USE globalData,  ONLY: basinID                ! basin ID
 USE dataTypes,   ONLY: remap                  ! remapping data type
 USE dataTypes,   ONLY: runoff                 ! runoff data type
 USE read_runoff, ONLY: read_runoff_metadata   ! read meta data from runoff data
 USE read_remap,  ONLY: get_remap_data         ! read remap data

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
                          ierr, cmessage)                       ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 !write(*,*) 'runoff_data_in%nSpace, nTime, trim(time_units) = ', runoff_data_in%nSpace(:), runoff_data_in%nTime, trim(time_units)

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

 else ! if runoff given in RN_HRU

   allocate(runoff_data_in%hru_ix(size(runoff_data_in%hru_id)), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data_in%hru_ix'; return; endif

   ! get indices of the HRU ids in the runoff file in the routing layer
   call get_qix(runoff_data_in%hru_id,  &    ! input: vector of ids in mapping file
                basinID,                &    ! input: vector of ids in the routing layer
                runoff_data_in%hru_ix,  &    ! output: indices of hru ids in routing layer
                ierr, cmessage)              ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! check
   if(count(runoff_data_in%hru_ix/=integerMissing)/=size(basinID))then
    message=trim(message)//'unable to identify all polygons in the mapping file'
    ierr=20; return
   endif

   ! check that the basins match
   do iHRU = 1, size(runoff_data_in%hru_ix)
     jHRU = runoff_data_in%hru_ix(iHRU)
     if (jHRU == integerMissing) cycle
     if( runoff_data_in%hru_id(iHRU) /= basinID(jHRU) )then
      message=trim(message)//'mismatch in HRU ids for basins in the routing layer'
      ierr=20; return
     endif
   end do

 endif

 END SUBROUTINE init_runoff


 ! *****
 ! private subroutine: get indices of mapping points within runoff file...
 ! ***********************************************************************
 SUBROUTINE get_qix(qid,qidMaster,qix,ierr,message)

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

 qix(1:size(qid)) = integerMissing
 nx=0
 jx=1
 ! loop through id vector
 do ix=1,size(qid)

  ! find match
  do ixMaster=jx,size(qidMaster) ! normally a very short loop

   ! keep track of trials
   nx=nx+1
   !write(*,*) 'qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) ) = ', qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) )

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
   write(iulog,*) trim(message)//'matching ids: ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) ) = ', &
                                                ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) )
  endif

 end do  ! looping through the vector

 write(iulog,*) 'nMissing = ', count(qix==integerMissing)

 ! check again
 do ix=1,size(qid)
  if(qix(ix) /= integerMissing)then
   if(qid(ix) /= qidMaster( qix(ix) ) )then
    message=trim(message)//'unable to find the match'
    ierr=20; return
   endif
  endif
 end do

 END SUBROUTINE get_qix

 ! *********************************************************************
 ! private subroutine: openmp domain decomposition for entire network
 ! *********************************************************************
 SUBROUTINE all_domain_omp_decomp(nRch_in,            & ! input: number of stream segments in whole domain
                                  structNTOPO,        & ! input: data structure for network toopology
                                  ierr,message)         ! output: error control

  ! Perform the entire network openmp domain decomposition
  ! This sub-routine is called only with single mpi proc use

  USE var_lookup,          ONLY: ixNTOPO                  ! index of variables for the network topology
  USE globalData,          ONLY: river_basin_main         ! OMP domain data structure for mainstem
  USE globalData,          ONLY: rch_per_proc             ! number of reaches assigned to each proc (size = num of procs+1)
  USE globalData,          ONLY: ixRch_order              ! global reach index in the order of proc assignment (size = total number of reaches in the entire network)
  USE domain_decomposition,ONLY: omp_domain_decomposition ! domain decomposition for omp

  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: nRch_in                  ! number of total reaches
  type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)           ! network topology
  ! Output error handling variables
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                   ! error message
  ! Local variables
  integer(i4b)                                :: iSeg                      ! reach and hru loop indices
  character(len=strLen)                       :: cmessage                  ! error message from subroutine

  ierr=0; message='all_domain_omp_decomp/'

  ! allocate local and global indices
  allocate(rch_per_proc(-1:0), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [rch_per_proc]'; return; endif

  allocate(ixRch_order(nRch_in), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for [ixRch_order]'; return; endif

  ! Count the number of reaches and hrus in each node
  rch_per_proc(-1) = nRch_in; rch_per_proc(0) = 0
  do iSeg =1,nRch_in
    ixRch_order(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
  end do

  ! OMP domain decomposition
  call omp_domain_decomposition(stream_order, rch_per_proc(-1), structNTOPO, river_basin_main, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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


END MODULE model_setup
