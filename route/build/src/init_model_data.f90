MODULE init_model_data

! data types
USE nrtype,    ONLY: i4b,dp,lgt,strLen
USE dataTypes, ONLY: var_ilength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_clength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_dlength         ! double precision type: var(:)%dat, or dat

USE var_lookup, ONLY: ixNTOPO            ! index of variables for the network topology
USE var_lookup, ONLY: ixHRU2SEG          ! index of variables for data structure
USE var_lookup, ONLY: ixPFAF             ! index of variables for the pfafstetter code

! Shared data
USE public_var, ONLY: iulog
USE public_var, ONLY: verySmall
USE public_var, ONLY: integerMissing
USE public_var, ONLY: realMissing
USE public_var, ONLY: charMissing

implicit none

integer(i4b),parameter  :: upstream_size=1
integer(i4b),parameter  :: stream_order=2

private
public :: get_mpi_omp
public :: init_model
public :: init_ntopo_data
public :: init_state_data
public :: update_time
public :: init_pio

CONTAINS

 ! *********************************************************************
 ! public subroutine: initialize pio system and decomposition
 ! *********************************************************************
 SUBROUTINE init_pio(pid, nNodes, ierr, message)

  USE pio_utils,       ONLY: pio_sys_init
  USE globalData,      ONLY: pioSystem
  USE globalData,      ONLY: mpicom_route
  USE globalData,      ONLY: pio_numiotasks
  USE globalData,      ONLY: pio_rearranger
  USE globalData,      ONLY: pio_root
  USE globalData,      ONLY: pio_stride
  USE globalData,      ONLY: runMode
  USE pio_decomp_data, ONLY: set_pio_decomp

  implicit none
  ! Argument variables
  integer(i4b), intent(in)    :: pid              ! proc id
  integer(i4b), intent(in)    :: nNodes           ! number of procs
  integer(i4b), intent(out)   :: ierr             ! error code
  character(*), intent(out)   :: message          ! error message
  ! local variables
  character(len=strLen)       :: cmessage         ! error message from subroutine

  ! pio system initialization
  if (trim(runMode)=='standalone') then

    allocate(pioSystem, stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//': pioSystem'; return; endif

    pio_numiotasks = nNodes/pio_stride
    call pio_sys_init(pid, mpicom_route,          & ! input: MPI related parameters
                      pio_stride, pio_numiotasks, & ! input: PIO related parameters
                      pio_rearranger, pio_root,   & ! input: PIO related parameters
                      pioSystem)                    ! output: PIO system descriptors
  end if

  ! pio domain decomposition
  call set_pio_decomp(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE init_pio

 ! *********************************************************************
 ! public subroutine: get mpi and omp info
 ! *********************************************************************
 SUBROUTINE get_mpi_omp(comm)

  ! Obtain mpi rank/ntasks and omp thread number

  USE public_var, ONLY: root         ! root proce id
  USE globalData, ONLY: nNodes       ! number of tasks
  USE globalData, ONLY: masterproc   ! root proc logical
  USE globalData, ONLY: multiProcs   ! mpi multi-procs logical (.true. -> use more than 1 processors)
  USE globalData, ONLY: pid          ! procs id (rank)
  USE globalData, ONLY: nThreads     ! number of OMP threads
  USE mpi_utils, ONLY: shr_mpi_commsize
  USE mpi_utils, ONLY: shr_mpi_commrank

  implicit none
  ! Argument variables
  integer(i4b),  intent(in)  :: comm      ! communicator
  ! local variables
  character(len=strLen)      :: message             ! error message
  integer(i4b)               :: omp_get_num_threads ! number of threads used for openMP

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

  USE public_var,          ONLY: ancil_dir
  USE public_var,          ONLY: input_dir
  USE public_var,          ONLY: param_nml
  USE public_var,          ONLY: gageMetaFile
  USE public_var,          ONLY: outputAtGage
  USE popMetadat_module,   ONLY: popMetadat        ! populate metadata
  USE read_control_module, ONLY: read_control      ! read the control file
  USE read_param_module,   ONLY: read_param        ! read the routing parameters
  USE process_gage_meta,   ONLY: read_gage_meta    ! process gauge metadata

  implicit none
  ! Argument variables
  character(*), intent(in)    :: cfile_name        ! name of the control file
  integer(i4b), intent(out)   :: ierr              ! error code
  character(*), intent(out)   :: message           ! error message
  ! Local variables
  character(len=strLen)       :: cmessage          ! error message of downwind routine

  ierr=0; message='init_model/'

  call popMetadat(ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call read_control(trim(cfile_name), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! read the routing parameter namelist
  call read_param(trim(input_dir)//trim(param_nml),ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! read gauge metadata if specified in control file
  if (outputAtGage .and. trim(gageMetaFile)/=charMissing) then
    call read_gage_meta(trim(ancil_dir)//trim(gageMetaFile),ierr,cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

  call init_route_method(ierr, cmessage)
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
  USE public_var,  ONLY: bypass_routing_option  ! cesm-coupling option
  USE public_var,  ONLY: is_lake_sim            ! logical if lakes are activated in simulation
  USE globalData,  ONLY: masterproc             ! root proc logical
  USE globalData,  ONLY: nHRU, nRch             ! number of HRUs and Reaches in the whole network
  USE globalData,  ONLY: nContribHRU            ! number of HRUs that are connected to any reaches
  USE globalData,  ONLY: basinID                ! HRU id vector
  USE globalData,  ONLY: reachID                ! reach ID vector
  USE globalData,  ONLY: runMode                ! mizuRoute run mode - standalone or ctsm-coupling
  ! external subroutines
  USE read_streamSeg,       ONLY: mod_meta_varFile         ! modify variable I/O options
  USE model_utils,          ONLY: model_finalize
  USE mpi_process,          ONLY: comm_ntopo_data          ! mpi routine: initialize river network data in slave procs (incl. river data transfer from root proc)
  USE mpi_utils,            ONLY: shr_mpi_initialized      ! If MPI is being used
  USE domain_decomposition, ONLY: mpi_domain_decomposition ! domain decomposition for mpi
  USE network_topo,         ONLY: lakeInlet                ! identify lake inlet reach
  USE network_topo,         ONLY: outletSegment            ! subroutine: find oultlet reach id, index  as a destination reach

  implicit none
  ! Argument variables
   integer(i4b),              intent(in)    :: pid              ! proc id
   integer(i4b),              intent(in)    :: nNodes           ! number of procs
   integer(i4b),              intent(in)    :: comm             ! communicator
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
  ! local variable
   type(var_dlength), allocatable           :: structHRU(:)     ! HRU properties
   type(var_dlength), allocatable           :: structSeg(:)     ! stream segment properties
   type(var_ilength), allocatable           :: structHRU2SEG(:) ! HRU-to-segment mapping
   type(var_ilength), allocatable           :: structNTOPO(:)   ! network topology
   type(var_clength), allocatable           :: structPFAF(:)    ! pfafstetter code
   integer(i4b)                             :: iHRU, iRch       ! loop index
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ierr=0; message='init_ntopo_data/'

   ! modify river and lake variable I/O options in meta
   call mod_meta_varFile(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! populate various river network data strucutures for each proc
   if (masterproc) then
     ! read the river network data and compute additonal network attributes (inncludes spatial decomposition)
     call init_ntopo(nHRU, nRch,                                                   & ! output: number of HRU and Reaches
                     structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                     ierr, cmessage)                                                 ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! ---- perform lake option is on
     if (is_lake_sim) then
       call lakeInlet(nRch, structNTOPO, ierr, cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end if

     ! ----  CESM-coupling only
     if (trim(runMode)=='cesm-coupling' .and. &
         trim(bypass_routing_option)=='direct_to_outlet') then
       call outletSegment(nRch, structNTOPO, ierr, message)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end if
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

   ! spatial domain decomposition and distribution for MPI parallelization
   if (masterproc) then
     call mpi_domain_decomposition(nNodes, nRch,               & ! input:
                                   structNTOPO, structHRU2SEG, & ! input:  input data structures
                                   nContribHRU,                & ! output: number of HRUs that are connected to any reaches
                                   ierr, cmessage)               ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   call comm_ntopo_data(pid, nNodes, comm,                                    & ! input: proc id, # of procs and commnicator
                        nRch, nHRU,                                           & ! input: number of reach and HRUs that contribut to any reaches
                        structHRU, structSEG, structHRU2SEG, structNTOPO,     & ! input: river network data structures for the entire network
                        ierr, cmessage)                                         ! output: error controls
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE init_ntopo_data


 ! *********************************************************************
 ! public subroutine: update time to next time
 ! *********************************************************************
 SUBROUTINE update_time(finished, ierr, message)

  USE public_var, ONLY: dt                ! time step [sec]
  USE globalData, ONLY: TSEC              ! beginning/ending of simulation time step [sec]
  USE globalData, ONLY: iTime             ! time index at simulation time step
  USE globalData, ONLY: timeVar           ! model time variables in time unit since reference time
  USE globalData, ONLY: endDatetime       ! model ending datetime
  USE globalData, ONLY: simDatetime       ! current model datetime
  USE write_simoutput_pio, ONLY: close_all

   implicit none
   ! Argument variables
   logical(lgt),              intent(out)   :: finished
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variables
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ierr=0; message='update_time/'

   if (simDatetime(1)>=endDatetime) then
     call close_all(ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     finished=.true.;return
   endif

   ! update model time step bound
   TSEC(1) = TSEC(2)
   TSEC(2) = TSEC(1) + dt

   ! update model time index
   iTime=iTime+1

   ! increment model datetime
   simDatetime(0) = simDatetime(1)
   simDatetime(1) = simDatetime(1)%add_sec(dt, ierr, cmessage)
   simDatetime(2) = simDatetime(1)%add_sec(dt, ierr, cmessage)

   ! increment time variable for history file - keep second now
   timeVar(1) = timeVar(2)
   timeVar(2) = timeVar(1) + dt

 END SUBROUTINE update_time


 ! *********************************************************************
 ! public subroutine: initialize channel state data
 ! *********************************************************************
 SUBROUTINE init_state_data(pid, nNodes, comm, ierr, message)

  ! external routines
  USE ascii_utils,  ONLY: lower             ! convert string to lower case
  USE read_restart, ONLY: read_state_nc     ! read netcdf state output file
  USE mpi_process,  ONLY: mpi_restart
  ! shared data
  USE public_var, ONLY: dt                     ! simulation time step (seconds)
  USE public_var, ONLY: restart_dir            ! directory containing output data
  USE public_var, ONLY: fname_state_in         ! name of state input file
  USE public_var, ONLY: accumRunoff            ! routing method ID
  USE public_var, ONLY: impulseResponseFunc    ! IRF routing ID = 1
  USE public_var, ONLY: kinematicWaveTracking  ! KWT routing ID = 2
  USE public_var, ONLY: kinematicWave          ! KW routing ID = 3
  USE public_var, ONLY: muskingumCunge         ! MC routing ID = 4
  USE public_var, ONLY: diffusiveWave          ! DW routing ID = 5
  USE public_var, ONLY: is_lake_sim            ! logical if lakes are activated in simulation
  USE globalData, ONLY: idxSUM, idxIRF, idxKWT, &
                        idxKW, idxMC, idxDW
  USE globalData, ONLY: nRoutes                ! number of available routing methods
  USE globalData, ONLY: routeMethods           ! ID of active routing method
  USE globalData, ONLY: onRoute                ! logical array for active routing method
  USE globalData, ONLY: masterproc             ! root proc logical
  USE globalData, ONLY: nRch_mainstem          ! number of mainstem reaches
  USE globalData, ONLY: nTribOutlet            ! number of mainstem reaches
  USE globalData, ONLY: rch_per_proc           ! number of tributary reaches
  USE globalData, ONLY: NETOPO_main            ! mainstem reach topology structure
  USE globalData, ONLY: NETOPO_trib            ! tributary reach topology structure
  USE globalData, ONLY: RCHFLX_trib            ! reach flux structure
  USE globalData, ONLY: RCHSTA_trib            ! reach flux structure
  USE globalData, ONLY: RCHFLX                 ! global Reach flux data structures
  USE globalData, ONLY: RCHSTA                 ! global Reach state data structures
  USE globalData, ONLY: nMolecule              ! computational molecule
  USE globalData, ONLY: TSEC                   ! begining/ending of simulation time step [sec]
  USE globalData, ONLY: initHvars              ! status of history variabiable data initialization
  USE globalData, ONLY: isColdStart            ! initial river state - cold start (T) or from restart file (F)
  USE write_simoutput_pio, ONLY: hVars         ! history variable data

  implicit none
  ! Argument variables:
  integer(i4b),        intent(in)  :: pid              ! proc id
  integer(i4b),        intent(in)  :: nNodes           ! number of procs
  integer(i4b),        intent(in)  :: comm             ! communicator
  integer(i4b),        intent(out) :: ierr             ! error code
  character(*),        intent(out) :: message          ! error message
  ! local variable
  real(dp)                         :: T0,T1            ! begining/ending of simulation time step [sec]
  integer(i4b)                     :: nRch_root        ! number of reaches in roor processors consisting (mainstem, halo, and tributary)
  integer(i4b)                     :: ix, iRoute       ! loop indices
  character(len=strLen)            :: cmessage         ! error message of downwind routine

  ierr=0; message='init_state_data/'

  do iRoute = 1, nRoutes
    if (routeMethods(iRoute)==kinematicWave) then
      nMolecule%KW_ROUTE = 2
    else if (routeMethods(iRoute)==muskingumCunge) then
      nMolecule%MC_ROUTE = 2
    else if (routeMethods(iRoute)==diffusiveWave) then
      nMolecule%DW_ROUTE = 20
    end if
  end do

  call init_pio(pid, nNodes, ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (isColdStart) then
    if (masterproc) then
      RCHFLX_trib(:)%BASIN_QI     = 0._dp
      RCHFLX_trib(:)%BASIN_QR(0)  = 0._dp
      RCHFLX_trib(:)%BASIN_QR(1)  = 0._dp
      nRch_root=nRch_mainstem+nTribOutlet+rch_per_proc(0)
      if (onRoute(accumRunoff)) then
        do ix = 1,nRch_root
          RCHFLX_trib(ix)%ROUTE(idxSUM)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(ix)%ROUTE(idxSUM)%REACH_Q        = 0._dp
        end do
      end if
      if (onRoute(impulseResponseFunc)) then
        do ix = 1,nRch_root
          RCHFLX_trib(ix)%ROUTE(idxIRF)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(ix)%ROUTE(idxIRF)%REACH_Q        = 0._dp
        end do
      end if
      if (onRoute(kinematicWaveTracking)) then
        do ix = 1, nRch_mainstem+nTribOutlet
          if (is_lake_sim .and. NETOPO_main(ix)%islake) then
            allocate(RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0:0),stat=ierr, errmsg=cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%LKW_ROUTE%KWAVE]'; return; endif
            RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%QF=-9999
            RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%TI=-9999
            RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%TR=-9999
            RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%RF=.False.
            RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%QM=-9999
          end if
        end do
        do ix = 1, rch_per_proc(0)
          if (is_lake_sim .and. NETOPO_trib(ix)%islake) then
            allocate(RCHSTA_trib(nRch_mainstem+nTribOutlet+ix)%LKW_ROUTE%KWAVE(0:0),stat=ierr, errmsg=cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%LKW_ROUTE%KWAVE]'; return; endif
            RCHSTA_trib(ix+nRch_mainstem+nTribOutlet)%LKW_ROUTE%KWAVE(0)%QF=-9999
            RCHSTA_trib(ix+nRch_mainstem+nTribOutlet)%LKW_ROUTE%KWAVE(0)%TI=-9999
            RCHSTA_trib(ix+nRch_mainstem+nTribOutlet)%LKW_ROUTE%KWAVE(0)%TR=-9999
            RCHSTA_trib(ix+nRch_mainstem+nTribOutlet)%LKW_ROUTE%KWAVE(0)%RF=.False.
            RCHSTA_trib(ix+nRch_mainstem+nTribOutlet)%LKW_ROUTE%KWAVE(0)%QM=-9999
          end if
        end do
        do ix = 1, nRch_root
          RCHFLX_trib(ix)%ROUTE(idxKWT)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(ix)%ROUTE(idxKWT)%REACH_Q        = 0._dp
        end do
      end if
      if (onRoute(kinematicWave)) then
        do ix = 1, nRch_root
          RCHFLX_trib(ix)%ROUTE(idxKW)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(ix)%ROUTE(idxKW)%REACH_Q        = 0._dp
          allocate(RCHSTA_trib(ix)%KW_ROUTE%molecule%Q(nMolecule%KW_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%KW_ROUTE%molecule%Q]'; return; endif
          RCHSTA_trib(ix)%KW_ROUTE%molecule%Q(:) = 0._dp
        end do
      end if
      if (onRoute(muskingumCunge)) then
        do ix = 1, nRch_root
          RCHFLX_trib(ix)%ROUTE(idxMC)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(ix)%ROUTE(idxMC)%REACH_Q        = 0._dp
          allocate(RCHSTA_trib(ix)%MC_ROUTE%molecule%Q(nMolecule%MC_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%MC_ROUTE%molecule%Q]'; return; endif
          RCHSTA_trib(ix)%MC_ROUTE%molecule%Q(:) = 0._dp
        end do
      end if
      if (onRoute(diffusiveWave)) then
        do ix = 1, nRch_root
          RCHFLX_trib(ix)%ROUTE(idxDW)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(ix)%ROUTE(idxDW)%REACH_Q        = 0._dp
          allocate(RCHSTA_trib(ix)%DW_ROUTE%molecule%Q(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%DW_ROUTE%molecule%Q]'; return; endif
          RCHSTA_trib(ix)%DW_ROUTE%molecule%Q(:) = 0._dp
        end do
      end if
    else
      if (rch_per_proc(pid) > 0) then
        RCHFLX_trib(:)%BASIN_QI     = 0._dp
        RCHFLX_trib(:)%BASIN_QR(0)  = 0._dp
        RCHFLX_trib(:)%BASIN_QR(1)  = 0._dp
        if (onRoute(accumRunoff)) then
          do ix = 1, size(RCHFLX_trib)
            RCHFLX_trib(ix)%ROUTE(idxSUM)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxSUM)%REACH_Q        = 0._dp
          end do
        end if
        if (onRoute(impulseResponseFunc)) then
          do ix = 1, size(RCHFLX_trib)
            RCHFLX_trib(ix)%ROUTE(idxIRF)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxIRF)%REACH_Q        = 0._dp
          end do
        end if
        if (onRoute(kinematicWaveTracking)) then
          do ix = 1, size(RCHSTA_trib)
            if (is_lake_sim .and. NETOPO_trib(ix)%islake) then
              allocate(RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0:0),stat=ierr, errmsg=cmessage)
              if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%LKW_ROUTE%KWAVE]'; return; endif
              RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%QF=-9999
              RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%TI=-9999
              RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%TR=-9999
              RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%RF=.False.
              RCHSTA_trib(ix)%LKW_ROUTE%KWAVE(0)%QM=-9999
            end if
          end do
          do ix = 1, size(RCHFLX_trib)
            RCHFLX_trib(ix)%ROUTE(idxKWT)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxKWT)%REACH_Q        = 0._dp
          end do
        end if
        if (onRoute(kinematicWave)) then
          do ix = 1, size(RCHSTA_trib)
            RCHFLX_trib(ix)%ROUTE(idxKW)%FLOOD_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxKW)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxKW)%REACH_Q        = 0._dp
            allocate(RCHSTA_trib(ix)%KW_ROUTE%molecule%Q(nMolecule%KW_ROUTE), stat=ierr, errmsg=cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%KW_ROUTE%molecule%Q]'; return; endif
            RCHSTA_trib(ix)%KW_ROUTE%molecule%Q(:) = 0._dp
          end do
        end if
        if (onRoute(muskingumCunge)) then
          do ix = 1, size(RCHSTA_trib)
            RCHFLX_trib(ix)%ROUTE(idxMC)%FLOOD_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxMC)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxMC)%REACH_Q        = 0._dp
            allocate(RCHSTA_trib(ix)%MC_ROUTE%molecule%Q(nMolecule%MC_ROUTE), stat=ierr, errmsg=cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%MC_ROUTE%molecule%Q]'; return; endif
            RCHSTA_trib(ix)%MC_ROUTE%molecule%Q(:) = 0._dp
          end do
        end if
        if (onRoute(diffusiveWave)) then
          do ix = 1, size(RCHSTA_trib)
            RCHFLX_trib(ix)%ROUTE(idxDW)%FLOOD_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxDW)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(ix)%ROUTE(idxDW)%REACH_Q        = 0._dp
            allocate(RCHSTA_trib(ix)%DW_ROUTE%molecule%Q(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%DW_ROUTE%molecule%Q]'; return; endif
            RCHSTA_trib(ix)%DW_ROUTE%molecule%Q(:) = 0._dp
          end do
        end if
      end if
    end if
    ! initialize time
    TSEC(1)=0._dp; TSEC(2)=dt
  else
    ! start with restart condition
    if (masterproc) then
      call read_state_nc(trim(restart_dir)//trim(fname_state_in), T0, T1, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! time bound [sec] is at previous time step, so need to add dt for curent time step
      TSEC(1)=T0+dt; TSEC(2)=T1+dt
    else ! if other processors, just allocate RCHSTA for dummy
      allocate(RCHSTA(1),RCHFLX(1), stat=ierr, errmsg=cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA,RCHFLX]'; return; endif
    end if

    call hVars%read_restart(trim(restart_dir)//trim(fname_state_in), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif

    initHvars = .true.

    call mpi_restart(pid, nNodes, comm, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
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
  USE ncio_utils,           ONLY: get_var_dims
  USE process_ntopo,        ONLY: augment_ntopo            ! compute all the additional network topology (only compute option = on)

  implicit none
  ! Argument variables
  integer(i4b)                  , intent(out) :: nHRU_out                 ! number of HRUs
  integer(i4b)                  , intent(out) :: nRch_out                 ! number of reaches
  type(var_dlength), allocatable, intent(out) :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(out) :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(out) :: structHRU2SEG(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable, intent(out) :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable, intent(out) :: structPFAF(:)            ! pfafstetter code
  integer(i4b)      , intent(out)             :: ierr                     ! error code
  character(*)      , intent(out)             :: message                  ! error message
  ! Local variables
  integer(i4b)                                :: tot_upstream             ! total number of all of the upstream stream segments for all stream segments
  integer(i4b)                                :: tot_upseg                ! total number of immediate upstream segments for all  stream segments
  integer(i4b)                                :: tot_hru                  ! total number of all the upstream hrus for all stream segments
  integer(i4b)                                :: tot_uh                   ! total number of unit hydrograph from all the stream segments
  integer(i4b),      allocatable              :: ixHRU_desired(:)         ! indices of desired hrus
  integer(i4b),      allocatable              :: ixSeg_desired(:)         ! indices of desired reaches
  integer(i4b)                                :: dummy(2)                 ! dummy variable to hold dimension length for 2D variables in netCDF
  integer(i4b)   , parameter                  :: maxUpstreamFile=90000000 ! 90 million: maximum number of upstream reaches to enable writing
  character(len=strLen)                       :: cmessage                 ! error message of downwind routine

  ierr=0; message='init_ntopo/'

  ! get the variable dimensions
  ! NOTE: need to update maxPfafLen to the exact character size for pfaf code in netCDF
  if (meta_PFAF(ixPFAF%code)%varFile) then
    call get_var_dims(trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
                      trim(meta_PFAF(ixPFAF%code)%varName), & ! input: pfaf code variable name in netcdf
                      ierr, cmessage,                       & ! output: error control
                      dlen=dummy)                             ! output optional: dimension length
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    maxPfafLen = dummy(1)
  end if

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
  ! public subroutine: initialize routing method object
  ! *********************************************************************
  SUBROUTINE init_route_method(ierr, message)
    !
    ! DESCRIPTION:
    ! Instantiate a collection of routing method objects
    !
    USE globalData,         ONLY: rch_routes            ! routing methods instantiated
    USE globalData,         ONLY: routeMethods          ! Active routing method
    USE public_var,         ONLY: accumRunoff           ! routing method ID
    USE public_var,         ONLY: impulseResponseFunc   ! routing method ID
    USE public_var,         ONLY: kinematicWaveTracking ! routing method ID
    USE public_var,         ONLY: kinematicWave         ! routing method ID
    USE public_var,         ONLY: muskingumCunge        ! routing method ID
    USE public_var,         ONLY: diffusiveWave         ! routing method ID
    USE accum_runoff_module,ONLY: accum_runoff_rch      ! routing routine: accumulation instantaneous runoff
    USE irf_route_module,   ONLY: irf_route_rch         ! routing routine: Impulse response function
    USE kwt_route_module,   ONLY: kwt_route_rch         ! routing routine: Lagrangian kinematic
    USE kw_route_module,    ONLY: kwe_route_rch         ! routing routine: kinematic
    USE mc_route_module,    ONLY: mc_route_rch          ! routing routine: muskingum
    USE dfw_route_module,   ONLY: dfw_route_rch         ! routing routine: diffusive

    implicit none
    ! Argument variables:
    integer(i4b),          intent(out)   :: ierr        ! error code
    character(*),          intent(out)   :: message     ! error message
    ! Local variables:
    character(len=strLen)                :: cmessage     ! error message from subroutine
    integer(i4b)                         :: ix

    ierr=0; message='init_route_method/'

    allocate(rch_routes(size(routeMethods)), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [rch_routes]'; return; endif

    do ix=1, size(routeMethods)
      select case (routeMethods(ix))
        case (accumRunoff)
          allocate(accum_runoff_rch :: rch_routes(ix)%rch_route)
        case (impulseResponseFunc)
          allocate(irf_route_rch :: rch_routes(ix)%rch_route)
        case (kinematicWaveTracking)
          allocate(kwt_route_rch :: rch_routes(ix)%rch_route)
        case (kinematicWave)
          allocate(kwe_route_rch :: rch_routes(ix)%rch_route)
        case (muskingumCunge)
          allocate(mc_route_rch  :: rch_routes(ix)%rch_route)
        case (diffusiveWave)
          allocate(dfw_route_rch :: rch_routes(ix)%rch_route)
        case default
          ierr=20; message=trim(message)//'no valid routing method'; return
      end select
    end do

  END SUBROUTINE init_route_method

END MODULE init_model_data
