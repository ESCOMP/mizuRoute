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

CONTAINS

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

  ! shared data used
  USE public_var, ONLY: ancil_dir
  USE public_var, ONLY: input_dir
  USE public_var, ONLY: param_nml
  USE public_var, ONLY: gageMetaFile
  USE public_var, ONLY: outputAtGage
  ! subroutines: populate metadata
  USE popMetadat_module, ONLY: popMetadat       ! populate metadata
  ! subroutines: model control
  USE read_control_module, ONLY: read_control     ! read the control file
  USE read_param_module,   ONLY: read_param       ! read the routing parameters
  USE process_gage_meta,   ONLY: read_gage_meta   ! process gauge metadata

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
  USE process_ntopo,        ONLY: put_data_struct          ! populate NETOPO and RPARAM data structure
  USE mpi_utils,            ONLY: shr_mpi_initialized      ! If MPI is being used
  USE domain_decomposition, ONLY: mpi_domain_decomposition ! domain decomposition for mpi
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
  USE public_var, ONLY: calendar          ! model calendar
  USE globalData, ONLY: TSEC              ! beginning/ending of simulation time step [sec]
  USE globalData, ONLY: iTime             ! time index at simulation time step
  USE globalData, ONLY: timeVar           ! model time variables in time unit since reference time
  USE public_var, ONLY: time_units        ! netcdf time units - t_unit since yyyy-mm-dd hh:mm:ss
  USE globalData, ONLY: endDatetime       ! model ending datetime
  USE globalData, ONLY: simDatetime       ! current model datetime
  USE write_simoutput_pio, ONLY: close_all

   implicit none
   ! Argument variables
   logical(lgt),              intent(out)   :: finished
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variables
   character(len=7)                         :: t_unit           ! time unit - sec, min, hr, day
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ierr=0; message='update_time/'

   if (simDatetime(1)==endDatetime) then
     call close_all()
     finished=.true.;return
   endif

   ! update model time step bound
   TSEC(0) = TSEC(0) + dt
   TSEC(1) = TSEC(0) + dt

   ! update model time index
   iTime=iTime+1

   ! increment model datetime
   simDatetime(0) = simDatetime(1)
   simDatetime(1) = simDatetime(1)%add_sec(dt, calendar, ierr, cmessage)
   simDatetime(2) = simDatetime(1)%add_sec(dt, calendar, ierr, cmessage)

   ! model time stamp variable for output - dt is in second
   t_unit = trim( time_units(1:index(time_units,' ')) )
   select case( trim(t_unit) )
     case('seconds','second','sec','s'); timeVar = timeVar+ dt
     case('minutes','minute','min');     timeVar = timeVar+ dt/60._dp
     case('hours','hour','hr','h');      timeVar = timeVar+ dt/3600._dp
     case('days','day','d');             timeVar = timeVar+ dt/86400._dp
     case default
       ierr=20; message=trim(message)//'<tunit>= '//trim(t_unit)//': <tunit> must be seconds, minutes, hours or days.'; return
   end select


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
  USE public_var, ONLY: impulseResponseFunc    ! IRF routing ID = 1
  USE public_var, ONLY: kinematicWaveTracking  ! KWT routing ID = 2
  USE public_var, ONLY: kinematicWave          ! KW routing ID = 3
  USE public_var, ONLY: muskingumCunge         ! MC routing ID = 4
  USE public_var, ONLY: diffusiveWave          ! DW routing ID = 5
  USE globalData, ONLY: idxIRF, idxKWT, &
                        idxKW, idxMC, idxDW
  USE globalData, ONLY: nRoutes                ! number of available routing methods
  USE globalData, ONLY: routeMethods           ! ID of active routing method
  USE globalData, ONLY: onRoute                ! logical array for active routing method
  USE globalData, ONLY: masterproc             ! root proc logical
  USE globalData, ONLY: nRch_mainstem          ! number of mainstem reaches
  USE globalData, ONLY: nTribOutlet            ! number of mainstem reaches
  USE globalData, ONLY: rch_per_proc           ! number of tributary reaches
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
  integer(i4b)                     :: iens             ! ensemble index (currently only 1)
  integer(i4b)                     :: ix, iRoute       ! loop indices
  character(len=strLen)            :: cmessage         ! error message of downwind routine

  ierr=0; message='init_state_data/'

  iens = 1_i4b

  do iRoute = 1, nRoutes
    if (routeMethods(iRoute)==kinematicWave) then
      nMolecule%KW_ROUTE = 2
    else if (routeMethods(iRoute)==muskingumCunge) then
      nMolecule%MC_ROUTE = 2
    else if (routeMethods(iRoute)==diffusiveWave) then
      nMolecule%DW_ROUTE = 20
    end if
  end do

  if (isColdStart) then
    if (masterproc) then
      RCHFLX_trib(:,:)%BASIN_QI     = 0._dp
      RCHFLX_trib(:,:)%BASIN_QR(0)  = 0._dp
      RCHFLX_trib(:,:)%BASIN_QR(1)  = 0._dp
      nRch_root=nRch_mainstem+nTribOutlet+rch_per_proc(0)
      if (onRoute(impulseResponseFunc)) then
        do ix = 1,nRch_root
          RCHFLX_trib(iens,ix)%ROUTE(idxIRF)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(iens,ix)%ROUTE(idxIRF)%REACH_Q        = 0._dp
        end do
      end if
      if (onRoute(kinematicWaveTracking)) then
        do ix = 1, nRch_root
          RCHFLX_trib(iens,ix)%ROUTE(idxKWT)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(iens,ix)%ROUTE(idxKWT)%REACH_Q        = 0._dp
        end do
      end if
      if (onRoute(kinematicWave)) then
        do ix = 1, nRch_root
          RCHFLX_trib(iens,ix)%ROUTE(idxKW)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(iens,ix)%ROUTE(idxKW)%REACH_Q        = 0._dp
          allocate(RCHSTA_trib(iens,ix)%KW_ROUTE%molecule%Q(nMolecule%KW_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%KW_ROUTE%molecule%Q]'; return; endif
          RCHSTA_trib(iens,ix)%KW_ROUTE%molecule%Q(:) = 0._dp
        end do
      end if
      if (onRoute(muskingumCunge)) then
        do ix = 1, nRch_root
          RCHFLX_trib(iens,ix)%ROUTE(idxMC)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(iens,ix)%ROUTE(idxMC)%REACH_Q        = 0._dp
          allocate(RCHSTA_trib(iens,ix)%MC_ROUTE%molecule%Q(nMolecule%MC_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%MC_ROUTE%molecule%Q]'; return; endif
          RCHSTA_trib(iens,ix)%MC_ROUTE%molecule%Q(:) = 0._dp
        end do
      end if
      if (onRoute(diffusiveWave)) then
        do ix = 1, nRch_root
          RCHFLX_trib(iens,ix)%ROUTE(idxDW)%REACH_VOL(0:1) = 0._dp
          RCHFLX_trib(iens,ix)%ROUTE(idxDW)%REACH_Q        = 0._dp
          allocate(RCHSTA_trib(iens,ix)%DW_ROUTE%molecule%Q(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%DW_ROUTE%molecule%Q]'; return; endif
          RCHSTA_trib(iens,ix)%DW_ROUTE%molecule%Q(:) = 0._dp
        end do
      end if
    else
      if (rch_per_proc(pid) > 0) then
        RCHFLX_trib(:,:)%BASIN_QI     = 0._dp
        RCHFLX_trib(:,:)%BASIN_QR(0)  = 0._dp
        RCHFLX_trib(:,:)%BASIN_QR(1)  = 0._dp
        if (onRoute(impulseResponseFunc)) then
          do ix = 1, size(RCHFLX_trib(1,:))
            RCHFLX_trib(iens,ix)%ROUTE(idxIRF)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(iens,ix)%ROUTE(idxIRF)%REACH_Q        = 0._dp
          end do
        end if
        if (onRoute(kinematicWaveTracking)) then
          do ix = 1, size(RCHFLX_trib(1,:))
            RCHFLX_trib(iens,ix)%ROUTE(idxKWT)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(iens,ix)%ROUTE(idxKWT)%REACH_Q        = 0._dp
          end do
        end if
        if (onRoute(kinematicWave)) then
          do ix = 1, size(RCHSTA_trib(iens,:))
            RCHFLX_trib(iens,ix)%ROUTE(idxKW)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(iens,ix)%ROUTE(idxKW)%REACH_Q        = 0._dp
            allocate(RCHSTA_trib(iens,ix)%KW_ROUTE%molecule%Q(nMolecule%KW_ROUTE), stat=ierr, errmsg=cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%KW_ROUTE%molecule%Q]'; return; endif
            RCHSTA_trib(iens,ix)%KW_ROUTE%molecule%Q(:) = 0._dp
          end do
        end if
        if (onRoute(muskingumCunge)) then
          do ix = 1, size(RCHSTA_trib(iens,:))
            RCHFLX_trib(iens,ix)%ROUTE(idxMC)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(iens,ix)%ROUTE(idxMC)%REACH_Q        = 0._dp
            allocate(RCHSTA_trib(iens,ix)%MC_ROUTE%molecule%Q(nMolecule%MC_ROUTE), stat=ierr, errmsg=cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%MC_ROUTE%molecule%Q]'; return; endif
            RCHSTA_trib(iens,ix)%MC_ROUTE%molecule%Q(:) = 0._dp
          end do
        end if
        if (onRoute(diffusiveWave)) then
          do ix = 1, size(RCHSTA_trib(iens,:))
            RCHFLX_trib(iens,ix)%ROUTE(idxDW)%REACH_VOL(0:1) = 0._dp
            RCHFLX_trib(iens,ix)%ROUTE(idxDW)%REACH_Q        = 0._dp
            allocate(RCHSTA_trib(iens,ix)%DW_ROUTE%molecule%Q(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage)
            if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA_trib%DW_ROUTE%molecule%Q]'; return; endif
            RCHSTA_trib(iens,ix)%DW_ROUTE%molecule%Q(:) = 0._dp
          end do
        end if
      end if
    end if
    ! initialize time
    TSEC(0)=0._dp; TSEC(1)=dt
  else
    ! start with restart condition
    if (masterproc) then
      call read_state_nc(trim(restart_dir)//trim(fname_state_in), T0, T1, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! time bound [sec] is at previous time step, so need to add dt for curent time step
      TSEC(0)=T0+dt; TSEC(1)=T1+dt
    else ! if other processors, just allocate RCHSTA for dummy
      allocate(RCHSTA(1,1),RCHFLX(1,1), stat=ierr, errmsg=cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHSTA,RCHFLX]'; return; endif
    end if

    call hVars%read_restart(trim(restart_dir)//trim(fname_state_in), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif

    initHvars = .true.

    call mpi_restart(pid, nNodes, comm, iens, ierr, cmessage)
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
  USE shr_sys_mod,          ONLY: shr_sys_system

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
  integer(i4b)   , parameter                  :: maxUpstreamFile=10000000 ! 10 million: maximum number of upstream reaches to enable writing
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

    call shr_sys_system('rm -f '//trim(ancil_dir)//trim(fname_ntopNew), ierr)
    if(ierr/=0)then; message=trim(message)//trim("Error in system call to remove fil"); return; endif

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

END MODULE init_model_data
