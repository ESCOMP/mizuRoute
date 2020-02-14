MODULE RtmMod

! !DESCRIPTION:
  USE shr_kind_mod  , ONLY : r8 => shr_kind_r8, CL => SHR_KIND_CL
  USE shr_sys_mod   , ONLY : shr_sys_flush
  USE RtmVar        , ONLY : nt_rtm, rtm_tracers, &
                             ice_runoff, do_rtm,  &
                             nsrContinue, nsrBranch, nsrStartup, nsrest, &
                             cfile_name, coupling_period, &
                             caseid, brnch_retain_casename, inst_suffix, &
                             barrier_timers
  USE RtmFileUtils,   ONLY : relavu, getavu, opnfil, getfil
  USE RtmTimeManager, ONLY : init_time
  USE RunoffMod     , ONLY : rtmCTL, RunoffInit
  USE public_var,     ONLY : dt                         ! routing time step
  USE public_var    , ONLY : iulog
  USE public_var    , ONLY : rpntfil
  USE globalData    , ONLY : iam        => pid
  USE globalData    , ONLY : npes       => nNodes
  USE globalData    , ONLY : mpicom_rof => mpicom_route
  USE globalData    , ONLY : masterproc
  USE perf_mod

! !PUBLIC TYPES:
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:
  public route_ini          ! Initialize mizuRoute
  public route_run          ! run routing

CONTAINS

! *********************************************************************!
! public subroutine: initialize
! *********************************************************************!
  SUBROUTINE route_ini(rtm_active,flood_active)
!
! !DESCRIPTION:
! Initialize mizuRoute network, decomp
! !CALLED FROM:
! subroutine initialize in module initializeMod
  USE public_var,  ONLY: isRestart
  USE globalData,  ONLY: ixHRU_order                ! global HRU index in the order of proc assignment (size = num of hrus contributing reach in entire network)
  USE globalData,  ONLY: hru_per_proc               ! number of hrus assigned to each proc (size = num of procs
  USE init_model_data, ONLY: init_ntopo_data
  USE init_model_data, ONLY: init_state_data

! !ARGUMENTS:
  implicit none
  logical, intent(out)       :: rtm_active
  logical, intent(out)       :: flood_active
  ! LOCAL VARIABLES:
  character(len=CL)          :: rtm_trstr                 ! tracer string
  integer                    :: ierr                      ! error code
  integer                    :: n, ix1, ix2, myid
  character(len= 7)          :: runtyp(4)                 ! run type
  character(len=CL)          :: cmessage
  character(len=*),parameter :: subname = '(route_ini) '
!-----------------------------------------------------------------------

  !-------------------------------------------------------
  ! mizuRoute setup
  !-------------------------------------------------------
  ! 1. populate meta data
  ! 2. read control file
  ! 3. read routing parameters

  call init_model(cfile_name, ierr, cmessage)
  if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

  ! if routing time step and coupling time step is different, match dt to coupling_period
  if (dt/=coupling_period) then
    write(iulog,*) 'WARNING: Adjust mizuRoute dt to coupling_period'
    dt = coupling_period
  endif

  ! Obtain restart file if appropriate
  if ((nsrest == nsrContinue) .or. &
      (nsrest == nsrBranch  )) then
     isRestart = .true.
     call RtmRestGetfile()
  else if (nsrest == nsrStartup ) then
     isRestart = .false.
  endif

  runtyp(:)               = 'missing'
  runtyp(nsrStartup  + 1) = 'initial'
  runtyp(nsrContinue + 1) = 'restart'
  runtyp(nsrBranch   + 1) = 'branch '

  if (masterproc) then
     write(iulog,*) 'define run:'
     write(iulog,*) '   run type              = ',runtyp(nsrest+1)
     write(iulog,*) '   coupling_period       = ',coupling_period
     write(iulog,*) '   delt_mizuRoute        = ',dt
  endif

    rtm_active = do_rtm

    if ( .not.do_rtm ) then
      if ( masterproc ) then
        write(iulog,*)'mizuRoute will not be active '
      endif
      RETURN
    end if

    if (coupling_period <= 0) then
       write(iulog,*) subname,' ERROR mizuRoute coupling_period invalid',coupling_period
       call shr_sys_abort( subname//' ERROR: coupling_period invalid' )
    endif

    if (dt <= 0) then
       write(iulog,*) subname,' ERROR mizuRoute dt invalid',dt
       call shr_sys_abort( subname//' ERROR: mizuRoute dt invalid' )
    endif

    !-------------------------------------------------------
    ! Initialize mizuRoute time
    !-------------------------------------------------------

    call init_time(ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    !-------------------------------------------------------
    ! Initialize rtm_trstr
    !-------------------------------------------------------

    rtm_trstr = trim(rtm_tracers(1))
    do n = 2,nt_rtm
       rtm_trstr = trim(rtm_trstr)//':'//trim(rtm_tracers(n))
    enddo
    if (masterproc) then
       write(iulog,*)'mizuRoute tracers = ',nt_rtm,trim(rtm_trstr)
    end if

    !-------------------------------------------------------
    ! Determine mosart ocn/land mask (global, all procs)
    !-------------------------------------------------------

    call init_ntopo_data(iam, npes, mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    !-------------------------------------------------------
    ! Allocate local flux variables
    !-------------------------------------------------------

    rtmCTL%begr  = 1
    if (masterproc) then
      rtmCTL%endr  = hru_per_proc(-1)+hru_per_proc(0)
      rtmCTL%lnumr = hru_per_proc(-1)+hru_per_proc(0)
    else
      rtmCTL%lnumr = hru_per_proc(iam)
      rtmCTL%endr  = hru_per_proc(iam)
    end if

    rtmCTL%numr   = sum(hru_per_proc)

    ix2=0
    do myid = 0, npes-1
      ix1 = ix2 + 1
      ix2 = ix1 + rtmCTL%lnumr - 1
      rtmCTL%gindex(rtmCTL%begr:rtmCTL%endr) = ixHRU_order(ix1:ix2)
    enddo

    call RunoffInit(rtmCTL%begr, rtmCTL%endr, rtmCTL%numr)

    !-------------------------------------------------------
    ! Read restart/initial info
    !-------------------------------------------------------

    call init_state_data(iam, npes, mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

  END SUBROUTINE route_ini


! *********************************************************************!
! public subroutine: run
! *********************************************************************!
  SUBROUTINE route_run(rstwr)

    USE dataTypes,         ONLY: strflx
    USE globalData,        ONLY: runoff_data      ! runoff data structure
    USE globalData,        ONLY: RCHFLX_trib      ! Reach flux data structures (per proc, tributary)
    USE globalData,        ONLY: RCHFLX_main      ! Reach flux data structures (master proc, mainstem)
    USE globalData,        ONLY: ixHRU_order      ! global HRU index in the order of proc assignment
    USE globalData,        ONLY: nRch_mainstem    ! number of mainstem reaches
    USE globalData,        ONLY: nHRU_mainstem    ! number of mainstem HRUs
    USE globalData,        ONLY: nContribHRU      ! number of reaches in the whoel river network
    USE globalData,        ONLY: hru_per_proc     ! number of hrus assigned to each proc (i.e., node)
    USE globalData,        ONLY: rch_per_proc     ! number of reaches assigned to each proc (i.e., node)
    USE mpi_routine,       ONLY: mpi_route        ! MPI routing call
    USE nr_utility_module, ONLY: arth
    USE init_model_data,   ONLY: update_time

! !DESCRIPTION:
! River routing model

! !ARGUMENTS:
    implicit none
    logical ,         intent(in) :: rstwr                  ! true => write restart file this step)
! !LOCAL VARIABLES:
    integer                      :: iens                   ! ensemble index (1 for now)
    integer                      :: ix                     ! loop index
    integer                      :: nr, ns, nt             ! indices
    integer                      :: yr, mon, day, ymd, tod ! time information
    integer                      :: nsub                   ! subcyling for cfl
    integer, parameter           :: gather=2
    integer , save               :: nsub_save              ! previous nsub
    real(r8), save               :: delt_save              ! previous delt
    logical,  save               :: first_call = .true.    ! first time flag (for backwards compatibility)
    real(r8)                     :: delt                   ! delt associated with subcycling
    real(r8)                     :: delt_coupling          ! real value of coupling_period
    real(r8),     allocatable    :: runoff_local(:)        !
    type(strflx), allocatable    :: RCHFLX_local(:)
    logical                      :: finished
    ! parameters used in negative runoff partitioning algorithm
    real(r8)                     :: irrig_volume           ! volume of irrigation demand during time step [m3]
    ! error handling
    character(len=*),parameter   :: subname = '(route_run) '
    character(len=CL)            :: cmessage
    integer                      :: ierr                   ! error code
!-----------------------------------------------------------------------

    call t_startf('mizuRoute_tot')
    call shr_sys_flush(iulog)

    delt_coupling = coupling_period*1.0_r8
    if (first_call) then
       delt_save    = dt
       if (masterproc) write(iulog,'(2a,g20.12)') trim(subname),' mizuRoute coupling period ',delt_coupling
    end if

    !-------------------------------------------------------
    ! Initialize mosart history handler and fields
    !-------------------------------------------------------
    call t_startf('mizuRoute_histinit')
    call prep_output(ierr, cmessage)
    call t_stopf('mizuRoute_histinit')

    rtmCTL%discharge = 0._r8
    rtmCTL%volr      = 0._r8
    rtmCTL%flood     = 0._r8

!       !-----------------------------------
!       ! Compute irrigation flux based on demand from clm
!       ! Must be calculated before volr is updated to be consistent with lnd
!       ! Just consider land points and only remove liquid water
!       !-----------------------------------
!
!       call t_startf('mizuRoute_irrig')
!       nt = 1
!       rtmCTL%qirrig_actual = 0._r8
!       do nr = rtmCTL%begr,rtmCTL%endr
!
!         ! calculate volume of irrigation flux during timestep
!         irrig_volume = -rtmCTL%qirrig(nr) * coupling_period
!
!         ! compare irrig_volume to main channel storage;
!         ! add overage to subsurface runoff
!         if(irrig_volume > RCHFLX_local(1,nr)%REACH_VOL(1)) then
!           rtmCTL%qsub(nr,nt) = rtmCTL%qsub(nr,nt) &
!                              + (RCHFLX_local(1,nr)%REACH_VOL(1) - irrig_volume) / coupling_period
!           irrig_volume = RCHFLX_local(1,nr)%REACH_VOL
!         endif
!
!   !scs: how to deal with sink points / river outlets?
!   !     if (rtmCTL%mask(nr) == 1) then
!
!         ! actual irrigation rate [m3/s]
!         ! i.e. the rate actually removed from the main channel
!         ! if irrig_volume is greater than TRunoff%wr
!         rtmCTL%qirrig_actual(nr) = - irrig_volume / coupling_period
!
!         ! remove irrigation from wr (main channel)
!         RCHFLX_local(1,nr)%REACH_VOL(1) = RCHFLX_local(1,nr)%REACH_VOL(1) - irrig_volume
!   !scs  endif
!       enddo
!       call t_stopf('mizuRoute_irrig')

    if (barrier_timers) then
       call t_startf('mizuRoute_SMdirect_barrier')
       call mpi_barrier(mpicom_rof,ierr)
       call t_stopf ('mizuRoute_SMdirect_barrier')
    endif

    ! Get total runoff for river routing
    ! Gather local total runoff (for each proc) into global total-runoff array (runoff_data%basinRunoff) in master
    ! This is temporal routine
    allocate(runoff_local(rtmCTL%lnumr), stat=ierr)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    do nr = rtmCTL%begr,rtmCTL%endr
      runoff_local(nr) = rtmCTL%qsur(nr,1)+rtmCTL%qsub(nr,1)+rtmCTL%qgwl(nr,1)
    end do

    call mpi_comm_single_flux(iam, npes, mpicom_rof,                       & ! input: mpi parameters
                              runoff_data%basinRunoff,                     & ! input: global array
                              runoff_local,                                & ! input: local array
                              rtmCTL%lnumr,                                & ! input: numpber of element for local
                              ixHRU_order(hru_per_proc(-1)+1:nContribHRU), & ! input: global index
                              arth(1,1,hru_per_proc(iam)),                 & ! input: local index
                              gather,                                      & ! input: type of communication
                              ierr, cmessage)
    if (masterproc) then
      do nr = 1,hru_per_proc(-1)
        runoff_data%basinRunoff(ixHRU_order(nr)) = runoff_local(nr)
      end do
    end if

    !-----------------------------------
    ! mizuRoute Subcycling
    !-----------------------------------
    ! subcycling needed when coupling frequency is greater than routing time steps
    call t_startf('mizuRoute_subcycling')

    nsub = coupling_period/dt
    if (nsub*dt < coupling_period) then
       nsub = nsub + 1
    end if
    delt = delt_coupling/float(nsub)
    if (delt /= delt_save) then
       if (masterproc) then
          write(iulog,'(2a,2g20.12,2i12)') trim(subname),' mizuRoute delt update from/to',delt_save,delt,nsub_save,nsub
       end if
    endif

    nsub_save = nsub
    delt_save = delt

    iens = 1
    do ns = 1,nsub

      call mpi_route(iam, npes, mpicom_rof, iens, ierr, cmessage)
      if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

      if (masterproc) then
        associate(nRch_main => rch_per_proc(-1), nRch_trib => rch_per_proc(0))
        allocate(RCHFLX_local(nRch_main+nRch_trib), stat=ierr)
        if (nRch_main/=0) then
          do ix = 1,nRch_main
           RCHFLX_local(ix) = RCHFLX_main(iens, ix)
          enddo
        end if
        RCHFLX_local(nRch_main+1:nRch_main+nRch_trib) = RCHFLX_trib(iens,:)
        end associate
      else
        allocate(RCHFLX_local(rch_per_proc(iam)), stat=ierr)
        RCHFLX_local(:) = RCHFLX_trib(iens, :)
      endif

      do nr = rtmCTL%begr,rtmCTL%endr
        rtmCTL%volr(nr)      = 0._r8
        rtmCTL%flood(nr)     = 0._r8
        do nt = 1,nt_rtm
          rtmCTL%discharge(nr,nt) = rtmCTL%discharge(nr,nt) + RCHFLX_local(nr)%REACH_Q/float(nsub)      ! Reach fluxes (ensembles, space [reaches]) for tributaries
        enddo
      enddo

    enddo
    call t_stopf('mizuRoute_subcycling')

    !-----------------------------------
    ! Write out mizuRoute history file
    !-----------------------------------
    call t_startf('mizuRoute_htapes')
    call output(ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif
    call t_stopf('mizuRoute_htapes')

    !-----------------------------------
    ! Write out mizuRoute restart file
    !-----------------------------------
    if (rstwr) then
      call t_startf('mizuRoute_rest')
      call output_state(ierr, cmessage)
      if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif
      call t_stopf('mizuRoute_rest')
    end if

    ! increment mizuRoute time step
    call update_time(finished, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    !-----------------------------------
    ! Done
    !-----------------------------------

    first_call = .false.

    call shr_sys_flush(iulog)
    call t_stopf('mizuRoute_tot')

  END SUBROUTINE route_run


  SUBROUTINE RtmRestGetfile()
    !---------------------------------------------------
    ! DESCRIPTION:
    ! Determine and obtain netcdf restart file
    USE public_var, ONLY: output_dir
    USE public_var, ONLY: fname_state_in

    ! ARGUMENTS:
    implicit none
    ! LOCAL VARIABLES:
    character(CL)      :: path           ! full pathname of netcdf restart file
    integer            :: status         ! return status
    integer            :: length         ! temporary
    character(len=256) :: ftest,ctest    ! temporaries
    !---------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file
    ! use "output_diro" for directory name for restart file
    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, output_dir, fname_state_in, 0 )
    end if

    ! Branch run:
    ! Restart file pathname is obtained from namelist "fname_state_in"
    if (nsrest==nsrBranch) then

       ! Check case name consistency (case name must be different
       ! for branch run, unless brnch_retain_casename is set)
       ctest = 'xx.'//trim(caseid)//'.mizuRoute'
       ftest = 'xx.'//trim(fname_state_in)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
          write(iulog,*) 'Must change case name on branch run if ',&
               'brnch_retain_casename namelist is not set'
          write(iulog,*) 'previous case filename= ',trim(fname_state_in),&
               ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
               ' ftest = ',trim(ftest)
          call shr_sys_abort()
       end if

    end if

  END SUBROUTINE RtmRestGetfile


  SUBROUTINE restFile_read_pfile( pnamer )

    ! !DESCRIPTION:
    ! Setup restart file and perform necessary consistency checks

    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: pnamer ! full path of restart file

    ! !LOCAL VARIABLES:
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
    !--------------------------------------------------------

    ! Obtain the restart file from the restart pointer file.
    ! For restart runs, the restart pointer file contains the full pathname
    ! of the restart file. For branch runs, the namelist variable
    ! [fname_state_in] contains the restart file name and restart file must be stored in output_dir.
    ! New history files are always created for branch runs.

    if (masterproc) then
       write(iulog,*) 'Reading restart pointer file....'
    endif

    nio = getavu()
    locfn = './'// trim(rpntfil)//trim(inst_suffix)
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    call relavu (nio)

    if (masterproc) then
       write(iulog,*) 'Reading restart data.....'
       write(iulog,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_read_pfile

END MODULE RtmMod
