MODULE RtmMod

  !DESCRIPTION:

  USE pio
  USE perf_mod
  USE shr_pio_mod,  ONLY: shr_pio_getiotype, shr_pio_getioformat, &
                          shr_pio_getrearranger, shr_pio_getioroot, shr_pio_getiosys
  USE shr_kind_mod, ONLY: r8 => shr_kind_r8, CL => SHR_KIND_CL
  USE shr_sys_mod,  ONLY: shr_sys_flush, shr_sys_abort
  USE RtmVar,       ONLY: nt_rof, rof_tracers, &
                          ice_runoff, do_rof, do_flood, &
                          nsrContinue, nsrBranch, nsrStartup, nsrest, &
                          cfile_name, coupling_period, &
                          caseid, brnch_retain_casename, inst_name, &
                          barrier_timers
  USE RunoffMod,    ONLY: rtmCTL, RunoffInit         ! rof related data objects
  USE public_var,   ONLY: dt                         ! routing time step
  USE public_var,   ONLY: iulog
  USE globalData,   ONLY: iam        => pid
  USE globalData,   ONLY: npes       => nNodes
  USE globalData,   ONLY: mpicom_rof => mpicom_route
  USE globalData,   ONLY: masterproc
  USE globalData,   ONLY: multiProcs

  implicit none

  private
  public route_ini          ! Initialize mizuRoute
  public route_run          ! run routing

CONTAINS

  ! *********************************************************************
  ! public subroutine: initialize
  ! *********************************************************************
  SUBROUTINE route_ini(rof_active,flood_active)

    ! DESCRIPTION: Initialize mizuRoute
    ! 1. initialize time
    ! 2. update mizuRoute PIO parameter from CIME
    ! 3. initialize river network topology, and routing parameters
    ! 4. Define Decomposed domain
    ! 5. For continue and/or Branch run-- Initialize restart and history files
    ! 6. initialize state variables

    USE globalData,          ONLY: runMode                     ! "ctsm-coupling" or "standalone" to differentiate some behaviours in mizuRoute
    USE globalData,          ONLY: ixHRU_order                 ! global HRU index in the order of proc assignment (size = num of hrus contributing reach in entire network)
    USE globalData,          ONLY: hru_per_proc                ! number of hrus assigned to each proc (size = num of procs
    USE globalData,          ONLY: nHRU_mainstem               ! number of mainstem HRUs
    USE globalData,          ONLY: pio_netcdf_format
    USE globalData,          ONLY: pio_typename
    USE globalData,          ONLY: pio_rearranger
    USE globalData,          ONLY: pio_root, pio_stride
    USE globalData,          ONLY: pioSystem
    USE init_model_data,     ONLY: init_ntopo_data, init_model !
    USE init_model_data,     ONLY: init_state_data
    USE RtmTimeManager,      ONLY: init_time
    USE mpi_process,         ONLY: pass_global_data
    USE write_simoutput_pio, ONLY: init_histFile

    implicit none
    !ARGUMENTS:
    logical, intent(out)       :: rof_active
    logical, intent(out)       :: flood_active
    ! LOCAL VARIABLES:
    character(len=CL)          :: rof_trstr                 ! tracer string
    integer                    :: ierr                      ! error code
    integer                    :: nt, ix, ix1, ix2
    character(len= 7)          :: runtyp(4)                 ! run type
    character(len=CL)          :: cmessage
    character(len=*),parameter :: subname = '(route_ini) '

    !-------------------------------------------------------
    ! 0. Run types etc.
    !-------------------------------------------------------
    runtyp(:)               = 'missing'
    runtyp(nsrStartup  + 1) = 'initial'
    runtyp(nsrContinue + 1) = 'restart'
    runtyp(nsrBranch   + 1) = 'branch '

    rof_active   = do_rof
    flood_active = do_flood

    if ( .not.do_rof ) then
      if ( masterproc ) then
        write(iulog,*)'mizuRoute will not be active '
      endif
      return
    end if

    ! Initialize tracers
    rof_trstr = trim(rof_tracers(1))
    do nt = 2,nt_rof
      rof_trstr = trim(rof_trstr)//':'//trim(rof_tracers(nt))
    enddo

    if (masterproc) then
      write(iulog,*)'mizuRoute tracers = ',nt_rof, trim(rof_trstr)
    end if

    !-------------------------------------------------------
    ! 1. Initialize time
    !-------------------------------------------------------
    ! If routing time step dt [sec] and coupling time step coupling_period [day] is different, match dt to coupling_period.
    if (dt/=coupling_period) then
      if (masterproc) then
        write(iulog,*) 'WARNING: Adjust mizuRoute dt to coupling_period'
      endif
      dt = coupling_period*60.0*60.0*24.0
    endif

    ! mizuRoute time initialize based on time from coupler
    call init_time(ierr, cmessage)
    if(ierr/=0) then; cmessage = trim(subname)//trim(cmessage); return; endif

    if (masterproc) then
      write(iulog,*) 'define run:'
      write(iulog,*) '   run type              = ',runtyp(nsrest+1)
      write(iulog,*) '   coupling_period       = ',coupling_period, '[day]'
      write(iulog,*) '   delt_mizuRoute        = ',dt, '[sec]'
      call shr_sys_flush(iulog)
    endif

    if (coupling_period <= 0) then
       write(iulog,*) subname,' ERROR mizuRoute coupling_period invalid',coupling_period
       call shr_sys_abort( subname//' ERROR: coupling_period invalid' )
    endif

    if (dt <= 0) then
      write(iulog,*) subname,' ERROR mizuRoute dt invalid',dt
      call shr_sys_abort( subname//' ERROR: mizuRoute dt invalid' )
    endif

    !-------------------------------------------------------
    ! 2. update mizuRoute PIO parameter from CIME
    !-------------------------------------------------------
    runMode='cesm-coupling'

    select case(shr_pio_getioformat(inst_name))
      case(PIO_64BIT_OFFSET); pio_netcdf_format = '64bit_offset'
      case(PIO_64BIT_DATA);   pio_netcdf_format = '64bit_data'
      case default; call shr_sys_abort(trim(subname)//'unexpected netcdf format index')
    end select

    select case(shr_pio_getiotype(inst_name))
      case(pio_iotype_netcdf);   pio_typename = 'netcdf'
      case(pio_iotype_pnetcdf);  pio_typename = 'pnetcdf'
      case(pio_iotype_netcdf4c); pio_typename = 'netcdf4c'
      case(pio_iotype_NETCDF4p); pio_typename = 'netcdf4p'
      case default; call shr_sys_abort(trim(subname)//'unexpected netcdf io type index')
    end select

    !pio_numiotasks    = shr_pio_(inst_name)    ! there is no function to extract pio_numiotasks in cime/src/drivers/nuops/nems/util/shr_pio_mod.F90
    pioSystem         = shr_pio_getiosys(inst_name)
    pio_rearranger    = shr_pio_getrearranger(inst_name)
    pio_root          = shr_pio_getioroot(inst_name)

    if (masterproc) then
      write(iulog,*) 'pio_netcdf_format = ', trim(pio_netcdf_format)
      write(iulog,*) 'pio_typename      = ', trim(pio_typename)
      write(iulog,*) 'pio_rearranger    = ', pio_rearranger
      write(iulog,*) 'pio_root          = ', pio_root
      write(iulog,*) 'pio_stride        = ', pio_stride
    end if

    !-------------------------------------------------------
    ! 3. Initialize river network topology and routing parameters
    !-------------------------------------------------------

    call init_ntopo_data(iam, npes, mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    call pass_global_data(mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    !-------------------------------------------------------
    ! 4. Define Decomposed domain
    !-------------------------------------------------------
    ! 1st and last indices and number of elements in local domain
    rtmCTL%begr  = 1

    ! total elements in whole domain
    rtmCTL%numr   = sum(hru_per_proc)

    if (npes==1) then
      rtmCTL%endr  = rtmCTL%numr
      rtmCTL%lnumr = rtmCTL%numr
      ix1 = 1
      ix2 = rtmCTL%lnumr
    else
      if (masterproc) then
        ! MAster proc uses the first two HRU sections
        rtmCTL%endr  = sum(hru_per_proc(-1:0))
        rtmCTL%lnumr = rtmCTL%endr
        ix1 = 1
        ix2 = rtmCTL%lnumr
      else
        rtmCTL%endr  = hru_per_proc(iam)
        rtmCTL%lnumr = hru_per_proc(iam)
        ix1 = sum(hru_per_proc(-1:iam-1)) + 1
        ix2 = ix1 + hru_per_proc(iam) - 1
      end if
    endif

    call RunoffInit(rtmCTL%begr, rtmCTL%endr, rtmCTL%numr)

    ! index wrt global domain
    rtmCTL%gindex(rtmCTL%begr:rtmCTL%endr) = ixHRU_order(ix1:ix2)

    if ( any(rtmCTL%gindex(rtmCTL%begr:rtmCTL%endr) < 1) )then
      call shr_sys_abort(trim(subname)//"bad gindex < 1")
    endif
    if ( any(rtmCTL%gindex(rtmCTL%begr:rtmCTL%endr) > rtmCTL%numr) )then
      call shr_sys_abort(trim(subname)//"bad gindex > max")
    endif

    !-------------------------------------------------------
    ! 5. For continue and/or Branch run-- Initialize restart and history files
    !-------------------------------------------------------
    ! Obtain last restart name, and obtain and open last history file depending on which run types-- contiuous or branch
    if ((nsrest == nsrContinue) .or. &
        (nsrest == nsrBranch)) then
      call RtmRestGetfile()
      if (nsrest == nsrContinue) then
        call init_histFile(ierr, cmessage)
        if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif
      end if
    endif

    !-------------------------------------------------------
    ! 6. initialize state variables
    !-------------------------------------------------------

    call init_state_data(iam, npes, mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

  END SUBROUTINE route_ini


  ! *********************************************************************!
  ! public subroutine: run
  ! *********************************************************************!
  SUBROUTINE route_run(rstwr)

    ! DESCRIPTION: Main routine for routing at all the river reaches per coupling period

    USE globalData,          ONLY: RCHFLX_trib      ! data structure: Reach flux variables (per proc, tributary)
    USE globalData,          ONLY: NETOPO_trib      ! data structure: River Network topology (other procs, tributary)
    USE globalData,          ONLY: NETOPO_main      ! data structure: River Network topology (main proc, mainstem)
    USE globalData,          ONLY: RPARAM_trib      ! data structure: River parameters (other procs, tributary)
    USE globalData,          ONLY: RPARAM_main      ! data structure: River parameters (main proc, mainstem)
    USE globalData,          ONLY: nHRU_mainstem    ! scalar data: number of mainstem HRUs
    USE globalData,          ONLY: nRch_mainstem    ! scalar data: number of mainstem reaches
    USE globalData,          ONLY: nTribOutlet      ! scalar data: number of tributaries flowing to mainstem
    USE globalData,          ONLY: hru_per_proc     ! array data: number of hrus assigned to each proc (i.e., node)
    USE globalData,          ONLY: rch_per_proc     ! array data: number of reaches assigned to each proc (i.e., node)
    USE globalData,          ONLY: basinRunoff_main ! array data: mainstem only HRU runoff
    USE globalData,          ONLY: basinRunoff_trib ! array data: tributary only HRU runoff
    USE globalData,          ONLY: flux_wm_main     ! array data: mainstem only irrigation demand (water abstract/injection)
    USE globalData,          ONLY: flux_wm_trib     ! array data: tributary only irrigation demand (water abstract/injection)
    USE write_simoutput_pio, ONLY: main_new_file    ! subroutine: create new history files if desired
    USE mpi_process,         ONLY: mpi_route        ! subroutine: MPI routing call
    USE write_simoutput_pio, ONLY: output           ! subroutine: write out history files
    USE write_restart_pio,   ONLY: restart_output   ! subroutine: write out restart file
    USE init_model_data,     ONLY: update_time      ! subroutine: increment time
    USE process_remap_module, ONLY: basin2reach     ! subroutine: mapping variables in hru domain to reach domian
    USE nr_utility_module,    ONLY: arth            ! utility: equivalent to python range function

    implicit none
    ! ARGUMENTS:
    logical ,         intent(in) :: rstwr                  ! true => write restart file this step)
    ! LOCAL VARIABLES:
    integer                      :: iens=1                 ! ensemble index (1 for now)
    integer                      :: ix                     ! loop index
    integer                      :: nr, ns, nt             ! indices
    integer                      :: lwr,upr                ! lower and upper bounds for array slicing
    integer                      :: nRch_trib              ! number of tributary reaches
    integer                      :: yr, mon, day, ymd, tod ! time information
    integer                      :: nsub                   ! subcyling for cfl
    integer , save               :: nsub_save              ! previous nsub
    real(r8), save               :: delt_save              ! previous delt
    logical,  save               :: first_call = .true.    ! first time flag (for backwards compatibility)
    real(r8)                     :: delt                   ! delt associated with subcycling
    real(r8)                     :: delt_coupling          ! real value of coupling_period [sec]
    real(r8), allocatable        :: array_temp(:)          ! temp array
    integer,  allocatable        :: ixRch(:)               ! temp array
    logical                      :: finished
    character(len=*),parameter   :: subname = '(route_run) '
    character(len=CL)            :: cmessage
    integer                      :: ierr                   ! error code

    call t_startf('mizuRoute_tot')
    call shr_sys_flush(iulog)

    delt_coupling = coupling_period*60.0*60.0*24.0
    if (first_call) then
      nsub_save = 1
      delt_save = dt
      if (masterproc) write(iulog,'(2a,g20.12,a)') trim(subname),' mizuRoute coupling period ',delt_coupling, '[sec]'
    end if

    !-------------------------------------------------------
    ! Initialize mizuRoute history handler and fields
    !-------------------------------------------------------
    call t_startf('mizuRoute_histinit')

    call main_new_file(ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    call t_stopf('mizuRoute_histinit')

    !----------------------------------------------------------
    ! Mapping CLM imported fluxes to mizuRoute HRUs or Reaches
    !-----------------------------------------------------------
    ! --- Mapping irrigation demand [mm/s] at HRU to river reach
    ! Must be calculated before volr is updated to be consistent with lnd
    call t_startf('mizuRoute_mapping_irrig')

    allocate(array_temp(rtmCTL%lnumr))
    array_temp = -1._r8*rtmCTL%qirrig
    if (multiProcs) then
      associate(nRch_trib => rch_per_proc(iam))   ! number of tributary reaches
      if (masterproc) then
        if (nRch_mainstem > 0) then ! mainstem
          if (.not. allocated(flux_wm_main)) then
            allocate(flux_wm_main(nRch_mainstem), stat=ierr)
            if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [flux_wm_main]'); endif
          end if
          allocate(ixRch(nRch_mainstem), stat=ierr)
          ixRch = arth(1,1,nRch_mainstem)
          call basin2reach(array_temp(1:nHRU_mainstem), NETOPO_main, RPARAM_main, flux_wm_main, ierr, cmessage, ixSubRch=ixRch)
          if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif
        end if
        if (nRch_trib > 0) then ! tributaries in main processor
          if (.not. allocated(flux_wm_trib)) then
            allocate(flux_wm_trib(nRch_trib), stat=ierr)
            if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [flux_wm_trib]'); endif
          end if
          call basin2reach(array_temp(nHRU_mainstem+1:rtmCTL%lnumr), NETOPO_trib, RPARAM_trib, flux_wm_trib, ierr, cmessage)
          if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif
        end if
      else ! other processors (tributary)
        if (.not. allocated(flux_wm_trib)) then
          allocate(flux_wm_trib(nRch_trib), stat=ierr)
          if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [flux_wm_trib]'); endif
        end if
        call basin2reach(array_temp, NETOPO_trib, RPARAM_trib, flux_wm_trib, ierr, cmessage)
        if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif
      end if
      end associate
    else ! if only single proc is used, all irrigation demand is stored in mainstem array
      if (.not. allocated(flux_wm_main)) then
        allocate(flux_wm_main(nRch_mainstem), stat=ierr)
        if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [flux_wm_main]'); endif
      end if
      call basin2reach(array_temp, NETOPO_main, RPARAM_main, flux_wm_main, ierr, cmessage)
      if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif
    end if

    call t_stopf('mizuRoute_mapping_irrig')

    ! --- Transfer total runoff [mm/s] at HRUs to mizuRoute array
    call t_startf('mizuRoute_mapping_runoff')

    if (npes==1) then
      ! if only single proc is used, all runoff is stored in mainstem runoff array
      if (.not. allocated(basinRunoff_main)) then
        allocate(basinRunoff_main(nHRU_mainstem), stat=ierr)
        if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [basinRunoff_main]'); endif
      end if
      do nr = 1,nHRU_mainstem
        basinRunoff_main(nr) = rtmCTL%qsur(nr,1)+rtmCTL%qsub(nr,1)!+rtmCTL%qgwl(nr,1)
      end do
    else
      if (masterproc) then
        associate(nHRU_trib => hru_per_proc(0))
        if (nHRU_mainstem > 0) then
          if (.not. allocated(basinRunoff_main)) then
            allocate(basinRunoff_main(nHRU_mainstem), stat=ierr)
            if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [basinRunoff_main]'); endif
          end if
          do nr = 1,nHRU_mainstem
            basinRunoff_main(nr) = rtmCTL%qsur(nr,1)+rtmCTL%qsub(nr,1)!+rtmCTL%qgwl(nr,1)
          end do
        end if
        if (.not. allocated(basinRunoff_trib)) then
          allocate(basinRunoff_trib(nHRU_trib), stat=ierr)
          if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [basinRunoff_trib]'); endif
        end if
        do nr = 1, nHRU_trib
          ix = nr + nHRU_mainstem
          basinRunoff_trib(nr) = rtmCTL%qsur(ix,1)+rtmCTL%qsub(ix,1)!+rtmCTL%qgwl(ix,1)
        end do
        end associate
      else
        if (.not. allocated(basinRunoff_trib)) then
          allocate(basinRunoff_trib(rtmCTL%lnumr), stat=ierr)
          if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [basinRunoff_trib]'); endif
        end if
        do nr = rtmCTL%begr,rtmCTL%endr
          basinRunoff_trib(nr) = rtmCTL%qsur(nr,1)+rtmCTL%qsub(nr,1)!+rtmCTL%qgwl(nr,1)
        end do
      end if
    end if

    call t_stopf('mizuRoute_mapping_runoff')

    if (barrier_timers) then
      call t_startf('mizuRoute_SMdirect_barrier')
      call mpi_barrier(mpicom_rof,ierr)
      call t_stopf ('mizuRoute_SMdirect_barrier')
    endif

    !-----------------------------------
    ! mizuRoute Subcycling
    !-----------------------------------
    ! subcycling needed when coupling frequency is greater than routing time steps
    call t_startf('mizuRoute_subcycling')

    nsub = delt_coupling/dt
    if (nsub*dt < delt_coupling) then
      nsub = nsub + 1
    end if
    delt = delt_coupling/float(nsub)
    if (delt /= delt_save) then
      if (masterproc) then
        write(iulog,'(a,2(a,i12,g20.12))') trim(subname),' mizuRoute sub-timestep, dt update from', nsub_save, delt_save, ' to ', nsub, delt
      end if
    endif

    nsub_save = nsub
    delt_save = delt

    do ns = 1,nsub
      if (masterproc) then
        if (allocated(basinRunoff_main)) then
          basinRunoff_main = basinRunoff_main/float(nsub)
        end if
        if (allocated(basinRunoff_trib)) then
          basinRunoff_trib = basinRunoff_trib/float(nsub)
        end if
      else
        basinRunoff_trib = basinRunoff_trib/float(nsub)
      end if

      call mpi_route(iam, npes, mpicom_rof, iens, ierr, cmessage, scatter_ro=.false.)
      if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif
    enddo

    call t_stopf('mizuRoute_subcycling')

    !-----------------------------------
    ! getting ready for exporting
    !-----------------------------------
    call t_startf('mizuRoute_prep_export')

    ! put reach flux variables to associated HRUs
    if (multiProcs) then
      if (masterproc) then
        nRch_trib=rch_per_proc(0)
        if (nRch_mainstem > 0) then
          lwr=1
          upr=nRch_mainstem
          call get_river_export_data(NETOPO_main, RCHFLX_trib(:,lwr:upr))
        end if
        if (nRch_trib > 0) then
          lwr = nRch_mainstem + nTribOutlet + 1
          upr = nRch_mainstem + nTribOutlet + nRch_trib
          call get_river_export_data(NETOPO_trib, RCHFLX_trib(:,lwr:upr))
        end if
      else ! other processors
        call get_river_export_data(NETOPO_trib, RCHFLX_trib)
      end if
    else ! using single processor
      call get_river_export_data(NETOPO_main, RCHFLX_trib(:,1:nRch_mainstem))
    end if

    call t_stopf('mizuRoute_prep_export')

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

      call restart_output(ierr, cmessage)
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

    CONTAINS

      SUBROUTINE get_river_export_data(NETOPO_in, RCHFLX_in)
        USE globalData, ONLY: idxIRF   ! routing method index
        USE dataTypes, ONLY: RCHTOPO   ! data structure - Network topology
        USE dataTypes, ONLY: STRFLX    ! data structure - fluxes in each reach
        implicit none
        ! ARGUMENTS:
        type(RCHTOPO), intent(in) :: NETOPO_in(:)
        type(STRFLX),  intent(in) :: RCHFLX_in(:,:)
        ! LOCAL VARIABLES:
        integer                   :: ix
        integer                   :: iens=1
        integer                   :: iRch, iHru
        integer                   :: nRch, nCatch
        nRch=size(NETOPO_in)
        do iRch =1, nRch
          nCatch = size(NETOPO_in(iRch)%HRUIX)
          do iHru = 1, nCatch
            ix = NETOPO_in(iRch)%HRUIX(iHru)
            rtmCTL%volr(ix)        = RCHFLX_in(iens, iRch)%ROUTE(idxIRF)%REACH_VOL(1)
            rtmCTL%flood(ix)       = 0._r8
            rtmCTL%discharge(ix,1) = RCHFLX_in(iens, iRch)%ROUTE(idxIRF)%REACH_Q
          end do
        end do
      END SUBROUTINE

  END SUBROUTINE route_run


  SUBROUTINE RtmRestGetfile()

    ! DESCRIPTION:
    ! Determine and obtain netcdf restart file
    USE RtmFileUtils, ONLY: getfil
    USE public_var,   ONLY: output_dir
    USE public_var,   ONLY: fname_state_in

    implicit none
    ! ARGUMENTS: None
    ! LOCAL VARIABLES:
    character(CL)      :: path           ! full pathname of netcdf restart file
    integer            :: status         ! return status
    integer            :: length         ! temporary
    character(len=256) :: ftest,ctest    ! temporaries

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
         write(iulog,*) 'Must change case name on branch run if ', &
                        'brnch_retain_casename namelist is not set'
         write(iulog,*) 'previous case filename= ',trim(fname_state_in), &
                        ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
                        ' ftest = ',trim(ftest)
         call shr_sys_abort()
       end if
    end if

  END SUBROUTINE RtmRestGetfile


  SUBROUTINE restFile_read_pfile( pnamer )

    ! DESCRIPTION:
    ! Setup restart file and perform necessary consistency checks
    ! Obtain the restart file from the restart pointer file.
    ! For restart runs, the restart pointer file contains the full pathname
    ! of the restart file. For branch runs, the namelist variable
    ! [fname_state_in] contains the restart file name and restart file must be stored in output_dir.
    ! New history files are always created for branch runs.

    USE RtmFileUtils, ONLY: relavu, getavu, opnfil
    USE RtmVar,       ONLY: inst_suffix
    USE public_var,   ONLY: rpntfil

    implicit none
    ! !ARGUMENTS:
    character(len=*), intent(out) :: pnamer ! full path of restart file
    ! !LOCAL VARIABLES:
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
    !--------------------------------------------------------

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

  END SUBROUTINE restFile_read_pfile

END MODULE RtmMod
