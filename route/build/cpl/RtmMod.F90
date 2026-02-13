MODULE RtmMod

  ! DESCRIPTION:

  USE pio
  USE perf_mod
  USE shr_pio_mod,  ONLY: shr_pio_getiotype, shr_pio_getioformat, &
                          shr_pio_getrearranger, shr_pio_getioroot, shr_pio_getiosys
  USE shr_kind_mod, ONLY: r8 => shr_kind_r8, CL => SHR_KIND_CL
  USE shr_sys_mod,  ONLY: shr_sys_flush, shr_sys_abort
  USE RtmVar,       ONLY: nt_rof, rof_tracers, &
                          ice_runoff, &
                          river_depth_minimum, &
                          nsrContinue, nsrBranch, nsrStartup, nsrest, &
                          cfile_name, coupling_period, nsub, &
                          caseid, brnch_retain_casename, inst_name, &
                          barrier_timers
  USE RunoffMod,    ONLY: ctl, RunoffInit         ! rof related data objects

  ! mizuRoute share routines
  USE public_var,   ONLY: secprday                   ! second per day
  USE public_var,   ONLY: dt                         ! routing time step
  USE public_var,   ONLY: iulog
  USE public_var,   ONLY: qgwl_runoff_option
  USE public_var,   ONLY: bypass_routing_option
  USE globalData,   ONLY: iam        => pid
  USE globalData,   ONLY: npes       => nNodes
  USE globalData,   ONLY: mpicom_rof => mpicom_route
  USE globalData,   ONLY: masterproc
  USE globalData,   ONLY: isColdStart                ! initial river state - cold start (T) or from restart file (F)
  USE mpi_utils,    ONLY: shr_mpi_barrier
  USE mpi_utils,    ONLY: shr_mpi_sparse_distribute

  implicit none
  logical, parameter :: verbose=.false.
  integer, parameter :: iRoute=1        ! index of routing method chosen. cesm-coupled run allows only one option (5 options available), this is always 1.

  private
  public route_ini          ! Initialize mizuRoute
  public route_run          ! run routing for one coupling step

CONTAINS

  ! *********************************************************************
  ! public subroutine: initialize river model
  ! *********************************************************************
  SUBROUTINE route_ini()

    ! DESCRIPTION: Initialize mizuRoute
    ! 1. initialize time
    ! 2. update mizuRoute PIO parameter from CIME
    ! 3. initialize river network topology, and routing parameters
    ! 4. Define Decomposed domain
    ! 5. For continue and/or Branch run-- Initialize restart and history files
    ! 6. initialize state variables

    ! mizuRoute share routines
    USE globalData,          ONLY: ixHRU_order                 ! global HRU index in the order of proc assignment (size = num of hrus contributing reach in entire network)
    USE globalData,          ONLY: hru_per_proc                ! number of hrus assigned to each proc (size = num of procs
    USE globalData,          ONLY: nRch_mainstem               ! scalar data: number of mainstem reaches
    USE globalData,          ONLY: nHRU_mainstem               ! scalar data: number of mainstem hrus
    USE globalData,          ONLY: nHRU_trib                   ! scalar data: number of tributary hrus
    USE globalData,          ONLY: nRch_trib                   ! scalar data: number of tributary reaches
    USE globalData,          ONLY: nTribOutlet                 ! scalar data: number of tributaries flowing to mainstem
    USE globalData,          ONLY: RCHFLX_trib                 ! data structure: Reach flux variables (per proc, tributary)
    USE globalData,          ONLY: NETOPO_trib                 ! data structure: River Network topology (other procs, tributary)
    USE globalData,          ONLY: NETOPO_main                 ! data structure: River Network topology (main proc, mainstem)
    USE globalData,          ONLY: RPARAM_trib                 ! data structure: River parameters (other procs, tributary)
    USE globalData,          ONLY: RPARAM_main                 ! data structure: River parameters (main proc, mainstem)
    USE globalData,          ONLY: pio_rearranger
    USE globalData,          ONLY: pio_root, pio_stride
    USE globalData,          ONLY: pioSystem
    USE public_var,          ONLY: pio_netcdf_format
    USE public_var,          ONLY: pio_typename
    USE init_model_data,     ONLY: init_ntopo_data             !
    USE init_model_data,     ONLY: init_state_data
    USE RtmTimeManager,      ONLY: init_time
    USE mpi_process,         ONLY: pass_global_data
    USE write_simoutput_pio, ONLY: init_histFile

    implicit none
    ! Argument variables:
    ! Local variables:
    character(len=CL)          :: rof_trstr                 ! tracer string
    integer                    :: ierr                      ! error code
    integer                    :: nt, ix, ix1, ix2
    integer                    :: lwr,upr                ! lower and upper bounds for array slicing
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

    ! Initialize tracers
    rof_trstr = trim(rof_tracers(1))
    do nt = 2,nt_rof
      rof_trstr = trim(rof_trstr)//':'//trim(rof_tracers(nt))
    enddo

    if (masterproc) then
      write(iulog,*)'mizuRoute tracers = ',nt_rof, trim(rof_trstr)
    end if

    !-------------------------------------------------------
    ! 1. Get ROF subcycling
    !-------------------------------------------------------
    if (abs(mod(coupling_period, dt)-0._r8)>0._r8) then
      if (masterproc) then
        write(iulog,'(2a,g20.12,a)') trim(subname),' mizuRoute coupling period ', coupling_period, '[sec]'
        write(iulog,'(2a,g20.12,a)') trim(subname),' mizuRoute simulation step <dt_qstim> ', dt, '[sec]'
        write(iulog, '(2A)') trim(subname),'WARNING: multiple of <dt_qsim> [sec] better to be coupling_period [sec]'
      end if
    end if

    nsub = coupling_period/dt
    if (nsub*dt < coupling_period) then
      nsub = nsub + 1
    end if

    if (masterproc) then
      write(iulog,'(2a,g20.12,1a,I3,1a)') trim(subname),' mizuRoute run at', dt, '[sec] step ', nsub, ' per a coupling'
    end if

    !-------------------------------------------------------
    ! 1. Initialize time
    !-------------------------------------------------------
    ! mizuRoute time initialize based on time from coupler
    call init_time(ierr, cmessage)
    if(ierr/=0) then; cmessage = trim(subname)//trim(cmessage); call shr_sys_flush(iulog); call shr_sys_abort( subname//cmessage ); endif

    if (masterproc) then
      write(iulog,*) 'define run:'
      write(iulog,*) '   run type              = ',runtyp(nsrest+1)
      write(iulog,*) '   coupling_frequency    = ',coupling_period, '[day]'
      write(iulog,*) '   mizuRoute timestep    = ',dt, '[sec]'
      call shr_sys_flush(iulog)
    endif

    if (coupling_period <= 0) then
       write(iulog,*) subname,' ERROR mizuRoute coupling_period invalid',coupling_period
       call shr_sys_flush(iulog)
       call shr_sys_abort( subname//' ERROR: coupling_period invalid' )
    endif

    if (dt <= 0) then
      write(iulog,*) subname,' ERROR mizuRoute dt invalid',dt
      call shr_sys_flush(iulog)
      call shr_sys_abort( subname//' ERROR: mizuRoute dt invalid' )
    endif

    !-------------------------------------------------------
    ! 2. update mizuRoute PIO parameter from CIME
    !-------------------------------------------------------
    select case(shr_pio_getioformat(inst_name))
      case(PIO_64BIT_OFFSET); pio_netcdf_format = '64bit_offset'
      case(PIO_64BIT_DATA);   pio_netcdf_format = '64bit_data'
      case default; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//'unexpected netcdf format index')
    end select

    select case(shr_pio_getiotype(inst_name))
      case(pio_iotype_netcdf);   pio_typename = 'netcdf'
      case(pio_iotype_pnetcdf);  pio_typename = 'pnetcdf'
      case(pio_iotype_netcdf4c); pio_typename = 'netcdf4c'
      case(pio_iotype_NETCDF4p); pio_typename = 'netcdf4p'
      case default; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//'unexpected netcdf io type index')
    end select

    !pio_numiotasks    = shr_pio_(inst_name)    ! there is no function to extract pio_numiotasks in cime/src/drivers/nuops/nems/util/shr_pio_mod.F90
    pioSystem         => shr_pio_getiosys(inst_name)
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
    if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    call pass_global_data(mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    !-------------------------------------------------------
    ! 4. Define Decomposed domain
    !-------------------------------------------------------
    ! 1st and last indices and number of elements in local domain
    ctl%begr  = 1

    ! total elements in whole domain
    ctl%numr   = sum(hru_per_proc)

    if (npes==1) then
      ctl%endr  = ctl%numr
      ctl%lnumr = ctl%numr
      ix1 = 1
      ix2 = ctl%lnumr
    else
      if (masterproc) then
        ! MAster proc uses the first two HRU sections
        ctl%endr  = sum(hru_per_proc(-1:0))
        ctl%lnumr = ctl%endr
        ix1 = 1
        ix2 = ctl%lnumr
      else
        ctl%endr  = hru_per_proc(iam)
        ctl%lnumr = hru_per_proc(iam)
        ix1 = sum(hru_per_proc(-1:iam-1)) + 1
        ix2 = ix1 + hru_per_proc(iam) - 1
      end if
    endif

    call RunoffInit(ctl%begr, ctl%endr, ctl%numr)

    ! index wrt global domain
    ctl%gindex(ctl%begr:ctl%endr) = ixHRU_order(ix1:ix2)

    ! additional river reach & catchment information
    if (masterproc) then
      if (nRch_mainstem > 0) then
        call get_hru_area(NETOPO_main, RPARAM_main, verbose=verbose)
      end if
      if (nRch_trib > 0) then
        call get_hru_area(NETOPO_trib, RPARAM_trib, offset=nHRU_mainstem, verbose=verbose)
      end if
    else ! other processors
      call get_hru_area(NETOPO_trib, RPARAM_trib, verbose=verbose)
    end if

    if ( any(ctl%gindex(ctl%begr:ctl%endr) < 1) )then
      call shr_sys_flush(iulog)
      call shr_sys_abort(trim(subname)//"bad gindex < 1")
    endif
    if ( any(ctl%gindex(ctl%begr:ctl%endr) > ctl%numr) )then
      call shr_sys_flush(iulog)
      call shr_sys_abort(trim(subname)//"bad gindex > max")
    endif

    !-------------------------------------------------------
    ! 5. For continue and/or Branch run-- Initialize restart and history files
    !-------------------------------------------------------
    ! Obtain last restart name, and obtain and open last history file depending on which run types-- contiuous or branch
    if (nsrest == nsrContinue) then
      call RtmRestGetfile()
      isColdStart=.false.
      call init_histFile(ierr, cmessage)
      if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif
    endif

    !-------------------------------------------------------
    ! 6. initialize state variables
    !-------------------------------------------------------

    call init_state_data(iam, npes, mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    ! put reach flux variables to associated HRUs
    if (.not. isColdStart) then
      if (masterproc) then
        if (nRch_mainstem > 0) then
          lwr=1
          upr=nRch_mainstem
          call get_river_export_data(NETOPO_main(lwr:upr), RCHFLX_trib(lwr:upr), discharge=.false.)
        end if
        if (nRch_trib > 0) then
          lwr = nRch_mainstem + nTribOutlet + 1
          upr = nRch_mainstem + nTribOutlet + nRch_trib
          call get_river_export_data(NETOPO_trib, RCHFLX_trib(lwr:upr), offset=nHRU_mainstem, discharge=.false.)
        end if
      else ! other processors
        call get_river_export_data(NETOPO_trib, RCHFLX_trib, discharge=.false.)
      end if
    end if

    !-------------------------------------------------------
    ! subroutines used only route_ini
    !-------------------------------------------------------
    CONTAINS

      SUBROUTINE get_hru_area(NETOPO_in, RPARAM_in, offset, verbose)

        ! Descriptions: Get HRU areas in ctl from mizuRoute NETOPO/RPARAM data structure
        ! Note: 1) mizuRoute holds contributory area [m2]-BASAREA, which CAN consists of multiple HRUs
        !       To get HRU area associated with each reach, need to compute based on areal weight
        !       2) offset (optional input): in main processor, index for ctl variable is based on 1 through nHRU_mainstem+nHRU_trib
        !       (combining mainstem and tributary hru in the order). Therefore, hru index based NETOPO for tributary in main processor
        !       need to be added by number of nHRU_mainstem. offset can be used for this.

        USE dataTypes, ONLY: RCHTOPO   ! data structure - Network topology
        USE dataTypes, ONLY: RCHPRP    ! data structure - Reach/hru physical parameters

        implicit none
        ! Arguments:
        type(RCHTOPO),     intent(in) :: NETOPO_in(:) ! network topology data structure
        type(RCHPRP),      intent(in) :: RPARAM_in(:) ! reach parameter data structure
        integer, optional, intent(in) :: offset       ! index offset to point correct location in ctl
        logical, optional, intent(in) :: verbose      ! verbose option
        ! Local variables:
        logical                   :: verb
        integer                   :: ix
        integer                   :: iRch, iHru
        integer                   :: nRch, nCatch
        real(r8)                  :: area_hru

        if (present(verbose)) then
          verb=verbose
        else
          verb=.false.
        end if

        nRch=size(NETOPO_in)
        do iRch =1, nRch
          nCatch = size(NETOPO_in(iRch)%HRUIX)
          do iHru = 1, nCatch
            ix = NETOPO_in(iRch)%HRUIX(iHru)
            if (present(offset)) ix = ix+offset
            ! BASAREA = total area [m2] of multiple HRUs contributing to a reach.
            ! To get HRU area [m2], mulitply HRUWGT (areal weight (HRU area/total contributory area)
            ctl%area(ix) = RPARAM_in(iRch)%BASAREA* NETOPO_in(iRch)%HRUWGT(iHru)
          end do
        end do

        if (verb) then
          do iRch =1, nRch
            nCatch = size(NETOPO_in(iRch)%HRUIX)
            do iHru = 1, nCatch
              ix = NETOPO_in(iRch)%HRUIX(iHru)
              if (present(offset)) ix = ix+offset
              write(iulog, '(a,1x,5(g20.12))') &
                    'reachID, hruID, basinArea [m2], weight[-], hruArea [m2]=', &
                    NETOPO_in(iRch)%REACHID, NETOPO_in(iRch)%HRUID(iHru), RPARAM_in(iRch)%BASAREA, &
                    NETOPO_in(iRch)%HRUWGT(iHru), ctl%area(ix)
            end do
          end do
        end if

      END SUBROUTINE get_hru_area

  END SUBROUTINE route_ini

  ! *********************************************************************!
  ! public subroutine: perform routing at one coupling step
  ! *********************************************************************!
  SUBROUTINE route_run(rstwr)

    ! DESCRIPTION: Main routine for putting irrigation, surface and subsurface runoff in HRUs and river reaches,
    !                               routing at all the river reaches per coupling period
    !                               preparing for exporting the discharge, river volume to the coupler
    ! NOTE: HRUs: Hydrologic Response Units, AKA catchments
    !       Coupler pass imported variables to HRUs first, then passed to corresponding river reaches

    USE globalData,          ONLY: RCHFLX_trib      ! data structure: Reach flux variables (per proc, tributary)
    USE globalData,          ONLY: NETOPO_trib      ! data structure: River Network topology (other procs, tributary)
    USE globalData,          ONLY: NETOPO_main      ! data structure: River Network topology (main proc, mainstem)
    USE globalData,          ONLY: RPARAM_trib      ! data structure: River parameters (other procs, tributary)
    USE globalData,          ONLY: RPARAM_main      ! data structure: River parameters (main proc, mainstem)
    USE globalData,          ONLY: nHRU_mainstem    ! scalar data: number of mainstem HRUs
    USE globalData,          ONLY: nRch_mainstem    ! scalar data: number of mainstem reaches
    USE globalData,          ONLY: nHRU_trib        ! scalar data: number of tributary HRUs
    USE globalData,          ONLY: nRch_trib        ! scalar data: number of tributary reaches
    USE globalData,          ONLY: nTribOutlet      ! scalar data: number of tributaries flowing to mainstem
    USE globalData,          ONLY: hru_per_proc     ! array data: number of hrus assigned to each proc (i.e., node)
    USE globalData,          ONLY: rch_per_proc     ! array data: number of reaches assigned to each proc (i.e., node)
    USE globalData,          ONLY: basinRunoff_main ! array data: mainstem only HRU runoff
    USE globalData,          ONLY: basinRunoff_trib ! array data: tributary only HRU runoff
    USE globalData,          ONLY: flux_wm_main     ! array data: mainstem only irrigation demand (water abstract/injection)
    USE globalData,          ONLY: flux_wm_trib     ! array data: tributary only irrigation demand (water abstract/injection)
    USE globalData,          ONLY: commRch          ! deriv.data: send-recieve river segment pair
    USE mpi_utils,           ONLY: shr_mpi_send     ! subroutine: point-to-point communication
    USE mpi_process,         ONLY: mpi_route        ! subroutine: MPI routing call
    USE write_simoutput_pio, ONLY: main_new_file    ! subroutine: create new history files if desired
    USE write_simoutput_pio, ONLY: output           ! subroutine: write out history files
    USE write_restart_pio,   ONLY: restart_output   ! subroutine: write out restart file
    USE init_model_data,     ONLY: update_time      ! subroutine: increment time
    USE process_remap_module,ONLY: basin2reach      ! subroutine: mapping variables in hru domain to reach domian
    USE nr_utils,            ONLY: arth             ! utility: equivalent to python range function

    implicit none
    ! Arguments:
    logical ,         intent(in) :: rstwr                  ! true => write restart file this step)
    ! Local variables:
    integer                      :: ix, nr, ns, nt         ! loop indices
    integer                      :: lwr,upr                ! lower and upper bounds for array slicing
    real(r8)                     :: qgwl_depth             ! depth of qgwl runoff during time step [mm]
    real(r8)                     :: irrig_depth            ! depth of irrigation demand during time step [mm]
    real(r8)                     :: river_depth            ! depth of river water during time step [mm]
    real(r8), allocatable        :: qSend(:)               ! array holding negative lateral flow to be sent to outlet
    logical                      :: finished               ! dummy arguments (not really used)
    character(len=CL)            :: cmessage               ! error message from subroutines
    integer                      :: ierr                   ! error code
    character(len=12),parameter  :: subname = '(route_run) '

    call t_startf('mizuRoute_tot')
    call shr_sys_flush(iulog)

    !-------------------------------------------------------
    ! Initialize mizuRoute history handler and fields
    !-------------------------------------------------------
    call t_startf('mizuRoute_histinit')

    call main_new_file(ierr, cmessage)
    if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    call t_stopf('mizuRoute_histinit')

    !----------------------------------------------------------
    ! Mapping CLM imported fluxes to mizuRoute HRUs or Reaches
    !-----------------------------------------------------------
    ! --- Mapping irrigation demand [mm/s] at HRU to river reach (qirrig is negative)
    ! qirrig is negative
    ! Must be calculated before volr is updated to be consistent with lnd
    call t_startf('mizuRoute_mapping_irrig')

    ctl%qirrig_actual = -1._r8* ctl%qirrig  ! actual water take [mm/s] - positive (take) or negative (inject)
    do nr = ctl%begr,ctl%endr
      ! calculate depth of irrigation [mm] during timestep
      irrig_depth = ctl%qirrig_actual(nr)* coupling_period
      river_depth = ctl%volr(nr)* 1000._r8 ! m to mm

      ! compare irrig_depth [mm] to previous channel storage [mm];
      ! add overage to subsurface runoff
      ! later check negative qsub is handle the same as qgwl
      if(irrig_depth > river_depth) then
        ctl%qsub(nr,1) = ctl%qsub(nr,1) + (river_depth-irrig_depth)/coupling_period
        irrig_depth = river_depth

        ! actual irrigation rate [mm/s]
        ! i.e. the rate actually removed from the river channel
        ctl%qirrig_actual(nr) = irrig_depth/coupling_period
      endif
    end do

    call t_stopf('mizuRoute_mapping_irrig')

    ! --- handling negative flow - qgwl and qsub after irrigation take
    !  1. direct_in_place:  place negative flow at that reach aside and coupler handle the flow
    !  2. direct_to_outlet: send negative flow to an outlet of the reach and subtract the flow from the outlet reach
    call t_startf('mizuRoute_bypass_route')

    select case(trim(bypass_routing_option))
      case('direct_in_place')
        ctl%direct = 0._r8
        do nr = ctl%begr,ctl%endr
          ! --- Transfer qgwl [mm/s] to ocean
          if (trim(qgwl_runoff_option) == 'all') then ! send all qgwl flow to ocean
            ctl%direct(nr,1) = ctl%qgwl(nr,1)
            ctl%qgwl(nr,1) = 0._r8
          else if (trim(qgwl_runoff_option) == 'negative') then ! send only negative qgwl flow to ocean
            if(ctl%qgwl(nr,1) < 0._r8) then
              ctl%direct(nr,1) = ctl%qgwl(nr,1)
              ctl%qgwl(nr,1) = 0._r8
            endif
          else if (trim(qgwl_runoff_option) == 'threshold') then
            ! if qgwl is negative, and adding it to the main channel
            ! would bring main channel storage below a threshold,
            ! send qgwl directly to ocean

            ! --- calculate depth of qgwl flux [mm] during timestep
            qgwl_depth = ctl%qgwl(nr,1)* coupling_period
            river_depth = ctl%volr(nr)* 1000._r8 ! convert m to mm

            if ((qgwl_depth + river_depth < river_depth_minimum) &
                 .and. (ctl%qgwl(nr,1) < 0._r8)) then
              ctl%direct(nr,1) = ctl%qgwl(nr,1)
              ctl%qgwl(nr,1) = 0._r8
            end if
          end if
          ! --- Transfer qsub to ocean [mm/s]
          if(ctl%qsub(nr,1) < 0._r8) then
            ctl%direct(nr,1) = ctl%direct(nr,1)+ ctl%qsub(nr,1)
            ctl%qsub(nr,1) = 0._r8
          endif
        end do
      case('direct_to_outlet')
        allocate(qSend(ctl%lnumr))
        qSend(:) = 0._r8     ! total negative q [mm/s] to be sent to the outlet (converted to +)
        do nr = ctl%begr,ctl%endr
          if(ctl%qgwl(nr,1) < 0._r8) then
            qSend(nr) = ctl%qgwl(nr,1)
            ctl%qgwl(nr,1) = 0._r8
          end if
          if(ctl%qsub(nr,1) < 0._r8) then
            qSend(nr) = qSend(nr) + ctl%qsub(nr,1)
            ctl%qsub(nr,1) = 0._r8
          end if
        end do

        ! Distribute "direct runoff to ocean" to targe reach (i.e., outlet of river network)
        call shr_mpi_sparse_distribute(qSend, commRch(:)%destTask, commRch(:)%destIndex, ctl%direct(:,1), fillvalue=0._r8)

        call shr_mpi_barrier(mpicom_rof, cmessage)
        if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif

      case default; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//'unexpected bypass_routing_option')
    end select

    call t_stopf('mizuRoute_bypass_route')

    ! --- Transfer total runoff [mm/s] at HRUs to mizuRoute array
    call t_startf('mizuRoute_mapping_runoff')

    ! Transfer actual irrigation rate [mm/s] to river segment
    if (masterproc) then
      if (nRch_mainstem > 0) then ! mainstem
        call basin2reach(ctl%qirrig_actual(1:nHRU_mainstem), NETOPO_main, RPARAM_main, flux_wm_main, &
                         ierr, cmessage, limitRunoff=.false.)
        if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif
      end if
      if (nRch_trib > 0) then ! tributaries in main processor
        call basin2reach(ctl%qirrig_actual(nHRU_mainstem+1:ctl%lnumr), NETOPO_trib, RPARAM_trib, flux_wm_trib, &
                         ierr, cmessage, limitRunoff=.false.)
        if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif
      end if
    else ! other processors (tributary)
      call basin2reach(ctl%qirrig_actual, NETOPO_trib, RPARAM_trib, flux_wm_trib, ierr, cmessage, limitRunoff=.false.)
      if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif
    end if

    if (masterproc) then
      if (nHRU_mainstem > 0) then
        do nr = 1,nHRU_mainstem
          basinRunoff_main(nr) = ctl%qsur(nr,1)+ctl%qsub(nr,1)+ctl%qgwl(nr,1)
        end do
      end if
      do nr = 1, nHRU_trib
        ix = nr + nHRU_mainstem
        basinRunoff_trib(nr) = ctl%qsur(ix,1)+ctl%qsub(ix,1)+ctl%qgwl(ix,1)
      end do
    else
      do nr = ctl%begr,ctl%endr
        basinRunoff_trib(nr) = ctl%qsur(nr,1)+ctl%qsub(nr,1)+ctl%qgwl(nr,1)
      end do
    end if

    call t_stopf('mizuRoute_mapping_runoff')

    if (barrier_timers) then
      call t_startf('mizuRoute_SMdirect_barrier')
      call shr_mpi_barrier(mpicom_rof, cmessage)
      if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//"mpi_barrier error"); endif
      call t_stopf ('mizuRoute_SMdirect_barrier')
    endif

    !-----------------------------------
    ! run mizuRoute
    !-----------------------------------
    ! subcycling needed when coupling frequency is greater than routing time steps
    do ns = 1,nsub
      call t_startf('mizuRoute_subcycling')

      call mpi_route(iam, npes, mpicom_rof, ierr, cmessage, scatter_ro=.false.)
      if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif

      call t_stopf('mizuRoute_subcycling')

      !-----------------------------------
      ! Write out mizuRoute history file
      !-----------------------------------
      call t_startf('mizuRoute_htapes')

      call output(ierr, cmessage)
      if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif

      call t_stopf('mizuRoute_htapes')

      !-----------------------------------
      ! Write out mizuRoute restart file
      !-----------------------------------
      if (rstwr) then
        call t_startf('mizuRoute_rest')

        call restart_output(ierr, cmessage)
        if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif

        call t_stopf('mizuRoute_rest')
      end if

      ! increment mizuRoute time step
      call update_time(finished, ierr, cmessage)
      if(ierr/=0)then; call shr_sys_flush(iulog); call shr_sys_abort(trim(subname)//trim(cmessage)); endif
    end do

    !-----------------------------------
    ! getting ready for exporting
    !-----------------------------------
    call t_startf('mizuRoute_prep_export')

    ! put reach flux variables to associated HRUs
    if (masterproc) then
      if (nRch_mainstem > 0) then
        lwr=1
        upr=nRch_mainstem
        call get_river_export_data(NETOPO_main(lwr:upr), RCHFLX_trib(lwr:upr))
      end if
      if (nRch_trib > 0) then
        lwr = nRch_mainstem + nTribOutlet + 1
        upr = nRch_mainstem + nTribOutlet + nRch_trib
        call get_river_export_data(NETOPO_trib, RCHFLX_trib(lwr:upr), offset=nHRU_mainstem)
      end if
    else ! other processors
      call get_river_export_data(NETOPO_trib, RCHFLX_trib)
    end if

    call t_stopf('mizuRoute_prep_export')

    !-----------------------------------
    ! Done
    !-----------------------------------
    call shr_sys_flush(iulog)
    call t_stopf('mizuRoute_tot')

  END SUBROUTINE route_run

  !-------------------------------------------------------
  ! subroutines used only route_run
  !-------------------------------------------------------
  SUBROUTINE get_river_export_data(NETOPO_in, RCHFLX_in, &
                                   volume, discharge, flood, &
                                   offset)

    ! Descriptions: get exporting variables
    !  - discharge [kg/m2/s]=[mm/s]
    !  - water volume per HUC area [m]
    !  - flood [mm/s]
    !  NOTE 1) offset (optional input): in main processor, index for ctl variable is based on 1 through nHRU_mainstem+nHRU_trib
    !       (combining mainstem and tributary hru in the order). Therefore, hru index based NETOPO for tributary in main processor
    !       need to be added by number of nHRU_mainstem. offset can be used for this.

    USE dataTypes, ONLY: RCHTOPO   ! data structure - Network topology
    USE dataTypes, ONLY: STRFLX    ! data structure - fluxes in each reach

    implicit none
    ! Arguments:
    type(RCHTOPO), intent(in)     :: NETOPO_in(:)   ! mizuRoute network topology data structure
    type(STRFLX),  intent(in)     :: RCHFLX_in(:)   ! mizuRoute reach flux data structure
    integer, optional, intent(in) :: offset         ! index offset to point correct location in ctl
    logical, optional, intent(in) :: volume         !
    logical, optional, intent(in) :: discharge      !
    logical, optional, intent(in) :: flood          !
    ! Local variables:
    integer                   :: ix
    integer                   :: iRch, iHru
    integer                   :: nRch, nCatch
    logical                   :: update_vol      !
    logical                   :: update_q        !
    logical                   :: update_fld      !

    update_vol = .true.
    update_q   = .true.
    update_fld = .true.
    if (present(volume)) then
      update_vol=volume
    end if
    if (present(discharge)) then
      update_q=discharge
    end if
    if (present(flood)) then
      update_fld=flood
    end if

    nRch=size(NETOPO_in)
    do iRch =1, nRch
      nCatch = size(NETOPO_in(iRch)%HRUIX)
      do iHru = 1, nCatch
        ix = NETOPO_in(iRch)%HRUIX(iHru)
        if (present(offset)) ix=ix+offset
        ! stream volume is split into HRUs based on HRU's areal weight for contributory area
        if (update_vol) then
          ctl%volr(ix) = RCHFLX_in(iRch)%ROUTE(iRoute)%REACH_VOL(1)*NETOPO_in(iRch)%HRUWGT(iHru)/ctl%area(ix)
        end if
        if (update_q) then
          if (NETOPO_in(iRch)%DREACHI==-1 .and. NETOPO_in(iRch)%DREACHK<=0) then ! if reach is the outlet
            ctl%discharge(ix,1) = RCHFLX_in(iRch)%ROUTE(iRoute)%REACH_Q* NETOPO_in(iRch)%HRUWGT(iHru)/ctl%area(ix)/0.001_r8
          end if
        end if
        if (update_fld) then
          ctl%flood(ix) = RCHFLX_in(iRch)%ROUTE(iRoute)%FLOOD_VOL(1)*NETOPO_in(iRch)%HRUWGT(iHru)/ctl%area(ix)/coupling_period
        end if
      end do
    end do
  END SUBROUTINE

  SUBROUTINE RtmRestGetfile()

    ! DESCRIPTION:
    ! Determine and obtain netcdf restart file

    USE RtmFileUtils, ONLY: getfil
    USE public_var,   ONLY: output_dir
    USE public_var,   ONLY: fname_state_in

    implicit none
    ! Arguments: None
    ! Local variables:
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
       ctest = 'xx.'//trim(caseid)//'.mizuroute'
       ftest = 'xx.'//trim(fname_state_in)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
         write(iulog,*) 'Must change case name on branch run if ', &
                        'brnch_retain_casename namelist is not set'
         write(iulog,*) 'previous case filename= ',trim(fname_state_in), &
                        ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
                        ' ftest = ',trim(ftest)
         call shr_sys_flush(iulog)
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

    USE RtmFileUtils, ONLY: opnfil
    USE RtmVar,       ONLY: inst_suffix
    USE public_var,   ONLY: rpntfil
    USE public_var,   ONLY: secprmin, secprhour
    USE globalData,   ONLY: simDatetime

    implicit none
    ! Arguments:
    character(len=*), intent(out) :: pnamer ! full path of restart file
    ! Local variables:
    integer :: i                    ! indices
    integer :: nio                  ! restart unit
    integer :: sec_in_day           ! time in second
    character(len=17) :: timestamp  ! datetime string in yyyy-mm-dd-sssss
    integer :: status               ! substring check status
    logical :: lexist               ! If local file exists
    character(len=256) :: locfn     ! Restart pointer file name
    !--------------------------------------------------------

    if (masterproc) then
      write(iulog,*) 'Reading restart pointer file....'
    endif

    ! construct rpointer file name with datetime - simDatetime has three datetime at previous(0), current(1) and next(2) time stamp
    sec_in_day = simDatetime(1)%hour()*nint(secprhour)+simDatetime(1)%minute()*nint(secprmin)+nint(simDatetime(1)%sec())
    write(timestamp,'(".",I4.4,"-",I2.2,"-",I2.2,"-",I5.5)') &
          simDatetime(1)%year(),simDatetime(1)%month(),simDatetime(1)%day(),sec_in_day
    locfn = './'// trim(rpntfil)//trim(inst_suffix)//timestamp
    inquire (file=locfn,exist=lexist)
    if (.not. lexist) then ! backward compatibility - rpointer file w/o datetime
      write(iulog,*) 'Could not find the rpointer file: ', trim(locfn)
      locfn = './'// trim(rpntfil)//trim(inst_suffix)
      write(iulog,*) 'Try file: ', trim(locfn)
    end if
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    close(nio)

    if (masterproc) then
      write(iulog,*) 'Reading restart data: ', trim(pnamer)
      write(iulog,'(72a1)') ("-",i=1,60)
    end if

  END SUBROUTINE restFile_read_pfile

END MODULE RtmMod
