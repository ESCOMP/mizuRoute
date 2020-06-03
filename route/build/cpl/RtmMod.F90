MODULE RtmMod

  !DESCRIPTION:
  USE pio
  USE perf_mod
  USE shr_pio_mod   , ONLY : shr_pio_getiotype, shr_pio_getioformat, &
                             shr_pio_getrearranger, shr_pio_getioroot
  USE shr_kind_mod  , ONLY : r8 => shr_kind_r8, CL => SHR_KIND_CL
  USE shr_sys_mod   , ONLY : shr_sys_flush, shr_sys_abort
  USE RtmVar        , ONLY : nt_rtm, rtm_tracers, &
                             ice_runoff, do_rtm, do_rtmflood, &
                             nsrContinue, nsrBranch, nsrStartup, nsrest, &
                             cfile_name, coupling_period, &
                             caseid, brnch_retain_casename, inst_suffix, inst_name, &
                             barrier_timers
  USE RtmFileUtils,   ONLY : relavu, getavu, opnfil, getfil
  USE RunoffMod     , ONLY : rtmCTL, RunoffInit
  USE public_var,     ONLY : dt                         ! routing time step
  USE public_var    , ONLY : iulog
  USE public_var    , ONLY : rpntfil
  USE globalData    , ONLY : iam        => pid
  USE globalData    , ONLY : npes       => nNodes
  USE globalData    , ONLY : mpicom_rof => mpicom_route
  USE globalData    , ONLY : masterproc
  USE globalData    , ONLY : pio_netcdf_format, pio_typename, pio_rearranger, &
                             pio_root, pio_stride

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

    ! !DESCRIPTION:
    ! Initialize mizuRoute network, decomp
    ! !CALLED FROM:
    ! subroutine initialize in module initializeMod
    USE globalData,      ONLY: ixHRU_order                 ! global HRU index in the order of proc assignment (size = num of hrus contributing reach in entire network)
    USE globalData,      ONLY: hru_per_proc                ! number of hrus assigned to each proc (size = num of procs
    USE globalData,      ONLY: nHRU_mainstem               ! number of mainstem HRUs
    USE init_model_data, ONLY: init_ntopo_data, init_model !
    USE init_model_data, ONLY: init_state_data
    USE mpi_routine,     ONLY: pass_global_data

    !ARGUMENTS:
    implicit none
    logical, intent(out)       :: rtm_active
    logical, intent(out)       :: flood_active
    ! LOCAL VARIABLES:
    character(len=CL)          :: rtm_trstr                 ! tracer string
    integer                    :: ierr                      ! error code
    integer                    :: nt, ix, ix1, ix2
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

!  call init_model(cfile_name, ierr, cmessage)
!  if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    ! If routing time step dt [sec] and coupling time step coupling_period [day] is different, match dt to coupling_period.
    if (dt/=coupling_period) then
      if (masterproc) then
        write(iulog,*) 'WARNING: Adjust mizuRoute dt to coupling_period'
      endif
      dt = coupling_period*60.0*60.0*24.0
    endif

    ! Obtain restart file if appropriate
    if ((nsrest == nsrContinue) .or. &
      (nsrest == nsrBranch  )) then
      call RtmRestGetfile()
    endif

    runtyp(:)               = 'missing'
    runtyp(nsrStartup  + 1) = 'initial'
    runtyp(nsrContinue + 1) = 'restart'
    runtyp(nsrBranch   + 1) = 'branch '

    if (masterproc) then
      write(iulog,*) 'define run:'
      write(iulog,*) '   run type              = ',runtyp(nsrest+1)
      write(iulog,*) '   coupling_period       = ',coupling_period, '[day]'
      write(iulog,*) '   delt_mizuRoute        = ',dt, '[sec]'
      call shr_sys_flush(iulog)
    endif

    rtm_active   = do_rtm
    flood_active = do_rtmflood

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
    ! Overwrite PIO parameter from CIME
    !-------------------------------------------------------
    select case(shr_pio_getioformat(inst_name))
      case(PIO_64BIT_OFFSET); pio_netcdf_format = '64bit_offset'
      case(PIO_64BIT_DATA);   pio_netcdf_format = '64bit_data'
      case default; call shr_sys_abort(trim(subname)//'unexpected netcdf format index')
    end select

    !select case(shr_pio_getiotype(inst_name))
    !  case(pio_iotype_netcdf);   pio_typename = 'netcdf'
    !  case(pio_iotype_pnetcdf);  pio_typename = 'pnetcdf'
    !  case(pio_iotype_netcdf4c); pio_typename = 'netcdf4c'
    !  case(pio_iotype_NETCDF4p); pio_typename = 'netcdf4p'
    !  case default; call shr_sys_abort(trim(subname)//'unexpected netcdf io type index')
    !end select

    !pio_numiotasks    = shr_pio_(inst_name)    ! there is no function to extract pio_numiotasks in cime/src/drivers/nuops/nems/util/shr_pio_mod.F90
     pio_rearranger    = shr_pio_getrearranger(inst_name)
     pio_root          = shr_pio_getioroot(inst_name)
    !pio_stride        = shr_pio_(inst_name)    ! there is no function to extract pio_stride

    write(iulog,*) 'pio_netcdf_format = ', trim(pio_netcdf_format)
    write(iulog,*) 'pio_typename      = ', trim(pio_typename)
    write(iulog,*) 'pio_rearranger    = ', pio_rearranger
    write(iulog,*) 'pio_root          = ', pio_root
    write(iulog,*) 'pio_stride        = ', pio_stride

    !-------------------------------------------------------
    ! Initialize rtm_trstr
    !-------------------------------------------------------

    rtm_trstr = trim(rtm_tracers(1))
    do nt = 2,nt_rtm
       rtm_trstr = trim(rtm_trstr)//':'//trim(rtm_tracers(nt))
    enddo

    if (masterproc) then
       write(iulog,*)'mizuRoute tracers = ',nt_rtm, trim(rtm_trstr)
    end if

    !-------------------------------------------------------
    ! Initialize river network connectivity (global, all procs)
    !-------------------------------------------------------

    call init_ntopo_data(iam, npes, mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    call pass_global_data(mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

    !-------------------------------------------------------
    ! Allocate local flux variables
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
    ! Read restart/initial info
    !-------------------------------------------------------

    call init_state_data(iam, npes, mpicom_rof, ierr, cmessage)
    if(ierr/=0)then; call shr_sys_abort(trim(subname)//trim(cmessage)); endif

  END SUBROUTINE route_ini


! *********************************************************************!
! public subroutine: run
! *********************************************************************!
  SUBROUTINE route_run(rstwr)

    USE dataTypes,           ONLY: strflx
    USE globalData,          ONLY: RCHFLX_trib      ! Reach flux data structures (per proc, tributary)
    USE globalData,          ONLY: RCHFLX_main      ! Reach flux data structures (master proc, mainstem)
    USE globalData,          ONLY: nHRU_mainstem    ! number of mainstem HRUs
    USE globalData,          ONLY: hru_per_proc     ! number of hrus assigned to each proc (i.e., node)
    USE globalData,          ONLY: rch_per_proc     ! number of reaches assigned to each proc (i.e., node)
    USE globalData,          ONLY: basinRunoff_main ! mainstem only HRU runoff
    USE globalData,          ONLY: basinRunoff_trib ! tributary only HRU runoff
    USE globalData,          ONLY: modJulday        ! julian day: current model time step
    USE globalData,          ONLY: restartJulday    ! julian dat: restart dropoff time
    USE write_simoutput_pio, ONLY: prep_output
    USE mpi_routine,         ONLY: mpi_route        ! MPI routing call
    USE write_simoutput_pio, ONLY: output
    USE write_restart_pio,   ONLY: output_state
    USE init_model_data,     ONLY: update_time

! !DESCRIPTION:
! Main routine for routing at all the river reaches per coupling period

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
    real(r8)                     :: delt_coupling          ! real value of coupling_period [sec]
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
    restartJulday = modJulday
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
!         irrig_volume = -rtmCTL%qirrig(nr) * delt_coupling
!
!         ! compare irrig_volume to main channel storage;
!         ! add overage to subsurface runoff
!         if(irrig_volume > RCHFLX_local(1,nr)%REACH_VOL(1)) then
!           rtmCTL%qsub(nr,nt) = rtmCTL%qsub(nr,nt) &
!                              + (RCHFLX_local(1,nr)%REACH_VOL(1) - irrig_volume) / delt_coupling
!           irrig_volume = RCHFLX_local(1,nr)%REACH_VOL
!         endif
!
!         !scs: how to deal with sink points / river outlets?
!         !if (rtmCTL%mask(nr) == 1) then
!
!         ! actual irrigation rate [m3/s]
!         ! i.e. the rate actually removed from the main channel
!         ! if irrig_volume is greater than TRunoff%wr
!         rtmCTL%qirrig_actual(nr) = - irrig_volume / delt_coupling
!
!         ! remove irrigation from wr (main channel)
!         RCHFLX_local(1,nr)%REACH_VOL(1) = RCHFLX_local(1,nr)%REACH_VOL(1) - irrig_volume
!         !scs  endif
!       enddo
!       call t_stopf('mizuRoute_irrig')

    ! Get total runoff for each catchment
    if (npes==1) then

      ! if only single proc is used, all runoff is stored in mainstem runoff array
      if (.not. allocated(basinRunoff_main)) then
        allocate(basinRunoff_main(nHRU_mainstem), stat=ierr)
        if(ierr/=0)then; call shr_sys_abort(trim(subname)//'problem allocating array for [basinRunoff_main]'); endif
      end if
      do nr = 1,nHRU_mainstem
        basinRunoff_main(nr) = rtmCTL%qsur(nr,1)+rtmCTL%qsub(nr,1)!+rtmCTL%qgwl(nr,1)
        write(iulog,'(a,2x,i8,2x,d21.14)')'nr, basinRunoff_main = ',nr, basinRunoff_main(nr)
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

    iens = 1
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

    if (npes==1) then
      allocate(RCHFLX_local(rch_per_proc(-1)), stat=ierr)
      RCHFLX_local(1:rch_per_proc(-1)) = RCHFLX_main(iens,1:rch_per_proc(-1))
    else
      if (masterproc) then
        associate(nRch_main => rch_per_proc(-1), nRch_trib => rch_per_proc(0))
        allocate(RCHFLX_local(nRch_main+nRch_trib), stat=ierr)
        if (nRch_main/=0) then
          RCHFLX_local(1:nRch_main) = RCHFLX_main(iens, 1:nRch_main)
        end if
        if (nRch_trib/=0) then
          RCHFLX_local(nRch_main+1:nRch_main+nRch_trib) = RCHFLX_trib(iens,1:nRch_trib)
        end if
        end associate
      else
        allocate(RCHFLX_local(rch_per_proc(iam)), stat=ierr)
        RCHFLX_local(1:rch_per_proc(iam)) = RCHFLX_trib(iens, 1:rch_per_proc(iam))
      endif
    endif

    do nr = rtmCTL%begr,rtmCTL%endr
      rtmCTL%volr(nr)      = 0._r8
      rtmCTL%flood(nr)     = 0._r8
      rtmCTL%discharge(nr,1) = rtmCTL%discharge(nr,1) + RCHFLX_local(nr)%REACH_Q
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
