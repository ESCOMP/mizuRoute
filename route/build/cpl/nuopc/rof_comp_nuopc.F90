module rof_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for mizuRoute
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model           , only : model_routine_SS           => SetServices
  use NUOPC_Model           , only : model_label_Advance        => label_Advance
  use NUOPC_Model           , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model           , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model           , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet
  use shr_kind_mod          , only : R8=>SHR_KIND_R8, CL=>SHR_KIND_CL
  use shr_sys_mod           , only : shr_sys_abort, shr_sys_flush
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_cal_mod           , only : shr_cal_ymd2date

  use public_var            , only : iulog
  use public_var            , only : calendar, simStart, simEnd, time_units
  use globalData            , only : masterproc  !create this  logical variable  in mizuRoute (masterproc=true => master task, false => other tasks
  use globalData            , only : iam        => pid
  use globalData            , only : npes       => nNodes
  use globalData            , only : mpicom_rof => mpicom_route
  use globalData            , only : nHRU
  use init_model_data       , only : get_mpi_omp, init_model
  use RunoffMod             , only : rtmCTL
  use RtmMod                , only : route_ini, route_run
  use RtmTimeManager        , only : init_time, shr_timeStr
  USE RtmVar                , ONLY : cfile_name
  use RtmVar                , only : inst_index, inst_suffix, inst_name, RtmVarSet
  use RtmVar                , only : nsrStartup, nsrContinue, nsrBranch
  use RtmVar                , only : coupling_period !day

  use perf_mod              , only : t_startf, t_stopf, t_barrierf
  use rof_import_export     , only : advertise_fields, realize_fields
  use rof_import_export     , only : import_fields, export_fields
  use rof_shr_methods       , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use rof_shr_methods       , only : set_component_logging, get_component_instance, log_clock_advance

  implicit none
  private ! except

  ! Module routines
  public  :: SetServices
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelSetRunClock
  private :: ModelAdvance
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(len=CL)       :: flds_scalar_name = ''
  integer                 :: flds_scalar_num = 0
  integer                 :: flds_scalar_index_nx = 0
  integer                 :: flds_scalar_index_ny = 0
  integer                 :: flds_scalar_index_nextsw_cday = 0._r8
!#ifdef NDEBUG
  logical, parameter      :: debug_write = .true.
!#else
  !logical, parameter      :: debug_write = .false.
!#endif

  character(*), parameter :: modName =  "(rof_comp_nuopc)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries

    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    integer           :: mpicom
    character(CL)     :: cvalue
    integer           :: shrlogunit
    integer           :: n
    character(len=CL) :: logmsg
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! initialize mizuRoute MPI communicator
    !----------------------------------------------------------------------------

    ! The following call initializees the module variable mpicom_rof in RtmSpmd
    mpicom_rof = mpicom
    call get_mpi_omp(mpicom)

!! WHAT IS THIS? ROFID need for mizuRoute??
!    ! Set ROFID - needed for the mizuRoute code that requires MCT
!    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    read(cvalue,*) ROFID  ! convert from string to integer
!! need for mizuRoute??

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    inst_name = "ROF"//trim(inst_suffix)

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    call set_component_logging(gcomp, masterproc, iulog, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! advertise fields
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxNextSwCday')
    endif

    call advertise_fields(gcomp, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)             :: Emesh
    type(ESMF_DistGrid)         :: DistGrid              ! esmf global index space descriptor
    type(ESMF_Time)             :: currTime              ! Current time
    type(ESMF_Time)             :: startTime             ! Start time
    type(ESMF_Time)             :: stopTime              ! Stop time
    type(ESMF_Time)             :: refTime               ! Ref time
    type(ESMF_TimeInterval)     :: timeStep              ! Model timestep
    type(ESMF_Calendar)         :: esmf_calendar         ! esmf calendar
    type(ESMF_CalKind_Flag)     :: esmf_caltype          ! esmf calendar type
    integer , allocatable       :: gindex(:)             ! global index space on my processor
    integer                     :: yy,mm,dd              ! Temporaries for time query
    integer                     :: ref_ymd               ! reference date (YYYYMMDD)
    integer                     :: ref_tod               ! reference time of day (sec)
    integer                     :: start_ymd             ! start date (YYYYMMDD)
    integer                     :: start_tod             ! start time of day (sec)
    integer                     :: stop_ymd              ! stop date (YYYYMMDD)
    integer                     :: stop_tod              ! stop time of day (sec)
    integer                     :: curr_ymd              ! Start date (YYYYMMDD)
    integer                     :: curr_tod              ! Start time of day (sec)
    logical                     :: flood_present         ! flag
    logical                     :: rof_prognostic        ! flag
    integer                     :: shrlogunit            ! original log unit
    integer                     :: lsize                 ! local size ofarrays
    integer                     :: n,ni                  ! indices
    integer                     :: lbnum                 ! input to memory diagnostic
    integer                     :: nsrest                ! restart type
    integer                     :: ierr                  ! error
    character(len=CL)           :: cmessage              ! error message
    character(len=CL)           :: simRef                ! date string defining the reference time
    character(len=CL)           :: username              ! user name
    character(len=CL)           :: caseid                ! case identifier name
    character(len=CL)           :: ctitle                ! case description title
    character(len=CL)           :: hostname              ! hostname of machine running on
    character(len=CL)           :: model_version         ! model version
    character(len=CL)           :: starttype             ! start-type (startup, continue, branch, hybrid)
    character(len=CL)           :: stdname, shortname    ! needed for advertise
    logical                     :: brnch_retain_casename ! flag if should retain the case name on a branch start type
    character(len=CL)           :: cvalue
    character(ESMF_MAXSTR)      :: convCIM, purpComp
    integer                     :: parametricDim         ! ESMF Mesh parametic dimension
    integer                     :: spatialDim            ! ESMF Mesh spatial dimension
    integer                     :: numOwnedNodes         ! ESMF Mesh number of owned nodes
    integer                     :: numOwnedElements      ! ESMF Mesh number of owned elements
    logical                     :: isMemFreed            ! ESMF Mesh if memory is freed or not
    integer                     :: dimCount              ! ESMF dist-grid dimension count
    integer                     :: tileCount             ! ESMF dist-grid tile count
    integer                     :: deCount               ! ESMF dist-grid processor element count
    integer                     :: localDeCount          ! ESMF dist-grid local processor element count
    logical                     :: regDecompFlag         ! ESMF dist-grid regular decomposition flag
    integer                     :: connectionCount       ! ESMF dist-grid connection count

    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if (masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_comp_nuopc_InitializeRealize:start::',lbnum)
    endif
#endif

    !----------------------
    ! Obtain attribute values
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    read(cvalue,*) caseid
    ctitle=trim(caseid)

    call NUOPC_CompAttributeGet(gcomp, name='brnch_retain_casename', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) brnch_retain_casename

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    call NUOPC_CompAttributeGet(gcomp, name='model_version', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) model_version

    call NUOPC_CompAttributeGet(gcomp, name='hostname', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) hostname

    call NUOPC_CompAttributeGet(gcomp, name='username', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) username

    ! Temporal fix
    call init_model(cfile_name, ierr, cmessage)
    if(ierr/=0) then; cmessage = trim(subname)//trim(cmessage); return; endif

    !----------------------
    ! Initialize time managers in mizuRoute
    !----------------------

    ! Get clock properties
    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get coupling interval in day
    call ESMF_TimeIntervalGet( timeStep, d=coupling_period, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get ymd integer and string for start time, end time and reference time
    ! for continuous run, start time == current time
    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    if (trim(starttype) == trim('continue') ) then
      call shr_timeStr( currTime, simStart )
    else
      call shr_timeStr( startTime, simStart )
    endif

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! adjust stopTime for mizuRoute clock
    stopTime = stopTime - timeStep
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)
    call shr_timeStr( stopTime, simEnd )

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_timeStr( refTime, simRef )

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = 'noleap'
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = 'gregorian'
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    time_units = 'days since '//trim(simRef)

    ! mizuRoute time initialize based on time from coupler
    call init_time(ierr, cmessage)
    if(ierr/=0) then; cmessage = trim(subname)//trim(cmessage); return; endif

    !----------------------
    ! Read namelist, grid and surface data
    !----------------------

    if (masterproc .and. debug_write) then
       write(iulog,*) "mizuRoute initialization"
       write(iulog,*) ' mizuRoute npes = ',npes
       write(iulog,*) ' mizuRoute iam  = ',iam
       write(iulog,*) ' mizuRoute cal  = ',trim(calendar)
       write(iulog,*) ' mizuRoute time = ',trim(time_units)
       write(iulog,*) ' inst_name = ',trim(inst_name)
       call shr_sys_flush(iulog)
    endif

    ! Initialize RtmVar module variables
    ! TODO: the following strings must not be hard-wired - must have module variables
    ! like seq_infodata_start_type_type - maybe another entry in seq_flds_mod?
    if (     trim(starttype) == trim('startup')) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim('continue') ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim('branch')) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( subname//' ERROR: unknown starttype' )
    end if

    call RtmVarSet(&
         caseid_in=caseid, &
         ctitle_in=ctitle,   &
         brnch_retain_casename_in=brnch_retain_casename, &
         nsrest_in=nsrest, &
         version_in=model_version,     &
         hostname_in=hostname, &
         username_in=username)

    !----------------------
    ! Initialize mizuRoute
    !----------------------

    ! - Read in control file
    ! - Initialize time manager
    ! - Initialize number of tracers (ice and liquid) -- NOT ACTIVATED
    ! - Read/process river network input data (global)
    ! - Allocate basins to pes
    ! - Count and distribute HRUs to rglo2gdc (determine rtmCTL%begr, rtmCTL%endr)
    ! - Initialize runoff datatype (rtmCTL)

    call route_ini(rtm_active=rof_prognostic, flood_active=flood_present)

    !--------------------------------
    ! generate the mesh and realize fields
    !--------------------------------

    ! determine global index array
    lsize = rtmCTL%endr - rtmCTL%begr + 1
    allocate(gindex(lsize))
    ni = 0
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       gindex(ni) = rtmCTL%gindex(n)
    end do

    if ( .not. allocated(gindex) )then
       cmessage = cmessage // " gindex is not allocated"
       return
    end if
    if ( debug_write ) then
       write(iulog,*) "iam, lsize = ", iam, lsize
       write(iulog,*) "iam, gindex(min,max,mid) = ", iam, minval(gindex), maxval(gindex), gindex(lsize/2)
       call shr_sys_flush(iulog)
    end if

    ! create distGrid from global index array
    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(gindex)
    if ( .not. ESMF_DistGridIsCreated( DistGrid, rc=rc ) )then
       cmessage = cmessage // " DistGrid is NOT created"
       return
    end if
    call ESMF_DistGridValidate( DistGrid, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get global dist grid attributes
    call ESMF_DistGridGet( DistGrid, dimCount=dimCount, tileCount=tileCount, deCount=deCount, localDeCount=localDeCount, &
     regDecompFlag=regDecompFlag, connectionCount=connectionCount, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (masterproc .and. debug_write) then
       write(iulog,*) "dimCount = ", dimCount
       write(iulog,*) "tileCount = ", tileCount
       write(iulog,*) "deCount = ", deCount
       write(iulog,*) "localDeCount = ", localDeCount
       write(iulog,*) "connectionCount = ", connectionCount
       write(iulog,*) "regDecompFlag = ", regDecompFlag
       call shr_sys_flush(iulog)
    end if

    ! read in the mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_rof', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (masterproc) then
       write(iulog,*)'mesh file for domain is ',trim(cvalue)
       call shr_sys_flush(iulog)
    end if

    ! Create the mesh in one step using the above distGrid
    EMesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, elementDistgrid=Distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshWrite(EMesh, filename='miz_mesh', rc=rc)

!    if ( debug_write )then
!       write(iulog,*) ' After ESMF_MeshCreate'
!       call shr_sys_flush(iulog)
!    end if
!    call ESMF_MeshGet( EMesh,  parametricDim=parametricDim, spatialDim=spatialDim, numOwnedNodes=numOwnedNodes, &
!                       numOwnedElements=numOwnedElements, isMemFreed=isMemFreed, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    if ( debug_write )then
!       write(iulog,*) ' parametricDim = ', parametricDim
!       write(iulog,*) ' spatialDim = ', spatialDim
!       write(iulog,*) ' numOwnedNodes = ', numOwnedNodes
!       write(iulog,*) ' numOwnedElements = ', numOwnedElements
!       write(iulog,*) ' isMemFreed = ', isMemFreed
!       call shr_sys_flush(iulog)
!    end if

    !--------------------------------
    ! realize actively coupled fields
    !--------------------------------

    call realize_fields(gcomp,  Emesh, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Create mizuRoute export state
    !--------------------------------

    call export_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set global grid size scalars in export state
    call State_SetScalar(dble(nHRU), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(1), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Reset shr logging
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug_write)then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "mizuRoute", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "mizuRoute River Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", "mizuRoute River Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2017", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "River", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Name", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", TBD, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
#endif

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_comp_nuopc_InitializeRealize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !------------------------
    ! Run mizuRoute
    !------------------------

    ! arguments:
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables:
    type(ESMF_Clock)  :: clock
    type(ESMF_Alarm)  :: alarm
    type(ESMF_Time)   :: currTime
    type(ESMF_Time)   :: nextTime
    type(ESMF_State)  :: importState
    type(ESMF_State)  :: exportState
    character(CL)     :: cvalue
    integer           :: shrlogunit    ! original log unit
    integer           :: ymd_sync, ymd ! current date (YYYYMMDD)
    integer           :: yr_sync, yr   ! current year
    integer           :: mon_sync, mon ! current month
    integer           :: day_sync, day ! current day
    integer           :: tod_sync, tod ! current time of day (sec)
    logical           :: rstwr         ! .true. ==> write restart file before returning
    logical           :: nlend         ! .true. ==> signaling last time-step
    integer           :: lbnum         ! input to memory diagnostic
    integer           :: g,i           ! indices
    character(len=32) :: rdate         ! date char string for restart file names
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','mizuRoute_comp_nuopc_ModelAdvance:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Unpack import state from mediator
    !--------------------------------

    call t_startf ('lc_mizuRoute_import')

    call import_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_stopf ('lc_mizuRoute_import')

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       rstwr = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       rstwr = .false.
    endif

    !--------------------------------
    ! Determine if time to stop
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       nlend = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       nlend = .false.
    endif

    !--------------------------------
    ! Run mizuRoute
    !--------------------------------

    ! Restart File - use nexttimestr rather than currtimestr here since that is the time at the end of
    ! the timestep and is preferred for restart file names

    call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(nexttime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync, mon_sync, day_sync, tod_sync

    ! Advance mizuRoute time step then run mizuRoute (MODIFIY THIS COMMENT FOR MIZUROUTE: export data is in rtmCTL and Trunoff data types)
    call route_run(rstwr)
    !call advance_timestep()

    !--------------------------------
    ! Pack export state to mediator
    !--------------------------------

    ! (MODIFIY THIS COMMENT FOR MIZUROUTE: input is rtmCTL%runoff, output is r2x)
    call t_startf ('lc_rof_export')

    call export_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_stopf ('lc_rof_export')

    !--------------------------------
    ! Check that internal clock is in sync with master clock
    !--------------------------------
! MIZUROUTE_TODO get ymd tod from mizuRoute
!    call get_curr_date( yr, mon, day, tod)
!    ymd = yr*10000 + mon*100 + day
!    tod = tod
!
!    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) ) then
!       write(iulog,*)' mizuRoute ymd=',ymd     ,'  mizuRoute tod= ',tod
!       write(iulog,*)'      sync ymd=',ymd_sync,'       sync tod= ',tod_sync
!       rc = ESMF_FAILURE
!       call ESMF_LogWrite(subname//" mizuRoute clock not in sync with Master Sync clock",ESMF_LOGMSG_ERROR)
!    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug_write)then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    end if

    if (masterproc) then
       call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing ROF from: ", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       call ESMF_ClockPrint(clock, options="stopTime", preString="--------------------------------> to: ", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    endif

    !--------------------------------
    ! Reset shr logging to my original values
    !--------------------------------

    call shr_file_setLogUnit (shrlogunit)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','mizuRoute_comp_nuopc_ModelAdvance:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(mizuRoute_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(mizuRoute_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (masterproc) then
       write(iulog,F91)
       write(iulog,F00) 'mizuRoute: end of main integration loop'
       write(iulog,F91)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================

end module rof_comp_nuopc
