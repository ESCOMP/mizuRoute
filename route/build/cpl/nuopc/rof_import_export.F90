module rof_import_export

  use ESMF            , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF            , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF            , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF            , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF            , only : operator(/=), operator(==)
  use NUOPC           , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model     , only : NUOPC_ModelGet
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use rof_shr_methods , only : chkerr

  use public_var      , only : iulog
  use globalData      , only : masterproc          !create this  logical variable  in mizuRoute (masterproc=true => master task, false => other tasks
  use globalData      , only : iTime               ! get step number in mizuRoute
  use RunoffMod       , only : rtmCTL
  use RtmVar          , only : nt_rtm, rtm_tracers

  implicit none
  private ! except

  public  :: advertise_fields
  public  :: realize_fields
  public  :: import_fields
  public  :: export_fields

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_getimport
  private :: state_setexport
  private :: state_getfldptr
  private :: check_for_nans

  type fld_list_type
     character(len=128) :: stdname
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToRof_num = 0
  integer                :: fldsFrRof_num = 0
  type (fld_list_type)   :: fldsToRof(fldsMax)
  type (fld_list_type)   :: fldsFrRof(fldsMax)

  logical     ,parameter :: debug_write = .true. ! internal debug level
  character(*),parameter :: F01 = "('(mizuRoute_import_export) ',a,i5,2x,i8,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: stdname
    character(ESMF_MAXSTR) :: cvalue
    integer                :: dbrc
    integer                :: n, num
    character(len=128)     :: fldname
    character(len=*), parameter :: subname='(rof_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Advertise export fields
    !--------------------------------
    call fldlist_add(fldsFrRof_num, fldsFrRof, trim(flds_scalar_name))
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofi')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_flood')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volr')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volrmch')

    do n = 1,fldsFrRof_num
       call NUOPC_Advertise(exportState, standardName=fldsFrRof(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    call fldlist_add(fldsToRof_num, fldsToRof, trim(flds_scalar_name))

    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsur')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofgwl')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsub')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofi')    ! NOT doing anything for now
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_irrig')   !

    do n = 1,fldsToRof_num
       call NUOPC_Advertise(importState, standardName=fldsToRof(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine advertise_fields

!===============================================================================

  subroutine realize_fields(gcomp, Emesh, flds_scalar_name, flds_scalar_num, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(in)    :: Emesh
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    character(len=*), parameter :: subname='(rof_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrRof, &
         numflds=fldsFrRof_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':mizuRouteExport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToRof, &
         numflds=fldsToRof_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':mizuRouteImport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine realize_fields

!===============================================================================

  subroutine import_fields( gcomp, rc )

    !---------------------------------------------------------------------------
    ! Obtain the runoff input from the mediator and convert from kg/m2s to m3/s
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State) :: importState
    integer          :: n,nt
    integer          :: begr, endr
    integer          :: nliq, nfrz
    integer          :: dbrc
    character(len=*), parameter :: subname='(rof_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set tracers
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') nliq = nt
       if (trim(rtm_tracers(nt)) == 'ICE') nfrz = nt
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(subname),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    begr = rtmCTL%begr
    endr = rtmCTL%endr

    ! determine output array and scale by unit convertsion
    ! NOTE: the call to state_getimport will not convert input unit kg/m2s (mm/s)
    call state_getimport(importState, 'Flrl_rofsur', begr, endr, output=rtmCTL%qsur(:,nliq), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofsub', begr, endr, output=rtmCTL%qsub(:,nliq), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofgwl', begr, endr, output=rtmCTL%qgwl(:,nliq), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofi', begr, endr, output=rtmCTL%qsur(:,nfrz), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_irrig', begr, endr, output=rtmCTL%qirrig(:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    rtmCTL%qsub(begr:endr, nfrz) = 0.0_r8
    rtmCTL%qgwl(begr:endr, nfrz) = 0.0_r8

  end subroutine import_fields

  !====================================================================================

  subroutine export_fields (gcomp, rc)

    !---------------------------------------------------------------------------
    ! Send the runoff model export state to the mediator and convert from m3/s to kg/m2s
    !---------------------------------------------------------------------------

    ! uses
    use RtmVar, only : ice_runoff

    ! input/output/variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State)  :: exportState
    integer           :: n,nt
    integer           :: begr,endr
    integer           :: nliq, nfrz
    real(r8), pointer :: rofl(:)
    real(r8), pointer :: rofi(:)
    real(r8), pointer :: flood(:)
    real(r8), pointer :: volr(:)
    real(r8), pointer :: volrmch(:)
    logical, save     :: first_time = .true.
    integer           :: dbrc
    character(len=*), parameter :: subname='(rof_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set tracers
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') nliq = nt
       if (trim(rtm_tracers(nt)) == 'ICE') nfrz = nt
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(subname),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    if (first_time) then
       if (masterproc) then
          if ( ice_runoff )then
             write(iulog,*)'Snow capping will flow out in frozen river runoff'
          else
             write(iulog,*)'Snow capping will flow out in liquid river runoff'
          endif
          call shr_sys_flush(iulog)
       endif
       first_time = .false.
    end if

    begr = rtmCTL%begr
    endr = rtmCTL%endr

    allocate(rofl(begr:endr))
    allocate(rofi(begr:endr))
    allocate(flood(begr:endr))
    allocate(volr(begr:endr))
    allocate(volrmch(begr:endr))

    if ( ice_runoff )then
       ! separate liquid and ice runoff
       call shr_sys_abort(trim(subname)//': ERROR ice_runoff can NOT be true')
       do n = begr,endr
          rofl(n) = rtmCTL%discharge(n,nliq)
          rofi(n) = rtmCTL%discharge(n,nfrz)
       end do
    else
       ! liquid and ice runoff added to liquid runoff, ice runoff is zero
       do n = begr,endr
          rofl(n) = (rtmCTL%discharge(n,nfrz) + rtmCTL%discharge(n,nliq))
          rofi(n) = 0._r8
       end do
    end if

    ! Flooding back to land, sign convention is positive in land->rof direction
    ! so if water is sent from rof to land, the flux must be negative.

    do n = begr, endr
       !flood(n)   = -rtmCTL%flood(n) / (rtmCTL%area(n)*0.001_r8)
       !flood(n)   = -rtmCTL%flood(n)
       flood(n)   = 0.0_r8
       !volr(n)    =  rtmCTL%volr(n)  / rtmCTL%area(n)
       !volr(n)    =  rtmCTL%volr(n)
       volr(n)    =  0.0_r8
       volrmch(n) =  0.0_r8
    end do

    call state_setexport(exportState, 'Forr_rofl', begr, endr, input=rofl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Forr_rofi', begr, endr, input=rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_flood', begr, endr, input=flood, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_volr', begr, endr, input=volr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_volrmch', begr, endr, input=volrmch, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(rofl, rofi, flood, volr, volrmch)

  end subroutine export_fields

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname

    ! local variables
    integer :: rc
    integer :: dbrc
    character(len=*), parameter :: subname='(rof_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

  end subroutine fldlist_add

  !===============================================================================

  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: dbrc
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(rof_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(rof_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================

  subroutine state_getimport(state, fldname, begr, endr, output, rc, unit_conversion)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    real(r8)            , intent(out)   :: output(begr:endr)
    integer             , intent(out)   :: rc
    real(r8),optional   , intent(in)    :: unit_conversion(begr:endr)

    ! local variables
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    real(r8)                    :: unitConversion(begr:endr)
    type(ESMF_StateItem_Flag)   :: itemFlag
    integer                     :: dbrc
    character(len=*), parameter :: subname='(rof_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), fldptr,  rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! determine runoff flux array and scale by unit convertsion. runoff flux from lnd model is mm/s
       ! mizuRoute expect runoff over area in length/time. length: mm or m, time: second, hour, second
       ! unit of runoff flux is specified in variable units_qsim in src/public_var.f90 and unit conversion
       ! from depth per time to m3/s in mizuRoute.

       unitConversion = 1.0_r8
       if (present(unit_conversion)) then
          unitConversion = unit_conversion
       end if

       do g = begr,endr
          output(g) = fldptr(g-begr+1) * unitConversion(g)
       end do

       ! write debug output if appropriate
       if (masterproc .and. debug_write .and. iTime < 5) then
          do g = begr,endr
             i = 1 + g - begr
             write(iulog,F01)'import: nstep, n, '//trim(fldname)//' = ',iTime, g, output(g)
          end do
          call shr_sys_flush(iulog)
       end if

       ! check for nans
       call check_for_nans(fldptr, trim(fldname), begr)
    end if

  end subroutine state_getimport

  !===============================================================================

  subroutine state_setexport(state, fldname, begr, endr, input, rc)

    use shr_const_mod, only : fillvalue=>SHR_CONST_SPVAL

    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    real(r8)            , intent(in)    :: input(begr:endr)
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    integer                     :: dbrc
    character(len=*), parameter :: subname='(rof_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), fldptr, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       fldptr(:) = 0._r8

       ! set fldptr values to input array
       do g = begr,endr
          fldptr(g-begr+1) = input(g)
       end do

       ! write debug output if appropriate
       if (masterproc .and. debug_write .and. iTime < 5) then
          do g = begr,endr
             i = 1 + g - begr
             write(iulog,F01)'export: nstep, n, '//trim(fldname)//' = ',iTime,i,input(g)
          end do
          call shr_sys_flush(iulog)
       end if

       ! check for nans
       call check_for_nans(fldptr, trim(fldname), begr)
    end if

  end subroutine state_setexport

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------
    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),  intent(in)    :: State
    character(len=*),  intent(in)    :: fldname
    real(R8), pointer, intent(out)   :: fldptr(:)
    integer,           intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: dbrc
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(rof_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif  ! status

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine state_getfldptr

  !===============================================================================

  subroutine check_for_nans(array, fname, begg)

    ! uses
    use shr_infnan_mod, only : isnan => shr_infnan_isnan

    ! input/output variables
    real(r8), pointer             :: array(:)
    character(len=*) , intent(in) :: fname
    integer          , intent(in) :: begg

    ! local variables
    integer :: i
    !-------------------------------------------------------------------------------

    ! Check if any input from mediator or output to mediator is NaN

    if (any(isnan(array))) then
       write(iulog,*) '# of NaNs = ', count(isnan(array))
       write(iulog,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array)
          if (isnan(array(i))) then
             write(iulog,*) "NaN found in field ", trim(fname), ' at gridcell index ',begg+i-1
          end if
       end do
       call shr_sys_abort(' ERROR: One or more of the output from mizuRoute to the coupler are NaN ' )
    end if
  end subroutine check_for_nans

end module rof_import_export
