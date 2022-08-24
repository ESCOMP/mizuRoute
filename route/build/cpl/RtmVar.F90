MODULE RtmVar

  ! Public variables used for only coupled version mizuRoute

  USE shr_kind_mod , ONLY: r8 => shr_kind_r8, CL => SHR_KIND_CL
  USE shr_sys_mod  , ONLY: shr_sys_abort
  USE globalData   , ONLY: masterproc

  implicit none

  private
  public rofVarSet          ! Initialize control variable

  save

  ! private variables
  integer, private, parameter          :: iundef = -9999999
  integer, private, parameter          :: rundef = -9999999._r8
  logical, private                     :: rofVar_isset = .false.

  ! public variables (saved)
  !TODO - nt_rof and rof_tracers need to be removed and set by access to the index array
  integer,           public, parameter :: nt_rof = 2                      ! number of tracers
  character(len=3),  public, parameter :: rof_tracers(nt_rof) = (/'LIQ','ICE'/)

  logical,           public            :: barrier_timers = .false.       ! barrier timers

  character(len=CL), public            :: caseid      = ' '              ! case id
  character(len=CL), public            :: ctitle      = ' '              ! case title
  integer,           public, parameter :: nsrStartup  = 0                ! Startup from initial conditions
  integer,           public, parameter :: nsrContinue = 1                ! Continue from restart files
  integer,           public, parameter :: nsrBranch   = 2                ! Branch from restart files
  integer,           public            :: nsrest = iundef                ! Type of run
  logical,           public            :: brnch_retain_casename = .false.! true => allow case name to remain the same for branch run
                                                                         ! by default this is not allowed
  character(len=CL), public            :: hostname = ' '                 ! Hostname of machine running on
  character(len=CL), public            :: username = ' '                 ! username of user running program
  character(len=CL), public            :: version  = 'v2.0'              ! version of program
  character(len=CL), public            :: conventions = 'CF-1.0'         ! dataset conventions
  character(len=CL), public            :: source   = 'mizuRoute'         ! description of this source
  character(len=CL), public            :: model_doi_url                  ! Web address of the Digital Object Identifier (DOI) for this model version

  ! Instance control
  integer,           public            :: inst_index
  character(len=16), public            :: inst_name
  character(len=16), public            :: inst_suffix

  ! rof control variables
  logical,           public            :: do_rof         = .true.                   !
  logical,           public            :: do_flood       = .false.                  !
  character(len=CL), public            :: nrevsn_rtm     = ' '                      ! restart data file name for branch run
  real(r8),          public            :: river_depth_minimum = 1.e-1               ! minimum river depth for water take [mm]
  integer,           public            :: coupling_period                           ! coupling period
  integer,           public            :: rtmhist_ndens  = 1                        ! namelist: output density of netcdf history files
  integer,           public            :: rtmhist_mfilt  = 30                       ! namelist: number of time samples per tape
  integer,           public            :: rtmhist_nhtfrq = 0                        ! namelist: history write freq(0=monthly)
  logical,           public            :: ice_runoff     = .false.                  ! true => runoff is split into liquid and ice,
  character(len=256),public            :: cfile_name     = 'mizuRoute.control'
  character(len=256),public            :: para_xxxx      = 'mizuRoute_in'

CONTAINS

  SUBROUTINE rofVarSet( caseid_in, ctitle_in, brnch_retain_casename_in,    &
                        nsrest_in, version_in, hostname_in, username_in,   &
                        model_doi_url_in )

    !-----------------------------------------------------------------------
    !  Set input control variables.
    !
    ! !ARGUMENTS:
    character(len=CL), optional, intent(in) :: caseid_in                ! case id
    character(len=CL), optional, intent(in) :: ctitle_in                ! case title
    integer          , optional, intent(in) :: nsrest_in                ! 0: initial run. 1: restart: 3: branch
    character(len=CL), optional, intent(in) :: version_in               ! model version
    character(len=CL), optional, intent(in) :: hostname_in              ! hostname running on
    character(len=CL), optional, intent(in) :: username_in              ! username running job
    character(len=CL), optional, intent(in) :: model_doi_url_in         ! web address of Digital Object Identifier (DOI) for model version
    logical          , optional, intent(in) :: brnch_retain_casename_in ! true => allow case name to
    !-----------------------------------------------------------------------

    if ( rofVar_isset )then
       call shr_sys_abort( 'rofVarSet ERROR:: control variables already set -- EXIT' )
    end if

    if (present(caseid_in)) caseid = caseid_in
    if (present(ctitle_in)) ctitle = ctitle_in
    if (present(nsrest_in)) nsrest = nsrest_in
    if (present(version_in)) version = version_in
    if (present(username_in)) username = username_in
    if (present(hostname_in)) hostname = hostname_in
    if (present(model_doi_url_in)) model_doi_url = model_doi_url_in
    if (present(brnch_retain_casename_in)) brnch_retain_casename = brnch_retain_casename_in

  END SUBROUTINE rofVarSet

!================================================================================

  SUBROUTINE RtmVarInit( )
    if (masterproc) then
       if (nsrest == iundef) then
          call shr_sys_abort( 'RtmVarInit ERROR:: must set nsrest' )
       end if
       if (nsrest == nsrBranch .and. nrevsn_rtm == ' ') then
          call shr_sys_abort( 'RtmVarInit ERROR: need to set restart data file name' )
       end if
       if (nsrest == nsrStartup ) then
          nrevsn_rtm = ' '
       end if
       if (nsrest == nsrContinue) then
          nrevsn_rtm = 'set by restart pointer file file'
       end if
       if (nsrest /= nsrStartup .and. nsrest /= nsrContinue .and. nsrest /= nsrBranch ) then
          call shr_sys_abort( 'RtmVarInit ERROR: nsrest NOT set to a valid value' )
       end if
    endif
    rofVar_isset = .true.
  END SUBROUTINE RtmVarInit

END MODULE RtmVar
