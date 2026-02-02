MODULE RtmVar

  ! Public variables used for only coupled version mizuRoute

  USE shr_kind_mod , ONLY: r8 => shr_kind_r8, CL => SHR_KIND_CL
  USE shr_sys_mod  , ONLY: shr_sys_abort, shr_sys_flush
  USE globalData   , ONLY: masterproc
  USE globalData   , ONLY: version
  USE public_var   , ONLY: iulog, integerMissing

  implicit none

  private

  save

  ! public variables (saved)
  !TODO - nt_rof and rof_tracers need to be removed and set by access to the index array
  integer,           public, parameter :: nt_rof = 2                      ! number of tracers
  character(len=3),  public, parameter :: rof_tracers(nt_rof) = (/'LIQ','ICE'/)

  logical,           public            :: barrier_timers = .false.       ! barrier timers

  ! Metadata variables used in history and restart files's global attributes
  character(len=CL), public            :: caseid      = ' '              ! case id
  character(len=CL), public            :: ctitle      = ' '              ! case title
  character(len=CL), public            :: hostname = ' '                 ! Hostname of machine running on
  character(len=CL), public            :: username = ' '                 ! username of user running program
  character(len=CL), public            :: conventions = 'CF-1.0'         ! dataset conventions
  character(len=CL), public            :: model_doi_url                  ! Web address of the Digital Object Identifier (DOI) for this model version
  character(len=CL), public            :: source   = 'mizuRoute'         ! description of this source

  ! Run startup
  integer,           public, parameter :: nsrStartup  = 0                ! Startup from initial conditions
  integer,           public, parameter :: nsrContinue = 1                ! Continue from restart files
  integer,           public, parameter :: nsrBranch   = 2                ! Branch from restart files
  integer,           public            :: nsrest = integerMissing        ! Type of run
  logical,           public            :: brnch_retain_casename = .false.! true => allow case name to remain the same for branch run

  ! Instance control
  integer,           public            :: inst_index
  character(len=16), public            :: inst_name
  character(len=16), public            :: inst_suffix

  ! rof control variables
  character(len=CL), public            :: nrevsn_rtm     = ' '                      ! restart data file name for branch run
  real(r8),          public            :: river_depth_minimum = 1.e-1               ! minimum river depth for water take [mm]
  integer,           public            :: nsub                                      ! number of subcyling for rof simulation per coupling
  real(r8),          public            :: coupling_period                           ! coupling period [sec]
  integer,           public            :: rtmhist_ndens  = 1                        ! namelist: output density of netcdf history files
  integer,           public            :: rtmhist_mfilt  = 30                       ! namelist: number of time samples per tape
  integer,           public            :: rtmhist_nhtfrq = 0                        ! namelist: history write freq(0=monthly)
  logical,           public            :: ice_runoff     = .false.                  ! true => runoff is split into liquid and ice,
  character(len=256),public            :: cfile_name     = 'mizuRoute.control'
  character(len=256),public            :: para_xxxx      = 'mizuRoute_in'

END MODULE RtmVar
