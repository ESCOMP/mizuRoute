! Unit tests verifying read_control_toml and read_control_legacy parity across all 226 keys
program test_read_control

  USE nrtype
  USE public_var
  USE globalData
  USE var_lookup
  USE dataTypes,         ONLY: var_info
  USE objTypes,          ONLY: var_info_new
  USE popMetadat_module, ONLY: popMetadat
  USE read_control_module, ONLY: read_control, read_control_toml, read_control_legacy

  implicit none

  type(var_info) :: default_meta_SEG(nVarsSEG)
  type(var_info) :: default_meta_HRU(nVarsHRU)
  type(var_info) :: default_meta_HRU2SEG(nVarsHRU2SEG)
  type(var_info) :: default_meta_NTOPO(nVarsNTOPO)
  type(var_info) :: default_meta_PFAF(nVarsPFAF)
  integer(i4b) :: init_err
  character(len=strLen) :: init_msg

  call popMetadat(init_err, init_msg)
  call assert_true(init_err == 0, 'popMetadat failed at startup')
  default_meta_SEG     = meta_SEG
  default_meta_HRU     = meta_HRU
  default_meta_HRU2SEG = meta_HRU2SEG
  default_meta_NTOPO   = meta_NTOPO
  default_meta_PFAF    = meta_PFAF

  call test_parity_all_226_keys()
  call test_output_frequency_types()
  call test_default_preservation()
  call test_error_handling()

  write(*, "(a)") "ALL read_control unit tests PASSED."

contains

  ! Reset all public_var and globalData control variables to their initial defaults
  subroutine reset_defaults()
    ancil_dir             = charMissing
    input_dir             = charMissing
    output_dir            = charMissing
    restart_dir           = charMissing
    case_name             = charMissing
    continue_run          = .false.
    simStart              = charMissing
    simEnd                = charMissing
    newFileFrequency      = 'yearly'
    outputFrequency       = '1'
    outputNameOption      = 'specific'
    nOutFreq              = integerMissing
    routOpt               = '1'
    doesBasinRoute        = 1
    floodplain            = .false.
    hw_drain_point        = 2
    tracer                = .false.
    is_lake_sim           = .false.
    lakeRegulate          = .true.
    LakeInputOption       = 0
    is_flux_wm            = .false.
    is_vol_wm             = .false.
    is_vol_wm_jumpstart   = .false.
    scale_factor_runoff   = realMissing
    offset_value_runoff   = realMissing
    scale_factor_Ep       = realMissing
    offset_value_Ep       = realMissing
    is_Ep_upward_negative = .false.
    scale_factor_prec     = realMissing
    offset_value_prec     = realMissing
    min_length_route      = 0.0_dp
    checkMassBalance      = .false.
    dt                    = realMissing
    fname_ntopOld         = charMissing
    ntopAugmentMode       = .false.
    fname_ntopNew         = charMissing
    dname_sseg            = 'seg'
    dname_nhru            = 'hru'
    fname_qsim            = charMissing
    vname_qsim            = 'runoff'
    vname_evapo           = 'evap'
    vname_precip          = 'precip'
    vname_solute          = charMissing
    vname_time            = 'time'
    vname_hruid           = 'hru'
    dname_time            = 'time'
    dname_hruid           = 'hru'
    dname_xlon            = 'lon'
    dname_ylat            = 'lat'
    units_qsim            = charMissing
    units_cc              = charMissing
    dt_ro                 = realMissing
    input_fillvalue       = realMissing
    ro_time_units         = charMissing
    ro_calendar           = charMissing
    ro_time_stamp         = 'start'
    runoffMin             = 0._dp
    fname_wm              = charMissing
    vname_flux_wm         = charMissing
    vname_vol_wm          = charMissing
    vname_time_wm         = 'time'
    vname_segid_wm        = 'seg'
    dname_time_wm         = 'time'
    dname_segid_wm        = 'seg'
    dt_wm                 = realMissing
    is_remap              = .false.
    fname_remap           = charMissing
    vname_hruid_in_remap  = 'polyid'
    vname_weight          = 'weight'
    vname_qhruid          = 'overlapPolyId'
    vname_num_qhru        = 'overlaps'
    vname_i_index         = 'i_index'
    vname_j_index         = 'j_index'
    dname_hru_remap       = 'polyid'
    dname_data_remap      = 'data'
    restart_write         = 'never'
    restart_date          = charMissing
    restart_month         = 1
    restart_day           = 1
    restart_hour          = 0
    fname_state_in        = charMissing
    param_nml             = ''
    gageMetaFile          = charMissing
    outputAtGage          = .false.
    fname_gageObs         = charMissing
    vname_gageFlow        = 'flow'
    vname_gageSite        = 'site'
    vname_gageTime        = 'time'
    dname_gageSite        = 'site'
    dname_gageTime        = 'time'
    strlen_gageSite       = 30
    pio_netcdf_format     = '64bit_offset'
    pio_typename          = 'pnetcdf'
    histTimeStamp_offset  = 0._dp
    outputInflow          = .false.
    qmodOption            = 0
    QerrTrend             = 1
    qBlendPeriod          = 10
    hydGeometryOption     = readFromFile
    topoNetworkOption     = compute
    computeReachList      = compute
    time_units            = charMissing
    calendar              = charMissing
    debug                 = .false.
    idSegOut              = integerMissing
    desireId              = integerMissing
    maxPfafLen            = 32
    pfafMissing           = '0'
    bypass_routing_option = 'direct_in_place'
    qgwl_runoff_option    = 'threshold'
    correct_area          = .false.
    ice_runoff            = .true.

    runMode               = 'standalone'
    isColdStart           = .true.
    masterproc            = .false.
    nRoutes               = 0
    if (allocated(routeMethods)) deallocate(routeMethods)
    onRoute               = .false.
    idxSUM                = integerMissing
    idxIRF                = integerMissing
    idxKWT                = integerMissing
    idxKW                 = integerMissing
    idxMC                 = integerMissing
    idxDW                 = integerMissing
    time_conv             = realMissing
    length_conv           = realMissing
    time_conv_solute      = realMissing
    mass_conv_solute      = realMissing

    meta_SEG     = default_meta_SEG
    meta_HRU     = default_meta_HRU
    meta_HRU2SEG = default_meta_HRU2SEG
    meta_NTOPO   = default_meta_NTOPO
    meta_PFAF    = default_meta_PFAF
  end subroutine reset_defaults

  ! Assert a boolean condition or stop with a non-zero exit code
  subroutine assert_true(cond, msg)
    logical(lgt), intent(in) :: cond
    character(*), intent(in) :: msg
    if (.not. cond) then
      write(*, '(a,a)') 'ASSERTION FAILED: ', trim(msg)
      error stop 1
    endif
  end subroutine assert_true

  ! Assert two strings are equal after trimming
  subroutine assert_char_eq(actual, expected, label)
    character(*), intent(in) :: actual, expected, label
    if (trim(actual) /= trim(expected)) then
      write(*, '(a,a,a,a,a,a)') 'CHAR MISMATCH [', trim(label), ']: actual="', &
                                trim(actual), '" expected="', trim(expected) // '"'
      error stop 1
    endif
  end subroutine assert_char_eq

  ! Assert two 32-bit integers are equal
  subroutine assert_int_eq(actual, expected, label)
    integer(i4b), intent(in) :: actual, expected
    character(*), intent(in) :: label
    if (actual /= expected) then
      write(*, '(a,a,a,i0,a,i0)') 'INT MISMATCH [', trim(label), ']: actual=', &
                                  actual, ' expected=', expected
      error stop 1
    endif
  end subroutine assert_int_eq

  ! Assert two double-precision reals are equal within tolerance
  subroutine assert_dp_eq(actual, expected, label)
    real(dp),     intent(in) :: actual, expected
    character(*), intent(in) :: label
    real(dp), parameter      :: tol = 1.0e-12_dp ! tight relative/absolute tolerance for double precision comparisons
    if (abs(actual - expected) > tol) then
      write(*, '(a,a,a,g0,a,g0)') 'DP MISMATCH [', trim(label), ']: actual=', &
                                  actual, ' expected=', expected
      error stop 1
    endif
  end subroutine assert_dp_eq

  ! Assert two logicals are equal
  subroutine assert_bool_eq(actual, expected, label)
    logical(lgt), intent(in) :: actual, expected
    character(*), intent(in) :: label
    if (actual .neqv. expected) then
      write(*, '(a,a,a,l1,a,l1)') 'BOOL MISMATCH [', trim(label), ']: actual=', &
                                  actual, ' expected=', expected
      error stop 1
    endif
  end subroutine assert_bool_eq

  ! Write one key-value pair to both the legacy .control unit and the .toml unit
  subroutine write_kv(u_ctl, u_toml, key, ctl_val, toml_val)
    integer(i4b), intent(in) :: u_ctl, u_toml
    character(*), intent(in) :: key, ctl_val, toml_val
    write(u_ctl,  '(a,a,a,a,a)') '<', trim(key), '>  ', trim(ctl_val), '  ! comment'
    write(u_toml, '(a,a,a)')     trim(key), ' = ', trim(toml_val)
  end subroutine write_kv

  ! Write a string key-value pair to both the legacy .control unit and the .toml unit
  subroutine write_kv_str(u_ctl, u_toml, key, val)
    integer(i4b), intent(in) :: u_ctl, u_toml
    character(*), intent(in) :: key, val
    call write_kv(u_ctl, u_toml, key, val, '"' // trim(val) // '"')
  end subroutine write_kv_str

  ! Test 1: Verify full 226-key parity between read_control_legacy and read_control_toml
  subroutine test_parity_all_226_keys()
    integer(i4b) :: u_ctl, u_toml, err, i
    character(len=strLen) :: message
    character(len=*), parameter :: ctl_file  = 'test_all_226.control'
    character(len=*), parameter :: toml_file = 'test_all_226.toml'

    ! Saved state from read_control_legacy
    character(len=strLen) :: s_ancil_dir, s_input_dir, s_output_dir, s_restart_dir
    character(len=strLen) :: s_case_name, s_simStart, s_simEnd, s_newFileFrequency
    character(len=strLen) :: s_outputFrequency, s_outputNameOption, s_routOpt
    character(len=strLen) :: s_fname_ntopOld, s_fname_ntopNew, s_dname_sseg, s_dname_nhru
    character(len=strLen) :: s_fname_qsim, s_vname_qsim, s_vname_evapo, s_vname_precip
    character(len=strLen) :: s_vname_solute, s_vname_time, s_vname_hruid
    character(len=strLen) :: s_dname_time, s_dname_hruid, s_dname_xlon, s_dname_ylat
    character(len=strLen) :: s_units_qsim, s_units_cc, s_ro_time_units, s_ro_calendar, s_ro_time_stamp
    character(len=strLen) :: s_fname_wm, s_vname_flux_wm, s_vname_vol_wm, s_vname_time_wm, s_vname_segid_wm
    character(len=strLen) :: s_dname_time_wm, s_dname_segid_wm, s_fname_remap
    character(len=strLen) :: s_vname_hruid_in_remap, s_vname_weight, s_vname_qhruid, s_vname_num_qhru
    character(len=strLen) :: s_vname_i_index, s_vname_j_index, s_dname_hru_remap, s_dname_data_remap
    character(len=strLen) :: s_restart_write, s_restart_date, s_fname_state_in, s_param_nml
    character(len=strLen) :: s_gageMetaFile, s_fname_gageObs, s_vname_gageFlow, s_vname_gageSite, s_vname_gageTime
    character(len=strLen) :: s_dname_gageSite, s_dname_gageTime, s_pio_netcdf_format, s_pio_typename
    character(len=strLen) :: s_time_units, s_pfafMissing, s_bypass_routing_option, s_qgwl_runoff_option
    integer(i4b)          :: s_nOutFreq, s_doesBasinRoute, s_hw_drain_point, s_LakeInputOption
    integer(i4b)          :: s_restart_month, s_restart_day, s_restart_hour, s_strlen_gageSite
    integer(i4b)          :: s_qmodOption, s_QerrTrend, s_qBlendPeriod, s_hydGeometryOption
    integer(i4b)          :: s_topoNetworkOption, s_computeReachList, s_idSegOut, s_desireId, s_maxPfafLen
    real(dp)              :: s_scale_factor_runoff, s_offset_value_runoff, s_scale_factor_Ep, s_offset_value_Ep
    real(dp)              :: s_scale_factor_prec, s_offset_value_prec, s_min_length_route, s_dt
    real(dp)              :: s_dt_ro, s_input_fillvalue, s_runoffMin, s_dt_wm, s_histTimeStamp_offset
    real(dp)              :: s_time_conv, s_length_conv, s_time_conv_solute, s_mass_conv_solute
    logical(lgt)          :: s_continue_run, s_floodplain, s_tracer, s_is_lake_sim, s_lakeRegulate
    logical(lgt)          :: s_is_flux_wm, s_is_vol_wm, s_is_vol_wm_jumpstart, s_is_Ep_upward_negative
    logical(lgt)          :: s_checkMassBalance, s_ntopAugmentMode, s_is_remap, s_outputAtGage
    logical(lgt)          :: s_outputInflow, s_debug, s_correct_area, s_ice_runoff, s_isColdStart
    integer(i4b)          :: s_nRoutes, s_idxSUM, s_idxIRF, s_idxKWT, s_idxKW, s_idxMC, s_idxDW
    integer(i4b), allocatable :: s_routeMethods(:)
    logical(lgt)          :: s_onRoute(0:nRouteMethods-1)
    type(var_info)        :: s_meta_HRU(nVarsHRU), s_meta_HRU2SEG(nVarsHRU2SEG)
    type(var_info)        :: s_meta_SEG(nVarsSEG), s_meta_NTOPO(nVarsNTOPO), s_meta_PFAF(nVarsPFAF)
    type(var_info_new)    :: s_meta_rflx(nVarsRFLX), s_meta_hflx(nVarsHFLX)

    ! Set up temporary control and TOML files containing all 226 parameters for parity verification
    open(newunit=u_ctl,  file=ctl_file,  status='replace', action='write')
    open(newunit=u_toml, file=toml_file, status='replace', action='write')

    ! 1-8: Directories & Run control
    call write_kv_str(u_ctl, u_toml, 'ancil_dir', '/path/ancil/')
    call write_kv_str(u_ctl, u_toml, 'input_dir', '/path/input/')
    call write_kv_str(u_ctl, u_toml, 'output_dir', '/path/output/')
    call write_kv_str(u_ctl, u_toml, 'restart_dir', '/path/restart/')
    call write_kv_str(u_ctl, u_toml, 'case_name', 'parity_case')
    call write_kv_str(u_ctl, u_toml, 'sim_start', '2000-01-01 00:00:00')
    call write_kv_str(u_ctl, u_toml, 'sim_end', '2000-12-31 00:00:00')
    call write_kv    (u_ctl, u_toml, 'continue_run', 'T', 'true')
    ! 9-28: Routing & physics options (route_opt as integer 123450 activates all 6 routing schemes 0..5)
    call write_kv    (u_ctl, u_toml, 'route_opt', '123450', '123450')
    call write_kv    (u_ctl, u_toml, 'doesBasinRoute', '1', '1')
    call write_kv    (u_ctl, u_toml, 'dt_qsim', '3600.0', '3600.0')
    call write_kv    (u_ctl, u_toml, 'floodplain', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'hw_drain_point', '1', '1')
    call write_kv    (u_ctl, u_toml, 'tracer', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'is_lake_sim', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'lakeRegulate', 'F', 'false')
    call write_kv    (u_ctl, u_toml, 'LakeInputOption', '2', '2')
    call write_kv    (u_ctl, u_toml, 'is_flux_wm', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'is_vol_wm', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'is_vol_wm_jumpstart', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'scale_factor_runoff', '1.25', '1.25')
    call write_kv    (u_ctl, u_toml, 'offset_value_runoff', '0.05', '0.05')
    call write_kv    (u_ctl, u_toml, 'scale_factor_Ep', '0.95', '0.95')
    call write_kv    (u_ctl, u_toml, 'offset_value_Ep', '-0.01', '-0.01')
    call write_kv    (u_ctl, u_toml, 'is_Ep_upward_negative', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'scale_factor_prec', '1.1', '1.1')
    call write_kv    (u_ctl, u_toml, 'offset_value_prec', '0.02', '0.02')
    call write_kv    (u_ctl, u_toml, 'min_length_route', '50.5', '50.5')
    ! 29-33: River network topology
    call write_kv_str(u_ctl, u_toml, 'fname_ntopOld', 'ntop_old.nc')
    call write_kv    (u_ctl, u_toml, 'ntopAugmentMode', 'T', 'true')
    call write_kv_str(u_ctl, u_toml, 'fname_ntopNew', 'ntop_new.nc')
    call write_kv_str(u_ctl, u_toml, 'dname_nhru', 'nhru_custom')
    call write_kv_str(u_ctl, u_toml, 'dname_sseg', 'sseg_custom')
    ! 34-52: Runoff, evaporation, precipitation
    call write_kv_str(u_ctl, u_toml, 'fname_qsim', 'qsim_custom.nc')
    call write_kv_str(u_ctl, u_toml, 'vname_qsim', 'qsim_var')
    call write_kv_str(u_ctl, u_toml, 'vname_evapo', 'evapo_var')
    call write_kv_str(u_ctl, u_toml, 'vname_precip', 'precip_var')
    call write_kv_str(u_ctl, u_toml, 'vname_solute', 'solute_var')
    call write_kv_str(u_ctl, u_toml, 'vname_time', 'time_var')
    call write_kv_str(u_ctl, u_toml, 'vname_hruid', 'hruid_var')
    call write_kv_str(u_ctl, u_toml, 'dname_time', 'time_dim')
    call write_kv_str(u_ctl, u_toml, 'dname_hruid', 'hruid_dim')
    call write_kv_str(u_ctl, u_toml, 'dname_xlon', 'xlon_dim')
    call write_kv_str(u_ctl, u_toml, 'dname_ylat', 'ylat_dim')
    call write_kv_str(u_ctl, u_toml, 'units_qsim', 'mm/h')
    call write_kv_str(u_ctl, u_toml, 'units_cc', 'g/h')
    call write_kv    (u_ctl, u_toml, 'dt_ro', '1800.0', '1800.0')
    call write_kv    (u_ctl, u_toml, 'input_fillvalue', '-8888.0', '-8888.0')
    call write_kv_str(u_ctl, u_toml, 'ro_calendar', 'noleap')
    call write_kv_str(u_ctl, u_toml, 'ro_time_units', 'hours since 2000-01-01')
    call write_kv_str(u_ctl, u_toml, 'ro_time_stamp', 'middle')
    call write_kv    (u_ctl, u_toml, 'runoffMin', '1.0e-9', '1.0e-9')
    ! 53-60: Water management
    call write_kv_str(u_ctl, u_toml, 'fname_wm', 'wm_custom.nc')
    call write_kv_str(u_ctl, u_toml, 'vname_flux_wm', 'flux_wm_var')
    call write_kv_str(u_ctl, u_toml, 'vname_vol_wm', 'vol_wm_var')
    call write_kv_str(u_ctl, u_toml, 'vname_time_wm', 'time_wm_var')
    call write_kv_str(u_ctl, u_toml, 'vname_segid_wm', 'segid_wm_var')
    call write_kv_str(u_ctl, u_toml, 'dname_time_wm', 'time_wm_dim')
    call write_kv_str(u_ctl, u_toml, 'dname_segid_wm', 'segid_wm_dim')
    call write_kv    (u_ctl, u_toml, 'dt_wm', '7200.0', '7200.0')
    ! 61-70: Remapping
    call write_kv    (u_ctl, u_toml, 'is_remap', 'T', 'true')
    call write_kv_str(u_ctl, u_toml, 'fname_remap', 'remap_custom.nc')
    call write_kv_str(u_ctl, u_toml, 'vname_hruid_in_remap', 'hruid_remap_var')
    call write_kv_str(u_ctl, u_toml, 'vname_weight', 'weight_var')
    call write_kv_str(u_ctl, u_toml, 'vname_qhruid', 'qhruid_var')
    call write_kv_str(u_ctl, u_toml, 'vname_num_qhru', 'num_qhru_var')
    call write_kv_str(u_ctl, u_toml, 'vname_i_index', 'i_idx_var')
    call write_kv_str(u_ctl, u_toml, 'vname_j_index', 'j_idx_var')
    call write_kv_str(u_ctl, u_toml, 'dname_hru_remap', 'hru_remap_dim')
    call write_kv_str(u_ctl, u_toml, 'dname_data_remap', 'data_remap_dim')
    ! 71-77: Restart & parameters
    call write_kv_str(u_ctl, u_toml, 'restart_write', 'specified')
    call write_kv_str(u_ctl, u_toml, 'restart_date', '2000-06-01 00:00:00')
    call write_kv    (u_ctl, u_toml, 'restart_month', '6', '6')
    call write_kv    (u_ctl, u_toml, 'restart_day', '15', '15')
    call write_kv    (u_ctl, u_toml, 'restart_hour', '12', '12')
    call write_kv_str(u_ctl, u_toml, 'fname_state_in', 'restart_in.nc')
    call write_kv_str(u_ctl, u_toml, 'param_nml', 'param_custom.nml')
    ! 78-83: User options
    call write_kv    (u_ctl, u_toml, 'qmodOption', '1', '1')
    call write_kv    (u_ctl, u_toml, 'qBlendPeriod', '24', '24')
    call write_kv    (u_ctl, u_toml, 'QerrTrend', '2', '2')
    call write_kv    (u_ctl, u_toml, 'hydGeometryOption', '1', '1')
    call write_kv    (u_ctl, u_toml, 'topoNetworkOption', '0', '0')
    call write_kv    (u_ctl, u_toml, 'computeReachList', '0', '0')
    ! 84-92: Gauge data
    call write_kv_str(u_ctl, u_toml, 'gageMetaFile', 'gage_meta.csv')
    call write_kv    (u_ctl, u_toml, 'outputAtGage', 'T', 'true')
    call write_kv_str(u_ctl, u_toml, 'fname_gageObs', 'gage_obs.nc')
    call write_kv_str(u_ctl, u_toml, 'vname_gageFlow', 'gage_flow_var')
    call write_kv_str(u_ctl, u_toml, 'vname_gageSite', 'gage_site_var')
    call write_kv_str(u_ctl, u_toml, 'vname_gageTime', 'gage_time_var')
    call write_kv_str(u_ctl, u_toml, 'dname_gageSite', 'gage_site_dim')
    call write_kv_str(u_ctl, u_toml, 'dname_gageTime', 'gage_time_dim')
    call write_kv    (u_ctl, u_toml, 'strlen_gageSite', '64', '64')
    ! 93-106: IO, Misc, Pfafcode, Output options
    call write_kv_str(u_ctl, u_toml, 'pio_netcdf_format', 'netcdf4')
    call write_kv_str(u_ctl, u_toml, 'pio_netcdf_type', 'netcdf4p')
    call write_kv    (u_ctl, u_toml, 'debug', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'seg_outlet', '-10', '-10')
    call write_kv    (u_ctl, u_toml, 'desireId', '42', '42')
    call write_kv    (u_ctl, u_toml, 'checkMassBalance', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'maxPfafLen', '48', '48')
    call write_kv_str(u_ctl, u_toml, 'pfafMissing', 'X')
    call write_kv_str(u_ctl, u_toml, 'time_units', 'seconds since 2000-01-01')
    call write_kv_str(u_ctl, u_toml, 'newFileFrequency', 'monthly')
    call write_kv_str(u_ctl, u_toml, 'outputFrequency', 'daily')
    call write_kv_str(u_ctl, u_toml, 'outputNameOption', 'generic')
    call write_kv    (u_ctl, u_toml, 'histTimeStamp_offset', '1800.0', '1800.0')
    call write_kv    (u_ctl, u_toml, 'outputInflow', 'T', 'true')
    ! 107-110: CESM coupling
    call write_kv_str(u_ctl, u_toml, 'qgwl_runoff_option', 'negative')
    call write_kv_str(u_ctl, u_toml, 'bypass_routing_option', 'direct_to_outlet')
    call write_kv    (u_ctl, u_toml, 'correct_area', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'ice_runoff', 'F', 'false')
    ! 111-203: All 93 varname_* keys (including the 57 lake/reservoir keys)
    call write_kv_str(u_ctl, u_toml, 'varname_area', 'v_area')
    call write_kv_str(u_ctl, u_toml, 'varname_HRUid', 'v_HRUid')
    call write_kv_str(u_ctl, u_toml, 'varname_HRUindex', 'v_HRUindex')
    call write_kv_str(u_ctl, u_toml, 'varname_hruSegId', 'v_hruSegId')
    call write_kv_str(u_ctl, u_toml, 'varname_hruSegIndex', 'v_hruSegIndex')
    call write_kv_str(u_ctl, u_toml, 'varname_length', 'v_length')
    call write_kv_str(u_ctl, u_toml, 'varname_slope', 'v_slope')
    call write_kv_str(u_ctl, u_toml, 'varname_width', 'v_width')
    call write_kv_str(u_ctl, u_toml, 'varname_depth', 'v_depth')
    call write_kv_str(u_ctl, u_toml, 'varname_sideSlope', 'v_sideSlope')
    call write_kv_str(u_ctl, u_toml, 'varname_man_n', 'v_man_n')
    call write_kv_str(u_ctl, u_toml, 'varname_floodplainSlope', 'v_floodplainSlope')
    call write_kv_str(u_ctl, u_toml, 'varname_hruArea', 'v_hruArea')
    call write_kv_str(u_ctl, u_toml, 'varname_weight', 'v_weight')
    call write_kv_str(u_ctl, u_toml, 'varname_timeDelayHist', 'v_timeDelayHist')
    call write_kv_str(u_ctl, u_toml, 'varname_upsArea', 'v_upsArea')
    call write_kv_str(u_ctl, u_toml, 'varname_basUnderLake', 'v_basUnderLake')
    call write_kv_str(u_ctl, u_toml, 'varname_rchUnderLake', 'v_rchUnderLake')
    call write_kv_str(u_ctl, u_toml, 'varname_minFlow', 'v_minFlow')
    call write_kv_str(u_ctl, u_toml, 'varname_D03_MaxStorage', 'v_D03_MaxStorage')
    call write_kv_str(u_ctl, u_toml, 'varname_D03_Coefficient', 'v_D03_Coefficient')
    call write_kv_str(u_ctl, u_toml, 'varname_D03_Power', 'v_D03_Power')
    call write_kv_str(u_ctl, u_toml, 'varname_D03_S0', 'v_D03_S0')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_E_emr', 'v_HYP_E_emr')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_E_lim', 'v_HYP_E_lim')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_E_min', 'v_HYP_E_min')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_E_zero', 'v_HYP_E_zero')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_Qrate_emr', 'v_HYP_Qrate_emr')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_Erate_emr', 'v_HYP_Erate_emr')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_Qrate_prim', 'v_HYP_Qrate_prim')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_Qrate_amp', 'v_HYP_Qrate_amp')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_Qrate_phs', 'v_HYP_Qrate_phs')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_prim_F', 'v_HYP_prim_F')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_A_avg', 'v_HYP_A_avg')
    call write_kv_str(u_ctl, u_toml, 'varname_HYP_Qsim_mode', 'v_HYP_Qsim_mode')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_Smax', 'v_H06_Smax')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_alpha', 'v_H06_alpha')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_envfact', 'v_H06_envfact')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_S_ini', 'v_H06_S_ini')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_c1', 'v_H06_c1')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_c2', 'v_H06_c2')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_exponent', 'v_H06_exponent')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_denominator', 'v_H06_denominator')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_c_compare', 'v_H06_c_compare')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_frac_Sdead', 'v_H06_frac_Sdead')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_E_rel_ini', 'v_H06_E_rel_ini')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Jan', 'v_H06_I_Jan')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Feb', 'v_H06_I_Feb')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Mar', 'v_H06_I_Mar')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Apr', 'v_H06_I_Apr')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_May', 'v_H06_I_May')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Jun', 'v_H06_I_Jun')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Jul', 'v_H06_I_Jul')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Aug', 'v_H06_I_Aug')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Sep', 'v_H06_I_Sep')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Oct', 'v_H06_I_Oct')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Nov', 'v_H06_I_Nov')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_Dec', 'v_H06_I_Dec')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Jan', 'v_H06_D_Jan')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Feb', 'v_H06_D_Feb')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Mar', 'v_H06_D_Mar')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Apr', 'v_H06_D_Apr')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_May', 'v_H06_D_May')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Jun', 'v_H06_D_Jun')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Jul', 'v_H06_D_Jul')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Aug', 'v_H06_D_Aug')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Sep', 'v_H06_D_Sep')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Oct', 'v_H06_D_Oct')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Nov', 'v_H06_D_Nov')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_Dec', 'v_H06_D_Dec')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_purpose', 'v_H06_purpose')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_mem_F', 'v_H06_I_mem_F')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_mem_F', 'v_H06_D_mem_F')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_I_mem_L', 'v_H06_I_mem_L')
    call write_kv_str(u_ctl, u_toml, 'varname_H06_D_mem_L', 'v_H06_D_mem_L')
    call write_kv_str(u_ctl, u_toml, 'varname_hruContribIx', 'v_hruContribIx')
    call write_kv_str(u_ctl, u_toml, 'varname_hruContribId', 'v_hruContribId')
    call write_kv_str(u_ctl, u_toml, 'varname_segId', 'v_segId')
    call write_kv_str(u_ctl, u_toml, 'varname_segIndex', 'v_segIndex')
    call write_kv_str(u_ctl, u_toml, 'varname_downSegId', 'v_downSegId')
    call write_kv_str(u_ctl, u_toml, 'varname_downSegIndex', 'v_downSegIndex')
    call write_kv_str(u_ctl, u_toml, 'varname_upSegIds', 'v_upSegIds')
    call write_kv_str(u_ctl, u_toml, 'varname_upSegIndices', 'v_upSegIndices')
    call write_kv_str(u_ctl, u_toml, 'varname_rchOrder', 'v_rchOrder')
    call write_kv_str(u_ctl, u_toml, 'varname_lakeId', 'v_lakeId')
    call write_kv_str(u_ctl, u_toml, 'varname_lakeIndex', 'v_lakeIndex')
    call write_kv_str(u_ctl, u_toml, 'varname_isLakeInlet', 'v_isLakeInlet')
    call write_kv_str(u_ctl, u_toml, 'varname_islake', 'v_islake')
    call write_kv_str(u_ctl, u_toml, 'varname_lakeModelType', 'v_lakeModelType')
    call write_kv_str(u_ctl, u_toml, 'varname_LakeTargVol', 'v_LakeTargVol')
    call write_kv_str(u_ctl, u_toml, 'varname_userTake', 'v_userTake')
    call write_kv_str(u_ctl, u_toml, 'varname_goodBasin', 'v_goodBasin')
    call write_kv_str(u_ctl, u_toml, 'varname_pfafCode', 'v_pfafCode')
    ! 204-226: All 23 output history flags
    call write_kv    (u_ctl, u_toml, 'basRunoff', 'F', 'false')
    call write_kv    (u_ctl, u_toml, 'instRunoff', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'dlayRunoff', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'sumUpstreamRunoff', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'KWTroutedRunoff', 'F', 'false')
    call write_kv    (u_ctl, u_toml, 'IRFroutedRunoff', 'F', 'false')
    call write_kv    (u_ctl, u_toml, 'KWroutedRunoff', 'F', 'false')
    call write_kv    (u_ctl, u_toml, 'DWroutedRunoff', 'F', 'false')
    call write_kv    (u_ctl, u_toml, 'MCroutedRunoff', 'F', 'false')
    call write_kv    (u_ctl, u_toml, 'IRFvolume', 'F', 'false')
    call write_kv    (u_ctl, u_toml, 'KWTvolume', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'KWvolume', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'MCvolume', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'DWvolume', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'KWfloodVolume', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'MCfloodVolume', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'DWfloodVolume', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'KWheight', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'MCheight', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'DWheight', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'localSolute', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'soluteFlux', 'T', 'true')
    call write_kv    (u_ctl, u_toml, 'soluteMass', 'T', 'true')

    close(u_ctl)
    close(u_toml)

    ! Parse legacy .control file and snapshot all variables
    call reset_defaults()
    call read_control(ctl_file, err, message)
    call assert_int_eq(err, 0, 'read_control legacy 226 keys err')

    s_ancil_dir             = ancil_dir
    s_input_dir             = input_dir
    s_output_dir            = output_dir
    s_restart_dir           = restart_dir
    s_case_name             = case_name
    s_simStart              = simStart
    s_simEnd                = simEnd
    s_continue_run          = continue_run
    s_routOpt               = routOpt
    s_doesBasinRoute        = doesBasinRoute
    s_dt                    = dt
    s_floodplain            = floodplain
    s_hw_drain_point        = hw_drain_point
    s_tracer                = tracer
    s_is_lake_sim           = is_lake_sim
    s_lakeRegulate          = lakeRegulate
    s_LakeInputOption       = LakeInputOption
    s_is_flux_wm            = is_flux_wm
    s_is_vol_wm             = is_vol_wm
    s_is_vol_wm_jumpstart   = is_vol_wm_jumpstart
    s_scale_factor_runoff   = scale_factor_runoff
    s_offset_value_runoff   = offset_value_runoff
    s_scale_factor_Ep       = scale_factor_Ep
    s_offset_value_Ep       = offset_value_Ep
    s_is_Ep_upward_negative = is_Ep_upward_negative
    s_scale_factor_prec     = scale_factor_prec
    s_offset_value_prec     = offset_value_prec
    s_min_length_route      = min_length_route
    s_fname_ntopOld         = fname_ntopOld
    s_ntopAugmentMode       = ntopAugmentMode
    s_fname_ntopNew         = fname_ntopNew
    s_dname_nhru            = dname_nhru
    s_dname_sseg            = dname_sseg
    s_fname_qsim            = fname_qsim
    s_vname_qsim            = vname_qsim
    s_vname_evapo           = vname_evapo
    s_vname_precip          = vname_precip
    s_vname_solute          = vname_solute
    s_vname_time            = vname_time
    s_vname_hruid           = vname_hruid
    s_dname_time            = dname_time
    s_dname_hruid           = dname_hruid
    s_dname_xlon            = dname_xlon
    s_dname_ylat            = dname_ylat
    s_units_qsim            = units_qsim
    s_units_cc              = units_cc
    s_dt_ro                 = dt_ro
    s_input_fillvalue       = input_fillvalue
    s_ro_calendar           = ro_calendar
    s_ro_time_units         = ro_time_units
    s_ro_time_stamp         = ro_time_stamp
    s_runoffMin             = runoffMin
    s_fname_wm              = fname_wm
    s_vname_flux_wm         = vname_flux_wm
    s_vname_vol_wm          = vname_vol_wm
    s_vname_time_wm         = vname_time_wm
    s_vname_segid_wm        = vname_segid_wm
    s_dname_time_wm         = dname_time_wm
    s_dname_segid_wm        = dname_segid_wm
    s_dt_wm                 = dt_wm
    s_is_remap              = is_remap
    s_fname_remap           = fname_remap
    s_vname_hruid_in_remap  = vname_hruid_in_remap
    s_vname_weight          = vname_weight
    s_vname_qhruid          = vname_qhruid
    s_vname_num_qhru        = vname_num_qhru
    s_vname_i_index         = vname_i_index
    s_vname_j_index         = vname_j_index
    s_dname_hru_remap       = dname_hru_remap
    s_dname_data_remap      = dname_data_remap
    s_restart_write         = restart_write
    s_restart_date          = restart_date
    s_restart_month         = restart_month
    s_restart_day           = restart_day
    s_restart_hour          = restart_hour
    s_fname_state_in        = fname_state_in
    s_param_nml             = param_nml
    s_qmodOption            = qmodOption
    s_qBlendPeriod          = qBlendPeriod
    s_QerrTrend             = QerrTrend
    s_hydGeometryOption     = hydGeometryOption
    s_topoNetworkOption     = topoNetworkOption
    s_computeReachList      = computeReachList
    s_gageMetaFile          = gageMetaFile
    s_outputAtGage          = outputAtGage
    s_fname_gageObs         = fname_gageObs
    s_vname_gageFlow        = vname_gageFlow
    s_vname_gageSite        = vname_gageSite
    s_vname_gageTime        = vname_gageTime
    s_dname_gageSite        = dname_gageSite
    s_dname_gageTime        = dname_gageTime
    s_strlen_gageSite       = strlen_gageSite
    s_pio_netcdf_format     = pio_netcdf_format
    s_pio_typename          = pio_typename
    s_debug                 = debug
    s_idSegOut              = idSegOut
    s_desireId              = desireId
    s_checkMassBalance      = checkMassBalance
    s_maxPfafLen            = maxPfafLen
    s_pfafMissing           = pfafMissing
    s_time_units            = time_units
    s_newFileFrequency      = newFileFrequency
    s_outputFrequency       = outputFrequency
    s_nOutFreq              = nOutFreq
    s_outputNameOption      = outputNameOption
    s_histTimeStamp_offset  = histTimeStamp_offset
    s_outputInflow          = outputInflow
    s_qgwl_runoff_option    = qgwl_runoff_option
    s_bypass_routing_option = bypass_routing_option
    s_correct_area          = correct_area
    s_ice_runoff            = ice_runoff
    s_isColdStart           = isColdStart
    s_nRoutes               = nRoutes
    s_routeMethods          = routeMethods
    s_onRoute               = onRoute
    s_idxSUM                = idxSUM
    s_idxIRF                = idxIRF
    s_idxKWT                = idxKWT
    s_idxKW                 = idxKW
    s_idxMC                 = idxMC
    s_idxDW                 = idxDW
    s_time_conv             = time_conv
    s_length_conv           = length_conv
    s_time_conv_solute      = time_conv_solute
    s_mass_conv_solute      = mass_conv_solute
    s_meta_HRU              = meta_HRU
    s_meta_HRU2SEG          = meta_HRU2SEG
    s_meta_SEG              = meta_SEG
    s_meta_NTOPO            = meta_NTOPO
    s_meta_PFAF             = meta_PFAF
    s_meta_rflx             = meta_rflx
    s_meta_hflx             = meta_hflx

    ! Parse TOML file from clean defaults and compare every single field
    call reset_defaults()
    call read_control(toml_file, err, message)
    call assert_int_eq(err, 0, 'read_control toml 226 keys err')

    call assert_char_eq(ancil_dir,             s_ancil_dir,             'ancil_dir')
    call assert_char_eq(input_dir,             s_input_dir,             'input_dir')
    call assert_char_eq(output_dir,            s_output_dir,            'output_dir')
    call assert_char_eq(restart_dir,           s_restart_dir,           'restart_dir')
    call assert_char_eq(case_name,             s_case_name,             'case_name')
    call assert_char_eq(simStart,              s_simStart,              'sim_start')
    call assert_char_eq(simEnd,                s_simEnd,                'sim_end')
    call assert_bool_eq(continue_run,          s_continue_run,          'continue_run')
    call assert_char_eq(routOpt,               s_routOpt,               'route_opt')
    call assert_int_eq (doesBasinRoute,        s_doesBasinRoute,        'doesBasinRoute')
    call assert_dp_eq  (dt,                    s_dt,                    'dt_qsim')
    call assert_bool_eq(floodplain,            s_floodplain,            'floodplain')
    call assert_int_eq (hw_drain_point,        s_hw_drain_point,        'hw_drain_point')
    call assert_bool_eq(tracer,                s_tracer,                'tracer')
    call assert_bool_eq(is_lake_sim,           s_is_lake_sim,           'is_lake_sim')
    call assert_bool_eq(lakeRegulate,          s_lakeRegulate,          'lakeRegulate')
    call assert_int_eq (LakeInputOption,       s_LakeInputOption,       'LakeInputOption')
    call assert_bool_eq(is_flux_wm,            s_is_flux_wm,            'is_flux_wm')
    call assert_bool_eq(is_vol_wm,             s_is_vol_wm,             'is_vol_wm')
    call assert_bool_eq(is_vol_wm_jumpstart,   s_is_vol_wm_jumpstart,   'is_vol_wm_jumpstart')
    call assert_dp_eq  (scale_factor_runoff,   s_scale_factor_runoff,   'scale_factor_runoff')
    call assert_dp_eq  (offset_value_runoff,   s_offset_value_runoff,   'offset_value_runoff')
    call assert_dp_eq  (scale_factor_Ep,       s_scale_factor_Ep,       'scale_factor_Ep')
    call assert_dp_eq  (offset_value_Ep,       s_offset_value_Ep,       'offset_value_Ep')
    call assert_bool_eq(is_Ep_upward_negative, s_is_Ep_upward_negative, 'is_Ep_upward_negative')
    call assert_dp_eq  (scale_factor_prec,     s_scale_factor_prec,     'scale_factor_prec')
    call assert_dp_eq  (offset_value_prec,     s_offset_value_prec,     'offset_value_prec')
    call assert_dp_eq  (min_length_route,      s_min_length_route,      'min_length_route')
    call assert_char_eq(fname_ntopOld,         s_fname_ntopOld,         'fname_ntopOld')
    call assert_bool_eq(ntopAugmentMode,       s_ntopAugmentMode,       'ntopAugmentMode')
    call assert_char_eq(fname_ntopNew,         s_fname_ntopNew,         'fname_ntopNew')
    call assert_char_eq(dname_nhru,            s_dname_nhru,            'dname_nhru')
    call assert_char_eq(dname_sseg,            s_dname_sseg,            'dname_sseg')
    call assert_char_eq(fname_qsim,            s_fname_qsim,            'fname_qsim')
    call assert_char_eq(vname_qsim,            s_vname_qsim,            'vname_qsim')
    call assert_char_eq(vname_evapo,           s_vname_evapo,           'vname_evapo')
    call assert_char_eq(vname_precip,          s_vname_precip,          'vname_precip')
    call assert_char_eq(vname_solute,          s_vname_solute,          'vname_solute')
    call assert_char_eq(vname_time,            s_vname_time,            'vname_time')
    call assert_char_eq(vname_hruid,           s_vname_hruid,           'vname_hruid')
    call assert_char_eq(dname_time,            s_dname_time,            'dname_time')
    call assert_char_eq(dname_hruid,           s_dname_hruid,           'dname_hruid')
    call assert_char_eq(dname_xlon,            s_dname_xlon,            'dname_xlon')
    call assert_char_eq(dname_ylat,            s_dname_ylat,            'dname_ylat')
    call assert_char_eq(units_qsim,            s_units_qsim,            'units_qsim')
    call assert_char_eq(units_cc,              s_units_cc,              'units_cc')
    call assert_dp_eq  (dt_ro,                 s_dt_ro,                 'dt_ro')
    call assert_dp_eq  (input_fillvalue,       s_input_fillvalue,       'input_fillvalue')
    call assert_char_eq(ro_calendar,           s_ro_calendar,           'ro_calendar')
    call assert_char_eq(ro_time_units,         s_ro_time_units,         'ro_time_units')
    call assert_char_eq(ro_time_stamp,         s_ro_time_stamp,         'ro_time_stamp')
    call assert_dp_eq  (runoffMin,             s_runoffMin,             'runoffMin')
    call assert_char_eq(fname_wm,              s_fname_wm,              'fname_wm')
    call assert_char_eq(vname_flux_wm,         s_vname_flux_wm,         'vname_flux_wm')
    call assert_char_eq(vname_vol_wm,          s_vname_vol_wm,          'vname_vol_wm')
    call assert_char_eq(vname_time_wm,         s_vname_time_wm,         'vname_time_wm')
    call assert_char_eq(vname_segid_wm,        s_vname_segid_wm,        'vname_segid_wm')
    call assert_char_eq(dname_time_wm,         s_dname_time_wm,         'dname_time_wm')
    call assert_char_eq(dname_segid_wm,        s_dname_segid_wm,        'dname_segid_wm')
    call assert_dp_eq  (dt_wm,                 s_dt_wm,                 'dt_wm')
    call assert_bool_eq(is_remap,              s_is_remap,              'is_remap')
    call assert_char_eq(fname_remap,           s_fname_remap,           'fname_remap')
    call assert_char_eq(vname_hruid_in_remap,  s_vname_hruid_in_remap,  'vname_hruid_in_remap')
    call assert_char_eq(vname_weight,          s_vname_weight,          'vname_weight')
    call assert_char_eq(vname_qhruid,          s_vname_qhruid,          'vname_qhruid')
    call assert_char_eq(vname_num_qhru,        s_vname_num_qhru,        'vname_num_qhru')
    call assert_char_eq(vname_i_index,         s_vname_i_index,         'vname_i_index')
    call assert_char_eq(vname_j_index,         s_vname_j_index,         'vname_j_index')
    call assert_char_eq(dname_hru_remap,       s_dname_hru_remap,       'dname_hru_remap')
    call assert_char_eq(dname_data_remap,      s_dname_data_remap,      'dname_data_remap')
    call assert_char_eq(restart_write,         s_restart_write,         'restart_write')
    call assert_char_eq(restart_date,          s_restart_date,          'restart_date')
    call assert_int_eq (restart_month,         s_restart_month,         'restart_month')
    call assert_int_eq (restart_day,           s_restart_day,           'restart_day')
    call assert_int_eq (restart_hour,          s_restart_hour,          'restart_hour')
    call assert_char_eq(fname_state_in,        s_fname_state_in,        'fname_state_in')
    call assert_char_eq(param_nml,             s_param_nml,             'param_nml')
    call assert_int_eq (qmodOption,            s_qmodOption,            'qmodOption')
    call assert_int_eq (qBlendPeriod,          s_qBlendPeriod,          'qBlendPeriod')
    call assert_int_eq (QerrTrend,             s_QerrTrend,             'QerrTrend')
    call assert_int_eq (hydGeometryOption,     s_hydGeometryOption,     'hydGeometryOption')
    call assert_int_eq (topoNetworkOption,     s_topoNetworkOption,     'topoNetworkOption')
    call assert_int_eq (computeReachList,      s_computeReachList,      'computeReachList')
    call assert_char_eq(gageMetaFile,          s_gageMetaFile,          'gageMetaFile')
    call assert_bool_eq(outputAtGage,          s_outputAtGage,          'outputAtGage')
    call assert_char_eq(fname_gageObs,         s_fname_gageObs,         'fname_gageObs')
    call assert_char_eq(vname_gageFlow,        s_vname_gageFlow,        'vname_gageFlow')
    call assert_char_eq(vname_gageSite,        s_vname_gageSite,        'vname_gageSite')
    call assert_char_eq(vname_gageTime,        s_vname_gageTime,        'vname_gageTime')
    call assert_char_eq(dname_gageSite,        s_dname_gageSite,        'dname_gageSite')
    call assert_char_eq(dname_gageTime,        s_dname_gageTime,        'dname_gageTime')
    call assert_int_eq (strlen_gageSite,       s_strlen_gageSite,       'strlen_gageSite')
    call assert_char_eq(pio_netcdf_format,     s_pio_netcdf_format,     'pio_netcdf_format')
    call assert_char_eq(pio_typename,          s_pio_typename,          'pio_netcdf_type')
    call assert_bool_eq(debug,                 s_debug,                 'debug')
    call assert_int_eq (idSegOut,              s_idSegOut,              'seg_outlet')
    call assert_int_eq (desireId,              s_desireId,              'desireId')
    call assert_bool_eq(checkMassBalance,      s_checkMassBalance,      'checkMassBalance')
    call assert_int_eq (maxPfafLen,            s_maxPfafLen,            'maxPfafLen')
    call assert_char_eq(pfafMissing,           s_pfafMissing,           'pfafMissing')
    call assert_char_eq(time_units,            s_time_units,            'time_units')
    call assert_char_eq(newFileFrequency,      s_newFileFrequency,      'newFileFrequency')
    call assert_char_eq(outputFrequency,       s_outputFrequency,       'outputFrequency')
    call assert_int_eq (nOutFreq,              s_nOutFreq,              'nOutFreq')
    call assert_char_eq(outputNameOption,      s_outputNameOption,      'outputNameOption')
    call assert_dp_eq  (histTimeStamp_offset,  s_histTimeStamp_offset,  'histTimeStamp_offset')
    call assert_bool_eq(outputInflow,          s_outputInflow,          'outputInflow')
    call assert_char_eq(qgwl_runoff_option,    s_qgwl_runoff_option,    'qgwl_runoff_option')
    call assert_char_eq(bypass_routing_option, s_bypass_routing_option, 'bypass_routing_option')
    call assert_bool_eq(correct_area,          s_correct_area,          'correct_area')
    call assert_bool_eq(ice_runoff,            s_ice_runoff,            'ice_runoff')

    ! Compare globalData derived state and metadata structures
    call assert_bool_eq(isColdStart,      s_isColdStart,      'isColdStart')
    call assert_int_eq (nRoutes,          s_nRoutes,          'nRoutes')
    call assert_int_eq (idxSUM,           s_idxSUM,           'idxSUM')
    call assert_int_eq (idxIRF,           s_idxIRF,           'idxIRF')
    call assert_int_eq (idxKWT,           s_idxKWT,           'idxKWT')
    call assert_int_eq (idxKW,            s_idxKW,            'idxKW')
    call assert_int_eq (idxMC,            s_idxMC,            'idxMC')
    call assert_int_eq (idxDW,            s_idxDW,            'idxDW')
    call assert_dp_eq  (time_conv,        s_time_conv,        'time_conv')
    call assert_dp_eq  (length_conv,      s_length_conv,      'length_conv')
    call assert_dp_eq  (time_conv_solute, s_time_conv_solute, 'time_conv_solute')
    call assert_dp_eq  (mass_conv_solute, s_mass_conv_solute, 'mass_conv_solute')

    do i = 1, nRoutes
      call assert_int_eq(routeMethods(i), s_routeMethods(i), 'routeMethods')
    end do
    do i = 0, nRouteMethods - 1
      call assert_bool_eq(onRoute(i), s_onRoute(i), 'onRoute')
    end do

    do i = 1, nVarsHRU
      call assert_char_eq(meta_HRU(i)%varName, s_meta_HRU(i)%varName, 'meta_HRU%varName')
    end do
    do i = 1, nVarsHRU2SEG
      call assert_char_eq(meta_HRU2SEG(i)%varName, s_meta_HRU2SEG(i)%varName, 'meta_HRU2SEG%varName')
    end do
    do i = 1, nVarsSEG
      call assert_char_eq(meta_SEG(i)%varName, s_meta_SEG(i)%varName, 'meta_SEG%varName')
    end do
    do i = 1, nVarsNTOPO
      call assert_char_eq(meta_NTOPO(i)%varName, s_meta_NTOPO(i)%varName, 'meta_NTOPO%varName')
    end do
    do i = 1, nVarsPFAF
      call assert_char_eq(meta_PFAF(i)%varName, s_meta_PFAF(i)%varName, 'meta_PFAF%varName')
    end do
    do i = 1, nVarsRFLX
      call assert_bool_eq(meta_rflx(i)%varFile, s_meta_rflx(i)%varFile, 'meta_rflx%varFile')
      call assert_char_eq(meta_rflx(i)%varName, s_meta_rflx(i)%varName, 'meta_rflx%varName')
    end do
    do i = 1, nVarsHFLX
      call assert_bool_eq(meta_hflx(i)%varFile, s_meta_hflx(i)%varFile, 'meta_hflx%varFile')
      call assert_char_eq(meta_hflx(i)%varUnit, s_meta_hflx(i)%varUnit, 'meta_hflx%varUnit')
    end do

    ! Spot-check that the newly added 57 lake/reservoir varname_* keys actually took the TOML values
    call assert_char_eq(meta_SEG(ixSEG%basUnderLake)%varName,   'v_basUnderLake',   'varname_basUnderLake')
    call assert_char_eq(meta_SEG(ixSEG%D03_MaxStorage)%varName, 'v_D03_MaxStorage', 'varname_D03_MaxStorage')
    call assert_char_eq(meta_SEG(ixSEG%D03_Power)%varName,      'v_D03_Power',      'varname_D03_Power')
    call assert_char_eq(meta_SEG(ixSEG%D03_S0)%varName,         'v_D03_S0',         'varname_D03_S0')
    call assert_char_eq(meta_SEG(ixSEG%HYP_Qsim_mode)%varName,  'v_HYP_Qsim_mode',  'varname_HYP_Qsim_mode')
    call assert_char_eq(meta_SEG(ixSEG%H06_Smax)%varName,       'v_H06_Smax',       'varname_H06_Smax')
    call assert_char_eq(meta_SEG(ixSEG%H06_D_mem_L)%varName,    'v_H06_D_mem_L',    'varname_H06_D_mem_L')

    open(newunit=u_ctl,  file=ctl_file,  status='old')
    close(u_ctl,  status='delete')
    open(newunit=u_toml, file=toml_file, status='old')
    close(u_toml, status='delete')

    write(*, '(a)') '  [PASS] test_parity_all_226_keys'
  end subroutine test_parity_all_226_keys

  ! Test 2: Verify outputFrequency works as both a TOML string and an unquoted TOML integer
  subroutine test_output_frequency_types()
    integer(i4b) :: u_toml, err
    character(len=strLen) :: message
    character(len=*), parameter :: toml_file = 'test_outfreq.toml'

    ! Case A: outputFrequency as a string ("monthly")
    ! Set up temporary TOML file with outputFrequency as a string keyword
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    write(u_toml, '(a)') 'newFileFrequency = "yearly"'
    write(u_toml, '(a)') 'outputFrequency = "monthly"'
    close(u_toml)

    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq (err,             0,              'outputFrequency string err')
    call assert_char_eq(outputFrequency, 'monthly',      'outputFrequency string value')
    call assert_int_eq (nOutFreq,        integerMissing, 'outputFrequency string nOutFreq')

    ! Case B: outputFrequency as an unquoted integer (3 timesteps of 3600s = 10800s, divides 86400s)
    ! Set up temporary TOML file with outputFrequency as an unquoted integer timestep count
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    write(u_toml, '(a)') 'outputFrequency = 3'
    close(u_toml)

    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq (err,             0,   'outputFrequency integer err')
    call assert_char_eq(outputFrequency, '3', 'outputFrequency integer char value')
    call assert_int_eq (nOutFreq,        3,   'outputFrequency integer nOutFreq')

    open(newunit=u_toml, file=toml_file, status='old')
    close(u_toml, status='delete')

    write(*, '(a)') '  [PASS] test_output_frequency_types'
  end subroutine test_output_frequency_types

  ! Test 3: Verify omitted optional keys retain their public_var.f90 and popMetadat defaults
  subroutine test_default_preservation()
    integer(i4b) :: u_toml, err
    character(len=strLen) :: message
    character(len=*), parameter :: toml_file = 'test_defaults.toml'

    ! Set up minimal temporary TOML file omitting optional keys to verify default preservation
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    close(u_toml)

    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, 0, 'default preservation err')

    ! Check representative character default: verifies intent(inout) dummy argument
    ! preserves pre-initialized character defaults when optional keys are omitted in TOML
    call assert_char_eq(dname_sseg, 'seg', 'default dname_sseg (intent(inout) character preservation)')

    ! Check representative integer, real, and logical defaults from public_var.f90
    call assert_int_eq (doesBasinRoute,   1,       'default doesBasinRoute (integer)')
    call assert_dp_eq  (min_length_route, 0.0_dp,  'default min_length_route (real)')
    call assert_bool_eq(lakeRegulate,     .true.,  'default lakeRegulate (logical)')

    ! Check representative popMetadat default is preserved when varname_* keys are omitted
    call assert_char_eq(meta_HRU(ixHRU%area)%varName, 'area', 'default meta_HRU area (popMetadat)')

    open(newunit=u_toml, file=toml_file, status='old')
    close(u_toml, status='delete')

    write(*, '(a)') '  [PASS] test_default_preservation'
  end subroutine test_default_preservation

  ! Test 4: Verify error handling (err = 20) for syntax errors, unknown keys, and type mismatches
  subroutine test_error_handling()
    integer(i4b)                :: u_toml, err
    character(len=strLen)       :: message
    character(len=*), parameter :: toml_file = 'test_errors.toml'
    integer(i4b), parameter     :: expected_parse_err = 20 ! mizuRoute error code for control parsing failures
    integer(i4b), parameter     :: expected_val_err   = 81 ! mizuRoute error code for control validation failures

    ! 4a: TOML syntax error (extra equals sign)
    ! Set up temporary TOML file with a syntax error (extra equals sign) to test parser error detection
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'dt_qsim = = 3600'
    close(u_toml)
    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, expected_parse_err, 'syntax error err=20')
    call assert_true(index(message, 'read_control_toml/') > 0, 'syntax error prefix')
    call assert_true(index(message, 'failed to parse TOML control file') > 0, 'syntax error message')

    ! 4b: Unknown key error
    ! Set up temporary TOML file with an unrecognized key to test unknown variable error detection
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'unknown_bogus_key = 123'
    close(u_toml)
    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, expected_parse_err, 'unknown key err=20')
    call assert_true(index(message, 'read_control_toml/') > 0, 'unknown key prefix')
    call assert_true(index(message, 'unexpected variable in TOML control file: unknown_bogus_key') > 0, 'unknown key message')

    ! 4c: Type mismatch - string for integer (route_opt = "5" instead of 5)
    ! Set up temporary TOML file with string route_opt to test integer type mismatch error detection
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    write(u_toml, '(a)') 'route_opt = "5"'
    close(u_toml)
    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, expected_parse_err, 'route_opt string type mismatch err=20')
    call assert_true(index(message, 'read_control_toml/') > 0, 'route_opt type mismatch prefix')
    call assert_true(index(message, 'problem with internal read of route_opt') > 0, 'route_opt type mismatch message')

    ! 4d: Type mismatch - pfafMissing must be a string, despite looking like an int (pfafMissing = 0),
    ! because Pfafstetter sentinels can be either 0-9 or non-numeric character sentinels (e.g. ' ', 'X', or '-')
    ! Set up temporary TOML file with integer pfafMissing to test string type mismatch error detection
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    write(u_toml, '(a)') 'pfafMissing = 0'
    close(u_toml)
    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, expected_parse_err, 'pfafMissing int type mismatch err=20')
    call assert_true(index(message, 'read_control_toml/') > 0, 'pfafMissing type mismatch prefix')
    call assert_true(index(message, 'problem with internal read of pfafMissing') > 0, 'pfafMissing type mismatch message')

    ! 4e: Type mismatch - string for boolean (is_lake_sim = "false" instead of false)
    ! Set up temporary TOML file with string boolean to test type mismatch detection
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    write(u_toml, '(a)') 'is_lake_sim = "false"'
    close(u_toml)
    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, expected_parse_err, 'is_lake_sim string type mismatch err=20')
    call assert_true(index(message, 'read_control_toml/') > 0, 'is_lake_sim type mismatch prefix')
    call assert_true(index(message, 'problem with internal read of is_lake_sim') > 0, 'is_lake_sim type mismatch message')

    ! Verify that unquoted false is accepted as valid boolean in TOML
    ! Set up temporary TOML file with unquoted false to validate correct TOML boolean parsing
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    write(u_toml, '(a)') 'is_lake_sim = false'
    close(u_toml)
    call reset_defaults()
    is_lake_sim = .true.
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, 0, 'is_lake_sim unquoted false err=0')
    call assert_bool_eq(is_lake_sim, .false., 'is_lake_sim unquoted false value')

    ! 4f: Type mismatch - boolean for outputFrequency (outputFrequency = true)
    ! Set up temporary TOML file with boolean outputFrequency to test type mismatch error detection
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    write(u_toml, '(a)') 'outputFrequency = true'
    close(u_toml)
    call reset_defaults()
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, expected_parse_err, 'outputFrequency bool type mismatch err=20')
    call assert_true(index(message, 'read_control_toml/') > 0, 'outputFrequency type mismatch prefix')
    call assert_true(index(message, 'problem with internal read of outputFrequency') > 0, 'outputFrequency type mismatch message')

    ! 4g: Validate child_message chaining and non-masterproc ro_time_stamp check
    ! Set up temporary TOML file with invalid ro_time_stamp to test validation error and error chaining
    open(newunit=u_toml, file=toml_file, status='replace', action='write')
    write(u_toml, '(a)') 'units_qsim = "mm/s"'
    write(u_toml, '(a)') 'dt_qsim = 3600'
    write(u_toml, '(a)') 'ro_time_stamp = "bad_stamp"'
    close(u_toml)
    call reset_defaults()
    masterproc = .false.
    call read_control_toml(toml_file, err, message)
    call assert_int_eq(err, expected_val_err, 'ro_time_stamp validation err=81 on non-masterproc')
    call assert_true(index(message, 'read_control_toml/validate_and_finalize_control/') > 0, 'child_message chaining prefix')

    open(newunit=u_toml, file=toml_file, status='old')
    close(u_toml, status='delete')

    write(*, '(a)') '  [PASS] test_error_handling'
  end subroutine test_error_handling

end program test_read_control
