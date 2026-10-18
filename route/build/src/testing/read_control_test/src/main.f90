! Unit tests verifying read_control_toml outputFrequency types, defaults, and error handling
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
