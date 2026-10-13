MODULE read_control_module

! reading mizuRoute control files and save them in variables defined in public_var.f90

USE nrtype
USE public_var
USE tomlf, ONLY: toml_table, toml_error, toml_key, toml_stat, &
                 toml_load, &
                 toml_get_value => get_value

implicit none

INTERFACE get_toml_val
   MODULE PROCEDURE get_toml_val_char
   MODULE PROCEDURE get_toml_val_int
   MODULE PROCEDURE get_toml_val_dp
   MODULE PROCEDURE get_toml_val_sp
   MODULE PROCEDURE get_toml_val_bool
END INTERFACE get_toml_val

private
public::read_control
public::read_control_toml
public::read_control_legacy

CONTAINS

 ! =======================================================================================================
 ! public subroutine: read the control file (generic dispatcher based on file suffix)
 ! =======================================================================================================
 SUBROUTINE read_control(ctl_fname, err, message)
   USE ascii_utils, ONLY: lower            ! convert string to lower case

   implicit none
   ! argument variables
   character(*), intent(in)          :: ctl_fname               ! name of the control file
   integer(i4b),intent(out)          :: err                     ! error code
   character(*),intent(out)          :: message                 ! error message
   ! Local variables
   character(len=strLen)             :: cmessage                ! error message of downwind routine
   logical(lgt)                      :: is_toml                 ! true if control file is in TOML format

   err=0; message='read_control/'

   is_toml = .false.
   if (len_trim(ctl_fname) >= 5) then
      if (lower(ctl_fname(len_trim(ctl_fname)-4:len_trim(ctl_fname))) == '_toml' .or. &
          lower(ctl_fname(len_trim(ctl_fname)-4:len_trim(ctl_fname))) == '.toml') then
         is_toml = .true.
      endif
   endif

   if (is_toml) then
      call read_control_toml(ctl_fname, err, cmessage)
   else
      call read_control_legacy(ctl_fname, err, cmessage)
   endif
   if (err /= 0) message = trim(message) // trim(cmessage)

 END SUBROUTINE read_control

 ! =======================================================================================================
 ! public subroutine: read TOML format control file
 ! =======================================================================================================
 SUBROUTINE read_control_toml(ctl_fname, err, message)
   ! global vars
   USE globalData
   USE var_lookup

   implicit none
   ! argument variables
   character(*), intent(in)          :: ctl_fname               ! name of the control file
   integer(i4b),intent(out)          :: err                     ! error code
   character(*),intent(out)          :: message                 ! error message
   ! Local variables
   type(toml_table), allocatable     :: table                   ! TOML root table structure
   type(toml_error), allocatable     :: error                   ! TOML parser error
   type(toml_key), allocatable       :: keys(:)                 ! list of keys in TOML table
   integer(i4b)                      :: iKey                    ! loop index

   err=0; message='read_control_toml/'

   if (masterproc) then
      write(iulog,'(2a)') new_line('a'), '---- read TOML control file --- '
   end if

   ! Open and parse TOML control file using toml-f library
   call toml_load(table, trim(ctl_fname), error=error)

   if (allocated(error)) then
      err = 20
      message = trim(message) // 'failed to parse TOML control file: ' // trim(error%message)
      return
   endif

   ! Validate keys in TOML control file against known control variables
   call table%get_keys(keys)
   if (allocated(keys)) then
      do iKey = 1, size(keys)
         if (.not. is_valid_control_key(keys(iKey)%key)) then
            err = 20
            message = trim(message) // 'unexpected variable in TOML control file: ' // trim(keys(iKey)%key)
            return
         endif
      end do
   endif

   ! Extract variables from TOML table using toml-f interface
   call get_toml_val(table, "ancil_dir", ancil_dir)
   call get_toml_val(table, "input_dir", input_dir)
   call get_toml_val(table, "output_dir", output_dir)
   call get_toml_val(table, "restart_dir", restart_dir)
   call get_toml_val(table, "case_name", case_name)
   call get_toml_val(table, "sim_start", simStart)
   call get_toml_val(table, "sim_end", simEnd)
   call get_toml_val(table, "continue_run", continue_run)
   call get_toml_val(table, "route_opt", routOpt)
   call get_toml_val(table, "doesBasinRoute", doesBasinRoute)
   call get_toml_val(table, "dt_qsim", dt)
   call get_toml_val(table, "floodplain", floodplain)
   call get_toml_val(table, "hw_drain_point", hw_drain_point)
   call get_toml_val(table, "tracer", tracer)
   call get_toml_val(table, "is_lake_sim", is_lake_sim)
   call get_toml_val(table, "lakeRegulate", lakeRegulate)
   call get_toml_val(table, "LakeInputOption", LakeInputOption)
   call get_toml_val(table, "is_flux_wm", is_flux_wm)
   call get_toml_val(table, "is_vol_wm", is_vol_wm)
   call get_toml_val(table, "is_vol_wm_jumpstart", is_vol_wm_jumpstart)
   call get_toml_val(table, "scale_factor_runoff", scale_factor_runoff)
   call get_toml_val(table, "offset_value_runoff", offset_value_runoff)
   call get_toml_val(table, "scale_factor_Ep", scale_factor_Ep)
   call get_toml_val(table, "offset_value_Ep", offset_value_Ep)
   call get_toml_val(table, "is_Ep_upward_negative", is_Ep_upward_negative)
   call get_toml_val(table, "scale_factor_prec", scale_factor_prec)
   call get_toml_val(table, "offset_value_prec", offset_value_prec)
   call get_toml_val(table, "min_length_route", min_length_route)
   call get_toml_val(table, "fname_ntopOld", fname_ntopOld)
   call get_toml_val(table, "ntopAugmentMode", ntopAugmentMode)
   call get_toml_val(table, "fname_ntopNew", fname_ntopNew)
   call get_toml_val(table, "dname_nhru", dname_nhru)
   call get_toml_val(table, "dname_sseg", dname_sseg)
   call get_toml_val(table, "fname_qsim", fname_qsim)
   call get_toml_val(table, "vname_qsim", vname_qsim)
   call get_toml_val(table, "vname_evapo", vname_evapo)
   call get_toml_val(table, "vname_precip", vname_precip)
   call get_toml_val(table, "vname_solute", vname_solute)
   call get_toml_val(table, "vname_time", vname_time)
   call get_toml_val(table, "vname_hruid", vname_hruid)
   call get_toml_val(table, "dname_time", dname_time)
   call get_toml_val(table, "dname_hruid", dname_hruid)
   call get_toml_val(table, "dname_xlon", dname_xlon)
   call get_toml_val(table, "dname_ylat", dname_ylat)
   call get_toml_val(table, "units_qsim", units_qsim)
   call get_toml_val(table, "units_cc", units_cc)
   call get_toml_val(table, "dt_ro", dt_ro)
   call get_toml_val(table, "input_fillvalue", input_fillvalue)
   call get_toml_val(table, "ro_calendar", ro_calendar)
   call get_toml_val(table, "ro_time_units", ro_time_units)
   call get_toml_val(table, "ro_time_stamp", ro_time_stamp)
   call get_toml_val(table, "runoffMin", runoffMin)
   call get_toml_val(table, "fname_wm", fname_wm)
   call get_toml_val(table, "vname_flux_wm", vname_flux_wm)
   call get_toml_val(table, "vname_vol_wm", vname_vol_wm)
   call get_toml_val(table, "vname_time_wm", vname_time_wm)
   call get_toml_val(table, "vname_segid_wm", vname_segid_wm)
   call get_toml_val(table, "dname_time_wm", dname_time_wm)
   call get_toml_val(table, "dname_segid_wm", dname_segid_wm)
   call get_toml_val(table, "dt_wm", dt_wm)
   call get_toml_val(table, "is_remap", is_remap)
   call get_toml_val(table, "fname_remap", fname_remap)
   call get_toml_val(table, "vname_hruid_in_remap", vname_hruid_in_remap)
   call get_toml_val(table, "vname_weight", vname_weight)
   call get_toml_val(table, "vname_qhruid", vname_qhruid)
   call get_toml_val(table, "vname_num_qhru", vname_num_qhru)
   call get_toml_val(table, "vname_i_index", vname_i_index)
   call get_toml_val(table, "vname_j_index", vname_j_index)
   call get_toml_val(table, "dname_hru_remap", dname_hru_remap)
   call get_toml_val(table, "dname_data_remap", dname_data_remap)
   call get_toml_val(table, "restart_write", restart_write)
   call get_toml_val(table, "restart_date", restart_date)
   call get_toml_val(table, "restart_month", restart_month)
   call get_toml_val(table, "restart_day", restart_day)
   call get_toml_val(table, "restart_hour", restart_hour)
   call get_toml_val(table, "fname_state_in", fname_state_in)
   call get_toml_val(table, "param_nml", param_nml)
   call get_toml_val(table, "qmodOption", qmodOption)
   call get_toml_val(table, "qBlendPeriod", qBlendPeriod)
   call get_toml_val(table, "QerrTrend", QerrTrend)
   call get_toml_val(table, "hydGeometryOption", hydGeometryOption)
   call get_toml_val(table, "topoNetworkOption", topoNetworkOption)
   call get_toml_val(table, "computeReachList", computeReachList)
   call get_toml_val(table, "gageMetaFile", gageMetaFile)
   call get_toml_val(table, "outputAtGage", outputAtGage)
   call get_toml_val(table, "fname_gageObs", fname_gageObs)
   call get_toml_val(table, "vname_gageFlow", vname_gageFlow)
   call get_toml_val(table, "vname_gageSite", vname_gageSite)
   call get_toml_val(table, "vname_gageTime", vname_gageTime)
   call get_toml_val(table, "dname_gageSite", dname_gageSite)
   call get_toml_val(table, "dname_gageTime", dname_gageTime)
   call get_toml_val(table, "strlen_gageSite", strlen_gageSite)
   call get_toml_val(table, "pio_netcdf_format", pio_netcdf_format)
   call get_toml_val(table, "pio_netcdf_type", pio_typename)
   call get_toml_val(table, "debug", debug)
   call get_toml_val(table, "seg_outlet", idSegOut)
   call get_toml_val(table, "desireId", desireId)
   call get_toml_val(table, "checkMassBalance", checkMassBalance)
   call get_toml_val(table, "maxPfafLen", maxPfafLen)
   call get_toml_val(table, "pfafMissing", pfafMissing)
   call get_toml_val(table, "time_units", time_units)
   call get_toml_val(table, "newFileFrequency", newFileFrequency)
   call get_toml_val(table, "outputFrequency", outputFrequency)
   call get_toml_val(table, "outputNameOption", outputNameOption)
   call get_toml_val(table, "histTimeStamp_offset", histTimeStamp_offset)
   call get_toml_val(table, "outputInflow", outputInflow)
   call get_toml_val(table, "qgwl_runoff_option", qgwl_runoff_option)
   call get_toml_val(table, "bypass_routing_option", bypass_routing_option)
   call get_toml_val(table, "correct_area", correct_area)
   call get_toml_val(table, "ice_runoff", ice_runoff)
   call get_toml_val(table, "varname_area", meta_HRU(ixHRU%area)%varName)
   call get_toml_val(table, "varname_HRUid", meta_HRU2SEG(ixHRU2SEG%HRUid)%varName)
   call get_toml_val(table, "varname_HRUindex", meta_HRU2SEG(ixHRU2SEG%HRUindex)%varName)
   call get_toml_val(table, "varname_hruSegId", meta_HRU2SEG(ixHRU2SEG%hruSegId)%varName)
   call get_toml_val(table, "varname_hruSegIndex", meta_HRU2SEG(ixHRU2SEG%hruSegIndex)%varName)
   call get_toml_val(table, "varname_length", meta_SEG(ixSEG%length)%varName)
   call get_toml_val(table, "varname_slope", meta_SEG(ixSEG%slope)%varName)
   call get_toml_val(table, "varname_width", meta_SEG(ixSEG%width)%varName)
   call get_toml_val(table, "varname_depth", meta_SEG(ixSEG%depth)%varName)
   call get_toml_val(table, "varname_sideSlope", meta_SEG(ixSEG%sideSlope)%varName)
   call get_toml_val(table, "varname_man_n", meta_SEG(ixSEG%man_n)%varName)
   call get_toml_val(table, "varname_floodplainSlope", meta_SEG(ixSEG%floodplainSlope)%varName)
   call get_toml_val(table, "varname_hruArea", meta_SEG(ixSEG%hruArea)%varName)
   call get_toml_val(table, "varname_weight", meta_SEG(ixSEG%weight)%varName)
   call get_toml_val(table, "varname_timeDelayHist", meta_SEG(ixSEG%timeDelayHist)%varName)
   call get_toml_val(table, "varname_upsArea", meta_SEG(ixSEG%upsArea)%varName)
   call get_toml_val(table, "varname_hruContribIx", meta_NTOPO(ixNTOPO%hruContribIx)%varName)
   call get_toml_val(table, "varname_hruContribId", meta_NTOPO(ixNTOPO%hruContribId)%varName)
   call get_toml_val(table, "varname_segId", meta_NTOPO(ixNTOPO%segId)%varName)
   call get_toml_val(table, "varname_segIndex", meta_NTOPO(ixNTOPO%segIndex)%varName)
   call get_toml_val(table, "varname_downSegId", meta_NTOPO(ixNTOPO%downSegId)%varName)
   call get_toml_val(table, "varname_downSegIndex", meta_NTOPO(ixNTOPO%downSegIndex)%varName)
   call get_toml_val(table, "varname_upSegIds", meta_NTOPO(ixNTOPO%upSegIds)%varName)
   call get_toml_val(table, "varname_upSegIndices", meta_NTOPO(ixNTOPO%upSegIndices)%varName)
   call get_toml_val(table, "varname_rchOrder", meta_NTOPO(ixNTOPO%rchOrder)%varName)
   call get_toml_val(table, "varname_lakeId", meta_NTOPO(ixNTOPO%lakeId)%varName)
   call get_toml_val(table, "varname_lakeIndex", meta_NTOPO(ixNTOPO%lakeIndex)%varName)
   call get_toml_val(table, "varname_isLakeInlet", meta_NTOPO(ixNTOPO%isLakeInlet)%varName)
   call get_toml_val(table, "varname_islake", meta_NTOPO(ixNTOPO%islake)%varName)
   call get_toml_val(table, "varname_lakeModelType", meta_NTOPO(ixNTOPO%lakeModelType)%varName)
   call get_toml_val(table, "varname_LakeTargVol", meta_NTOPO(ixNTOPO%LakeTargVol)%varName)
   call get_toml_val(table, "varname_userTake", meta_NTOPO(ixNTOPO%userTake)%varName)
   call get_toml_val(table, "varname_goodBasin", meta_NTOPO(ixNTOPO%goodBasin)%varName)
   call get_toml_val(table, "varname_pfafCode", meta_PFAF(ixPFAF%code)%varName)
   call get_toml_val(table, "varname_D03_Coefficient", meta_SEG(ixSEG%D03_Coefficient)%varName)
   call get_toml_val(table, "varname_H06_denominator", meta_SEG(ixSEG%H06_denominator)%varName)

   ! Output history flags
   call get_toml_val(table, "basRunoff", meta_hflx(ixHFLX%basRunoff)%varFile)
   call get_toml_val(table, "instRunoff", meta_rflx(ixRFLX%instRunoff)%varFile)
   call get_toml_val(table, "dlayRunoff", meta_rflx(ixRFLX%dlayRunoff)%varFile)
   call get_toml_val(table, "sumUpstreamRunoff", meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile)
   call get_toml_val(table, "KWTroutedRunoff", meta_rflx(ixRFLX%KWTroutedRunoff)%varFile)
   call get_toml_val(table, "IRFroutedRunoff", meta_rflx(ixRFLX%IRFroutedRunoff)%varFile)
   call get_toml_val(table, "KWroutedRunoff", meta_rflx(ixRFLX%KWroutedRunoff)%varFile)
   call get_toml_val(table, "DWroutedRunoff", meta_rflx(ixRFLX%DWroutedRunoff)%varFile)
   call get_toml_val(table, "MCroutedRunoff", meta_rflx(ixRFLX%MCroutedRunoff)%varFile)
   call get_toml_val(table, "IRFvolume", meta_rflx(ixRFLX%IRFvolume)%varFile)
   call get_toml_val(table, "KWTvolume", meta_rflx(ixRFLX%KWTvolume)%varFile)
   call get_toml_val(table, "KWvolume", meta_rflx(ixRFLX%KWvolume)%varFile)
   call get_toml_val(table, "MCvolume", meta_rflx(ixRFLX%MCvolume)%varFile)
   call get_toml_val(table, "DWvolume", meta_rflx(ixRFLX%DWvolume)%varFile)
   call get_toml_val(table, "KWfloodVolume", meta_rflx(ixRFLX%KWheight)%varFile)
   call get_toml_val(table, "KWheight", meta_rflx(ixRFLX%KWheight)%varFile)
   call get_toml_val(table, "MCfloodVolume", meta_rflx(ixRFLX%MCheight)%varFile)
   call get_toml_val(table, "MCheight", meta_rflx(ixRFLX%MCheight)%varFile)
   call get_toml_val(table, "DWfloodVolume", meta_rflx(ixRFLX%DWheight)%varFile)
   call get_toml_val(table, "DWheight", meta_rflx(ixRFLX%DWheight)%varFile)
   call get_toml_val(table, "localSolute", meta_rflx(ixRFLX%localSolute)%varFile)
   call get_toml_val(table, "soluteFlux", meta_rflx(ixRFLX%DWsoluteFlux)%varFile)
   call get_toml_val(table, "soluteMass", meta_rflx(ixRFLX%DWsoluteMass)%varFile)

   call validate_and_finalize_control(err, message)

 END SUBROUTINE read_control_toml

 ! =======================================================================================================
 ! public subroutine: read legacy format control file
 ! =======================================================================================================
 SUBROUTINE read_control_legacy(ctl_fname, err, message)
   ! global vars
   USE globalData
   USE var_lookup
   USE ascii_utils, ONLY: file_open        ! open file (performs a few checks as well)
   USE ascii_utils, ONLY: get_vlines       ! get a list of character strings from non-comment lines

   implicit none
   ! argument variables
   character(*), intent(in)          :: ctl_fname               ! name of the control file
   integer(i4b),intent(out)          :: err                     ! error code
   character(*),intent(out)          :: message                 ! error message
   ! Local variables
   character(len=strLen),allocatable :: cLines(:)               ! vector of character strings
   character(len=strLen)             :: cName,cData             ! name and data from cLines(iLine)
   integer(i4b)                      :: ibeg_name               ! start index of variable name in string cLines(iLine)
   integer(i4b)                      :: iend_name               ! end index of variable name in string cLines(iLine)
   integer(i4b)                      :: iend_data               ! end index of data in string cLines(iLine)
   integer(i4b)                      :: iLine                   ! index of line in cLines
   integer(i4b)                      :: iunit                   ! file unit
   integer(i4b)                      :: io_error                ! error in I/O
   character(len=strLen)             :: cmessage                ! error message from subroutine

   err=0; message='read_control_legacy/'

   if (masterproc) then
      write(iulog,'(2a)') new_line('a'), '---- read control file --- '
   end if

   call file_open(trim(ctl_fname),iunit,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage);return;endif

   call get_vlines(iunit,cLines,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage);return;endif

   close(iunit)

 ! loop through the non-comment lines in the input file, and extract the name and the information
 do iLine=1,size(cLines)

   ! initialize io_error
   io_error=0

   ! identify start and end of the name and the data
   ibeg_name = index(cLines(iLine),'<'); if(ibeg_name==0) err=20
   iend_name = index(cLines(iLine),'>'); if(iend_name==0) err=20
   iend_data = index(cLines(iLine),'!'); if(iend_data==0) err=20
   if(err/=0)then; message=trim(message)//'problem disentangling cLines(iLine) [string='//trim(cLines(iLine))//']';return;endif

   ! extract name of the information, and the information itself
   cName = adjustl(cLines(iLine)(ibeg_name:iend_name))
   cData = adjustl(cLines(iLine)(iend_name+1:iend_data-1))
   if (masterproc) then
     write(iulog,'(1x,a,a,a)') trim(cName), ' --> ', trim(cData)
   endif

   if (index(cData, achar(9)) > 0) then
     err=10; message=trim(message)//'Control variable:'//trim(cData)//' includes TABs. Use spaces'
     return
   end if

   ! populate variables
   select case(trim(cName))

   ! DIRECTORIES
   case('<ancil_dir>');            ancil_dir   = trim(cData)                           ! directory containing ancillary data (network, mapping, namelist)
   case('<input_dir>');            input_dir   = trim(cData)                           ! directory containing input forcing netCDF, e.g. runoff
   case('<output_dir>');           output_dir  = trim(cData)                           ! directory for routed flow output (netCDF)
   case('<restart_dir>');          restart_dir = trim(cData)                           ! directory for restart output (netCDF)
   ! RUN CONTROL
   case('<case_name>');            case_name   = trim(cData)                           ! name of simulation. used as head of model output and restart file
   case('<sim_start>');            simStart    = trim(cData)                           ! date string defining the start of the simulation
   case('<sim_end>');              simEnd      = trim(cData)                           ! date string defining the end of the simulation
   case('<continue_run>');         read(cData,*,iostat=io_error) continue_run          ! logical; T-> append output in existing history files. F-> write output in new history file
   case('<route_opt>');            routOpt     = trim(cData)                           ! routing scheme options  0-> accumRunoff, 1->IRF, 2->KWT, 3-> KW, 4->MC, 5->DW
   case('<doesBasinRoute>');       read(cData,*,iostat=io_error) doesBasinRoute        ! basin routing options   0-> no, 1->IRF, otherwise error
   case('<dt_qsim>');              read(cData,*,iostat=io_error) dt                    ! time interval of the simulation [sec] (To-do: change dt to dt_sim)
   case('<floodplain>');           read(cData,*,iostat=io_error) floodplain            ! logical: floodwater is computed, otherwise, channel is unlimited bank depth
   case('<hw_drain_point>');       read(cData,*,iostat=io_error) hw_drain_point        ! integer: how to add inst. runoff in reach for headwater HRUs. 1->top of reach, 2->bottom of reach (default)
   case('<tracer>');               read(cData,*,iostat=io_error) tracer                ! logical: tracer is activated to compute a solute or constituent transport.
   case('<is_lake_sim>');          read(cData,*,iostat=io_error) is_lake_sim           ! logical; lakes are simulated
   case('<lakeRegulate>');         read(cData,*,iostat=io_error) lakeRegulate          ! logical: F -> turn all the lakes into natural (lakeType=1) regardless of lakeModelType defined individually. T (default)
   case('<LakeInputOption>');      read(cData,*,iostat=io_error) LakeInputOption       ! fluxes for lake simulation; 0->evaporation+precipitation (default), 1->runoff, 2->evaporation+precipitation+runoff
   case('<is_flux_wm>');           read(cData,*,iostat=io_error) is_flux_wm            ! logical; provided fluxes to or from seg/lakes should be considered
   case('<is_vol_wm>');            read(cData,*,iostat=io_error) is_vol_wm             ! logical; provided target volume for managed lakes are considered
   case('<is_vol_wm_jumpstart>');  read(cData,*,iostat=io_error) is_vol_wm_jumpstart   ! logical; jump to the first time step target volume is set to true
   case('<scale_factor_runoff>');  read(cData,*,iostat=io_error) scale_factor_runoff   ! float; factor to scale the runoff values
   case('<offset_value_runoff>');  read(cData,*,iostat=io_error) offset_value_runoff   ! float; offset for runoff values
   case('<scale_factor_Ep>');      read(cData,*,iostat=io_error) scale_factor_Ep       ! float; factor to scale the evaporation values
   case('<offset_value_Ep>');      read(cData,*,iostat=io_error) offset_value_Ep       ! float; offset for evaporation values
   case('<is_Ep_upward_negative>'); read(cData,*,iostat=io_error) is_Ep_upward_negative ! logical; flip evaporation in case upward direction is negative in input values convention
   case('<scale_factor_prec>');    read(cData,*,iostat=io_error) scale_factor_prec     ! float; factor to scale the precipitation values
   case('<offset_value_prec>');    read(cData,*,iostat=io_error) offset_value_prec     ! float; offset for precipitation values
   case('<min_length_route>');     read(cData,*,iostat=io_error) min_length_route      ! float; minimum reach length for routing to be performed. pass-through is performed for length less than this threshold
   ! RIVER NETWORK TOPOLOGY
   case('<fname_ntopOld>');        fname_ntopOld = trim(cData)                         ! name of file containing stream network topology information
   case('<ntopAugmentMode>');      read(cData,*,iostat=io_error) ntopAugmentMode       ! option for river network augmentation mode. terminate the program after writing augmented ntopo.
   case('<fname_ntopNew>');        fname_ntopNew = trim(cData)                         ! name of file containing stream segment information
   case('<dname_nhru>');           dname_nhru    = trim(cData)                         ! dimension name of the HRUs
   case('<dname_sseg>');           dname_sseg    = trim(cData)                         ! dimension name of the stream segments
   ! RUNOFF, EVAPORATION AND PRECIPITATION FILE
   case('<fname_qsim>');           fname_qsim   = trim(cData)                          ! name of text file listing netcdf names. netCDF include runoff, evaporation and precipitation varialbes
   case('<vname_qsim>');           vname_qsim   = trim(cData)                          ! name of runoff variable
   case('<vname_evapo>');          vname_evapo  = trim(cData)                          ! name of actual evapoartion variable
   case('<vname_precip>');         vname_precip = trim(cData)                          ! name of precipitation variable
   case('<vname_solute>');         vname_solute = trim(cData)                          ! name of solute mass flux variable
   case('<vname_time>');           vname_time   = trim(cData)                          ! name of time variable in the runoff file
   case('<vname_hruid>');          vname_hruid  = trim(cData)                          ! name of the HRU id
   case('<dname_time>');           dname_time   = trim(cData)                          ! name of time variable in the runoff file
   case('<dname_hruid>');          dname_hruid  = trim(cData)                          ! name of the HRU id dimension
   case('<dname_xlon>');           dname_xlon   = trim(cData)                          ! name of x (j,lon) dimension
   case('<dname_ylat>');           dname_ylat   = trim(cData)                          ! name of y (i,lat) dimension
   case('<units_qsim>');           units_qsim   = trim(cData)                          ! units of runoff
   case('<units_cc>');             units_cc     = trim(cData)                          ! units of concentration
   case('<dt_ro>');                read(cData,*,iostat=io_error) dt_ro                 ! time interval of the runoff data [sec]
   case('<input_fillvalue>');      read(cData,*,iostat=io_error) input_fillvalue       ! fillvalue used for input variable
   case('<ro_calendar>');          ro_calendar  = trim(cData)                          ! name of calendar used in runoff input netcdfs
   case('<ro_time_units>');        ro_time_units = trim(cData)                         ! time units used in runoff input netcdfs
   case('<ro_time_stamp>');        ro_time_stamp = trim(cData)                         ! time stamp used input - start, middle, or end, otherwise error
   case('<runoffMin>');            read(cData,*,iostat=io_error) runoffMin             ! minimum runoff volume [m3/s] from HRU.
   ! Water-management input netCDF - water abstraction/infjection or lake target volume
   case('<fname_wm>');             fname_wm        = trim(cData)                       ! name of text file containing ordered nc file names
   case('<vname_flux_wm>');        vname_flux_wm   = trim(cData)                       ! name of varibale for fluxes to and from seg (reachs/lakes)
   case('<vname_vol_wm>');         vname_vol_wm    = trim(cData)                       ! name of varibale for target volume for managed lakes
   case('<vname_time_wm>');        vname_time_wm   = trim(cData)                       ! name of time variable
   case('<vname_segid_wm>');       vname_segid_wm  = trim(cData)                       ! name of the segid varibale in nc files
   case('<dname_time_wm>');        dname_time_wm   = trim(cData)                       ! name of time dimension
   case('<dname_segid_wm>');       dname_segid_wm  = trim(cData)                       ! name of the routing HRUs dimension
   case('<dt_wm>');                read(cData,*,iostat=io_error) dt_wm                 ! time interval of the water-management data [sec]
   ! RUNOFF REMAPPING
   case('<is_remap>');             read(cData,*,iostat=io_error) is_remap              ! logical case runnoff needs to be mapped to river network HRU
   case('<fname_remap>');          fname_remap          = trim(cData)                  ! name of runoff mapping netCDF
   case('<vname_hruid_in_remap>'); vname_hruid_in_remap = trim(cData)                  ! name of variable containing ID of river network HRU
   case('<vname_weight>');         vname_weight         = trim(cData)                  ! name of variable contating areal weights of runoff HRUs within each river network HRU
   case('<vname_qhruid>');         vname_qhruid         = trim(cData)                  ! name of variable containing ID of runoff HRU
   case('<vname_num_qhru>');       vname_num_qhru       = trim(cData)                  ! name of variable containing numbers of runoff HRUs within each river network HRU
   case('<vname_i_index>');        vname_i_index        = trim(cData)                  ! name of variable containing index of xlon dimension in runoff grid (if runoff file is grid)
   case('<vname_j_index>');        vname_j_index        = trim(cData)                  ! name of variable containing index of ylat dimension in runoff grid (if runoff file is grid)
   case('<dname_hru_remap>');      dname_hru_remap      = trim(cData)                  ! name of dimension of river network HRU ID
   case('<dname_data_remap>');     dname_data_remap     = trim(cData)                  ! name of dimension of runoff HRU overlapping with river network HRU
   ! RESTART
   case('<restart_write>');        restart_write        = trim(cData)                  ! restart write option (case-insensitive): never, last, specified, yearly, monthly, or daily
   case('<restart_date>');         restart_date         = trim(cData)                  ! specified restart date, yyyy-mm-dd (hh:mm:ss) for Specified option
   case('<restart_month>');        read(cData,*,iostat=io_error) restart_month         ! restart periodic month
   case('<restart_day>');          read(cData,*,iostat=io_error) restart_day           ! restart periodic day
   case('<restart_hour>');         read(cData,*,iostat=io_error) restart_hour          ! restart periodic hour
   case('<fname_state_in>');       fname_state_in       = trim(cData)                  ! filename for the channel states
   ! SPATIAL CONSTANT PARAMETERS
   case('<param_nml>');            param_nml            = trim(cData)                  ! name of namelist including routing parameter value
   ! USER OPTIONS: Define options to include/skip calculations
   case('<qmodOption>');           read(cData,*,iostat=io_error) qmodOption            ! option for streamflow modification (DA): 0->no DA, 1->direct insertion
   case('<qBlendPeriod>');         read(cData,*,iostat=io_error) qBlendPeriod          ! number of time steps for which streamflow modification is performed through blending observation
   case('<QerrTrend>');            read(cData,*,iostat=io_error) QerrTrend             ! temporal dischargge error decreasing trend. 1->constant, 2->linear, 3->logistic, 4->exponential
   case('<hydGeometryOption>');    read(cData,*,iostat=io_error) hydGeometryOption     ! option for hydraulic geometry calculations (0=read from file, 1=compute)
   case('<topoNetworkOption>');    read(cData,*,iostat=io_error) topoNetworkOption     ! option for network topology calculations (0=read from file, 1=compute)
   case('<computeReachList>');     read(cData,*,iostat=io_error) computeReachList      ! option to compute list of upstream reaches (0=do not compute, 1=compute)
   ! GAUGE DATA
   case('<gageMetaFile>');         gageMetaFile = trim(cData)                          ! name of csv file containing gauge metadata (gauge id, reach id, gauge lat/lon)
   case('<outputAtGage>');         read(cData,*,iostat=io_error) outputAtGage          ! logical; T-> history file output at only gauge points
   case('<fname_gageObs>');        fname_gageObs = trim(cData)                         ! netcdf name of gauge observed data
   case('<vname_gageFlow>');       vname_gageFlow  = trim(cData)                       ! varialbe name for gauge site flow data
   case('<vname_gageSite>');       vname_gageSite  = trim(cData)                       ! variable name for site name data
   case('<vname_gageTime>');       vname_gageTime  = trim(cData)                       ! variable name for time data
   case('<dname_gageSite>');       dname_gageSite  = trim(cData)                       ! dimension name for gauge site
   case('<dname_gageTime>');       dname_gageTime  = trim(cData)                       ! dimension name for time
   case('<strlen_gageSite>');      read(cData,*,iostat=io_error) strlen_gageSite       ! site name max character length
   ! IO
   case('<pio_netcdf_format>');    read(cData,*,iostat=io_error) pio_netcdf_format     ! netCDF format: 64bit_offset (default) for PIO use or netCDF-4
   case('<pio_netcdf_type>');      read(cData,*,iostat=io_error) pio_typename          ! netCDF type: pnetcdf (default), netcdf, netcdf4c, or netcdf4p
   ! MISCELLANEOUS
   case('<debug>');                read(cData,*,iostat=io_error) debug                 ! print out detailed information throught the probram
   case('<seg_outlet>'   );        read(cData,*,iostat=io_error) idSegOut              ! desired outlet reach id (if -9999 --> route over the entire network)
   case('<desireId>'   );          read(cData,*,iostat=io_error) desireId              ! turn off checks or speficy reach ID if necessary to print on screen
   case('<checkMassBalance>');     read(cData,*,iostat=io_error) checkMassBalance      ! turn on/off domain wide water balance printing on screen
   ! PFAFCODE
   case('<maxPfafLen>');           read(cData,*,iostat=io_error) maxPfafLen            ! maximum digit of pfafstetter code (default 32)
   case('<pfafMissing>');          pfafMissing = trim(cData)                           ! missing pfafcode (e.g., reach without any upstream area)
   ! OUTPUT OPTIONS
   case('<time_units>');           time_units = trim(cData)                            ! time units used in history file output. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted
   case('<newFileFrequency>');     newFileFrequency = trim(cData)                      ! frequency for new history files (daily, monthly, yearly, single)
   case('<outputFrequency>');      outputFrequency  = trim(cData)                      ! output frequency (integer for multiple of simulation time step or daily, monthly or yearly)
   case('<outputNameOption>');     outputNameOption = trim(cData)                      ! option for routing method dependent output names (e.g., routedRunoff) - generic or specific (default)
   case('<histTimeStamp_offset>'); read(cData,*,iostat=io_error) histTimeStamp_offset  ! time stamp offset [second] from a start of time step
   case('<basRunoff>');            read(cData,*,iostat=io_error) meta_hflx(ixHFLX%basRunoff        )%varFile  ! default: true
   case('<instRunoff>');           read(cData,*,iostat=io_error) meta_rflx(ixRFLX%instRunoff       )%varFile  ! default: false
   case('<dlayRunoff>');           read(cData,*,iostat=io_error) meta_rflx(ixRFLX%dlayRunoff       )%varFile  ! default: false
   case('<sumUpstreamRunoff>');    read(cData,*,iostat=io_error) meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile  ! default: false
   case('<KWTroutedRunoff>');      read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWTroutedRunoff  )%varFile  ! default: true (turned off if inactive)
   case('<IRFroutedRunoff>');      read(cData,*,iostat=io_error) meta_rflx(ixRFLX%IRFroutedRunoff  )%varFile  ! default: true (turned off if inactive)
   case('<KWroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWroutedRunoff   )%varFile  ! default: true (turned off if inactive)
   case('<DWroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%DWroutedRunoff   )%varFile  ! default: true (turned off if inactive)
   case('<MCroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%MCroutedRunoff   )%varFile  ! default: true (turned off if inactive)
   case('<IRFvolume>');            read(cData,*,iostat=io_error) meta_rflx(ixRFLX%IRFvolume        )%varFile  ! default: true (turned off if inactive)
   case('<KWTvolume>');            read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWTvolume        )%varFile  ! default: true (turned off if inactive)
   case('<KWvolume>');             read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWvolume         )%varFile  ! default: true (turned off if inactive)
   case('<MCvolume>');             read(cData,*,iostat=io_error) meta_rflx(ixRFLX%MCvolume         )%varFile  ! default: true (turned off if inactive)
   case('<DWvolume>');             read(cData,*,iostat=io_error) meta_rflx(ixRFLX%DWvolume         )%varFile  ! default: true (turned off if inactive)
   case('<KWfloodVolume>');        read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWheight         )%varFile  ! default: true (turned off if floodplain is inactive)
   case('<MCfloodVolume>');        read(cData,*,iostat=io_error) meta_rflx(ixRFLX%MCheight         )%varFile  ! default: true (turned off if floodplain is inactive)
   case('<DWfloodVolume>');        read(cData,*,iostat=io_error) meta_rflx(ixRFLX%DWheight         )%varFile  ! default: true (turned off if floodplain is inactive)
   case('<KWheight>');             read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWheight         )%varFile  ! default: true (turned off if floodplain is inactive)
   case('<MCheight>');             read(cData,*,iostat=io_error) meta_rflx(ixRFLX%MCheight         )%varFile  ! default: true (turned off if floodplain is inactive)
   case('<DWheight>');             read(cData,*,iostat=io_error) meta_rflx(ixRFLX%DWheight         )%varFile  ! default: true (turned off if floodplain is inactive)
   case('<localSolute>');          read(cData,*,iostat=io_error) meta_rflx(ixRFLX%localSolute      )%varFile  ! default: false (turned off if tracer is inactive)
   case('<soluteFlux>');           read(cData,*,iostat=io_error) meta_rflx(ixRFLX%DWsoluteFlux     )%varFile  ! default: false (turned off if tracer is inactive)
   case('<soluteMass>');           read(cData,*,iostat=io_error) meta_rflx(ixRFLX%DWsoluteMass     )%varFile  ! default: false (turned off if tracer is inactive)
   case('<outputInflow>');         read(cData,*,iostat=io_error) outputInflow

   ! VARIABLE NAMES for data (overwrite default name in popMeta.f90)
   ! HRU structure
   case('<varname_area>'         ); meta_HRU    (ixHRU%area            )%varName = trim(cData) ! HRU area
   ! Mapping from HRUs to stream segments
   case('<varname_HRUid>'        ); meta_HRU2SEG(ixHRU2SEG%HRUid       )%varName = trim(cData) ! HRU id
   case('<varname_HRUindex>'     ); meta_HRU2SEG(ixHRU2SEG%HRUindex    )%varName = trim(cData) ! HRU index
   case('<varname_hruSegId>'     ); meta_HRU2SEG(ixHRU2SEG%hruSegId    )%varName = trim(cData) ! the stream segment id below each HRU
   case('<varname_hruSegIndex>'  ); meta_HRU2SEG(ixHRU2SEG%hruSegIndex )%varName = trim(cData) ! the stream segment index below each HRU
   ! reach properties
   case('<varname_length>'         ); meta_SEG    (ixSEG%length          )%varName =trim(cData)   ! length of segment  (m)
   case('<varname_slope>'          ); meta_SEG    (ixSEG%slope           )%varName =trim(cData)   ! slope of segment   (-)
   case('<varname_width>'          ); meta_SEG    (ixSEG%width           )%varName =trim(cData)   ! width of segment   (m)
   case('<varname_depth>'          ); meta_SEG    (ixSEG%depth           )%varName =trim(cData)   ! bankful depth of segment   (m)
   case('<varname_sideSlope>'      ); meta_SEG    (ixSEG%sideSlope       )%varName =trim(cData)   ! River channel side slope (-): v:h=1:sideSlope
   case('<varname_man_n>'          ); meta_SEG    (ixSEG%man_n           )%varName =trim(cData)   ! Manning's n        (weird units)
   case('<varname_floodplainSlope>'); meta_SEG    (ixSEG%floodplainSlope )%varName =trim(cData)   ! Floodplain slope (-): v:h=1:floodplainSlope
   case('<varname_hruArea>'        ); meta_SEG    (ixSEG%hruArea         )%varName =trim(cData)   ! local basin area (m2)
   case('<varname_weight>'         ); meta_SEG    (ixSEG%weight          )%varName =trim(cData)   ! HRU weight
   case('<varname_timeDelayHist>'  ); meta_SEG    (ixSEG%timeDelayHist   )%varName =trim(cData)   ! time delay histogram for each reach (s)
   case('<varname_upsArea>'        ); meta_SEG    (ixSEG%upsArea         )%varName =trim(cData)   ! area above the top of the reach -- zero if headwater (m2)
   case('<varname_basUnderLake>'   ); meta_SEG    (ixSEG%basUnderLake    )%varName =trim(cData)   ! Area of basin under lake  (m2)
   case('<varname_rchUnderLake>'   ); meta_SEG    (ixSEG%rchUnderLake    )%varName =trim(cData)   ! Length of reach under lake (m)
   case('<varname_minFlow>'        ); meta_SEG    (ixSEG%minFlow         )%varName =trim(cData)   ! minimum environmental flow
   case('<varname_D03_MaxStorage>' ); meta_SEG    (ixSEG%D03_MaxStorage  )%varName =trim(cData)   ! Doll 2006; maximum active storage for Doll 2003 formulation
   case('<varname_D03_Coefficient>'); meta_SEG    (ixSEG%D03_Coefficient )%varName =trim(cData)   ! Doll 2006; coefficient for Doll 2003 formulation (day-1)
   case('<varname_D03_Power>'      ); meta_SEG    (ixSEG%D03_Power       )%varName =trim(cData)   ! Doll 2006; power for Doll 2003 formulation
   case('<varname_D03_S0>'         ); meta_SEG    (ixSEG%D03_S0          )%varName =trim(cData)   ! Doll 2006; additional parameter to represent inactive storage

   case('<varname_HYP_E_emr>'      ); meta_SEG    (ixSEG%HYP_E_emr       )%varName =trim(cData)   ! HYPE; elevation of emergency spillway [m]
   case('<varname_HYP_E_lim>'      ); meta_SEG    (ixSEG%HYP_E_lim       )%varName =trim(cData)   ! HYPE; elevation below which primary spillway flow is restrcited [m]
   case('<varname_HYP_E_min>'      ); meta_SEG    (ixSEG%HYP_E_min       )%varName =trim(cData)   ! HYPE; elevation below which outflow is zero [m]
   case('<varname_HYP_E_zero>'     ); meta_SEG    (ixSEG%HYP_E_zero      )%varName =trim(cData)   ! HYPE; elevation at which lake/reservoir storage is zero [m]
   case('<varname_HYP_Qrate_emr>'  ); meta_SEG    (ixSEG%HYP_Qrate_emr   )%varName =trim(cData)   ! HYPE; emergency rate of flow for each unit of elevation above HYP_E_emr [m3/s]
   case('<varname_HYP_Erate_emr>'  ); meta_SEG    (ixSEG%HYP_Erate_emr   )%varName =trim(cData)   ! HYPE; power for the rate of flow for each unit of elevation above HYP_E_emr [-]
   case('<varname_HYP_Qrate_prim>' ); meta_SEG    (ixSEG%HYP_Qrate_prim  )%varName =trim(cData)   ! HYPE; the average yearly or long term output from primary spillway [m3/s]
   case('<varname_HYP_Qrate_amp>'  ); meta_SEG    (ixSEG%HYP_Qrate_amp   )%varName =trim(cData)   ! HYPE; amplitude of the Qrate_main [-]
   case('<varname_HYP_Qrate_phs>'  ); meta_SEG    (ixSEG%HYP_Qrate_phs   )%varName =trim(cData)   ! HYPE; phase of the Qrate_main based on the day of the year [-]; default 100
   case('<varname_HYP_prim_F>'     ); meta_SEG    (ixSEG%HYP_prim_F      )%varName =trim(cData)   ! HYPE; if the reservoir has a primary spillway then set to 1 otherwise 0
   case('<varname_HYP_A_avg>'      ); meta_SEG    (ixSEG%HYP_A_avg       )%varName =trim(cData)   ! HYPE; average area for the lake; this might not be used if bathymetry is provided [m]
   case('<varname_HYP_Qsim_mode>'  ); meta_SEG    (ixSEG%HYP_Qsim_mode   )%varName =trim(cData)   ! HYPE; the outflow is sum of emergency and primary spillways if 1, otherwise the maximum

   case('<varname_H06_Smax>'       ); meta_SEG    (ixSEG%H06_Smax        )%varName =trim(cData)   ! Hanasaki 2006; maximume reservoir storage [m3]
   case('<varname_H06_alpha>'      ); meta_SEG    (ixSEG%H06_alpha       )%varName =trim(cData)   ! Hanasaki 2006; fraction of active storage compared to total storage [-]
   case('<varname_H06_envfact>'    ); meta_SEG    (ixSEG%H06_envfact     )%varName =trim(cData)   ! Hanasaki 2006; fraction of inflow that can be used to meet demand [-]
   case('<varname_H06_S_ini>'      ); meta_SEG    (ixSEG%H06_S_ini       )%varName =trim(cData)   ! Hanasaki 2006; initial storage used for initial estimation of release coefficient [m3]
   case('<varname_H06_c1>'         ); meta_SEG    (ixSEG%H06_c1          )%varName =trim(cData)   ! Hanasaki 2006; coefficient 1 for target release for irrigation reseroir [-]
   case('<varname_H06_c2>'         ); meta_SEG    (ixSEG%H06_c2          )%varName =trim(cData)   ! Hanasaki 2006; coefficient 2 for target release for irrigation reseroir [-]
   case('<varname_H06_exponent>'   ); meta_SEG    (ixSEG%H06_exponent    )%varName =trim(cData)   ! Hanasaki 2006; Exponenet of actual release for "within-a-year" reservoir [-]
   case('<varname_H06_denominator>'); meta_SEG    (ixSEG%H06_denominator )%varName =trim(cData)   ! Hanasaki 2006; Denominator of actual release for "within-a-year" reservoir [-]
   case('<varname_H06_c_compare>'  ); meta_SEG    (ixSEG%H06_c_compare   )%varName =trim(cData)   ! Hanasaki 2006; Criterion for distinguish of "within-a-year" or "multi-year" reservoir [-]
   case('<varname_H06_frac_Sdead>' ); meta_SEG    (ixSEG%H06_frac_Sdead  )%varName =trim(cData)   ! Hanasaki 2006; Fraction of dead storage to maximume storage [-]
   case('<varname_H06_E_rel_ini>'  ); meta_SEG    (ixSEG%H06_E_rel_ini   )%varName =trim(cData)   ! Hanasaki 2006; Initial release coefficient [-]
   case('<varname_H06_I_Jan>'      ); meta_SEG    (ixSEG%H06_I_Jan       )%varName =trim(cData)   ! Hanasaki 2006; Average January   inflow [m3/s]
   case('<varname_H06_I_Feb>'      ); meta_SEG    (ixSEG%H06_I_Feb       )%varName =trim(cData)   ! Hanasaki 2006; Average Februrary inflow [m3/s]
   case('<varname_H06_I_Mar>'      ); meta_SEG    (ixSEG%H06_I_Mar       )%varName =trim(cData)   ! Hanasaki 2006; Average March     inflow [m3/s]
   case('<varname_H06_I_Apr>'      ); meta_SEG    (ixSEG%H06_I_Apr       )%varName =trim(cData)   ! Hanasaki 2006; Average April     inflow [m3/s]
   case('<varname_H06_I_May>'      ); meta_SEG    (ixSEG%H06_I_May       )%varName =trim(cData)   ! Hanasaki 2006; Average May       inflow [m3/s]
   case('<varname_H06_I_Jun>'      ); meta_SEG    (ixSEG%H06_I_Jun       )%varName =trim(cData)   ! Hanasaki 2006; Average June      inflow [m3/s]
   case('<varname_H06_I_Jul>'      ); meta_SEG    (ixSEG%H06_I_Jul       )%varName =trim(cData)   ! Hanasaki 2006; Average July      inflow [m3/s]
   case('<varname_H06_I_Aug>'      ); meta_SEG    (ixSEG%H06_I_Aug       )%varName =trim(cData)   ! Hanasaki 2006; Average August    inflow [m3/s]
   case('<varname_H06_I_Sep>'      ); meta_SEG    (ixSEG%H06_I_Sep       )%varName =trim(cData)   ! Hanasaki 2006; Average September inflow [m3/s]
   case('<varname_H06_I_Oct>'      ); meta_SEG    (ixSEG%H06_I_Oct       )%varName =trim(cData)   ! Hanasaki 2006; Average October   inflow [m3/s]
   case('<varname_H06_I_Nov>'      ); meta_SEG    (ixSEG%H06_I_Nov       )%varName =trim(cData)   ! Hanasaki 2006; Average November  inflow [m3/s]
   case('<varname_H06_I_Dec>'      ); meta_SEG    (ixSEG%H06_I_Dec       )%varName =trim(cData)   ! Hanasaki 2006; Average December  inflow [m3/s]
   case('<varname_H06_D_Jan>'      ); meta_SEG    (ixSEG%H06_D_Jan       )%varName =trim(cData)   ! Hanasaki 2006; Average January   demand [m3/s]
   case('<varname_H06_D_Feb>'      ); meta_SEG    (ixSEG%H06_D_Feb       )%varName =trim(cData)   ! Hanasaki 2006; Average Februrary demand [m3/s]
   case('<varname_H06_D_Mar>'      ); meta_SEG    (ixSEG%H06_D_Mar       )%varName =trim(cData)   ! Hanasaki 2006; Average March     demand [m3/s]
   case('<varname_H06_D_Apr>'      ); meta_SEG    (ixSEG%H06_D_Apr       )%varName =trim(cData)   ! Hanasaki 2006; Average April     demand [m3/s]
   case('<varname_H06_D_May>'      ); meta_SEG    (ixSEG%H06_D_May       )%varName =trim(cData)   ! Hanasaki 2006; Average May       demand [m3/s]
   case('<varname_H06_D_Jun>'      ); meta_SEG    (ixSEG%H06_D_Jun       )%varName =trim(cData)   ! Hanasaki 2006; Average June      demand [m3/s]
   case('<varname_H06_D_Jul>'      ); meta_SEG    (ixSEG%H06_D_Jul       )%varName =trim(cData)   ! Hanasaki 2006; Average July      demand [m3/s]
   case('<varname_H06_D_Aug>'      ); meta_SEG    (ixSEG%H06_D_Aug       )%varName =trim(cData)   ! Hanasaki 2006; Average August    demand [m3/s]
   case('<varname_H06_D_Sep>'      ); meta_SEG    (ixSEG%H06_D_Sep       )%varName =trim(cData)   ! Hanasaki 2006; Average September demand [m3/s]
   case('<varname_H06_D_Oct>'      ); meta_SEG    (ixSEG%H06_D_Oct       )%varName =trim(cData)   ! Hanasaki 2006; Average October   demand [m3/s]
   case('<varname_H06_D_Nov>'      ); meta_SEG    (ixSEG%H06_D_Nov       )%varName =trim(cData)   ! Hanasaki 2006; Average November  demand [m3/s]
   case('<varname_H06_D_Dec>'      ); meta_SEG    (ixSEG%H06_D_Dec       )%varName =trim(cData)   ! Hanasaki 2006; Average December  demand [m3/s]
   case('<varname_H06_purpose>'    ); meta_SEG    (ixSEG%H06_purpose     )%varName =trim(cData)   ! Hanasaki 2006; reservoir purpose; (0= non-irrigation, 1=irrigation) [-]
   case('<varname_H06_I_mem_F>'    ); meta_SEG    (ixSEG%H06_I_mem_F     )%varName =trim(cData)   ! Hanasaki 2006; Flag to transition to modelled inflow [-]
   case('<varname_H06_D_mem_F>'    ); meta_SEG    (ixSEG%H06_D_mem_F     )%varName =trim(cData)   ! Hanasaki 2006; Flag to transition to modelled/provided demand [-]
   case('<varname_H06_I_mem_L>'    ); meta_SEG    (ixSEG%H06_I_mem_L     )%varName =trim(cData)   ! Hanasaki 2006; Memory length in years for inflow [year]
   case('<varname_H06_D_mem_L>'    ); meta_SEG    (ixSEG%H06_D_mem_L     )%varName =trim(cData)   ! Hanasaki 2006; Memory length in years for demand [year]

   ! network topology
   case('<varname_hruContribIx>' ); meta_NTOPO  (ixNTOPO%hruContribIx  )%varName =trim(cData)  ! indices of the vector of HRUs that contribute flow to each segment
   case('<varname_hruContribId>' ); meta_NTOPO  (ixNTOPO%hruContribId  )%varName =trim(cData)  ! ids of the vector of HRUs that contribute flow to each segment
   case('<varname_segId>'        ); meta_NTOPO  (ixNTOPO%segId         )%varName =trim(cData)  ! unique id of each stream segment
   case('<varname_segIndex>'     ); meta_NTOPO  (ixNTOPO%segIndex      )%varName =trim(cData)  ! index of each stream segment (1, 2, 3, ..., n)
   case('<varname_downSegId>'    ); meta_NTOPO  (ixNTOPO%downSegId     )%varName =trim(cData)  ! unique id of the next downstream segment
   case('<varname_downSegIndex>' ); meta_NTOPO  (ixNTOPO%downSegIndex  )%varName =trim(cData)  ! index of downstream reach index
   case('<varname_upSegIds>'     ); meta_NTOPO  (ixNTOPO%upSegIds      )%varName =trim(cData)  ! ids for the vector of immediate upstream stream segments
   case('<varname_upSegIndices>' ); meta_NTOPO  (ixNTOPO%upSegIndices  )%varName =trim(cData)  ! indices for the vector of immediate upstream stream segments
   case('<varname_rchOrder>'     ); meta_NTOPO  (ixNTOPO%rchOrder      )%varName =trim(cData)  ! order that stream segments are processed`
   case('<varname_lakeId>'       ); meta_NTOPO  (ixNTOPO%lakeId        )%varName =trim(cData)  ! unique id of each lake in the river network
   case('<varname_lakeIndex>'    ); meta_NTOPO  (ixNTOPO%lakeIndex     )%varName =trim(cData)  ! index of each lake in the river network
   case('<varname_isLakeInlet>'  ); meta_NTOPO  (ixNTOPO%isLakeInlet   )%varName =trim(cData)  ! flag to define if reach is a lake inlet (1=inlet, 0 otherwise)
   case('<varname_islake>'       ); meta_NTOPO  (ixNTOPO%islake        )%varName =trim(cData)  ! flag to define a lake (1=lake, 0=reach)
   case('<varname_lakeModelType>'); meta_NTOPO  (ixNTOPO%lakeModelType )%varName =trim(cData)  ! defines the lake model type (1=Doll, 2=Hanasaki, 3=etc)
   case('<varname_LakeTargVol>'  ); meta_NTOPO  (ixNTOPO%LakeTargVol   )%varName =trim(cData)  ! flag to follow the provided target volume (1=yes, 0=no)
   case('<varname_userTake>'     ); meta_NTOPO  (ixNTOPO%userTake      )%varName =trim(cData)  ! flag to define if user takes water from reach (1=extract, 0 otherwise)
   case('<varname_goodBasin>'    ); meta_NTOPO  (ixNTOPO%goodBasin     )%varName =trim(cData)  ! flag to define a good basin (1=good, 0=bad)
   ! pfafstetter code
   case('<varname_pfafCode>'     ); meta_PFAF   (ixPFAF%code           )%varName =trim(cData)  ! pfafstetter code

   ! CESM coupling variables (not used for stand-alone)
   case('<qgwl_runoff_option>'   ); qgwl_runoff_option    = trim(cData)                        ! handling negative qgwl runoff: all, negative, threshold
   case('<bypass_routing_option>'); bypass_routing_option = trim(cData)                        ! routing bypass option: direct_in_place, direct_to_outlet, none
   case('<correct_area>');          read(cData,*,iostat=io_error) correct_area                 ! logical flag: area correction between model and coupler areas.
   case('<ice_runoff>');            read(cData,*,iostat=io_error) ice_runoff                   ! logical flag: true => ice sent to coupler separately from liquid, otherwise combined with liquid

   ! if not in list then keep going
   case default
    message=trim(message)//'unexpected text in control file provided: '//trim(cName) &
                         //' (note strings in control file must match the variable names in public_var.f90)'
    err=20; return
  end select

  ! check I/O error
  if(io_error/=0)then
    message=trim(message)//'problem with internal read of '//trim(cName); err=20; return
  endif

 end do  ! looping through lines in the control file

 call validate_and_finalize_control(err, message)

 END SUBROUTINE read_control_legacy

 ! =======================================================================================================
 ! private subroutine: common validation and post-processing for control variables
 ! =======================================================================================================
 SUBROUTINE validate_and_finalize_control(err, message)
   ! global vars
   USE globalData, ONLY: time_conv,length_conv   ! conversion factors
   USE globalData, ONLY: time_conv_solute        ! time conversion factor for solute mass
   USE globalData, ONLY: mass_conv_solute        ! mass conversion factors
   USE globalData, ONLY: masterproc              ! procs id and number of procs
   ! metadata structures
   USE globalData, ONLY: meta_HRU                ! HRU properties
   USE globalData, ONLY: meta_HRU2SEG            ! HRU-to-segment mapping
   USE globalData, ONLY: meta_SEG                ! stream segment properties
   USE globalData, ONLY: meta_NTOPO              ! network topology
   USE globalData, ONLY: meta_PFAF               ! pfafstetter code
   USE globalData, ONLY: meta_rflx               ! river flux variables
   USE globalData, ONLY: meta_hflx               ! river flux variables
   USE globalData, ONLY: isColdStart             ! initial river state - cold start (T) or from restart file (F)
   USE globalData, ONLY: nRoutes                 ! number of active routing methods
   USE globalData, ONLY: routeMethods            ! active routing method index and id
   USE globalData, ONLY: onRoute                 ! logical to indicate actiive routing method(s)
   USE globalData, ONLY: idxSUM,idxIRF,idxKWT, &
                         idxKW,idxMC,idxDW
   USE globalData, ONLY: runMode                 ! mizuRoute run mode: standalone, cesm-coupling
   ! index of named variables in each structure
   USE var_lookup, ONLY: ixHRU
   USE var_lookup, ONLY: ixHRU2SEG
   USE var_lookup, ONLY: ixSEG
   USE var_lookup, ONLY: ixNTOPO
   USE var_lookup, ONLY: ixPFAF
   USE var_lookup, ONLY: ixRFLX
   USE var_lookup, ONLY: ixHFLX
   USE nr_utils,    ONLY: char2int         ! convert integer number to a array containing individual digits
   USE ascii_utils, ONLY: lower           ! convert string to lower case

   implicit none
   integer(i4b), intent(out)          :: err                     ! error code
   character(*), intent(out)          :: message                 ! error message
   character(len=strLen)             :: cLength,cTime           ! length and time units
   character(len=strLen)             :: cMass                   ! mass units needed only when tracer is on
   integer(i4b)                      :: ipos                    ! index of character string
   logical(lgt)                      :: isGeneric, onlyOneRouting
   integer(i4b)                      :: iRoute                  ! loop index

   err=0; message='validate_and_finalize_control/'

   ! ---------- Perform minor processing and checking control variables ----------------------------------------

 ! ---------- directory option  ---------------------------------------------------------------------
 if (trim(restart_dir)==charMissing) then
   restart_dir = output_dir
 endif

 ! ---------- restart option  ---------------------------------------------------------------------
 if (trim(runMode)=='standalone' .or. .not. continue_run) then
   if (trim(fname_state_in)==charMissing .or. lower(trim(fname_state_in))=='none' &
      .or. lower(trim(fname_state_in))=='coldstart') then
     isColdStart=.true.
   else
     isColdStart=.false.
   end if
 end if

 ! ---------- control river network writing option  ---------------------------------------------------------------------
 ! option 1- river network subset mode (idSegOut>0):  Write the network variables read from file over only upstream network specified idSegOut
 ! option 2- river network augment mode: Write full network variables over the entire network
 ! River network subset mode turnes off augmentation mode.
 if (idSegOut>0) then
   ntopAugmentMode = .false.
 endif

 ! ---------- time variables  --------------------------------------------------------------------------------------------
 if (masterproc) then
   if (trim(runMode)=='standalone') then
     write(iulog,'(2a)') new_line('a'), '---- calendar --- '
     write(iulog,'(a)') '  calendar used for simulation and history output is the same as runoff input'
     if (trim(ro_calendar)/=charMissing) then
       write(iulog,'(a)') '  calendar used in runoff input is provided in control file: '//trim(ro_calendar)
       write(iulog,'(a)') '  However, this will be overwritten the one read from '//trim(fname_qsim)
     else
       write(iulog,'(a)') '  calendar used in runoff input will be read from '//trim(fname_qsim)
     end if
     write(iulog,'(2a)') new_line('a'), '---- time units --- '
     if (trim(time_units)/=charMissing) then
       write(iulog,'(a)') '  time_unit for history files is provided in control file: '//trim(time_units)
     else
       write(iulog,'(a)') '  time_unit for history files will be the same as runoff input'
       if (trim(ro_time_units)/=charMissing) then
         write(iulog,'(a)') '  time_unit in runoff input is provided in control file: '//trim(ro_time_units)
       else
         write(iulog,'(a)') '  time_unit used in runoff input will be read from '//trim(fname_qsim)
       end if
     end if
   end if
 end if

 ! ---------- runoff unit conversion --------------------------------------------------------------------------------------------
 if (masterproc) then
   write(iulog,'(2a)') new_line('a'), '---- runoff unit --- '
   write(iulog,'(a)') '  runoff unit is provided as: '//trim(units_qsim)
 end if
 ! find the position of the "/" character
 ipos = index(trim(units_qsim),'/')
 if(ipos==0)then
   message=trim(message)//'expect the character "/" exists in the units string [units='//trim(units_qsim)//']'
   err=80; return
 endif

 ! set runoff depth units correctly
 meta_hflx(ixHFLX%basRunoff)%varUnit = units_qsim

 ! get the length and time units
 cLength = units_qsim(1:ipos-1)
 cTime   = units_qsim(ipos+1:len_trim(units_qsim))

 ! get the conversion factor for length
 select case(trim(cLength))
   case('m');  length_conv = 1._dp
   case('mm'); length_conv = 1._dp/1000._dp
   case default
     message=trim(message)//'expect the length units of runoff to be "m" or "mm" [units='//trim(cLength)//']'
     err=81; return
 end select

 ! get the conversion factor for time
 select case(trim(cTime))
   case('d','day');          time_conv = 1._dp/secprday
   case('h','hr','hour');    time_conv = 1._dp/secprhour
   case('s','sec','second'); time_conv = 1._dp
   case default
     message=trim(message)//'expect the time units of runoff to be "day"("d"), "hour"("h") or "second"("s") &
             & [time units = '//trim(cTime)//']'
     err=81; return
 end select

 if (tracer) then
   ! find the position of the "/" character
   ipos = index(trim(units_cc),'/')
   if(ipos==0)then
     message=trim(message)//'expect the character "/" exists in the mass flux units string [units='//trim(units_qsim)//']'
     err=80; return
   endif

   ! get the mass and time units
   cMass = units_cc(1:ipos-1)
   cTime = units_cc(ipos+1:len_trim(units_cc))

   ! get the conversion factor for length
   select case(trim(cMass))
     case('mg');  mass_conv_solute = 1._dp
     case('g');   mass_conv_solute = 1000._dp
     case('kg');  mass_conv_solute = 1000000._dp
     case default
       message=trim(message)//'expect the mass units of solute to be "mg" [units='//trim(cMass)//']'
       err=81; return
   end select

   ! get the conversion factor for time
   select case(trim(cTime))
     case('d','day');          time_conv_solute = 1._dp/secprday
     case('h','hr','hour');    time_conv_solute = 1._dp/secprhour
     case('s','sec','second'); time_conv_solute = 1._dp
     case default
       message=trim(message)//'expect the time units of mass flux to be "day"("d"), "hour"("h") or "second"("s") &
             & [time units = '//trim(cTime)//']'
       err=81; return
   end select
 end if

 ! ---------- I/O time stamp -------
 if (masterproc) then
   write(iulog,'(2a)') new_line('a'), '---- input time stamp --- '
   write(iulog,'(2A)')      '  Input time stamp <ro_time_stamp>:  ', trim(ro_time_stamp)
   if (trim(ro_time_stamp)=='start' .or. trim(ro_time_stamp)=='end' .or. trim(ro_time_stamp)=='middle') then
     write(iulog,'(2A)')      '  The same time stamp is used for history output'
   else
     write(message, '(2A)') trim(message), 'ERROR: Input time stamp <ro_time_stamp> must be start, end, or middle'
     err=81; return
   end if
 end if

 ! ---------- simulation time step, output frequency, file frequency -------
 if (masterproc) then
   write(iulog,'(2a)') new_line('a'), '---- output/simulation time steps --- '
   write(iulog,'(A,F10.1)') '  simulation time step <dt_qsim>:             ', dt
   write(iulog,'(2A)')      '  history file freqeuncy <newFileFrequency>:  ', trim(newFileFrequency)
   write(iulog,'(2A)')      '  history output freqeuncy <outputFrequency>: ', trim(outputFrequency)
 end if

 ! 1. Process history output frequency
 select case(trim(outputFrequency))
   case('daily', 'monthly', 'yearly') ! do nothing
   case default
     read(outputFrequency,'(I5)',iostat=err) nOutFreq
     if (err/=0) then
       message=trim(message)//'<outputFrequency> is invalid: must be "daily", "monthly", "yearly" &
             & or positive integer (number of time steps)'
       return
     end if
     if (nOutFreq<0) then
       message=trim(message)//'<outputFrequency> is invalid: must be positive integer'; return
     end if
 end select

 ! 2. Check simulation time step
 ! 2.1 must be less than one day
 if (dt>86400._dp) then
   write(message, '(2A)') trim(message), '<dt_qsim> must be less than one-day (86400 sec)'
   err=81; return
 end if
 ! 2.2. multiple of simulation time step must be one day
 if (mod(86400._dp, dt)>0._dp) then
   write(message, '(2A)') trim(message), 'multiple of <dt_qsim> [sec] must be 86400 [sec] (one day)'
   err=81; return
 end if

 ! 3. Check output frequency against simulation time step if outputFrequency is numeric
 if (nOutFreq/=integerMissing) then
   if (mod(86400._dp, real(nOutFreq,kind=dp)*dt)>0._dp) then
     write(message, '(2A)') trim(message), 'multiple of <outputFrequency> x <dt_qsim> [sec] must be 86400 [sec] (one day)'
     err=81; return
   end if
 end if

 ! 4. Check new history frequency against output frequency
 if (trim(newFileFrequency)=='daily') then
   if (trim(outputFrequency)=='monthly' .or. trim(outputFrequency)=='yearly') then
     write(message, '(2A)') trim(message), 'you cannot output monthly or yearly output in daily file'
     err=81; return
   end if
 else if (trim(newFileFrequency)=='monthly') then
   if (trim(outputFrequency)=='yearly') then
     write(message, '(2A)') trim(message), 'you cannot output yearly output in monthly file'
     err=81; return
   end if
 end if

 ! 5. routing options
 ! ---- routing methods
 ! Assign index for each active routing method
 ! Make sure to turn off write option for routines not used
 if (trim(routOpt)=='0')then; write(iulog,'(a)') 'WARNING: routOpt=0 is accumRunoff option now.'; endif
 call char2int(trim(routOpt), routeMethods, invalid_value=0)
 nRoutes = size(routeMethods)
 onRoute = .false.
 do iRoute = 1, nRoutes
   select case(routeMethods(iRoute))
     case(accumRunoff);           idxSUM = iRoute; onRoute(accumRunoff)=.true.
     case(kinematicWaveTracking); idxKWT = iRoute; onRoute(kinematicWaveTracking)=.true.
     case(impulseResponseFunc);   idxIRF = iRoute; onRoute(impulseResponseFunc)=.true.
     case(muskingumCunge);        idxMC  = iRoute; onRoute(muskingumCunge)=.true.
     case(kinematicWave);         idxKW  = iRoute; onRoute(kinematicWave)=.true.
     case(diffusiveWave);         idxDW  = iRoute; onRoute(diffusiveWave)=.true.
     case default
       message=trim(message)//'routOpt may include invalid digits; expect digits 1-5 in routOpt'; err=81; return
   end select
 end do

 if (masterproc) then
   write(iulog,'(2a)') new_line('a'), '---- mis. routing options --- '
   if (min_length_route>0._dp)then
     write(iulog,'(a,F6.1,a)') '  pass-through is activated for <', min_length_route, ' m reaches only for IRF, KWE, MC, DW'
   end if
   if (hw_drain_point==1)then
     write(iulog,'(a)') '  Lateral flow drains at the top of headwater reaches only for IRF, KWE, MC, DW'
   end if
 end if

 ! ---- history Output variables
 if (outputInflow) then ! if outputInflow is active, turn on all the routing methdos for now
   meta_rflx(ixRFLX%KWTinflow)%varFile=.true.
   meta_rflx(ixRFLX%IRFinflow)%varFile=.true.
   meta_rflx(ixRFLX%MCinflow)%varFile=.true.
   meta_rflx(ixRFLX%KWinflow)%varFile=.true.
   meta_rflx(ixRFLX%DWinflow)%varFile=.true.
 end if
 if (.not. floodplain) then ! if floodplain is inactive, turn off history file writing
   meta_rflx(ixRFLX%MCfloodVolume)%varFile=.false.
   meta_rflx(ixRFLX%KWfloodVolume)%varFile=.false.
   meta_rflx(ixRFLX%DWfloodVolume)%varFile=.false.
   meta_rflx(ixRFLX%MCheight)%varFile=.false.
   meta_rflx(ixRFLX%KWheight)%varFile=.false.
   meta_rflx(ixRFLX%DWheight)%varFile=.false.
 end if
 ! Make sure turned off if the corresponding routing is inactive
 do iRoute = 0, nRouteMethods-1
   select case(iRoute)
     case(accumRunoff)
        if (.not. onRoute(iRoute)) then
          meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile=.false.
        end if
     case(kinematicWaveTracking)
       if (.not. onRoute(iRoute)) then
         meta_rflx(ixRFLX%KWTroutedRunoff)%varFile=.false.
         meta_rflx(ixRFLX%KWTvolume)%varFile=.false.
         meta_rflx(ixRFLX%KWTinflow)%varFile=.false.
       end if
     case(impulseResponseFunc)
        if (.not. onRoute(iRoute)) then
          meta_rflx(ixRFLX%IRFroutedRunoff)%varFile=.false.
          meta_rflx(ixRFLX%IRFvolume)%varFile=.false.
          meta_rflx(ixRFLX%IRFinflow)%varFile=.false.
       end if
     case(muskingumCunge)
       if (.not. onRoute(iRoute)) then
         meta_rflx(ixRFLX%MCroutedRunoff)%varFile=.false.
         meta_rflx(ixRFLX%MCvolume)%varFile=.false.
         meta_rflx(ixRFLX%MCfloodVolume)%varFile=.false.
         meta_rflx(ixRFLX%MCheight)%varFile=.false.
         meta_rflx(ixRFLX%MCinflow)%varFile=.false.
       end if
     case(kinematicWave)
       if (.not. onRoute(iRoute)) then
         meta_rflx(ixRFLX%KWroutedRunoff)%varFile=.false.
         meta_rflx(ixRFLX%KWvolume)%varFile=.false.
         meta_rflx(ixRFLX%KWfloodVolume)%varFile=.false.
         meta_rflx(ixRFLX%KWheight)%varFile=.false.
         meta_rflx(ixRFLX%KWinflow)%varFile=.false.
       end if
     case(diffusiveWave)
       if (.not. onRoute(iRoute)) then
         meta_rflx(ixRFLX%DWroutedRunoff)%varFile=.false.
         meta_rflx(ixRFLX%DWvolume)%varFile=.false.
         meta_rflx(ixRFLX%DWfloodVolume)%varFile=.false.
         meta_rflx(ixRFLX%DWheight)%varFile=.false.
         meta_rflx(ixRFLX%DWinflow)%varFile=.false.
       end if
     case default; message=trim(message)//'expect digits from 0 and 5'; err=81; return
   end select
 end do

 ! Control routing method dependent variable name - routedRunoff, volume, elevation etc.
 ! use generic name if outputNameOption is set to 'generic' AND only one routing method is activated w or w/o accumRunoff
 isGeneric = trim(lower(outputNameOption))=='generic'
 onlyOneRouting = ((nRoutes==2 .and. any(routeMethods==accumRunoff)) .or. nRoutes==1)
 if (onlyOneRouting .and. isGeneric) then
   do iRoute = 1, nRoutes
     select case(routeMethods(iRoute))
       case(accumRunoff)
         ! nothing to do
       case(kinematicWaveTracking)
         meta_rflx(ixRFLX%KWTroutedRunoff)%varName='RoutedRunoff'
         meta_rflx(ixRFLX%KWTvolume)%varName='volume'
         meta_rflx(ixRFLX%KWTinflow)%varName='inflow'
       case(impulseResponseFunc)
         meta_rflx(ixRFLX%IRFroutedRunoff)%varName='RoutedRunoff'
         meta_rflx(ixRFLX%IRFroutedRunoff)%varName='volume'
         meta_rflx(ixRFLX%IRFinflow)%varName='inflow'
       case(muskingumCunge)
         meta_rflx(ixRFLX%MCroutedRunoff)%varName='RoutedRunoff'
         meta_rflx(ixRFLX%MCvolume)%varName='volume'
         meta_rflx(ixRFLX%MCinflow)%varName='inflow'
       case(kinematicWave)
         meta_rflx(ixRFLX%KWroutedRunoff)%varName='RoutedRunoff'
         meta_rflx(ixRFLX%KWvolume)%varName='volume'
         meta_rflx(ixRFLX%KWinflow)%varName='inflow'
       case(diffusiveWave)
         meta_rflx(ixRFLX%DWroutedRunoff)%varName='RoutedRunoff'
         meta_rflx(ixRFLX%DWvolume)%varName='volume'
         meta_rflx(ixRFLX%DWinflow)%varName='inflow'
       case default
         message=trim(message)//'routOpt may include invalid digits; expect digits 0-5 in routOpt'; err=81; return
     end select
   end do
 end if

 ! basin runoff routing option
 if (doesBasinRoute==0) meta_rflx(ixRFLX%instRunoff)%varFile = .false.

 ! tracer mass and flux
 if (.not. tracer) then
   meta_rflx(ixRFLX%DWsoluteFlux)%varFile = .false.
   meta_rflx(ixRFLX%DWsoluteMass)%varFile = .false.
 end if
 ! Temporary tracer output option control as of Apr 18, 2025:
 ! Curently tracer is computed using discharge based on any routing methods, but only tracer with DW is properly output in history file
 ! if routing methods are not including DW and tracer is on, execution is aborted and instruct user to include 5(DW) in route_opt
 if (all(routeMethods/=diffusiveWave)) then
   if (tracer) then
     message=trim(message)//'Tracer is on w/o DW routig method. Please include 5 (==DW method) in route_opt when tracer is on.'
     err=20; return
   end if
   meta_rflx(ixRFLX%DWsoluteFlux)%varFile = .false.
   meta_rflx(ixRFLX%DWsoluteMass)%varFile = .false.
 endif

  END SUBROUTINE validate_and_finalize_control

 ! =======================================================================================================
 ! private function: validate control variable key name
 ! =======================================================================================================
 PURE LOGICAL FUNCTION is_valid_control_key(key)
   character(*), intent(in) :: key
   select case(trim(key))
   case('ancil_dir', 'input_dir', 'output_dir', 'restart_dir', &
        'case_name', 'sim_start', 'sim_end', 'continue_run', 'route_opt', &
        'doesBasinRoute', 'dt_qsim', 'floodplain', 'hw_drain_point', 'tracer', &
        'is_lake_sim', 'lakeRegulate', 'LakeInputOption', 'is_flux_wm', 'is_vol_wm', &
        'is_vol_wm_jumpstart', 'scale_factor_runoff', 'offset_value_runoff', &
        'scale_factor_Ep', 'offset_value_Ep', 'is_Ep_upward_negative', &
        'scale_factor_prec', 'offset_value_prec', 'min_length_route', &
        'fname_ntopOld', 'ntopAugmentMode', 'fname_ntopNew', 'dname_nhru', 'dname_sseg', &
        'fname_qsim', 'vname_qsim', 'vname_evapo', 'vname_precip', 'vname_solute', &
        'vname_time', 'vname_hruid', 'dname_time', 'dname_hruid', 'dname_xlon', 'dname_ylat', &
        'units_qsim', 'units_cc', 'dt_ro', 'input_fillvalue', 'ro_calendar', &
        'ro_time_units', 'ro_time_stamp', 'runoffMin', 'fname_wm', 'vname_flux_wm', &
        'vname_vol_wm', 'vname_time_wm', 'vname_segid_wm', 'dname_time_wm', &
        'dname_segid_wm', 'dt_wm', 'is_remap', 'fname_remap', 'vname_hruid_in_remap', &
        'vname_weight', 'vname_qhruid', 'vname_num_qhru', 'vname_i_index', 'vname_j_index', &
        'dname_hru_remap', 'dname_data_remap', 'restart_write', 'restart_date', &
        'restart_month', 'restart_day', 'restart_hour', 'fname_state_in', 'param_nml', &
        'qmodOption', 'qBlendPeriod', 'QerrTrend', 'hydGeometryOption', 'topoNetworkOption', &
        'computeReachList', 'gageMetaFile', 'outputAtGage', 'fname_gageObs', 'vname_gageFlow', &
        'vname_gageSite', 'vname_gageTime', 'dname_gageSite', 'dname_gageTime', 'strlen_gageSite', &
        'pio_netcdf_format', 'pio_netcdf_type', 'debug', 'seg_outlet', 'desireId', &
        'checkMassBalance', 'maxPfafLen', 'pfafMissing', 'time_units', 'newFileFrequency', &
        'outputFrequency', 'outputNameOption', 'histTimeStamp_offset', 'outputInflow', &
        'qgwl_runoff_option', 'bypass_routing_option', 'correct_area', 'ice_runoff', &
        'varname_area', 'varname_HRUid', 'varname_HRUindex', 'varname_hruSegId', &
        'varname_hruSegIndex', 'varname_length', 'varname_slope', 'varname_width', &
        'varname_depth', 'varname_sideSlope', 'varname_man_n', 'varname_floodplainSlope', &
        'varname_hruArea', 'varname_weight', 'varname_timeDelayHist', 'varname_upsArea', &
        'varname_basUnderLake', 'varname_rchUnderLake', 'varname_minFlow', 'varname_D03_MaxStorage', &
        'varname_D03_Coefficient', 'varname_D03_Power', 'varname_D03_S0', 'varname_HYP_E_emr', &
        'varname_HYP_E_lim', 'varname_HYP_E_min', 'varname_HYP_E_zero', 'varname_HYP_Qrate_emr', &
        'varname_HYP_Erate_emr', 'varname_HYP_Qrate_prim', 'varname_HYP_Qrate_amp', &
        'varname_HYP_Qrate_phs', 'varname_HYP_prim_F', 'varname_HYP_A_avg', 'varname_HYP_Qsim_mode', &
        'varname_H06_Smax', 'varname_H06_alpha', 'varname_H06_envfact', 'varname_H06_S_ini', &
        'varname_H06_c1', 'varname_H06_c2', 'varname_H06_exponent', 'varname_H06_denominator', &
        'varname_H06_c_compare', 'varname_H06_frac_Sdead', 'varname_H06_E_rel_ini', &
        'varname_H06_I_Jan', 'varname_H06_I_Feb', 'varname_H06_I_Mar', 'varname_H06_I_Apr', &
        'varname_H06_I_May', 'varname_H06_I_Jun', 'varname_H06_I_Jul', 'varname_H06_I_Aug', &
        'varname_H06_I_Sep', 'varname_H06_I_Oct', 'varname_H06_I_Nov', 'varname_H06_I_Dec', &
        'varname_H06_D_Jan', 'varname_H06_D_Feb', 'varname_H06_D_Mar', 'varname_H06_D_Apr', &
        'varname_H06_D_May', 'varname_H06_D_Jun', 'varname_H06_D_Jul', 'varname_H06_D_Aug', &
        'varname_H06_D_Sep', 'varname_H06_D_Oct', 'varname_H06_D_Nov', 'varname_H06_D_Dec', &
        'varname_H06_purpose', 'varname_H06_I_mem_F', 'varname_H06_D_mem_F', 'varname_H06_I_mem_L', &
        'varname_H06_D_mem_L', 'varname_hruContribIx', 'varname_hruContribId', 'varname_segId', &
        'varname_segIndex', 'varname_downSegId', 'varname_downSegIndex', 'varname_upSegIds', &
        'varname_upSegIndices', 'varname_rchOrder', 'varname_lakeId', 'varname_lakeIndex', &
        'varname_isLakeInlet', 'varname_islake', 'varname_lakeModelType', 'varname_LakeTargVol', &
        'varname_userTake', 'varname_goodBasin', 'varname_pfafCode', 'basRunoff', 'instRunoff', &
        'dlayRunoff', 'sumUpstreamRunoff', 'KWTroutedRunoff', 'IRFroutedRunoff', 'KWroutedRunoff', &
        'DWroutedRunoff', 'MCroutedRunoff', 'IRFvolume', 'KWTvolume', 'KWvolume', 'MCvolume', &
        'DWvolume', 'KWfloodVolume', 'KWheight', 'MCfloodVolume', 'MCheight', 'DWfloodVolume', &
        'DWheight', 'localSolute', 'soluteFlux', 'soluteMass', 'KWTinflow', 'IRFinflow', &
        'KWinflow', 'MCinflow', 'DWinflow')
      is_valid_control_key = .true.
   case default
      is_valid_control_key = .false.
   end select
 END FUNCTION is_valid_control_key

 SUBROUTINE get_toml_val_char(table, key, var)
    type(toml_table), intent(inout) :: table
    character(*), intent(in) :: key
    character(*), intent(out) :: var
    character(len=:), allocatable :: str_val
    integer :: stat
    call toml_get_value(table, key, str_val, stat=stat)
    if (stat == toml_stat%success .and. allocated(str_val)) then
       var = str_val
    endif
 END SUBROUTINE get_toml_val_char

 SUBROUTINE get_toml_val_int(table, key, var)
    type(toml_table), intent(inout) :: table
    character(*), intent(in) :: key
    integer(i4b), intent(inout) :: var
    integer :: int_val
    integer :: stat
    call toml_get_value(table, key, int_val, stat=stat)
    if (stat == toml_stat%success) then
       var = int(int_val, i4b)
    endif
 END SUBROUTINE get_toml_val_int

 SUBROUTINE get_toml_val_dp(table, key, var)
    type(toml_table), intent(inout) :: table
    character(*), intent(in) :: key
    real(dp), intent(inout) :: var
    real(dp) :: real_val
    integer :: stat
    call toml_get_value(table, key, real_val, stat=stat)
    if (stat == toml_stat%success) then
       var = real_val
    endif
 END SUBROUTINE get_toml_val_dp

 SUBROUTINE get_toml_val_sp(table, key, var)
    type(toml_table), intent(inout) :: table
    character(*), intent(in) :: key
    real(sp), intent(inout) :: var
    real(sp) :: real_val
    integer :: stat
    call toml_get_value(table, key, real_val, stat=stat)
    if (stat == toml_stat%success) then
       var = real_val
    endif
 END SUBROUTINE get_toml_val_sp

 SUBROUTINE get_toml_val_bool(table, key, var)
    type(toml_table), intent(inout) :: table
    character(*), intent(in) :: key
    logical(lgt), intent(inout) :: var
    logical :: bool_val
    integer :: stat
    call toml_get_value(table, key, bool_val, stat=stat)
    if (stat == toml_stat%success) then
       var = bool_val
    endif
 END SUBROUTINE get_toml_val_bool

END MODULE read_control_module
