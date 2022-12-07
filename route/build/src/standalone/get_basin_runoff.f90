MODULE get_runoff

USE nrtype

implicit none

private
public::get_hru_runoff

CONTAINS

 ! *********************************************************************
 ! public subroutine: get forcing data at hru and current model time step
 ! *********************************************************************
 SUBROUTINE get_hru_runoff(ierr, message)

 ! get forcing data at current step at all the hrus

  USE public_var,           ONLY: input_dir         ! directory containing input data
  USE public_var,           ONLY: vname_qsim        ! varibale runoff in netCDF file
  USE public_var,           ONLY: vname_evapo       ! varibale actual evaporation in netCDF file
  USE public_var,           ONLY: vname_precip      ! varibale precipitation in netCDF file
  USE public_var,           ONLY: vname_flux_wm     ! varibale precipitation in netCDF file
  USE public_var,           ONLY: vname_vol_wm      ! varibale precipitation in netCDF file
  USE public_var,           ONLY: is_remap          ! logical runnoff needs to be mapped to river network HRU
  USE public_var,           ONLY: is_lake_sim       ! logical lake should be simulated
  USE public_var,           ONLY: is_flux_wm        ! logical water management components fluxes should be read
  USE public_var,           ONLY: is_vol_wm         ! logical water management components target volume should be read
  USE public_var,           ONLY: suppress_runoff   ! logical suppress the read runoff to zero (0)
  USE public_var,           ONLY: suppress_P_Ep     ! logical suppress the read precipitation/evaporation to zero (0)
  USE globalData,           ONLY: iTime             ! simulation time step index
  USE globalData,           ONLY: runoff_data       ! data structure to hru runoff data
  USE globalData,           ONLY: wm_data           ! data strcuture for water management
  USE globalData,           ONLY: remap_data        ! data structure to remap data
  USE globalData,           ONLY: tmap_sim_ro       ! time-step mapping betweein simulation and runoff
  USE globalData,           ONLY: tmap_sim_wm       ! time-step mapping betweein simulation and water management
  USE globalData,           ONLY: inFileInfo_ro     ! metadata for input files for runoff, evapo and precip
  USE globalData,           ONLY: inFileInfo_wm     ! metadata for water-management input files
  USE read_runoff,          ONLY: read_forcing_data ! read forcing variable into data data strucuture
  USE process_remap_module, ONLY: remap_runoff      ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE process_remap_module, ONLY: sort_flux         ! mapping runoff, fluxes based on order of HRUs, Reaches in the network

  implicit none
  ! Argument variables
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  logical(lgt)                  :: remove_negatives   ! flag to replace the negative values to zeros
  character(len=strLen)         :: cmessage           ! error message from subroutine

  ierr=0; message='get_hru_runoff/'

  if (suppress_runoff) then  ! not reading runoff data
    runoff_data%basinRunoff = 0._dp
  else
    ! get the simulated runoff for the current time step - runoff_data%sim(:) or %sim2D(:,:)
    call read_forcing_data(input_dir,              & ! input: directory
                          inFileInfo_ro,          & ! input: meta for forcing input files
                          vname_qsim,             & ! input: varname
                          tmap_sim_ro(iTime),     & ! input: ro-sim time mapping at current simulation step
                          runoff_data,            & ! inout: forcing data structure
                          ierr, cmessage)           ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Get river network HRU runoff into runoff_data data structure
    if (is_remap) then ! remap LSM simulated flux to the HRUs in the river network
      call remap_runoff(runoff_data, remap_data, runoff_data%basinRunoff, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    else ! runoff is already remapped to river network HRUs
      remove_negatives = .true.
      call sort_flux (runoff_data%hru_id,         &
                      runoff_data%hru_ix,         &
                      runoff_data%sim,            &
                      remove_negatives,           &
                      runoff_data%basinRunoff,    &
                      ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end if
  end if

  ! Optional: lake module on -> read actual evaporation and preciptation
  if (is_lake_sim) then

    if (suppress_P_Ep) then
      ! suppressing the read runoff or precipitation and evaporation
      runoff_data%basinPrecip = 0._dp
      runoff_data%basinEvapo  = 0._dp
    else
      ! get the actual evaporation - runoff_data%sim(:) or %sim2D(:,:)
      call read_forcing_data(input_dir,              & ! input: directory
                             inFileInfo_ro,          & ! input: meta for forcing input files
                             vname_evapo,            & ! input: varname
                             tmap_sim_ro(iTime),     & ! input: ro-sim time mapping at current simulation step
                             runoff_data,            & ! inout: forcing data structure
                             ierr, cmessage)           ! output: error control
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! Get river network HRU runoff into runoff_data data structure
      if (is_remap) then ! remap LSM simulated flux to the HRUs in the river network
        call remap_runoff(runoff_data, remap_data, runoff_data%basinEvapo, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      else ! runoff is already remapped to river network HRUs
        remove_negatives = .true.
        call sort_flux  (runoff_data%hru_id,        &
                         runoff_data%hru_ix,        &
                         runoff_data%sim,           &
                         remove_negatives,          &
                         runoff_data%basinEvapo,    &
                         ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

      ! get the precepitation - runoff_data%sim(:) or %sim2D(:,:)
      call read_forcing_data(input_dir,              & ! input: directory
                             inFileInfo_ro,          & ! input: meta for forcing input files
                             vname_precip,           & ! input: varname
                             tmap_sim_ro(iTime),     & ! input: ro-sim time mapping at current simulation step
                             runoff_data,            & ! inout: forcing data structure
                             ierr, cmessage)           ! output: error control
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! Get river network HRU runoff into runoff_data data structure
      if (is_remap) then ! remap LSM simulated flux to the HRUs in the river network
        call remap_runoff(runoff_data, remap_data, runoff_data%basinPrecip, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      else ! runoff is already remapped to river network HRUs
        remove_negatives = .true.
        call sort_flux  (runoff_data%hru_id,        &
                         runoff_data%hru_ix,        &
                         runoff_data%sim,           &
                         remove_negatives,          &
                         runoff_data%basinPrecip,   &
                         ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if ! is_remap
    end if ! suppress_P_Ep

    ! Optional: target volume based water release on -> read target lake volume
    if (is_vol_wm) then

      ! get the added or subtracted discharge from river segments
      call read_forcing_data(input_dir,          & ! input: input directory
                             inFileInfo_wm,      & ! input: meta for water-management input files
                             vname_vol_wm,       & ! input: varname
                             tmap_sim_wm(iTime), & ! input: wm-sim time mapping at current simulation step
                             wm_data,            & ! inout: water management data structure
                             ierr, cmessage)       ! output: error control
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! sorting wm_data%sim(:) into wm_data%vol_wm
      remove_negatives = .true.
      call sort_flux  (wm_data%seg_id,        &
                       wm_data%seg_ix,        &
                       wm_data%sim,           &
                       remove_negatives,      &
                       wm_data%vol_wm,        &
                       ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end if
  end if ! is_lake_sim

  ! Optional: reading water abstraction and subtraction
  if (is_flux_wm) then
    ! get the added or subtracted discharge from river segments
    call read_forcing_data(input_dir,             & ! input: input directory
                           inFileInfo_wm,         & ! input: meta for water-management input files
                           vname_flux_wm,         & ! input: varname
                           tmap_sim_wm(iTime),    & ! input: wm-sim time mapping at current simulation step
                           wm_data,               & ! inout: water management data structure
                           ierr, cmessage)          ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! sorting wm_data%sim(:) into wm_data%flux_wm
    remove_negatives = .false.
    call sort_flux  (wm_data%seg_id,        &
                     wm_data%seg_ix,        &
                     wm_data%sim,           &
                     remove_negatives,      &
                     wm_data%flux_wm,       &
                     ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

 END SUBROUTINE get_hru_runoff

END MODULE get_runoff
