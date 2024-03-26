MODULE get_runoff

USE nrtype
USE datetime_data,  ONLY: datetime       ! data type for datetime
USE dataTypes,      ONLY: map_time       ! data type for time-step mapping between two time series
USE dataTypes,      ONLY: inFileInfo     ! data type for storing the infromation of the nc files and its attributes
USE nr_utils,       ONLY: arth

implicit none

private
public::get_hru_runoff

CONTAINS

 ! *********************************************************************
 ! public subroutine: get forcing data at hru and current model time step
 ! *********************************************************************
 SUBROUTINE get_hru_runoff(ierr, message)

 ! get forcing data at current step at all the hrus

  USE public_var,           ONLY: input_dir           ! directory containing input data
  USE public_var,           ONLY: vname_qsim          ! varibale runoff in netCDF file
  USE public_var,           ONLY: vname_evapo         ! varibale actual evaporation in netCDF file
  USE public_var,           ONLY: vname_precip        ! varibale precipitation in netCDF file
  USE public_var,           ONLY: vname_flux_wm       ! varibale precipitation in netCDF file
  USE public_var,           ONLY: vname_vol_wm        ! varibale precipitation in netCDF file
  USE public_var,           ONLY: is_remap            ! logical runnoff needs to be mapped to river network HRU
  USE public_var,           ONLY: is_lake_sim         ! logical lake should be simulated
  USE public_var,           ONLY: is_flux_wm          ! logical water management components fluxes should be read
  USE public_var,           ONLY: is_vol_wm           ! logical water management components target volume should be read
  USE public_var,           ONLY: scale_factor_runoff ! float; factor to scale the runoff values
  USE public_var,           ONLY: offset_value_runoff ! float; offset for runoff values
  USE public_var,           ONLY: scale_factor_Ep     ! float; factor to scale the evaporation values
  USE public_var,           ONLY: offset_value_Ep     ! float; offset for evaporation values
  USE public_var,           ONLY: scale_factor_prec   ! float; factor to scale the precipitation values
  USE public_var,           ONLY: offset_value_prec   ! float; offset for precipitation values
  USE public_var,           ONLY: is_Ep_upward_negative ! logical; flag to flip evaporation values if convention is negative for upward fluxes
  USE public_var,           ONLY: dt_ro               ! forcing (ro,p,evap) input time step [sec]
  USE public_var,           ONLY: dt_wm               ! water-management input time step [sec]
  USE public_var,           ONLY: verySmall           ! very small values
  USE public_var,           ONLY: realMissing         ! real missing value
  USE globalData,           ONLY: iTime               ! simulation time step index
  USE globalData,           ONLY: runoff_data         ! data structure to hru runoff data
  USE globalData,           ONLY: wm_data             ! data strcuture for water management
  USE globalData,           ONLY: remap_data          ! data structure to remap data
  USE globalData,           ONLY: inFileInfo_ro       ! metadata for input files for runoff, evapo and precip
  USE globalData,           ONLY: inFileInfo_wm       ! metadata for water-management input files
  USE globalData,           ONLY: begDatetime         ! simulation begin datetime data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,           ONLY: roBegDatetime       ! forcing data start datetime data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,           ONLY: wmBegDatetime       ! water-managment data start datetime data (yyyy:mm:dd:hh:mm:sec)
  USE read_runoff,          ONLY: read_forcing_data   ! read forcing variable into data data strucuture
  USE process_remap_module, ONLY: remap_runoff        ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE process_remap_module, ONLY: sort_flux           ! mapping runoff, fluxes based on order of HRUs, Reaches in the network

  implicit none
  ! Argument variables
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  type(map_time)                :: tmap_sim_ro        ! time-steps mapping data
  type(map_time)                :: tmap_sim_wm        ! time-steps mapping data
  logical(lgt)                  :: remove_negatives   ! flag to replace the negative values to zeros
  character(len=strLen)         :: cmessage           ! error message from subroutine

  ierr=0; message='get_hru_runoff/'

  if ((abs(scale_factor_runoff)<verySmall) .and. ((abs(offset_value_runoff)<verySmall) .or. (abs(offset_value_runoff-realmissing)<verySmall))) then  ! not reading runoff data
    runoff_data%basinRunoff = 0._dp ! replacing with zeros
  else
    ! time step mapping from runoff time step to simulation time step
    call timeMap_sim_forc(tmap_sim_ro, begDatetime, roBegDatetime, dt_ro, iTime, inFileInfo_ro, ierr, cmessage)
    if(ierr/=0) then; message=trim(message)//trim(cmessage); return; endif

    ! get the simulated runoff for the current time step - runoff_data%sim(:) or %sim2D(:,:)
    call read_forcing_data(input_dir,             & ! input: directory
                          inFileInfo_ro,          & ! input: meta for forcing input files
                          vname_qsim,             & ! input: varname
                          tmap_sim_ro,            & ! input: ro-sim time mapping at current simulation step
                          runoff_data,            & ! inout: forcing data structure
                          ierr, cmessage)           ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! scale and offset
    call scale_forcing (runoff_data, scale_factor_runoff, offset_value_runoff)

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

  if ((is_flux_wm).or.(is_vol_wm.and.is_lake_sim)) then
    call timeMap_sim_forc(tmap_sim_wm, begDatetime, wmBegDatetime, dt_wm, iTime, inFileInfo_wm, ierr, cmessage)
    if(ierr/=0) then; message=trim(message)//trim(cmessage); return; endif
  end if

  ! Optional: lake module on -> read actual evaporation and preciptation
  if (is_lake_sim) then

    if ((abs(scale_factor_Ep)<verySmall) .and. ((abs(offset_value_Ep)<verySmall) .or. (abs(offset_value_Ep-realmissing)<verySmall))) then  ! not reading evaporation
      runoff_data%basinEvapo  = 0._dp ! set evaporation to zero
    else
      ! get the actual evaporation - runoff_data%sim(:) or %sim2D(:,:)
      call read_forcing_data(input_dir,              & ! input: directory
                             inFileInfo_ro,          & ! input: meta for forcing input files
                             vname_evapo,            & ! input: varname
                             tmap_sim_ro,            & ! input: ro-sim time mapping at current simulation step
                             runoff_data,            & ! inout: forcing data structure
                             ierr, cmessage)           ! output: error control
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! if evaporation are provided in negative values (convention for upward) then flip them
      if (is_Ep_upward_negative) then
        call scale_forcing (runoff_data,  -1._dp, 0._dp)
      end if

      ! scale and offset
      call scale_forcing (runoff_data,  scale_factor_Ep, offset_value_Ep)

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
      end if ! is_remap
    end if ! supress evaporation

    if ((abs(scale_factor_prec)<verySmall) .and. ((abs(offset_value_prec)<verySmall) .or. (abs(offset_value_prec-realmissing)<verySmall))) then  ! not reading precipitation
      runoff_data%basinPrecip  = 0._dp ! set precipitation to zero
    else
      ! get the precepitation - runoff_data%sim(:) or %sim2D(:,:)
      call read_forcing_data(input_dir,              & ! input: directory
                             inFileInfo_ro,          & ! input: meta for forcing input files
                             vname_precip,           & ! input: varname
                             tmap_sim_ro,            & ! input: ro-sim time mapping at current simulation step
                             runoff_data,            & ! inout: forcing data structure
                             ierr, cmessage)           ! output: error control
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! scale and offset
      call scale_forcing (runoff_data, scale_factor_prec, offset_value_prec)

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
    end if ! supress precipitation

    ! Optional: target volume based water release on -> read target lake volume
    if (is_vol_wm) then

      ! get the added or subtracted discharge from river segments
      call read_forcing_data(input_dir,          & ! input: input directory
                             inFileInfo_wm,      & ! input: meta for water-management input files
                             vname_vol_wm,       & ! input: varname
                             tmap_sim_wm,        & ! input: wm-sim time mapping at current simulation step
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
                           tmap_sim_wm,           & ! input: wm-sim time mapping at current simulation step
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

 ! *********************************************************************
 ! private subroutine: initialize simulation time
 ! *********************************************************************
 SUBROUTINE timeMap_sim_forc(tmap_sim_forc,     & ! inout: time-steps mapping data
                             startSim, startRo, & ! input: starting datetime for simulation and input
                             dt_in,             & ! input: input data time step [sec]
                             ixTime,            & ! input: simulation time index
                             inputFileInfo,     & ! input: input file metadata
                             ierr, message)

   USE public_var,        ONLY: verySmall      ! smallest real values
   USE public_var,        ONLY: secprday       ! day to second conversion factor
   USE public_var,        ONLY: dt             ! simulation time step

   implicit none
   ! Argument variables
   type(map_time),              intent(out)   :: tmap_sim_forc    ! time-steps mapping data
   type(datetime),              intent(in)    :: startSim         ! simulation start datetime
   type(datetime),              intent(in)    :: startRo          ! runoff data start datetime
   real(dp),                    intent(in)    :: dt_in            ! time step [sec]
   integer(i4b),                intent(in)    :: ixTime           ! number of simulation time steps
   type(inFileInfo),            intent(in)    :: inputFileInfo(:) ! forcing input meta data structure
   integer(i4b),                intent(out)   :: ierr             ! error code
   character(*),                intent(out)   :: message          ! error message
   ! Local variables
   character(len=strLen)                      :: cmessage         ! error message from subroutine
   real(dp), allocatable                      :: frcLapse(:)      !
   real(dp)                                   :: simLapse(2)      ! time-bounds of current simulation time step from input start datetime [sec]
   real(dp)                                   :: startRoSec       ! starting runoff time relative to starting simulation time [sec]
   real(dp)                                   :: juldayRo         ! starting julian day in runoff time step [day]
   real(dp)                                   :: juldaySim        ! starting julian day in simulation time step [day]
   integer(i4b)                               :: nRo              ! number of runoff data time steps
   integer(i4b)                               :: ctr              ! counter
   integer(i4b)                               :: nRoSub           !
   integer(i4b)                               :: iFile            ! loop index of input file
   integer(i4b)                               :: iRo              ! loop index of runoff time step
   integer(i4b)                               :: idxFront         ! index of r of which top is within ith model layer (i=1..nLyr)
   integer(i4b)                               :: idxEnd           ! index of the lowest soil layer of which bottom is within ith model layer (i=1..nLyr)

   ierr=0; message='timeMap_sim_forc/'

   call startRo%julianday(juldayRo, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   call startSim%julianday(juldaySim,  ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   nRo = sum(inputFileInfo(:)%nTime)
   allocate(frcLapse(nRo+1), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': frcLapse'; return; endif

   startRoSec = (juldaySim - juldayRo)*secprday ! day->sec

   simLapse(1) = startRoSec+dt*(ixTime-1)
   simLapse(2) = simLapse(1) + dt
   frcLapse = arth(0._dp,      dt_in, nRo+1)

   !-- Find index of runoff time period of which end is within simulation period
   ! condition: from beginning to end of runoff time, first index of time period whose end exceeds the beginning of simulation time step
     do iRo = 1,nRo
       if (simLapse(1)<frcLapse(iRo+1)) then
         idxFront = iRo; exit
       end if
     end do
   !-- Find index of runoff time period of which beginning is within simulation period
   ! condition: from beginning to end of runoff time, first index of time period whose end exceeds the beginning of simulation time step
     do iRo = 1,nRo
       if ( simLapse(2)<frcLapse(iRo+1) .or. abs(simLapse(2)-frcLapse(iRo+1))<verySmall ) then
         idxEnd = iRo; exit
       end if
     end do

   ! Error check
   if (idxFront-idxEnd>0)then;ierr=30;message=trim(message)//'index of idxFront lower than idxEnd';return;endif

   !-- Compute weight of soil layer contributing to each model layer and populate lyrmap variable
   ctr = 1
   nRoSub = idxEnd-idxFront + 1

   allocate(tmap_sim_forc%iTime(nRoSub), tmap_sim_forc%iFile(nRoSub), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//'tmap_sim_forc%iTime or iFile'; return; endif

   if (idxFront == idxEnd)then ! if simulation period is completely within runoff period - one runoff time period per simulation period
     do iFile=1,size(inputFileInfo)
       if ( idxFront >= inputFileInfo(iFile)%iTimebound(1) .and. &
            idxFront <= inputFileInfo(iFile)%iTimebound(2)) then
         tmap_sim_forc%iFile(ctr) = iFile
         tmap_sim_forc%iTime(ctr) = idxFront - inputFileInfo(iFile)%iTimebound(1) + 1
         exit
       end if
     end do
   else
     allocate(tmap_sim_forc%frac(nRoSub), stat=ierr, errmsg=cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage)//'tmap_sim_forc%frac'; return; endif

     ! loop frm the ealiest runoff time index to the latest, and compute fraction
     do iRo=idxFront, idxEnd
       do iFile=1,size(inputFileInfo)
         if ( iRo >= inputFileInfo(iFile)%iTimebound(1) .and. &
              iRo <= inputFileInfo(iFile)%iTimebound(2)) then
           tmap_sim_forc%iFile(ctr) = iFile
           tmap_sim_forc%iTime(ctr) = iRo- inputFileInfo(iFile)%iTimebound(1) + 1
           exit
         end if
       end do

       if (iRo == idxFront)then      ! front side of simulation time step
         tmap_sim_forc%frac(ctr)   = (frcLapse(iRo+1) - simLapse(1))/dt
       else if ( iRo == idxEnd ) then  ! end side of simulation time step
         tmap_sim_forc%frac(ctr)   = (simLapse(2)-frcLapse(iRo))/dt
       else                                    ! for soil layers that completely in model layer
         tmap_sim_forc%frac(ctr)   = (frcLapse(iRo+1)-frcLapse(iRo))/dt
       endif
       ctr = ctr+1
     end do
   end if

 END SUBROUTINE timeMap_sim_forc


 ! *********************************************************************
 ! private subroutine: scale and offset of forcing input
 ! *********************************************************************
 SUBROUTINE scale_forcing (forcing_data_in,      & ! inout: array of forcing, runoff, evaporation or precipitation
                           scale,                & ! in: scale factor
                           offset)                 ! in: offset values

  USE public_var,     ONLY: verySmall          ! directory containing input data
  USE public_var,     ONLY: realMissing        ! real missing value
  USE dataTypes,      ONLY: inputData          ! input data class keeps runoff, evaporation, and precipitation

  implicit none
  ! input, output
  class(inputData), intent(inout)   :: forcing_data_in    ! forcing data structure
  REAL(dp), INTENT(in)              :: scale, offset

  ! local varibales
  real(dp)                          :: scale_local
  real(dp)                          :: offset_local

  ! assign the inputs to local variables
  scale_local   = scale
  offset_local  = offset

  ! one of the values of scale and offset are provided and is not real missing
  if ( (verySmall < abs(scale_local-realMissing)) .or. (verySmall < abs(offset_local-realMissing)) ) then

    ! scale is not provided and is set to one
    if (abs(scale_local-realMissing)<verySmall) then
      scale_local = 1.0_dp
    end if

    ! offset is not provided and is set to zero
    if (abs(offset_local-realMissing)<verySmall) then
      offset_local = 0.0_dp
    end if

    if (allocated(forcing_data_in%sim)) then
      where (verySmall < abs(forcing_data_in%sim - realMissing))
        forcing_data_in%sim = scale_local * forcing_data_in%sim + offset_local
      end where
    end if

    if (allocated(forcing_data_in%sim2d)) then
      where (verySmall < abs(forcing_data_in%sim2d - realMissing))
        forcing_data_in%sim2d = scale_local * forcing_data_in%sim2d + offset_local
      end where
    end if

  end if

END SUBROUTINE scale_forcing


END MODULE get_runoff
