module get_runoff

USE nrtype

implicit none

private

public::get_hru_runoff

contains

 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 subroutine get_hru_runoff(ierr, message)     ! output: error control

 ! populate runoff_data with runoff values at LSM domain and at iTime step

  ! shared data
  USE public_var,  only:input_dir               ! directory containing input data
  USE public_var,  only:fname_qsim              ! simulated runoff netCDF name
  USE public_var,  only:fname_wm                ! abstraction and injection netCDF name
  USE public_var,  only:vname_qsim              ! varibale runoff in netCDF file
  USE public_var,  only:vname_evapo             ! varibale actual evaporation in netCDF file
  USE public_var,  only:vname_precip            ! varibale precipitation in netCDF file
  USE public_var,  only:vname_AbsInj            ! varibale precipitation in netCDF file
  USE public_var,  only:vname_TargVol           ! varibale precipitation in netCDF file
  USE public_var,  only:is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE public_var,  only:is_lake_sim             ! logical whether or not lake should be simulated
  USE public_var,  only:is_AbsInj               ! logical whether or not abstraction, injection are active
  USE public_var,  only:is_TargVol              ! logical whether or not target volume is provided for the lakes
  USE globalData,  only:iTime_local             ! iTime index for the given netcdf file
  USE globalData,  only:iTime_local_wm          ! iTime index for the given netcdf file
  USE globalData,  only:nHRU                    ! number of routing sub-basin
  USE globalData,  only:runoff_data             ! data structure to hru runoff data
  USE globalData,  only:remap_data              ! data structure to remap data
  ! subroutines
  USE read_runoff, only:read_runoff_data        ! read runoff value into runoff_data data strucuture
  USE remapping,   only:remap_runoff            ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE remapping,   only:sort_runoff             ! mapping HM runoff to river network HRU runoff (HM_HRU == RN_HRU)

  implicit none
  ! input variables: none
  ! output variables
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  character(len=strLen)         :: cmessage           ! error message from subroutine

  ! initialize error control
  ierr=0; message='get_hru_runoff/'

  print*, "inside get basin runoff"
  print*, iTime_local
  print*, iTime_local_wm
  print*, fname_qsim
  print*, fname_wm

  ! get the simulated runoff for the current time step - runoff_data%sim(:) or %sim2D(:,:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        trim(vname_qsim),                  & ! input: varname
                        iTime_local,                       & ! input: time index
                        runoff_data,                       & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize runoff_data%basinRunoff
  if ( allocated(runoff_data%basinRunoff) ) then
    deallocate(runoff_data%basinRunoff, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinRunoff(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! Get river network HRU runoff into runoff_data data structure
  if (is_remap) then ! remap LSM simulated flux to the HRUs in the river network
   call remap_runoff(runoff_data, remap_data, runoff_data%basinRunoff, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  else ! runoff is already remapped to river network HRUs
   call sort_runoff(runoff_data, runoff_data%basinRunoff, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

  if (is_lake_sim) then ! if is_lake_sim if true then read actual evaporation and preciptation

   ! get the actual evaporation - runoff_data%sim(:) or %sim2D(:,:)
   call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                         trim(vname_evapo),                 & ! input: varname
                         iTime_local,                       & ! input: time index
                         runoff_data,                       & ! inout: runoff data structure
                         ierr, cmessage)                      ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! initialize runoff_data%basinEvapo
   if ( allocated(runoff_data%basinEvapo) ) then
     deallocate(runoff_data%basinEvapo, stat=ierr)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if
   allocate(runoff_data%basinEvapo(nHRU), stat=ierr)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! Get river network HRU runoff into runoff_data data structure
   if (is_remap) then ! remap LSM simulated flux to the HRUs in the river network
    call remap_runoff(runoff_data, remap_data, runoff_data%basinEvapo, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else ! runoff is already remapped to river network HRUs
    call sort_runoff(runoff_data, runoff_data%basinEvapo, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if


   ! get the precepitation - runoff_data%sim(:) or %sim2D(:,:)
   call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                         trim(vname_precip),                & ! input: varname
                         iTime_local,                       & ! input: time index
                         runoff_data,                       & ! inout: runoff data structure
                         ierr, cmessage)                      ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! initialize runoff_data%basinPrecip
   if ( allocated(runoff_data%basinPrecip) ) then
    deallocate(runoff_data%basinPrecip, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if
   allocate(runoff_data%basinPrecip(nHRU), stat=ierr)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! Get river network HRU runoff into runoff_data data structure
   if (is_remap) then ! remap LSM simulated flux to the HRUs in the river network
    call remap_runoff(runoff_data, remap_data, runoff_data%basinPrecip, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else ! runoff is already remapped to river network HRUs
    call sort_runoff(runoff_data, runoff_data%basinPrecip, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if
  end if

  stop

 end subroutine get_hru_runoff

end module get_runoff
