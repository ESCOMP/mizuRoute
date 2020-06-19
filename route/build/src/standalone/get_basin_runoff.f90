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

  ! shared data strcuture
  USE DataTypes,   only:runoff_temp             ! temporary data strcture to pass the read variable to runoff_data
  ! shared data
  USE public_var,  only:input_dir               ! directory containing input data
  USE public_var,  only:fname_qsim              ! simulated runoff netCDF name
  USE public_var,  only:vname_qsim              ! varibale runoff in netCDF file
  USE public_var,  only:is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE globalData,  only:iTime_local             ! iTime index for the given netcdf file
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
  type(runoff_temp)             :: runoff_data_temp   ! temporary runoff data to pass to the model
  character(len=strLen)         :: cmessage           ! error message from subroutine

  ! initialize error control
  ierr=0; message='get_hru_runoff/'

  ! populate the information from runoff_data to runoff_data_temp assuming that that runoff
  runoff_data_temp%nSpace(:) = runoff_data%nSpace(:) ! pass the dimenstion of varibale to the data, this can be better written if evapo and rain are per lake and not model grid
  runoff_data_temp%hru_id = runoff_data%hru_id ! should runoff_data_temp%hru_id be allocated first?
  !print*, "runoff_data%hru_id", runoff_data%hru_id
  !print*, "runoff_data_temp%hru_id", runoff_data_temp%hru_id
  !runoff_data_temp%hru_ix(:) = runoff_data%hru_ix(:) !

  ! get the first simulated flux, runoff, for the current time step - runoff_data_temp%sim(:) or %sim2D(:,:)

  ! get the simulated runoff for the current time step - runoff_data_temp%sim(:) or %sim2D(:,:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        trim(vname_qsim),                  & ! input: varname
                        iTime_local,                       & ! input: time index
                        runoff_data_temp,                  & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize runoff_data_temp%basinSim
  if ( allocated(runoff_data_temp%basinSim) ) then
    deallocate(runoff_data_temp%basinSim, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data_temp%basinSim(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! Get river network HRU runoff into runoff_data data structure
  if (is_remap) then ! remap LSM simulated flux to the HRUs in the river network
   call remap_runoff(runoff_data_temp, remap_data, runoff_data_temp%basinSim, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  else ! runoff is already remapped to river network HRUs
   call sort_runoff(runoff_data_temp, runoff_data_temp%basinSim, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

  ! pass that to the model simulation
  if ( allocated(runoff_data%basinRunoff) ) then
    deallocate(runoff_data%basinRunoff, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinRunoff(nHRU), stat=ierr)
  runoff_data%basinRunoff = runoff_data_temp%basinSim ! pass simulated runoff to the runoff_data

 end subroutine get_hru_runoff

end module get_runoff