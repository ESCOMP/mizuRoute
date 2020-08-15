MODULE get_runoff

USE nrtype

implicit none

private

public::get_hru_runoff

CONTAINS

 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 SUBROUTINE get_hru_runoff(ierr, message)     ! output: error control

 ! populate runoff_data with runoff values at LSM domain and at iTime step

  ! shared data
  USE public_var,  ONLY:input_dir               ! directory containing input data
  USE public_var,  ONLY:fname_qsim              ! simulated runoff netCDF name
  USE public_var,  ONLY:vname_qsim              ! varibale runoff in netCDF file
  USE public_var,  ONLY:vname_evapo             ! varibale actual evaporation in netCDF file
  USE public_var,  ONLY:vname_precip            ! varibale precipitation in netCDF file
  USE public_var,  ONLY:is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE public_var,  ONLY:is_lake_sim             ! logical whether or not lake should be simulated
  USE public_var,  ONLY:is_wm_sim               ! logical whether or not water management components should be read, abstraction, injection and target volume
  USE globalData,  ONLY:iTime_local             ! iTime index for the given netcdf file
  USE globalData,  ONLY:nHRU                    ! number of routing sub-basin
  USE globalData,  ONLY:runoff_data             ! data structure to hru runoff data
  USE globalData,  ONLY:remap_data              ! data structure to remap data
  ! subroutines
  USE read_runoff, ONLY:read_runoff_data        ! read runoff value into runoff_data data strucuture
  USE remapping,   ONLY:remap_runoff            ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE remapping,   ONLY:sort_flux               ! mapping HM runoff to river network HRU runoff (HM_HRU == RN_HRU)

  implicit none
  ! input variables: none
  ! output variables
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  character(len=strLen)         :: cmessage           ! error message from subroutine

  ierr=0; message='get_hru_runoff/'

  call infile_name(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the simulated runoff for the current time step - runoff_data%sim(:) or %sim2D(:,:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        trim(vname_qsim),                  & ! input: varname
                        iTime_local,                       & ! input: time index
                        runoff_data,                       & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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
   call sort_flux  (runoff_data%hru_id,         &
                    runoff_data%hru_ix,         &
                    runoff_data%sim,            &
                    runoff_data%basinRunoff,    &
                    ierr, cmessage)
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
    call sort_flux  (runoff_data%hru_id,        &
                     runoff_data%hru_ix,        &
                     runoff_data%sim,           &
                     runoff_data%basinEvapo,    &
                     ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   ! get the precepitation - runoff_data%sim(:) or %sim2D(:,:)
   call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                         trim(vname_precip),                & ! input: varname
                         iTime_local,                       & ! input: time index
                         runoff_data,                       & ! inout: runoff data structure
                         ierr, cmessage)                      ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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
    call sort_flux  (runoff_data%hru_id,        &
                     runoff_data%hru_ix,        &
                     runoff_data%sim,           &
                     runoff_data%basinPrecip,   &
                     ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if
  end if

 END SUBROUTINE get_hru_runoff

 ! *********************************************************************
 ! private subroutine: get the name of input file based current time
 ! *********************************************************************
 SUBROUTINE infile_name(ierr, message)  ! output

  ! Shared data
  USE public_var, ONLY: fname_qsim      ! simulated runoff netCDF name
  USE globalData, ONLY: iTime           ! time index at simulation time step
  USE globalData, ONLY: infileinfo_data ! the information of the input files
  USE globalData, ONLY: iTime_local     ! iTime index for the given netcdf file

  implicit none

  ! output:
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  integer(i4b)                             :: ix
  !character(len=strLen)                    :: cmessage         ! error message of downwind routine

  ierr=0; message='infile_name/'

  ! fast forward time to time index at simStart and save iTime and modJulday
  ixloop: do ix = 1, size(infileinfo_data) !loop over number of file
   if ((iTime >= infileinfo_data(ix)%iTimebound(1)).and.(iTime <= infileinfo_data(ix)%iTimebound(2))) then
    iTime_local = iTime - infileinfo_data(ix)%iTimebound(1) + 1
    fname_qsim = trim(infileinfo_data(ix)%infilename)
    exit ixloop
   endif
  enddo ixloop

 END SUBROUTINE infile_name

END MODULE get_runoff
