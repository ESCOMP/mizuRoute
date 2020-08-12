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
  USE public_var,  ONLY:fname_wm                ! flux and vol netCDF name
  USE public_var,  ONLY:vname_qsim              ! varibale runoff in netCDF file
  USE public_var,  ONLY:vname_evapo             ! varibale actual evaporation in netCDF file
  USE public_var,  ONLY:vname_precip            ! varibale precipitation in netCDF file
  USE public_var,  ONLY:vname_flux_wm           ! varibale precipitation in netCDF file
  USE public_var,  ONLY:vname_vol_wm            ! varibale precipitation in netCDF file
  USE public_var,  ONLY:is_remap                ! logical runnoff needs to be mapped to river network HRU
  USE public_var,  ONLY:is_lake_sim             ! logical lake should be simulated
  USE public_var,  ONLY:is_flux_wm              ! logical water management components fluxes should be read
  USE public_var,  ONLY:is_vol_wm               ! logical water management components target volume should be read
  USE globalData,  ONLY:basinID                 ! basin ID
  USE globalData,  ONLY:reachID                 ! reach ID
  USE globalData,  ONLY:iTime_local             ! iTime index for the given netcdf file
  USE globalData,  ONLY:iTime_local_wm          ! iTime index for the given netcdf file
  USE globalData,  ONLY:nHRU                    ! number of routing sub-basin
  USE globalData,  ONLY:nRch                    ! number of routing seg (reaches and lakes)
  USE globalData,  ONLY:runoff_data             ! data structure to hru runoff data
  USE globalData,  ONLY:wm_data                 ! data structure to hru runoff data
  USE globalData,  ONLY:remap_data              ! data structure to remap data
  ! subroutines
  USE read_runoff, ONLY:read_runoff_data        ! read runoff value into runoff_data data strucuture
  USE remapping,   ONLY:remap_runoff            ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE remapping,   ONLY:sort_flux               ! mapping runoff, fluxes based on order of HRUs, Reaches in the network

  implicit none
  ! input variables: none
  ! output variables
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  character(len=strLen)         :: cmessage           ! error message from subroutine

  ierr=0; message='get_hru_runoff/'

  call infile_name(ierr, cmessage) ! read the infile name for given iTime
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the simulated runoff for the current time step - runoff_data%sim(:) or %sim2D(:,:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        trim(vname_qsim),                  & ! input: varname
                        iTime_local,                       & ! input: time index
                        runoff_data%nSpace,                & ! inout: runoff data structure
                        runoff_data%sim,                & ! inout: runoff data structure
                        runoff_data%sim2D,                & ! inout: runoff data structure
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
                         runoff_data%nSpace,                & ! inout: runoff data structure
                        runoff_data%sim,                & ! inout: runoff data structure
                        runoff_data%sim2D,                & ! inout: runoff data structure
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
                         runoff_data%nSpace,                & ! inout: runoff data structure
                        runoff_data%sim,                & ! inout: runoff data structure
                        runoff_data%sim2D,                & ! inout: runoff data structure
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
  USE public_var, ONLY: fname_qsim             ! simulated runoff netCDF name
  USE public_var, ONLY: fname_wm               ! water management netCDF name
  USE public_var, ONLY: is_flux_wm             ! logical whether or not abstraction or injection should be read
  USE public_var, ONLY: is_vol_wm              ! logical whether or not target volume should be read
  USE globalData, ONLY: iTime                  ! time index at simulation time step
  USE globalData, ONLY: infileinfo_data        ! the information of the input files for runoff, evapo and precip
  USE globalData, ONLY: infileinfo_data_wm     ! the information of the input files
  USE globalData, ONLY: iTime_local            ! iTime index for the given netcdf file
  USE globalData, ONLY: iTime_local_wm         ! iTime index for the given netcdf file

  implicit none

  ! output:
  integer(i4b),              intent(out)    :: ierr             ! error code
  character(*),              intent(out)    :: message          ! error message
  ! local variable
  integer(i4b)                              :: ix
  integer(i4b)                              :: counter          ! to check if both iTime_local are read properly
  !character(len=strLen)                    :: cmessage         ! error message

  ! initialize error control
  ierr=0; message='infile_name/'; counter = 0 ;

  ! fast forward time to time index at simStart and save iTime and modJulday
  ixloop: do ix = 1, size(infileinfo_data) !loop over number of file
   if ((iTime >= infileinfo_data(ix)%iTimebound(1)).and.(iTime <= infileinfo_data(ix)%iTimebound(2))) then
    iTime_local = iTime - infileinfo_data(ix)%iTimebound(1) + 1
    fname_qsim = trim(infileinfo_data(ix)%infilename)
    counter = counter + 1
    exit ixloop
   endif
  enddo ixloop

  ! fast forward time to time index at simStart and save iTime and modJulday for water management nc file
  if ((is_flux_wm).or.(is_vol_wm)) then
    iyloop: do ix = 1, size(infileinfo_data_wm) !loop over number of file
     if ((iTime >= infileinfo_data_wm(ix)%iTimebound(1)).and.(iTime <= infileinfo_data_wm(ix)%iTimebound(2))) then
      iTime_local_wm = iTime - infileinfo_data_wm(ix)%iTimebound(1) + 1
      fname_wm = trim(infileinfo_data_wm(ix)%infilename)
      counter = counter + 1
      exit iyloop
     endif
    enddo iyloop
  endif

  ! check if the two files are identified in case is flux and vol flags are set to true
  if ((counter /= 2).and.((is_flux_wm).or.(is_vol_wm))) then
    ierr=20; message=trim(message)//'iTime local is out of bound for the netcdf file inputs based on given simulation time'; print*, ierr ; print*, message ; return ;
  endif

  !print*, counter
  !print*, ierr
  !print*, message
  !print*, infileinfo_data_wm(1)%iTimebound(1)
  !print*, infileinfo_data_wm(1)%iTimebound(2)
  !print*, infileinfo_data(1)%iTimebound(1)
  !print*, infileinfo_data(1)%iTimebound(2)

 END SUBROUTINE infile_name

END MODULE get_runoff
