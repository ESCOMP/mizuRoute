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
  USE public_var,  only:is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE globalData,  only:iTime
  USE globalData,  only:nHRU
  USE globalData,  only:runoff_data             ! data structure to hru runoff data
  USE globalData,  only:remap_data              ! data structure to remap data
  ! subroutines
  USE read_runoff, only:read_runoff_data        ! read runoff value into runoff_data data strucuture
  USE remapping,   only:remap_runoff            ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE remapping,   only:sort_runoff             ! mapping HM runoff to river network HRU runoff (HM_HRU == RN_HRU)

  implicit none
  ! input variables: none
  ! output variables
  integer(i4b), intent(out)     :: ierr                                 ! error code
  character(*), intent(out)     :: message                              ! error message
  ! local variables
  real(dp)    , allocatable     :: basinRunoff(:)                       ! basin runoff (m/s)
  real(dp)    , allocatable     :: basinEvapo(:)                        ! basin runoff (m/s)
  real(dp)    , allocatable     :: basinPrecip(:)                       ! basin runoff (m/s)
  character(len=strLen)         :: runoff = 'runoff'                    ! flag in case the flux is runoff
  character(len=strLen)         :: evaporation = 'evaporation'          ! flag in case the flux is evaporation
  character(len=strLen)         :: precipitation = 'precipitation'      ! flag in case the flux is precipitation
  character(len=strLen)         :: cmessage                             ! error message from subroutine

  ! initialize error control
  ierr=0; message='get_hru_runoff/'

  ! get the simulated runoff for the current time step - runoff_data%qsim(:), %qsim2D(:,:), easim(:), easim2d(:,:), precip(:) and precip2d(:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        iTime,                             & ! input: time index
                        runoff_data,                       & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control

  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate basinRunoff (local array)
  allocate(basinRunoff(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize runoff_data%basinRunoff
  if ( allocated(runoff_data%basinRunoff) ) then
    deallocate(runoff_data%basinRunoff, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinRunoff(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate basinEvapo (local array)
  allocate(basinEvapo(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize runoff_data%basinEvapo
  if ( allocated(runoff_data%basinEvapo) ) then
    deallocate(runoff_data%basinEvapo, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinEvapo(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate basinPrecip (local array)
  allocate(basinPrecip(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize runoff_data%basinPrecip
  if ( allocated(runoff_data%basinPrecip) ) then
    deallocate(runoff_data%basinPrecip, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinPrecip(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif


  ! Get river network HRU runoff into runoff_data data structure
  if (is_remap) then ! remap LSM simulated runoff to the HRUs in the river network

   call remap_runoff(runoff_data, remap_data, runoff,        basinRunoff, ierr, cmessage) ! if flage runoff then it uses qsim or qsim2d
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   call remap_runoff(runoff_data, remap_data, evaporation,   basinEvapo, ierr, cmessage) ! if flag evaporation then it uses easim or easim2d
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   call remap_runoff(runoff_data, remap_data, Precipitation, basinPrecip, ierr, cmessage) ! if flag precipitation then it uses precip or precip2d
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  else ! runoff is already remapped to river network HRUs

   call sort_runoff(runoff_data, runoff, basinRunoff, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call sort_runoff(runoff_data, evaporation, basinEvapo,  ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call sort_runoff(runoff_data, precipitation, basinPrecip, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  end if

  runoff_data%basinRunoff = basinRunoff
  runoff_data%basinEvapo = basinEvapo
  runoff_data%basinPrecip = basinPrecip

 end subroutine get_hru_runoff

end module get_runoff
