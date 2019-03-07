module get_hru_runoff

USE nrtype

implicit none

private

public::read_runoff_data

contains

 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 subroutine get_hru_runoff(iTime,          &  ! input: time index
                           ierr, message)     ! output: error control

  USE public_var,  only:input_dir               ! directory containing input data
  USE public_var,  only:fname_qsim              ! simulated runoff netCDF name
  USE public_var,  only:is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE globalData,  only:runoff_data             ! data structure to hru runoff data
  USE globalData,  only:remap_data              ! data structure to remap data

  implicit none
  ! input variables
  integer(i4b), intent(in)      :: iTime              ! index of time element
  ! output variables
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  character(len=strLen)         :: cmessage           ! error message from subroutine

  ! initialize error control
  ierr=0; message='get_runoff/'

  ! get the simulated runoff for the current time step - runoff_data%qsim(:) or %qsim2D(:,:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        iTime,                             & ! input: time index
                        runoff_data,                       & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! map simulated runoff to the basins in the river network
  if (is_remap) then
    call remap_runoff(runoff_data, remap_data, basinRunoff, ierr, cmessage)
    if(ierr/=0) call handle_err(ierr,cmessage)
  else
    basinRunoff=runoff_data%qsim
  end if

 end subroutine get_hru_runoff


end module get_hru_runoff
