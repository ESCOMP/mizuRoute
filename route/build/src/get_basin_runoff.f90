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
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  character(len=strLen)         :: cmessage           ! error message from subroutine
  ! timing
  integer*8                     :: startTime,endTime,cr ! star and end time stamp, rate
  real(dp)                      :: elapsedTime          ! elapsed time for the process

  ! initialize error control
  ierr=0; message='get_hru_runoff/'
  call system_clock(count_rate=cr)

  call system_clock(startTime)
  ! get the simulated runoff for the current time step - runoff_data%qsim(:) or %qsim2D(:,:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        iTime,                             & ! input: time index
                        runoff_data,                       & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  call system_clock(endTime)
  elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
!  write(*,"(A,1PG15.7,A)") '  elapsed-time [runoff_input/read] = ', elapsedTime, ' s'


  ! initialize runoff_data%basinRunoff
  if ( allocated(runoff_data%basinRunoff) ) then
    deallocate(runoff_data%basinRunoff, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinRunoff(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call system_clock(startTime)
  ! Get river network HRU runoff into runoff_data data structure
  if (is_remap) then ! remap LSM simulated runoff to the HRUs in the river network

   call remap_runoff(runoff_data, remap_data, runoff_data%basinRunoff, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  else ! runoff is already remapped to river network HRUs

   call sort_runoff(runoff_data, runoff_data%basinRunoff, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  end if
  call system_clock(endTime)
  elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
!  write(*,"(A,1PG15.7,A)") '  elapsed-time [runoff_input/remap] = ', elapsedTime, ' s'


 end subroutine get_hru_runoff

end module get_runoff
