MODULE get_runoff

USE nrtype

implicit none

private
public::get_hru_runoff

CONTAINS

 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 SUBROUTINE get_hru_runoff(ierr, message)

 ! populate runoff_data with runoff values at LSM domain and at iTime step

  USE public_var,  only:input_dir               ! directory containing input data
  USE public_var,  only:fname_qsim              ! simulated runoff netCDF name
  USE public_var,  only:is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE public_var,  ONLY: qmodOption             !
  USE public_var,  ONLY: integerMissing         !
  USE globalData,  only:iTime
  USE globalData,  only:nHRU
  USE globalData,  only:runoff_data             ! data structure to hru runoff data
  USE globalData,  only:remap_data              ! data structure to remap data
  USE globalData,  ONLY: modTime
  USE globalData,  ONLY: gage_obs_data
  USE globalData,  ONLY: rch_qtake_data
  USE globalData,  ONLY: RCHFLX
  USE read_runoff, only:read_runoff_data        ! read runoff value into runoff_data data strucuture
  USE remapping,   only:remap_runoff            ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE remapping,   only:sort_runoff             ! mapping HM runoff to river network HRU runoff (HM_HRU == RN_HRU)

  implicit none
  ! argument variables
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  character(len=strLen)         :: cmessage           ! error message from subroutine
  integer(i4b)                  :: iens=1
  integer(i4b)                  :: ix, jx
  integer(i4b), allocatable     :: reach_ix(:)
  integer(i4b), parameter       :: no_mod=0
  integer(i4b), parameter       :: direct_insert=1
  integer(i4b), parameter       :: qtake=2
  ! timing
!  integer*8                     :: startTime,endTime,cr ! star and end time stamp, rate
!  real(dp)                      :: elapsedTime          ! elapsed time for the process

  ierr=0; message='get_hru_runoff/'

!  call system_clock(count_rate=cr)

!  call system_clock(startTime)
  ! get the simulated runoff for the current time step - runoff_data%qsim(:) or %qsim2D(:,:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        iTime,                             & ! input: time index
                        runoff_data,                       & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!  call system_clock(endTime)
!  elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
!  write(*,"(A,1PG15.7,A)") '  elapsed-time [runoff_input/read] = ', elapsedTime, ' s'

  ! initialize runoff_data%basinRunoff
  if ( allocated(runoff_data%basinRunoff) ) then
    deallocate(runoff_data%basinRunoff, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinRunoff(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

!  call system_clock(startTime)
  ! Get river network HRU runoff into runoff_data data structure
  if (is_remap) then ! remap LSM simulated runoff to the HRUs in the river network
    call remap_runoff(runoff_data, remap_data, runoff_data%basinRunoff, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  else ! runoff is already remapped to river network HRUs
    call sort_runoff(runoff_data, runoff_data%basinRunoff, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
!  call system_clock(endTime)
!  elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
!  write(*,"(A,1PG15.7,A)") '  elapsed-time [runoff_input/remap] = ', elapsedTime, ' s'

  ! initialize TAKE for water take or direct insersion
  RCHFLX(:,:)%TAKE = 0.0_dp
  select case(qmodOption)
    case(no_mod) ! do nothing
    case(direct_insert)
      ! read gage observation [m3/s] at current time
      jx = gage_obs_data%time_ix(modTime(1))

      if (jx/=integerMissing) then
        call gage_obs_data%read_obs(ierr, cmessage, index_time=jx)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! put qmod at right reach
        reach_ix = gage_obs_data%link_ix()
        do ix=1,size(reach_ix)
          if (reach_ix(ix)==integerMissing) cycle
          RCHFLX(iens,reach_ix(ix))%TAKE = gage_obs_data%get_obs(tix=1, six=ix)
        end do
      end if
    case(qtake)
      ! read reach water take [m3/s] at current time
      jx = rch_qtake_data%time_ix(modTime(1))

      if (jx/=integerMissing) then
        call rch_qtake_data%read_obs(ierr, cmessage, index_time=jx)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! put qmod at right reach
        reach_ix = rch_qtake_data%link_ix()
        do ix=1,size(reach_ix)
          if (reach_ix(ix)==integerMissing) cycle
          RCHFLX(iens,reach_ix(ix))%TAKE = rch_qtake_data%get_obs(tix=1, six=ix)
        end do
      end if
    case default
      ierr=1; message=trim(message)//"Error: qmodOption invalid"; return
  end select

 END SUBROUTINE get_hru_runoff

END MODULE get_runoff
