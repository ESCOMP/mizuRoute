!! ======================================================================================================
!! mizuRoute stand-alone driver
!!
!! ======================================================================================================
program route_runoff

! ******
! provide access to external data, subroutines
! ****************************************************
! variable types
USE nrtype                                       ! variable types, etc.
USE globalData,          only : nThreads         ! a number of threads
! subroutines: model set up
USE model_setup,         only : init_model       ! model setupt - reading control file, populate metadata, read parameter file
USE model_setup,         only : init_data        ! initialize river reach data
USE model_setup,         only : update_time      ! Update simulation time information at each time step
! subroutines: routing
USE main_route_module,   only : main_route       !
! subroutines: model I/O
USE get_runoff        ,  only : get_hru_runoff   !
USE write_simoutput,     only : prep_output      !
USE write_simoutput,     only : output           !
USE write_restart,       only : output_state     ! write netcdf state output file

implicit none

! ******
! define variables
! ************************
character(len=strLen)         :: cfile_name          ! name of the control file
integer(i4b)                  :: ierr                ! error code
character(len=strLen)         :: cmessage            ! error message of downwind routine
integer(i4b)                  :: iens = 1
logical(lgt)                  :: finished=.false.
integer(i4b)                  :: omp_get_num_threads ! number of threads used for openMP
!Timing
integer*8                     :: cr, startTime, endTime
real(dp)                      :: elapsedTime

! ******
! Initialize the system_clock
CALL system_clock(count_rate=cr)

! ******
! get command-line argument defining the full path to the control file
! ***********************************
 call getarg(1,cfile_name)
 if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')

!  Get number of threads
nThreads = 1
!$OMP PARALLEL
!$ nThreads = omp_get_num_threads()
!$OMP END PARALLEL

! *****
! *** model setup
!    - read control files and namelist
!    - broadcast to all processors
! ************************
call init_model(cfile_name, ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! *****
! *** data initialization
!    - river topology, properties, river network domain decomposition
!    - runoff data (datetime, domain)
!    - runoff remapping data
!    - channel states
! ***********************************
call init_data(ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! ***********************************
! start of time-stepping simulation
! ***********************************
do while (.not.finished)

  ! prepare simulation output netCDF
  call prep_output(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

  ! Get river network hru runoff at current time step
call system_clock(startTime)
  call get_hru_runoff(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
write(*,"(A,1PG15.7,A)") '   elapsed-time [read_ro] = ', elapsedTime, ' s'

call system_clock(startTime)
  call main_route(iens, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
write(*,"(A,1PG15.7,A)") '   elapsed-time [routing] = ', elapsedTime, ' s'

call system_clock(startTime)
  call output(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
write(*,"(A,1PG15.7,A)") '   elapsed-time [output] = ', elapsedTime, ' s'

  call update_time(finished, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

end do  ! looping through time

! write state netCDF
!call write_state_nc(trim(output_dir)//trim(fname_state_out), routOpt, runoff_data%time, 1, T0, T1, reachID, ierr, cmessage)
call output_state(ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

stop

contains

 subroutine handle_err(err,message)
 ! handle error codes
 implicit none
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 if(err/=0)then
  print*,'FATAL ERROR: '//trim(message)
  call flush(6)
  stop
 endif
 end subroutine handle_err

end program route_runoff
