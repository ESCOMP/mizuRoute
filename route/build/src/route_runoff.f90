!! ======================================================================================================
!! mizuRoute stand-alone driver
!!
!! ======================================================================================================
program route_runoff

! ******
! provide access to desired data types / structures...
! ****************************************************
! variable types
USE nrtype                                     ! variable types, etc.
! shared data
USE public_var, only : iulog                   ! i/o logical unit number
USE globalData, only : pid, nNodes             ! procs id and number of procs
USE globalData, only : mpicom_route            ! communicator

! ******
! provide access to desired subroutines...
! ****************************************
USE mpi_mod,             only : shr_mpi_finalize
! subroutines: model set up
USE model_setup,         only : init_mpi         ! initialize MPI for this program
USE model_setup,         only : init_model       ! model setupt - reading control file, populate metadata, read parameter file
USE model_setup,         only : init_data        ! initialize river reach data
USE model_setup,         only : update_time      ! Update simulation time information at each time step
! subroutines: routing
USE mpi_routine,         only : mpi_route        ! Distribute runoff to proc, route them, and gather,
! subroutines: model I/O
USE get_runoff        ,  only : get_hru_runoff   !
USE write_simoutput_pio, only : prep_output      !
USE write_simoutput_pio, only : output           !
USE write_restart_pio,   only : output_state     ! write netcdf state output file

implicit none

! ******
! define variables
! ************************
character(len=strLen)         :: cfile_name          ! name of the control file
integer(i4b)                  :: ierr                ! error code
character(len=strLen)         :: cmessage            ! error message of downwind routine
integer(i4b)                  :: iens = 1
logical(lgt)                  :: finished=.false.
!Timing
integer*8                     :: cr, startTime, endTime
real(dp)                      :: elapsedTime

! ******
! Initialize the system_clock
! ***********************************
call system_clock(count_rate=cr)

! ******
! get command-line argument defining the full path to the control file
! ***********************************
 call getarg(1,cfile_name)
 if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')

! *****
! *** MPI initialization ....
! ***********************************
call init_mpi()

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
call init_data(pid, nNodes, mpicom_route, ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! ***********************************
! start of time-stepping simulation
! ***********************************
do while (.not.finished)

  ! prepare simulation output netCDF
  call prep_output(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

  ! Get river network hru runoff at current time step
  if(pid==0)then
call system_clock(startTime)
    call get_hru_runoff(ierr, cmessage)
    if(ierr/=0) call handle_err(ierr, cmessage)
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
write(iulog,"(A,1PG15.7,A)") '   elapsed-time [read_ro] = ', elapsedTime, ' s'
  endif

  ! process routing at each proc
call system_clock(startTime)
  call mpi_route(pid, nNodes, mpicom_route, iens, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
write(iulog,"(A,1PG15.7,A)") '   elapsed-time [routing] = ', elapsedTime, ' s'

call system_clock(startTime)
  call output(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
write(iulog,"(A,1PG15.7,A)") '   elapsed-time [output] = ', elapsedTime, ' s'

  call update_time(finished, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

end do  ! looping through time

call output_state(ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

call shr_mpi_finalize(mpicom_route)

stop

contains

 subroutine handle_err(err,message)
 ! handle error codes
 implicit none
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 if(err/=0)then
  write(iulog,*) 'FATAL ERROR: '//trim(message)
  call flush(6)
  stop
 endif
 end subroutine handle_err

end program route_runoff
