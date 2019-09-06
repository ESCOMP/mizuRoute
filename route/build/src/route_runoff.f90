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
USE globalData, only : pid, nNodes, nThreads   ! procs id and number of procs and threads

! ******
! provide access to desired subroutines...
! ****************************************
! Library
USE mpi                                          ! MPI
! subroutines: model set up
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
integer(i4b)                  :: omp_get_num_threads ! number of threads used for openMP
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
call MPI_INIT(ierr); call mpi_handle_err(ierr,0)
!  Get the number of processes
call MPI_COMM_SIZE(MPI_COMM_WORLD, nNodes, ierr)
!  Get the individual process ID
call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
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
call init_data(pid, nNodes, ierr, cmessage)
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
write(*,"(A,1PG15.7,A)") '   elapsed-time [read_ro] = ', elapsedTime, ' s'
  endif

  ! process routing at each proc
call system_clock(startTime)
  call mpi_route(pid, nNodes, iens, ierr, cmessage)
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
 call output_state(ierr, cmessage)
 if(ierr/=0) call handle_err(ierr, cmessage)

!  Shut down MPI
call MPI_FINALIZE(ierr)

stop

contains

 subroutine mpi_handle_err(ierr,pid)
 ! handle MPI error codes
 implicit none
 ! arguments
 integer(i4b),intent(in):: ierr   ! error code
 integer(i4b),intent(in):: pid    ! process ID
 ! local variables
 integer(i4b)           :: jerr   ! error code with message string call
 integer(i4b)           :: errLen ! length of error message
 character(len=strLen)  :: errMsg ! error message

 ! check errors
 if(ierr/=0)then

  ! get error string
  call MPI_Error_String(ierr, errMsg, errLen, jerr)
  if(jerr==0)errMsg='problemIdentifyingErrorMessage'
  if(errLen>strLen)errMsg='errorMessageLengthTooLong'

  ! include process ID
  write(*,'(a,1x,i4)') 'FATAL ERROR (MPI): '//trim(errMsg)//' for process ID ', pid

  ! finalize MPI
  call MPI_FINALIZE(jerr)
  call flush(6)
  stop

 endif

 end subroutine mpi_handle_err

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
