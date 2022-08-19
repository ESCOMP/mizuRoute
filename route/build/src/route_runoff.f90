!! ======================================================================================================
!! mizuRoute stand-alone driver
!!
!! ======================================================================================================
PROGRAM route_runoff

! ******
! provide access to external data, subroutines
! ****************************************************
! numeric types
USE nrtype                                       ! variable types, etc.
! subroutines: model set up
USE model_setup,         ONLY: init_model       ! model setupt - reading control file, populate metadata, read parameter file
USE model_setup,         ONLY: init_data        ! initialize river reach data
USE model_setup,         ONLY: update_time      ! Update simulation time information at each time step
! subroutines: routing
USE main_route_module,   ONLY: main_route       ! main routing routine
! subroutines: model I/O
USE get_runoff        ,  ONLY: get_hru_runoff   !
USE write_simoutput,     ONLY: prep_output      !
USE write_simoutput,     ONLY: output           !
USE write_restart,       ONLY: main_restart     ! write netcdf restart file
USE model_finalize,      ONLY: finalize
USE model_finalize,      ONLY: handle_err

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
! system_clock rate
call system_clock(count_rate=cr)

! ******
! get command-line argument defining the full path to the control file
! ***********************************
 call getarg(1,cfile_name)
 if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')

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

  call prep_output(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

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

  call main_restart(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

  call update_time(finished, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

end do

call finalize()

END PROGRAM route_runoff
