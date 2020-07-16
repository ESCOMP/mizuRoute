!! ======================================================================================================
!! mizuRoute stand-alone driver
!!
!! ======================================================================================================
PROGRAM route_runoff

! ******
! provide access to desired data types / structures...
! ****************************************************
! variable types
USE nrtype                                     ! variable types, etc.
! shared data
USE globalData, ONLY: pid, nNodes             ! procs id and number of procs
USE globalData, ONLY: mpicom_route            ! communicator

! ******
! provide access to desired subroutines...
! ****************************************
! subroutines: model set up
USE perf_mod,            ONLY: t_initf          ! initialize timing routines (GPTL library)
USE perf_mod,            ONLY: t_startf,t_stopf ! timing start/stop (GPTL library)
USE model_setup,         ONLY: init_mpi         ! initialize MPI for this program
USE model_setup,         ONLY: init_data        ! initialize river reach data
USE model_setup,         ONLY: infile_name      ! updating the file name and iTime_local based on iTime
USE init_model_data,     ONLY: init_model       ! model setupt - reading control file, populate metadata, read parameter file
USE init_model_data,     ONLY: update_time      ! Update simulation time information at each time step
! subroutines: model finalize
USE model_utils,         ONLY: model_finalize
USE model_utils,         ONLY: handle_err
! subroutines: routing
USE mpi_routine,         ONLY: mpi_route        ! Distribute runoff to proc, route them, and gather,
! subroutines: model I/O
USE get_runoff,          ONLY: get_hru_runoff   !
USE write_simoutput_pio, ONLY: prep_output      !
USE write_simoutput_pio, ONLY: output           !
USE write_restart_pio,   ONLY: output_state     ! write netcdf state output file

implicit none

! ******
! define variables
! ************************
character(len=strLen)         :: cfile_name          ! name of the control file
integer(i4b)                  :: ierr                ! error code
character(len=strLen)         :: cmessage            ! error message of downwind routine
integer(i4b)                  :: iens = 1
logical(lgt)                  :: finished=.false.

! ******
! get command-line argument defining the full path to the control file
! ***********************************
 call getarg(1,cfile_name)
 if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')

! *****
! *** MPI/timer initialization ....
! ***********************************
call init_mpi()
call t_initf('dummy_empty_namelist', LogPrint=.true., mpicom=mpicom_route, mastertask=.true.)

! *****
! *** model setup
!    - read control files and namelist
! ************************
call init_model(cfile_name, ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! *****
! *** data initialization
!    - river topology, properties, river network domain decomposition
!    - runoff data (datetime, domain)
!    - runoff remapping data
!    - channel states
!    - broadcast to all processors
! ***********************************
call init_data(pid, nNodes, mpicom_route, ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! ***********************************
! start of time-stepping simulation
! ***********************************
do while (.not.finished)

  ! update the name of the file and iTime_local
  call infile_name(ierr, cmessage)         ! output

  call prep_output(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

  if(pid==0)then
    call t_startf ('input')
    call get_hru_runoff(ierr, cmessage)
    if(ierr/=0) call handle_err(ierr, cmessage)
    call t_stopf ('input')
  endif

  call t_startf ('route-total')
  call mpi_route(pid, nNodes, mpicom_route, iens, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
  call t_stopf ('route-total')

  call t_startf ('output')
  call output(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
  call t_stopf ('output')

  call output_state(ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

  call update_time(finished, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)

end do

call model_finalize(mpicom_route)

END PROGRAM route_runoff
