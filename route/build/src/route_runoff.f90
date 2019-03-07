!! ======================================================================================================
!! mizuRoute stand-alone driver
!
!
!! ======================================================================================================
program route_runoff

! ******
! provide access to desired data types / structures...
! ****************************************************
! variable types
USE nrtype                                    ! variable types, etc.

! global data
USE public_var

!USE globalData, only : NETOPO                 ! river network data (tmp)
USE globalData, only : TSEC                    ! time bounds [sec]
USE globalData, only : runoff_data             ! hru runoff data structure
USE globalData, only : nHRU, nRch              ! number of HRUs and reaches in the whoel river network

! ******
! provide access to desired subroutines...
! ****************************************
! Library
USE mpi                                          ! MPI

! subroutines: model set up
USE model_setup,         only : model_setup      ! model setupt - reading control file, populate metadata, read parameter file
USE data_initialization, only : init_data        ! initialize river segment data

! subroutines: netcdf output
USE write_simoutput,     only : prep_output
USE write_simoutput,     only : output
USE write_restart,       only : write_state_nc   ! write netcdf state output file

! subroutines: routing per proc
USE mpi_routine,         only : comm_runoff_data

! ******
! define variables
! ************************
implicit none

! error control
integer(i4b)                  :: ierr                ! error code
character(len=strLen)         :: cmessage            ! error message of downwind routine

! index of looping variables
integer(i4b)                  :: ix                  ! general loop index
integer(i4b)                  :: iTime               ! index for time

! ancillary data on model input
integer(i4b),allocatable      :: ixSubHRU(:)         ! global HRU index in the order of domains
integer(i4b),allocatable      :: ixSubSEG(:)         ! global reach index in the order of domains
integer(i4b),allocatable      :: hru_per_proc(:)     ! number of hru assigned to each proc (i.e., node)
integer(i4b),allocatable      :: seg_per_proc(:)     ! number of reaches assigned to each proc (i.e., node)

! time variables
real(dp)                      :: modJulday           ! julian day: model simulation

! routing variables
real(dp)    , allocatable     :: basinRunoff(:)      ! basin runoff (m/s)
real(dp)    , allocatable     :: reachRunoff(:)      ! reach runoff (m/s)

! MPI variables
integer(i4b)                  :: pid                 ! process id
integer(i4b)                  :: nNodes              ! number of nodes
! ======================================================================================================
! ======================================================================================================

! *****
! *** MPI initialization ....
! ***********************************
call MPI_INIT(ierr)
!  Get the number of processes
call MPI_COMM_SIZE(MPI_COMM_WORLD, nNodes, ierr)
!  Get the individual process ID
call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)

! *****
! *** model setup
! ************************
call model_setup(ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! *****
! *** data initialization
!     - river topology, properties (time static)
!     - runoff input
!     - runoff remapping
!     - channel states
! ***********************************
call init_data(nNodes,       &  ! input:  number of procs
               pid,          &  ! input:  proc id
               ixSubHRU,     &  ! output: sorted HRU index array
               ixSubSEG,     &  ! output: sorted reach Index array
               hru_per_proc, &  ! output: number of hrus per proc
               seg_per_proc, &  ! output: number of reaches per proc
               ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
!print*, 'PAUSE: after getting network topology'; read(*,*)

!if (pid==7)then
!  print*, 'pid,reachID,downreachID,order'
!  do ix=1,size(NETOPO)
!   write(*,"(I1,A,I9,A,I9,A,I5)") pid,',',NETOPO(ix)%REACHID,',',NETOPO(ix)%DREACHK,',',NETOPO(ix)%RHORDER
!  enddo
!endif

if (pid==0) then
 ! allocate space for runoff vectors
 allocate(basinRunoff(nHRU), reachRunoff(nRch), stat=ierr)
 if(ierr/=0) call handle_err(ierr, 'unable to allocate space for runoff vectors')
end if

! start of time-stepping simulation code
do iTime=1,runoff_data%nTime

  if (pid==0) then

    ! prepare simulation output netCDF
    call prep_output(iTime,          & ! input:  time loop index
                     modJulday,      & ! out:    current simulation julian day
                     ierr, cmessage)   ! output: error control
    call handle_err(ierr,cmessage)

    ! check we are within the desired time interval
    if(modJulday < startJulday .or. modJulday > endJulday) cycle

    call get_hru_runoff(iTime, ierr, cmessage)
    call handle_err(ierr, cmessage)

  end if
 ! print*, 'PAUSE: after getting simulated runoff'; read(*,*)

  ! process routing at each proc
  call comm_runoff_data(pid,           &  ! input: proc id
                        nNodes,        &  ! input: number of procs
                        ixSubHRU,      &  ! input: global HRU index in the order of domains
                        hru_per_proc,  &  ! input: number of hrus assigned to each proc
                        seg_per_proc,  &  ! input: number of hrus assigned to each proc
                        basinRunoff,   &  ! input: runoff data structures
                        ierr, cmessage)   ! output: error control
  call handle_err(ierr,cmessage)

  ! *****
  ! * Output
  ! ************************
  call output(ierr, cmessage)
  call handle_err(ierr,cmessage)

  ! increment time bounds
  TSEC(0) = TSEC(0) + dt
  TSEC(1) = TSEC(0) + dt

end do  ! looping through time

!! Write states to netCDF
!call write_state_nc(trim(output_dir)//trim(fname_state_out), routOpt, runoff_data%time, 1, TSEC(0), TSEC(1), reachID, ierr, cmessage)
!if(ierr/=0) call handle_err(ierr, cmessage)

!  Shut down MPI
call MPI_FINALIZE(ierr)

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

end
