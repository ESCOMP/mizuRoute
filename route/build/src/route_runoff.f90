program route_runoff

! ******
! provide access to desired data types / structures...
! ****************************************************
! variable types
USE nrtype                                    ! variable types, etc.

! data structures
USE dataTypes,  only : remap                  ! remapping data type
USE dataTypes,  only : runoff                 ! runoff data type

! global data
USE public_var
!USE globalData, only : NETOPO                 ! river network data (tmp)
USE globalData, only : TSEC                    ! time bounds [sec]

! ******
! provide access to desired subroutines...
! ****************************************
! Library
USE mpi                                          ! MPI

USE model_setup,         only : model_setup      ! model setupt - reading control file, populate metadata, read parameter file

! subroutines: model set up
USE data_initialization, only : init_data        ! initialize river segment data

! subroutines: netcdf output
USE write_simoutput,     only : prep_output
USE write_restart,       only : write_state_nc   ! write netcdf state output file
USE write_netcdf,        only : write_nc         ! write a variable to the NetCDF file

! subroutines: get runoff for each basin in the routing layer
USE read_runoff,         only : get_runoff       ! read simulated runoff data
USE remapping,           only : remap_runoff     ! remap runoff from input polygons to routing basins

! subroutines: routing per proc
USE mpi_routine,         only : comm_runoff_data

! ******
! define variables
! ************************
implicit none

! model control
integer(i4b),parameter        :: nEns=1              ! number of ensemble members

! index of looping variables
integer(i4b)                  :: iens                ! ensemble member
integer(i4b)                  :: ix                  ! general loop index
integer(i4b)                  :: iTime,jTime         ! index for time

! error control
integer(i4b)                  :: ierr                ! error code
character(len=strLen)         :: cmessage            ! error message of downwind routine

! define desired reaches
integer(i4b)                  :: nHRU                ! number of HRUs
integer(i4b)                  :: nRch                ! number of desired reaches

! ancillary data on model input
type(remap)                   :: remap_data          ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
type(runoff)                  :: runoff_data         ! runoff for one time step for all HRUs
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
               nEns,         &  ! input:  number of ensembles
               ixSubHRU,     &  ! output: sorted HRU index array
               ixSubSEG,     &  ! output: sorted reach Index array
               hru_per_proc, &  ! output: number of hrus per proc
               seg_per_proc, &  ! output: number of reaches per proc
               remap_data,   &  ! output: remapping data strucutures
               runoff_data,  &  ! output: runoff data strucutures
               ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
!print*, 'PAUSE: after getting network topology'; read(*,*)

!if (pid==7)then
!  print*, 'pid,reachID,downreachID,order'
!  do ix=1,size(NETOPO)
!   write(*,"(I1,A,I9,A,I9,A,I5)") pid,',',NETOPO(ix)%REACHID,',',NETOPO(ix)%DREACHK,',',NETOPO(ix)%RHORDER
!  enddo
!endif

! *****
! *** Allocate space...
! *********************
! if the master node
if (pid==0) then

 ! allocate space for runoff vectors
 nHRU = size(ixSubHRU); nRch=size(ixSubSEG)
 allocate(basinRunoff(nHRU), reachRunoff(nRch), stat=ierr)
 if(ierr/=0) call handle_err(ierr, 'unable to allocate space for runoff vectors')

 ! define ensemble member (temporarily)
 iens=1

end if

! start of time-stepping simulation code
do iTime=1,runoff_data%nTime

  if (pid==0) then

    ! prepare simulation output netCDF
    call prep_output(nEns,           & ! input:  number of ensembles
                     nHRU,           & ! input:  number of HRUs
                     nRch,           & ! input:  number of reaches
                     iTime,          & ! input:  time loop index
                     modJulday,      & ! out:    current simulation julian day
                     jTime,          & ! output: time step in output netCDF
                     ierr, cmessage)   ! output: error control
    call handle_err(ierr,cmessage)

    ! check we are within the desired time interval
    if(modJulday < startJulday .or. modJulday > endJulday) cycle

    ! get the simulated runoff for the current time step
    call get_runoff(trim(input_dir)//trim(fname_qsim), & ! input: filename
                    iTime,                             & ! input: time index
                    runoff_data,                       & ! inout: runoff data structure
                    ierr, cmessage)                      ! output: error control
    call handle_err(ierr, cmessage)

    ! map simulated runoff to the basins in the river network
    if (is_remap) then
     call remap_runoff(runoff_data, remap_data, basinRunoff, ierr, cmessage)
     if(ierr/=0) call handle_err(ierr,cmessage)
    else
     basinRunoff=runoff_data%qsim
    end if

  end if

 ! print*, 'PAUSE: after getting simulated runoff'; read(*,*)

  ! process routing at each proc
  call comm_runoff_data(pid,           &  ! input: proc id
                        nNodes,        &  ! input: number of procs
                        nEns,          &  ! input: number of ensembles
                        nHRU,          &  ! input: number of HRUs in whole domain
                        ixSubHRU,      &  ! input: global HRU index in the order of domains
                        hru_per_proc,  &  ! input: number of hrus assigned to each proc
                        seg_per_proc,  &  ! input: number of hrus assigned to each proc
                        basinRunoff,   &  ! input: runoff data structures
                        ierr, cmessage)   ! output: error control
  call handle_err(ierr,cmessage)

! ! gether fluxes information from each proc and store it in global data structure

! ! *****
! ! * Output
! ! ************************
! ! write time -- note time is just carried across from the input
! call write_nc(trim(fileout), 'time', (/runoff_data%time/), (/jTime/), (/1/), ierr, cmessage)
! call handle_err(ierr,cmessage)
! ! write the basin runoff to the netcdf file
! call write_nc(trim(fileout), 'basRunoff', basinRunoff, (/1,jTime/), (/nHRU,1/), ierr, cmessage)
! call handle_err(ierr,cmessage)
! if (doesBasinRoute == 1) then
!  ! write instataneous local runoff in each stream segment (m3/s)
!  call write_nc(trim(fileout), 'instRunoff', RCHFLX(iens,:)%BASIN_QI, (/1,jTime/), (/nRch,1/), ierr, cmessage)
!  call handle_err(ierr,cmessage)
! endif
! ! write routed local runoff in each stream segment (m3/s)
! call write_nc(trim(fileout), 'dlayRunoff', RCHFLX(iens,:)%BASIN_QR(1), (/1,jTime/), (/nRch,1/), ierr, cmessage)
! call handle_err(ierr,cmessage)
! ! write accumulated runoff (m3/s)
! call write_nc(trim(fileout), 'sumUpstreamRunoff', RCHFLX(iens,:)%UPSTREAM_QI, (/1,jTime/), (/nRch,1/), ierr, cmessage)
! if (ierr/=0) call handle_err(ierr,cmessage)
! if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
!  ! write routed runoff (m3/s)
!  call write_nc(trim(fileout), 'KWTroutedRunoff', RCHFLX(iens,:)%REACH_Q, (/1,jTime/), (/nRch,1/), ierr, cmessage)
!  call handle_err(ierr,cmessage)
! endif
! if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
!  ! write routed runoff (m3/s)
!  call write_nc(trim(fileout), 'IRFroutedRunoff', RCHFLX(iens,:)%REACH_Q_IRF, (/1,jTime/), (/nRch,1/), ierr, cmessage)
!  call handle_err(ierr,cmessage)
! endif
!
 ! increment time bounds
 TSEC(0) = TSEC(0) + dt
 TSEC(1) = TSEC(0) + dt

! !print*, 'PAUSE: after routing'; read(*,*)

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
