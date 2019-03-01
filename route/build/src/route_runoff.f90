program route_runoff

! ******
! provide access to desired data types / structures...
! ****************************************************
! Library
USE mpi                                       ! MPI

! variable types
USE nrtype                                    ! variable types, etc.
USE dataTypes, only : var_ilength             ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength             ! double precision type: var(:)%dat
USE dataTypes, only : var_clength             ! character type: var(:)%dat

! data structures
USE dataTypes,  only : remap                  ! remapping data type
USE dataTypes,  only : runoff                 ! runoff data type
USE dataTypes,  only : time                   ! time data type
USE dataTypes,  only : basin                  ! river basin data type

! global data
USE public_var
USE globalData, only : NETOPO                 ! river network data (tmp)
USE globalData, only : RCHFLX                 ! reach flux structure

! metadata
USE var_lookup, only : ixHRU                  ! index of variables for data structure
USE var_lookup, only : ixHRU2SEG              ! index of variables for data structure
USE var_lookup, only : ixSEG                  ! index of variables for data structure
USE var_lookup, only : ixNTOPO                ! index of variables for data structure
USE var_lookup, only : ixPFAF                 ! index of variables for data structure

! ******
! provide access to desired subroutines...
! ****************************************

! subroutines: populate metadata
USE popMetadat_module, only : popMetadat      ! populate metadata

! subroutines: model control
USE read_control_module, only : read_control  ! read the control file
USE read_param_module,   only : read_param    ! read the routing parameters

! subroutines netcdf input
USE read_netcdf, only:get_nc                  ! netcdf input

! subroutines: netcdf output
USE write_simoutput, only : defineFile        ! define netcdf output file
USE write_restart,   only : define_state_nc,& ! define netcdf state output file
                            write_state_nc    ! write netcdf state output file
USE read_restart,    only : read_state_nc     ! read netcdf state output file
USE write_netcdf,    only : write_nc          ! write a variable to the NetCDF file

! subroutines: model set up
USE data_initialization, only : init_data ! initialize river segment data

! subroutines: routing per proc
USE mpi_routine,         only : comm_runoff_data

! subroutines: model time info
USE time_utils_module,   only : compCalday        ! compute calendar day
USE time_utils_module,   only : compCalday_noleap ! compute calendar day
USE process_time_module, only : process_time  ! process time information

! subroutines: get runoff for each basin in the routing layer
USE read_runoff, only : get_runoff            ! read simulated runoff data
USE remapping,   only : remap_runoff          ! remap runoff from input polygons to routing basins

! ******
! define variables
! ************************
implicit none

! model control
integer(i4b),parameter        :: nEns=1              ! number of ensemble members
character(len=strLen)         :: fileout             ! name of the output file
logical(lgt)                  :: defNewOutputFile    ! flag to define new output file

! index of looping variables
integer(i4b)                  :: iens                ! ensemble member
integer(i4b)                  :: iHRU                ! index for HRU
integer(i4b)                  :: iRch                ! index for the stream segment
integer(i4b)                  :: iTime,jTime         ! index for time
integer(i4b)                  :: iRoute              ! index in routing vector

! error control
integer(i4b)                  :: ierr                ! error code
character(len=strLen)         :: cmessage            ! error message of downwind routine

!  network topology data structures
type(var_ilength),allocatable :: structHRU2seg(:)    ! HRU-to-segment mapping
type(var_ilength),allocatable :: structNTOPO(:)      ! network topology

! read control file
character(len=strLen)         :: cfile_name          ! name of the control file

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

! time structures
type(time)                    :: modTime             ! model time
type(time)                    :: prevTime            ! previous model time
real(dp)                      :: startJulday         ! julian day: start
real(dp)                      :: endJulday           ! julian day: end
real(dp)                      :: refJulday           ! julian day: reference
real(dp)                      :: modJulday           ! julian day: model simulation
real(dp)                      :: convTime2Days       ! conversion factor to convert time to units of days
real(dp)    , allocatable     :: timeVec(:)          ! time vector

! routing variables
real(dp)                      :: T0,T1               ! entry/exit time for the reach
integer(i4b), allocatable     :: basinID(:)          ! basin ID
integer(i4b), allocatable     :: reachID(:)          ! reach ID
real(dp)    , allocatable     :: basinRunoff(:)      ! basin runoff (m/s)
real(dp)    , allocatable     :: reachRunoff(:)      ! reach runoff (m/s)
integer(i4b)                  :: ixDesire            ! desired reach index

! MPI variables
integer(i4b)                  :: pid                 ! process id
integer(i4b)                  :: nNodes              ! number of nodes
! ======================================================================================================
! ======================================================================================================
! ======================================================================================================
! ======================================================================================================
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

! if the master node
if (pid==0) then

 ! populate the metadata files
 call popMetadat(ierr,cmessage)
 if(ierr/=0) call handle_err(ierr, cmessage)

 ! get command-line argument defining the full path to the control file
 call getarg(1,cfile_name)
 if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')

 ! read the control file
 call read_control(trim(cfile_name), ierr, cmessage)
 if(ierr/=0) call handle_err(ierr, cmessage)

 ! read the routing parameter namelist
 call read_param(trim(ancil_dir)//trim(param_nml),ierr,cmessage)
 if(ierr/=0) call handle_err(ierr, cmessage)

endif  ! if the master node

! *****
! *** data initialization
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
!  do iRch=1,size(NETOPO)
!   write(*,"(I1,A,I9,A,I9,A,I5)") pid,',',NETOPO(iRch)%REACHID,',',NETOPO(iRch)%DREACHK,',',NETOPO(iRch)%RHORDER
!  enddo
!endif

! *****
! *** time prep
! ***********************************
! if the master node
if (pid==0) then

! call prep_time(&
!                ierr, cmessage)
! call handle_err(cmessage)

 ! allocate space for the time data
 allocate(timeVec(runoff_data%nTime), stat=ierr)
 if(ierr/=0) call handle_err(ierr, 'unable to allocate space for timeVec')

 ! get the time data
 call get_nc(trim(input_dir)//trim(fname_qsim), vname_time, timeVec, 1, runoff_data%nTime, ierr, cmessage)
 if(ierr/=0) call handle_err(ierr, cmessage)

 ! get the time multiplier needed to convert time to units of days
 select case( trim( time_units(1:index(time_units,' ')) ) )
  case('seconds'); convTime2Days=86400._dp
  case('hours');   convTime2Days=24._dp
  case('days');    convTime2Days=1._dp
  case default;    call handle_err(20, 'unable to identify time units')
 end select

 ! extract time information from the control information
 call process_time(time_units,    calendar, refJulday,   ierr, cmessage); if(ierr/=0) call handle_err(ierr, trim(cmessage)//' [refJulday]')
 call process_time(trim(simStart),calendar, startJulday, ierr, cmessage); if(ierr/=0) call handle_err(ierr, trim(cmessage)//' [startJulday]')
 call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage); if(ierr/=0) call handle_err(ierr, trim(cmessage)//' [endJulday]')

 ! check that the dates are aligned
 if(endJulday<startJulday) call handle_err(20,'simulation end is before simulation start')

 ! initialize previous model time
 prevTime = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

end if

! *****
! *** Initialize state
! *************************************
! if the master node
if (pid==0) then

 ! read restart file and initialize states
 if (isRestart)then
  call read_state_nc(trim(output_dir)//trim(fname_state_in), routOpt, T0, T1, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
 else
  ! Cold start .......
  ! initialize flux structures
  RCHFLX(:,:)%BASIN_QI = 0._dp
  forall(iRoute=0:1) RCHFLX(:,:)%BASIN_QR(iRoute) = 0._dp

  ! initialize time
  T0 = 0._dp
  T1 = dt
 endif

endif

! *****
! * Restart file output
! ************************
! if the master node
if (pid==0) then
 ! Define state netCDF
 call define_state_nc(trim(output_dir)//trim(fname_state_out), time_units, routOpt, ierr, cmessage)
 if(ierr/=0) call handle_err(ierr, cmessage)

! *****
! *** Allocate space...
! *********************

 ! allocate space for runoff vectors
 nHRU = size(ixSubHRU); nRch=size(ixSubSEG)
 allocate(basinID(nHRU), reachID(nRch), basinRunoff(nHRU), reachRunoff(nRch), stat=ierr)
 if(ierr/=0) call handle_err(ierr, 'unable to allocate space for runoff vectors')

 ! define ensemble member (temporarily)
 iens=1

end if
! *****
! *** Route runoff...
! *******************

call MPI_FINALIZE(ierr)
stop

! start of time-stepping simulation code

! loop through time
do iTime=1,runoff_data%nTime

 if (pid==0) then

  ! get the julian day of the model simulation
  modJulday = refJulday + timeVec(iTime)/convTime2Days

  ! check we are within the desired time interval
  if(modJulday < startJulday .or. modJulday > endJulday) cycle

  ! get the time
  select case(trim(calendar))
   case('noleap')
    call compCalday_noleap(modJulday,modTime%iy,modTime%im,modTime%id,modTime%ih,modTime%imin,modTime%dsec,ierr,cmessage)
   case ('standard','gregorian','proleptic_gregorian')
    call compCalday(modJulday,modTime%iy,modTime%im,modTime%id,modTime%ih,modTime%imin,modTime%dsec,ierr,cmessage)
   case default; call handle_err(20, 'calendar name: '//trim(calendar)//' invalid')
  end select
  call handle_err(ierr, cmessage)

  ! print progress
  !print*, modTime%iy,modTime%im,modTime%id,modTime%ih,modTime%imin

  ! *****
  ! *** Define model output file...
  ! *******************************

  ! check need for the new file
  select case(newFileFrequency)
   case(annual); defNewOutputFile=(modTime%iy/=prevTime%iy)
   case(month);  defNewOutputFile=(modTime%im/=prevTime%im)
   case(day);    defNewOutputFile=(modTime%id/=prevTime%id)
   case default; call handle_err(20,'unable to identify the option to define new output files')
  end select

  ! define new file
  if(defNewOutputFile)then

   ! initialize time
   jTime=1

   ! update filename
   write(fileout,'(a,3(i0,a))') trim(output_dir)//trim(fname_output)//'_', modTime%iy, '-', modTime%im, '-', modTime%id, '.nc'

   ! define output file
   call defineFile(trim(fileout),                         &  ! input: file name
                   nEns,                                  &  ! input: number of ensembles
                   nHRU,                                  &  ! input: number of HRUs
                   nRch,                                  &  ! input: number of stream segments
                   time_units,                            &  ! input: time units
                   calendar,                              &  ! input: calendar
                   ierr,cmessage)                            ! output: error control
   if(ierr/=0) call handle_err(ierr, cmessage)

   ! define basin ID
   forall(iHRU=1:nHRU) basinID(iHRU) = structHRU2seg(iHRU)%var(ixHRU2seg%hruId)%dat(1)
   call write_nc(trim(fileout), 'basinID', basinID, (/1/), (/nHRU/), ierr, cmessage)
   call handle_err(ierr,cmessage)

   ! define reach ID
   forall(iRch=1:nRch) reachID(iRch) = structNTOPO(iRch)%var(ixNTOPO%segId)%dat(1)
   call write_nc(trim(fileout), 'reachID', reachID, (/1/), (/nRch/), ierr, cmessage)
   call handle_err(ierr,cmessage)

  ! no new file requested: increment time
  else
   jTime = jTime+1
  endif

 endif

 ! *****
 ! * Get the simulated runoff for the current time step...
 ! *******************************************************
 if (pid==0) then

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

  !print*, 'PAUSE: after getting simulated runoff'; read(*,*)

 endif

 ! process routing at each proc
 call comm_runoff_data(pid,           &  ! input: proc id
                       nNodes,        &  ! input: number of procs
                       nEns,          &  ! input: number of ensembles
                       nHRU,          &  ! input: number of HRUs in whole domain
                       ixSubHRU,      &  ! input: global HRU index in the order of domains
                       T0,T1,         &  ! input: begining and ending of time step (sec)
                       hru_per_proc,  &  ! input: number of hrus assigned to each proc
                       seg_per_proc,  &  ! input: number of hrus assigned to each proc
                       basinRunoff,   &  ! input: runoff data structures
                       ierr, cmessage)   ! output: error control
 call handle_err(ierr,cmessage)

 ! gether fluxes information from each proc and store it in global data structure


 ! *****
 ! * Output
 ! ************************
 ! write time -- note time is just carried across from the input
 call write_nc(trim(fileout), 'time', (/runoff_data%time/), (/jTime/), (/1/), ierr, cmessage)
 call handle_err(ierr,cmessage)
 ! write the basin runoff to the netcdf file
 call write_nc(trim(fileout), 'basRunoff', basinRunoff, (/1,jTime/), (/nHRU,1/), ierr, cmessage)
 call handle_err(ierr,cmessage)
 if (doesBasinRoute == 1) then
  ! write instataneous local runoff in each stream segment (m3/s)
  call write_nc(trim(fileout), 'instRunoff', RCHFLX(iens,:)%BASIN_QI, (/1,jTime/), (/nRch,1/), ierr, cmessage)
  call handle_err(ierr,cmessage)
 endif
 ! write routed local runoff in each stream segment (m3/s)
 call write_nc(trim(fileout), 'dlayRunoff', RCHFLX(iens,:)%BASIN_QR(1), (/1,jTime/), (/nRch,1/), ierr, cmessage)
 call handle_err(ierr,cmessage)
 ! write accumulated runoff (m3/s)
 call write_nc(trim(fileout), 'sumUpstreamRunoff', RCHFLX(iens,:)%UPSTREAM_QI, (/1,jTime/), (/nRch,1/), ierr, cmessage)
 if (ierr/=0) call handle_err(ierr,cmessage)
 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
  ! write routed runoff (m3/s)
  call write_nc(trim(fileout), 'KWTroutedRunoff', RCHFLX(iens,:)%REACH_Q, (/1,jTime/), (/nRch,1/), ierr, cmessage)
  call handle_err(ierr,cmessage)
 endif
 if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
  ! write routed runoff (m3/s)
  call write_nc(trim(fileout), 'IRFroutedRunoff', RCHFLX(iens,:)%REACH_Q_IRF, (/1,jTime/), (/nRch,1/), ierr, cmessage)
  call handle_err(ierr,cmessage)
 endif

 ! increment time bounds
 T0 = T0 + dt
 T1 = T0 + dt

 ! update previous time
 prevTime = modTime

 !print*, 'PAUSE: after routing'; read(*,*)

end do  ! looping through time

! Write states to netCDF
call write_state_nc(trim(output_dir)//trim(fname_state_out), routOpt, runoff_data%time, 1, T0, T1, reachID, ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

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
