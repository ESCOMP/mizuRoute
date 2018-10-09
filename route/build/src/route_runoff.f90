program route_runoff

! ******
! provide access to desired data types / structures...
! ****************************************************

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
USE globalData, only : RCHFLX                 ! reach flux structure
USE globalData, only : KROUTE                 ! routing states

! metadata
USE var_lookup, only : ixHRU     , nVarsHRU      ! index of variables for data structure
USE var_lookup, only : ixHRU2SEG , nVarsHRU2SEG  ! index of variables for data structure
USE var_lookup, only : ixSEG     , nVarsSEG      ! index of variables for data structure
USE var_lookup, only : ixNTOPO   , nVarsNTOPO    ! index of variables for data structure
USE var_lookup, only : ixPFAF    , nVarsPFAF     ! index of variables for data structure

! ******
! provide access to desired subroutines...
! ****************************************

! subroutines: utility
USE nr_utility_module, only: findIndex        ! find index within a vector

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
USE process_ntopo, only : ntopo               ! process the network topology
USE getAncillary_module, only : getAncillary  ! get ancillary data

! Subroutines : pfafstetter code
USE pfafstetter_module,  only: classify_river_basin ! group river network segments into computational groups

! subroutines: model time info
USE time_utils_module,   only : compCalday        ! compute calendar day
USE time_utils_module,   only : compCalday_noleap ! compute calendar day
USE process_time_module, only : process_time  ! process time information

! subroutines: basin routing
USE basinUH_module, only : IRF_route_basin    ! perform UH convolution for basin routing

! subroutines: get runoff for each basin in the routing layer
USE read_runoff, only : get_runoff            ! read simulated runoff data
USE remapping,   only : remap_runoff          ! remap runoff from input polygons to routing basins
USE remapping,   only : basin2reach           ! remap runoff from routing basins to routing reaches

! subroutines: river routing
USE accum_runoff_module, only : accum_runoff  ! upstream flow accumulation
USE kwt_route_module,    only : kwt_route     ! kinematic wave routing method
USE irf_route_module,    only : irf_route     ! river network unit hydrograph method

! subroutines: river network unit hydrograph routing
USE irf_route_module, only : make_uh          ! reach unit hydrograph
USE irf_route_module, only : irf_route        ! river network unit hydrograph method

! ******
! define variables
! ************************
implicit none

! desired routing ids
integer(i4b), parameter       :: desireId=integerMissing  ! turn off checks or speficy reach ID if necessary to print on screen

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
type(var_dlength),allocatable :: structHRU(:)        ! HRU properties
type(var_dlength),allocatable :: structSeg(:)        ! stream segment properties
type(var_ilength),allocatable :: structHRU2seg(:)    ! HRU-to-segment mapping
type(var_ilength),allocatable :: structNTOPO(:)      ! network topology
type(var_clength),allocatable :: structPFAF(:)       ! pfafstetter code

! read control file
character(len=strLen)         :: cfile_name          ! name of the control file

! define desired reaches
integer(i4b)                  :: nHRU                ! number of HRUs
integer(i4b)                  :: nRch                ! number of desired reaches

! ancillary data on model input
type(remap)                   :: remap_data          ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
type(runoff)                  :: runoff_data         ! runoff for one time step for all HRUs
type(basin),     allocatable  :: river_basin(:)     !
integer(i4b)                  :: nSpatial(1:2)       ! number of spatial elements
integer(i4b)                  :: nTime               ! number of time steps
character(len=strLen)         :: time_units          ! time units
character(len=strLen)         :: calendar            ! calendar

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
integer(i4b), parameter       :: doesBasinRoute=1    ! integer to specify runoff input type (0 -> runoff input is already delayed, 1 -> runoff input is instantaneous)
integer(i4b)                  :: ixDesire            ! desired reach index
integer(i4b)                  :: ixOutlet            ! outlet reach index

! ======================================================================================================
! ======================================================================================================
! ======================================================================================================
! ======================================================================================================
! ======================================================================================================
! ======================================================================================================

! start of model/network configuration code

! *****
! *** Populate metadata...
! ************************

! populate the metadata files
call popMetadat(ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! *****
! *** Read control files...
! *************************

! get command-line argument defining the full path to the control file
call getarg(1,cfile_name)
if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')

! read the control file
call read_control(trim(cfile_name), ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! read the routing parameter namelist
call read_param(trim(ancil_dir)//trim(param_nml),ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! *****
! *** Process the network topology...
! ***********************************

! get the network topology
call ntopo(&
           ! output: model control
           nHRU,             & ! number of HRUs
           nRch,             & ! number of stream segments
           ! output: populate data structures
           structHRU,        & ! ancillary data for HRUs
           structSeg,        & ! ancillary data for stream segments
           structHRU2seg,    & ! ancillary data for mapping hru2basin
           structNTOPO,      & ! ancillary data for network topology
           structPFAF,       & ! ancillary data for pfafstetter
           ! output: error control
           ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
!print*, 'PAUSE: after getting network topology'; read(*,*)

! *****
! *** Get ancillary data for routing...
! *************************************

! get ancillary data for routing
call getAncillary(&
                  ! data structures
                  nHRU,            & ! input:  number of HRUs in the routing layer
                  structHRU2seg,   & ! input:  ancillary data for mapping hru2basin
                  is_remap,        & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                  remap_data,      & ! output: data structure to remap data
                  runoff_data,     & ! output: data structure for runoff
                  ! dimensions
                  nSpatial,        & ! output: number of spatial elements in runoff data
                  nTime,           & ! output: number of time steps
                  time_units,      & ! output: time units
                  calendar,        & ! output: calendar
                  ! error control
                  ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! ----------  pfafstetter code process to group segments -------------------------------------------------------
call classify_river_basin(nRch, structPFAF, structNTOPO, river_basin, ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! allocate space for the time data
allocate(timeVec(ntime), stat=ierr)
if(ierr/=0) call handle_err(ierr, 'unable to allocate space for timeVec')

! get the time data
call get_nc(trim(input_dir)//trim(fname_qsim), vname_time, timeVec, 1, nTime, ierr, cmessage)
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

! *****
! *** Initialize state
! *************************************

! allocate space for the routing structures
allocate(RCHFLX(nEns,nRch), KROUTE(nEns,nRch), stat=ierr)
if(ierr/=0) call handle_err(ierr, 'unable to allocate space for routing structures')

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

! *****
! * Restart file output
! ************************

! Define state netCDF
call define_state_nc(trim(output_dir)//trim(fname_state_out), time_units, routOpt, ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! ======================================================================================================
! ======================================================================================================
! ======================================================================================================
! ======================================================================================================
! ======================================================================================================
! ======================================================================================================

! start of time-stepping simulation code

! *****
! *** Allocate space...
! *********************

! allocate space for runoff vectors
allocate(basinID(nHRU), reachID(nRch), basinRunoff(nHRU), reachRunoff(nRch), stat=ierr)
if(ierr/=0) call handle_err(ierr, 'unable to allocate space for runoff vectors')

! *****
! *** define indices...
! *********************

! find index of desired reach
ixDesire = findIndex(reachID,desireId,integerMissing)

! find index of desired reach
ixOutlet = findIndex(reachID,idSegOut,integerMissing)

! define ensemble member
iens=1

! **
! **
! **

! *****
! *** Route runoff...
! *******************

! loop through time
do iTime=1,nTime

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
 ! * Get the simulated runoff for the current time step...
 ! *******************************************************

 ! get the simulated runoff for the current time step
 call get_runoff(trim(input_dir)//trim(fname_qsim), & ! input: filename
                 iTime,                             & ! input: time index
                 nSpatial,                          & ! input: runoff data spatial size
                 runoff_data,                       & ! output: runoff data structure
                 ierr, cmessage)                      ! output: error control
 call handle_err(ierr, cmessage)

 ! map simulated runoff to the basins in the river network
 if (is_remap) then
  call remap_runoff(runoff_data, remap_data, structHRU2seg, nSpatial, basinRunoff, ierr, cmessage)
  if(ierr/=0) call handle_err(ierr,cmessage)
 else
  basinRunoff=runoff_data%qsim
 end if

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
  jtime=1

  ! update filename
  write(fileout,'(a,3(i0,a))') trim(output_dir)//'runoff_', modTime%iy, '-', modTime%im, '-', modTime%id, '.nc'

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
  jtime = jtime+1
 endif

 ! write time -- note time is just carried across from the input
 call write_nc(trim(fileout), 'time', (/runoff_data%time/), (/jTime/), (/1/), ierr, cmessage)
 call handle_err(ierr,cmessage)

 ! write the basin runoff to the netcdf file
 call write_nc(trim(fileout), 'basRunoff', basinRunoff, (/1,jTime/), (/nHRU,1/), ierr, cmessage)
 call handle_err(ierr,cmessage)

 !print*, 'PAUSE: after getting simulated runoff'; read(*,*)

 ! *****
 ! * Map the basin runoff to the stream network...
 ! ***********************************************

 ! map the basin runoff to the stream network...
 call basin2reach(&
                  ! input
                  basinRunoff,       & ! intent(in):  basin runoff (m/s)
                  structNTOPO,       & ! intent(in):  Network topology structure
                  structSEG,         & ! intent(in):  Network attributes structure
                  ! output
                  reachRunoff,       & ! intent(out): reach runoff (m3/s)
                  ierr, cmessage)      ! intent(out): error control
 if(ierr/=0) call handle_err(ierr,cmessage)

 ! perform within-basin routing
 if (doesBasinRoute == 1) then

  ! instantaneous runoff volume (m3/s) to data structure
  RCHFLX(iens,:)%BASIN_QI = reachRunoff(:)

  ! write instataneous local runoff in each stream segment (m3/s)
  call write_nc(trim(fileout), 'instRunoff', RCHFLX(iens,:)%BASIN_QI, (/1,jTime/), (/nRch,1/), ierr, cmessage)
  call handle_err(ierr,cmessage)

  ! perform Basin routing
  call IRF_route_basin(iens, nRch, ierr, cmessage)
  call handle_err(ierr,cmessage)

 ! no within-basin routing (assume handled by the hydrology model)
 else
  RCHFLX(iens,:)%BASIN_QR(0) = RCHFLX(iens,:)%BASIN_QR(1)       ! streamflow from previous step
  RCHFLX(iens,:)%BASIN_QR(1) = reachRunoff(:)                   ! streamflow (m3/s)
 end if

 ! write routed local runoff in each stream segment (m3/s)
 call write_nc(trim(fileout), 'dlayRunoff', RCHFLX(iens,:)%BASIN_QR(1), (/1,jTime/), (/nRch,1/), ierr, cmessage)
 call handle_err(ierr,cmessage)

 !print*, 'PAUSE: after getting reach runoff'; read(*,*)

 ! *****
 ! * Perform the routing...
 ! ************************

 ! perform upstream flow accumulation
 call accum_runoff(iens,          &    ! input: ensemble index
                   nRch,          &    ! input: number of reaches in the river network
                   ixDesire,      &    ! input: index of verbose reach
                   ierr, cmessage)     ! output: error controls
 call handle_err(ierr,cmessage)

 ! write accumulated runoff (m3/s)
 call write_nc(trim(fileout), 'sumUpstreamRunoff', RCHFLX(iens,:)%UPSTREAM_QI, (/1,jTime/), (/nRch,1/), ierr, cmessage)
 if (ierr/=0) call handle_err(ierr,cmessage)

 ! perform KWT routing
 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then

  ! kinematic wave routing
  call kwt_route(iens,                 & ! input: ensemble index
                 nRch,                 & ! input: number of reach in the river network
                 ixDesire,             & ! input: index of the desired reach
                 ixOutlet,             & ! input: index of the outlet reach
                 T0,T1,                & ! input: start and end of the time step
                 ierr,cmessage)          ! output: error control
  call handle_err(ierr,cmessage)

  ! write routed runoff (m3/s)
  call write_nc(trim(fileout), 'KWTroutedRunoff', RCHFLX(iens,:)%REACH_Q, (/1,jTime/), (/nRch,1/), ierr, cmessage)
  call handle_err(ierr,cmessage)

 endif

 ! perform IRF routing
 if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then

  ! IRF routing
  call irf_route(iens,                 & ! input: ensemble index
                 nRch,                 & ! input: number of reach in the river network
                 river_basin,          & ! input: river basin data type
                 ixDesire,             & ! input: index of the desired reach
                 ierr,cmessage)          ! output: error control
  call handle_err(ierr,cmessage)

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
