program route_runoff

! ******
! provide access to desired data types / structures...
! ****************************************************

! variable types
USE nrtype                                    ! variable types, etc.
USE dataTypes, only : var_ilength             ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength             ! double precision type: var(:)%dat

! data structures
USE dataTypes,  only : remap                  ! remapping data type
USE dataTypes,  only : runoff                 ! runoff data type
USE dataTypes,  only : time                   ! time data type

! global data
USE public_var
USE globalData, only : NETOPO                 ! network topology structure
USE globalData, only : RPARAM                 ! reach parameter structure
USE globalData, only : RCHFLX                 ! reach flux structure
USE globalData, only : KROUTE                 ! routing states

! metadata
USE var_lookup, only : ixHRU     , nVarsHRU      ! index of variables for data structure
USE var_lookup, only : ixHRU2SEG , nVarsHRU2SEG  ! index of variables for data structure
USE var_lookup, only : ixSEG     , nVarsSEG      ! index of variables for data structure
USE var_lookup, only : ixNTOPO   , nVarsNTOPO    ! index of variables for data structure

! ******
! provide access to desired subroutines...
! ****************************************

! subroutines: utility
USE nr_utility_module, only: findIndex        ! find index within a vector

! subroutines: populate metadata
USE popMetadat_module, only : popMetadat      ! populate metadata

! subroutines: model control
USE read_control_module, only : read_control  ! read the control file
USE ascii_util_module, only : file_open       ! open file (performs a few checks as well)

! subroutines netcdf input
USE read_netcdf, only:get_nc                  ! netcdf input

! subroutines: netcdf output
USE write_simoutput, only : defineFile        ! define netcdf output file
USE write_netcdf,    only : write_nc          ! write a variable to the NetCDF file

! subroutines: model set up
USE process_ntopo, only : ntopo               ! process the network topology
USE getAncillary_module, only : getAncillary  ! get ancillary data

! subroutines: model time info
USE time_utils_module, only : extractTime     ! get time from units string
USE time_utils_module, only : compJulday      ! compute julian day
USE time_utils_module, only : compCalday      ! compute calendar day

! subroutines: unit hydrographs
USE basinUH_module, only : basinUH            ! basin unit hydrograph
USE irf_route, only : make_uh                 ! network unit hydrograph

! subroutines: get runoff for each basin in the routing layer
USE read_runoff, only : get_runoff            ! read simulated runoff data
USE remapping,   only : remap_runoff          ! remap runoff from input polygons to routing basins
USE remapping,   only : basin2reach           ! remap runoff from routing basins to routing reaches

! subroutines: routing
USE kwt_route,   only : QROUTE_RCH            ! kinematic wave routing method

! ******
! define variables
! ************************
implicit none

! index for printing (set to negative to supress printing
integer(i4b),parameter        :: ixPrint = -9999     ! index for printing

! model control
integer(i4b),parameter        :: nEns=1              ! number of ensemble members
character(len=strLen)         :: fileout             ! name of the output file
logical(lgt)                  :: defNewOutputFile    ! flag to define new output file

! index of looping variables
integer(i4b)                  :: iens                ! ensemble member
integer(i4b)                  :: iHRU                ! index for HRU
integer(i4b)                  :: iRch,jRch           ! index for the stream segment
integer(i4b)                  :: itime               ! index for time
integer(i4b)                  :: jtime               ! index for time
integer(i4b)                  :: iRoute              ! index in routing vector

! error control
integer(i4b)                  :: ierr                ! error code
character(len=strLen)         :: cmessage            ! error message of downwind routine

!  network topology data structures
type(var_dlength),allocatable :: structHRU(:)        ! HRU properties
type(var_dlength),allocatable :: structSeg(:)        ! stream segment properties
type(var_ilength),allocatable :: structHRU2seg(:)    ! HRU-to-segment mapping
type(var_ilength),allocatable :: structNTOPO(:)      ! network topology

! read control file
character(len=strLen)         :: cfile_name          ! name of the control file
integer(i4b)                  :: iunit               ! file unit

! define desired reaches
integer(i4b)                  :: nHRU                ! number of HRUs
integer(i4b)                  :: nRch                ! number of desired reaches

! ancillary data on model input
type(remap)                   :: remap_data          ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
type(runoff)                  :: runoff_data         ! runoff for one time step for all HRUs
integer(i4b)                  :: nSpatial            ! number of spatial elements
integer(i4b)                  :: nTime               ! number of time steps
character(len=strLen)         :: time_units          ! time units

! time structures
type(time)                    :: startTime           ! start time
type(time)                    :: endTime             ! end time
type(time)                    :: refTime             ! reference time
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
integer(i4b), parameter       :: lakeFlag=0          ! no lakes
integer(i4b)                  :: ixDesire            ! desired reach index
integer(i4b)                  :: ixOutlet            ! outlet reach index

! desired routing ids
integer(i4b), parameter       :: desireId=integerMissing  ! turn off checks

! namelist parameters
real(dp)                      :: fshape              ! shape parameter in time delay histogram (=gamma distribution) [-]
real(dp)                      :: tscale              ! scaling factor for the time delay histogram [sec]
real(dp)                      :: velo                ! velocity [m/s] for Saint-Venant equation   added by NM
real(dp)                      :: diff                ! diffusivity [m2/s] for Saint-Venant equation   added by NM
real(dp)                      :: mann_n              ! manning's roughness coefficient [unitless]  added by NM
real(dp)                      :: wscale              ! scaling factor for river width [-] added by NM
namelist /HSLOPE/fshape,tscale  ! route simulated runoff through the local basin
namelist /IRF_UH/velo,diff      ! route delayed runoff through river network with St.Venant UH
namelist /KWT/mann_n,wscale     ! route kinematic waves through the river network

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

! read the name list
call file_open(trim(param_nml),iunit,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
read(iunit, nml=HSLOPE)
read(iunit, nml=IRF_UH)
read(iunit, nml=KWT)
close(iunit)

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
           ! output: error control
           ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
!print*, 'PAUSE: after getting network topology'; read(*,*)

! specify some additional routing parameters (temporary "fix")
! NOTE: include here because using namelist parameters
if(hydGeometryOption==compute)then
 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
  RPARAM(:)%R_WIDTH = wscale * sqrt(RPARAM(:)%TOTAREA)  ! channel width (m)
  RPARAM(:)%R_MAN_N = mann_n                            ! Manning's "n" paramater (unitless)
 end if
endif  ! computing network topology

! *****
! *** Get ancillary data for routing...
! *************************************

! compute the time-delay histogram (to route runoff within basins)
! NOTE: allocates and populates global data FRAC_FUTURE
call basinUH(dt, fshape, tscale, ierr, cmessage)
call handle_err(ierr, cmessage)

! For IRF routing scheme: Compute unit hydrograph for each segment
! NOTE: include here because using namelist parameters
if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
 call make_uh(nRch, dt, velo, diff, ierr, cmessage)
 call handle_err(ierr, cmessage)
end if

! get ancillary data for routing
call getAncillary(&
                  ! data structures
                  nHRU,            & ! input:  number of HRUs in the routing layer
                  structHRU2seg,   & ! input:  ancillary data for mapping hru2basin
                  remap_data,      & ! output: data structure to remap data
                  runoff_data,     & ! output: data structure for runoff
                  ! dimensions
                  nSpatial,        & ! output: number of spatial elements in runoff data
                  nTime,           & ! output: number of time steps
                  time_units,      & ! output: time units
                  ! error control
                  ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! allocate space for the time data
allocate(timeVec(ntime), stat=ierr); call handle_err(ierr, 'problem allocating timeVec')

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

! * extract reference time from the time string

! extract reference time from the units string
call extractTime(time_units,refTime%iy,refTime%im,refTime%id,refTime%ih,refTime%imin,refTime%dsec,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! extract start time from the simStart string
call extractTime(trim(simStart),startTime%iy,startTime%im,startTime%id,startTime%ih,startTime%imin,startTime%dsec,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! extract end time from the simStart string
call extractTime(trim(simEnd),endTime%iy,endTime%im,endTime%id,endTime%ih,endTime%imin,endTime%dsec,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! * compute the julian day from the time structures

! calculate the reference julian day
call compjulday(refTime%iy,refTime%im,refTime%id,refTime%ih,refTime%imin,refTime%dsec,refJulday,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! calculate the julian day for the start time
call compjulday(startTime%iy,startTime%im,startTime%id,startTime%ih,startTime%imin,startTime%dsec,startJulday,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! calculate the julian day for the end time
call compjulday(endTime%iy,endTime%im,endTime%id,endTime%ih,endTime%imin,endTime%dsec,endJulday,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! check that the dates are aligned
if(endJulday<startJulday) call handle_err(20,'simulation end is before simulation start')

! initialize previous model time
prevTime = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

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

! allocate space for the routing structures
allocate(RCHFLX(nens,nRch), KROUTE(nens,nRch), stat=ierr)
if(ierr/=0) call handle_err(ierr, 'unable to allocate space for routing structures')

! *****
! *** iniitalize vectors...
! *************************

! initialize flux structures
RCHFLX(:,:)%BASIN_QI = 0._dp
forall(iRoute=0:1) RCHFLX(:,:)%BASIN_QR(iRoute) = 0._dp

! find index of desired reach
ixDesire = findIndex(reachID,desireId,integerMissing)

! find index of desired reach
ixOutlet = findIndex(reachID,idSegOut,integerMissing)

! define ensemble member
iens=1

! initialize time
T0 = 0._dp
T1 = dt

! **
! **
! **

! loop through time
do iTime=1,nTime

 ! get the julian day of the model simulation
 modJulday = refJulday + timeVec(iTime)/convTime2Days

 ! check we are within the desired time interval
 if(modJulday < startJulday .or. modJulday > endJulday) cycle

 ! get the time
 call compCalday(modJulday,modTime%iy,modTime%im,modTime%id,modTime%ih,modTime%imin,modTime%dsec,ierr,cmessage)
 call handle_err(ierr, cmessage)

 ! print progress
 print*, modTime%iy,modTime%im,modTime%id,modTime%ih,modTime%imin

 ! *****
 ! * Get the simulated runoff for the current time step...
 ! *******************************************************

 ! get the simulated runoff for the current time step
 call get_runoff(trim(input_dir)//trim(fname_qsim), & ! input: filename
                 iTime,                             & ! input: time index
                 nSpatial,                          & ! input:number of HRUs
                 runoff_data%time,                  & ! output: time
                 runoff_data%qSim,                  & ! output: runoff data
                 ierr, cmessage)                      ! output: error control
 call handle_err(ierr, cmessage)

 ! map simulated runoff to the basins in the river network
 if (is_remap) then
  call remap_runoff(runoff_data, remap_data, structHRU2seg, basinRunoff, ierr, cmessage)
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
                  nHRU,                                  &  ! input: number of HRUs
                  nRch,                                  &  ! input: number of stream segments
                  time_units,                            &  ! input: time units
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
                  reachRunoff,       & ! intent(out): reach runoff (m/s)
                  ierr, cmessage)      ! intent(out): error control
 if(ierr/=0) call handle_err(ierr,cmessage)

 ! convert runoff to m3/s
 ! NOTE: Use BASIN_QR here because input runoff is already routed
 RCHFLX(iens,:)%BASIN_QR(0) = RCHFLX(iens,:)%BASIN_QR(1)       ! streamflow from previous step
 RCHFLX(iens,:)%BASIN_QR(1) = reachRunoff(:)*RPARAM(:)%BASAREA ! streamflow (m3/s)

 ! ensure that routed streamflow is non-zero
 do iRch=1,nRch
  if(RCHFLX(iens,iRch)%BASIN_QR(1) < runoffMin) RCHFLX(iens,iRch)%BASIN_QR(1)=runoffMin
 end do

 ! write routed local runoff in each stream segment (m3/s)
 call write_nc(trim(fileout), 'dlayRunoff', RCHFLX(iens,:)%BASIN_QR(1), (/1,jTime/), (/nRch,1/), ierr, cmessage)
 call handle_err(ierr,cmessage)

 !print*, 'PAUSE: after getting reach runoff'; read(*,*)

 ! *****
 ! * Perform the routing...
 ! ************************

 ! route streamflow through the river network
 do iRch=1,nRch

  ! identify reach to process
  jRch = NETOPO(iRch)%RHORDER

  ! check
  if(reachId(jRch) == desireId)then
   print*, 'reachRunoff(jRch), RPARAM(jRch)%BASAREA, RCHFLX(iens,jRch)%BASIN_QR(1) = ', &
            reachRunoff(jRch), RPARAM(jRch)%BASAREA, RCHFLX(iens,jRch)%BASIN_QR(1)
  endif

  ! route kinematic waves through the river network
  call QROUTE_RCH(iens,jrch,           & ! input: array indices
                  ixDesire,            & ! input: index of the desired reach
                  ixOutlet,            & ! input: index of the outlet reach
                  T0,T1,               & ! input: start and end of the time step
                  MAXQPAR,             & ! input: maximum number of particle in a reach
                  LAKEFLAG,            & ! input: flag if lakes are to be processed
                  ierr,cmessage)         ! output: error control
  if (ierr/=0) call handle_err(ierr,cmessage)

 end do  ! (looping through stream segments)

 ! write routed runoff (m3/s)
 call write_nc(trim(fileout), 'KWTroutedRunoff', RCHFLX(iens,:)%REACH_Q, (/1,jTime/), (/nRch,1/), ierr, cmessage)
 call handle_err(ierr,cmessage)

 ! increment time bounds
 T0 = T0 + dt
 T1 = T0 + dt

 ! update previous time
 prevTime = modTime

 !print*, 'PAUSE: after routing'; read(*,*)

end do  ! looping through time

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
