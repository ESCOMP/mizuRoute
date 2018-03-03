program route_runoff

! ******
! provide access to desired modules
! ************************

! variable types
USE nrtype                                    ! variable types, etc.
USE dataTypes, only : var_ilength             ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength             ! double precision type: var(:)%dat

! data structures
USE dataTypes,  only : remap                  ! remapping data type
USE dataTypes,  only : runoff                 ! runoff data type

! global data
USE public_var
USE globalData, only : RPARAM                 ! reach parameter structure

!USE globalData, only:time_conv,length_conv    ! conversion factors

! general data structures
!USE globalData, only:remap_data
!USE globalData, only:runoff_data





! metadata
USE var_lookup, only : ixHRU     , nVarsHRU      ! index of variables for data structure
USE var_lookup, only : ixHRU2SEG , nVarsHRU2SEG  ! index of variables for data structure
USE var_lookup, only : ixSEG     , nVarsSEG      ! index of variables for data structure
USE var_lookup, only : ixNTOPO   , nVarsNTOPO    ! index of variables for data structure

! subroutines: populate metadata
USE popMetadat_module, only : popMetadat      ! populate metadata

! subroutines: model control
USE read_control_module, only : read_control  ! read the control file
USE ascii_util_module, only : file_open       ! open file (performs a few checks as well)

! subroutines: netcdf output
USE write_simoutput, only : defineFile        ! define netcdf output file
USE write_netcdf,    only : write_nc          ! write a variable to the NetCDF file

! subroutines: model set up
USE process_ntopo, only : ntopo               ! process the network topology
USE getAncillary_module, only : getAncillary  ! get ancillary data

! subroutines: unit hydrographs
USE basinUH_module, only : basinUH            ! basin unit hydrograph
USE irf_route, only : make_uh                 ! network unit hydrograph

! subroutines: get runoff for each basin in the routing layer
USE read_runoff, only : get_runoff            ! read simulated runoff data
USE remapping,   only : remap_runoff          ! remap runoff

! ******
! define variables
! ************************
implicit none

! index for printing (set to negative to supress printing
integer(i4b),parameter        :: ixPrint = -9999     ! index for printing

! model control
integer(i4b),parameter        :: nEns=1              ! number of ensemble members

! index of looping variables
integer(i4b)                  :: iens                ! ensemble member
integer(i4b)                  :: iHRU                ! HRU index
integer(i4b)                  :: itime               ! time

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

! routing variables
real(dp)                      :: T0,T1               ! entry/exit time for the reach
integer(i4b), allocatable     :: basinID(:)          ! basin ID
real(dp)    , allocatable     :: basinRunoff(:)      ! basin runoff (m/s)
integer(i4b), parameter       :: lakeFlag=0          ! no lakes

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
           structNTOPO,      & ! ancillary data for network toopology
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

! allocate space
allocate(basinID(nHRU), basinRunoff(nHRU), stat=ierr)
if(ierr/=0) call handle_err(ierr, 'unable to allocate space for basinRunoff')

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

! *****
! *** Define model output file...
! *******************************

! define output file
call defineFile(trim(output_dir)//trim(fname_output),  &  ! input: file name
                nHRU,                                  &  ! input: number of HRUs
                nRch,                                  &  ! input: number of stream segments
                time_units,                            &  ! input: time units
                ierr,cmessage)                            ! output: error control
if(ierr/=0) call handle_err(ierr, cmessage)

! define basin ID
forall(iHRU=1:nHRU) basinID(iHRU) = structHRU2seg(iHRU)%var(ixHRU2seg%hruId)%dat(1)
call write_nc(trim(output_dir)//trim(fname_output), 'basinID', basinID, (/1/), (/nHRU/), ierr, cmessage)
call handle_err(ierr,cmessage)

! *****
! *** Route runoff...
! *******************

! define time
T0 = 0._dp
T1 = dt

! define ensemble member
iens=1

! loop through time
do iTime=1,nTime

 ! *****
 ! * Get the simulated runoff for the current time step...
 ! *******************************************************

 ! get the simulated runoff for the current time step
 call get_runoff(trim(input_dir)//trim(fname_qsim), & ! input: filename
                 iTime,                             & ! input: time index
                 nSpatial,                          &  ! input:number of HRUs
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

 ! ensure that simulated runoff is non-zero
 where(basinRunoff < runoffMin) basinRunoff=runoffMin

 ! write time -- note time is just carried across from the input
 call write_nc(trim(output_dir)//trim(fname_output), 'time', (/runoff_data%time/), (/iTime/), (/1/), ierr, cmessage)
 call handle_err(ierr,cmessage)

 ! write the basin runoff to the netcdf file
 call write_nc(trim(output_dir)//trim(fname_output), 'basRunoff', basinRunoff, (/1,iTime/), (/nHRU,1/), ierr, cmessage)
 call handle_err(ierr,cmessage)

 print*, 'PAUSE: after getting simulated runoff'; read(*,*)

 ! *****
 ! * Map the basin runoff to the stream network...
 ! ***********************************************

 ! add



 ! *****
 ! * Perform the routing...
 ! ************************

 ! add





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
