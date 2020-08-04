module globalData
  ! This module includes global data structures

  use public_var, only : integerMissing

  ! data types
  use nrtype

  ! metadata types
  use dataTypes,  only : struct_info   ! metadata type
  use dataTypes,  only : dim_info      ! metadata type
  use dataTypes,  only : var_info      ! metadata type
  use objTypes,   only : var_info_new  ! metadata type

  ! parameter structures
  USE dataTypes,  only : RCHPRP        ! Reach parameters (properties)
  USE dataTypes,  only : RCHTOPO       ! Network topology

  ! routing structures
  USE dataTypes,  only : KREACH        ! Collection of flow particles in each reach
  USE dataTypes,  only : STRFLX        ! fluxes in each reach

  ! lake structures
  USE dataTypes,  only : LAKPRP        ! lake properties
  USE dataTypes,  only : LAKTOPO       ! lake topology
  USE dataTypes,  only : LKFLX         ! lake fluxes

  ! remapping structures
  use dataTypes,  only : remap         ! remapping data type
  use dataTypes,  only : runoff        ! runoff data type

  ! basin data structure
  use dataTypes,  only : subbasin_omp  ! mainstem+tributary data structures

  ! time data structure
  use dataTypes,  only : time         ! time data

  ! time data structure
  use dataTypes,  only : nc           ! netCDF data

  ! data size
  USE var_lookup, only : nStructures   ! number of variables for data structure
  USE var_lookup, only : nDimensions   ! number of variables for data structure
  USE var_lookup, only : nStateDims    ! number of variables for data structure
  USE var_lookup, only : nQdims        ! number of variables for data structure
  USE var_lookup, only : nVarsHRU      ! number of variables for data structure
  USE var_lookup, only : nVarsHRU2SEG  ! number of variables for data structure
  USE var_lookup, only : nVarsSEG      ! number of variables for data structure
  USE var_lookup, only : nVarsNTOPO    ! number of variables for data structure
  USE var_lookup, only : nVarsPFAF     ! number of variables for data structure
  USE var_lookup, only : nVarsRFLX     ! number of variables for data structure
  USE var_lookup, only : nVarsIRFbas   ! number of variables for data structure
  USE var_lookup, only : nVarsIRF      ! number of variables for data structure
  USE var_lookup, only : nVarsKWT      ! number of variables for data structure

  implicit none

  save

  ! ---------- MPI/OMP variables -------------------------------------------------------------------
  integer(i4b)                  :: nThreads            ! number of threads

  ! ---------- conversion factors -------------------------------------------------------------------

  real(dp)                       , public :: convTime2Days              ! conversion factor to convert time to units of days
  real(dp)                       , public :: time_conv                  ! time conversion factor -- used to convert to mm/s
  real(dp)                       , public :: length_conv                ! length conversion factor -- used to convert to mm/s

  ! ---------- routing parameter names -------------------------------------------------------------------

  real(dp)                       , public :: fshape                     ! shape parameter in time delay histogram (=gamma distribution) [-]
  real(dp)                       , public :: tscale                     ! scaling factor for the time delay histogram [sec]
  real(dp)                       , public :: velo                       ! velocity [m/s] for Saint-Venant equation   added by NM
  real(dp)                       , public :: diff                       ! diffusivity [m2/s] for Saint-Venant equation   added by NM
  real(dp)                       , public :: mann_n                     ! manning's roughness coefficient [unitless]  added by NM
  real(dp)                       , public :: wscale                     ! scaling factor for river width [-] added by NM

  ! ---------- general structure information --------------------------------------------------------

  type(struct_info)              , public :: meta_struct(nStructures)   ! metadata on the data structures
  type(dim_info)                 , public :: meta_dims(nDimensions)     ! metadata on the dimensions for network topology
  type(dim_info)                 , public :: meta_stateDims(nStateDims) ! metadata on the dimensions for state variables
  type(dim_info)                 , public :: meta_qDims(nQdims)         ! metadata on the dimensions for state variables

  ! ---------- metadata structures ------------------------------------------------------------------

  ! define vectors of metadata
  type(var_info)                 , public :: meta_HRU    (nVarsHRU    ) ! HRU properties
  type(var_info)                 , public :: meta_HRU2SEG(nVarsHRU2SEG) ! HRU-to-segment mapping
  type(var_info)                 , public :: meta_SEG    (nVarsSEG    ) ! stream segment properties
  type(var_info)                 , public :: meta_NTOPO  (nVarsNTOPO  ) ! network topology
  type(var_info)                 , public :: meta_PFAF   (nVarsPFAF   ) ! pfafstetter code
  type(var_info_new)             , public :: meta_rflx   (nVarsRFLX   ) ! reach flux variables
  type(var_info_new)             , public :: meta_irf_bas(nVarsIRFbas ) ! basin IRF routing fluxes/states
  type(var_info_new)             , public :: meta_kwt    (nVarsKWT    ) ! KWT routing fluxes/states
  type(var_info_new)             , public :: meta_irf    (nVarsIRF    ) ! IRF routing fluxes/states

  ! ---------- data structures ----------------------------------------------------------------------
  integer(i4b)                   , public :: nEns=1               ! number of ensemble
  ! number of spatial elements
  integer(i4b)                   , public :: nHRU                 ! number of HRUs in the whole river network
  integer(i4b)                   , public :: nRch                 ! number of reaches in the whole river network

  ! basin and reach IDs (to be removed)
  integer(i8b)    , allocatable  , public :: basinID(:)           ! HRU id
  integer(i4b)    , allocatable  , public :: reachID(:)           ! reach id

  ! DataTime data/variables
  integer(i4b)                   , public :: iTime                ! time index at simulation time step
  real(dp)                       , public :: startJulday          ! julian day: start of routing simulation
  real(dp)                       , public :: endJulday            ! julian day: end of routing simulation
  real(dp)                       , public :: refJulday            ! julian day: reference
  real(dp)                       , public :: modJulday            ! julian day: simulation time step
  real(dp)        , allocatable  , public :: roJulday(:)          ! julian day: runoff input time
  real(dp)        , allocatable  , public :: timeVar(:)           ! time variables (unit given by time variable)
  real(dp)                       , public :: TSEC(0:1)            ! begning and end of time step (sec)
  type(time)                     , public :: modTime(0:1)         ! previous and current model time (yyyy:mm:dd:hh:mm:ss)
  type(time)                     , public :: restCal              ! desired restart date/time (yyyy:mm:dd:hh:mm:ss)
  type(time)                     , public :: dropCal              ! restart dropoff date/time (yyyy:mm:dd:hh:mm:ss)
  logical(lgt)                   , public :: restartAlarm         ! alarm to triger restart file writing

  ! simulation output netcdf
  type(nc)                       , public :: simout_nc

  ! river topology and parameter structures
  type(RCHPRP)    , allocatable  , public :: RPARAM(:)            ! Reach Parameters for whole domain
  type(RCHTOPO)   , allocatable  , public :: NETOPO(:)            ! River Network topology for whole domain

  ! time delay histogram
  REAL(DP)        , ALLOCATABLE  , public :: FRAC_FUTURE(:)       ! fraction of runoff in future time steps

  ! routing data structures
  TYPE(KREACH)    , allocatable  , public :: KROUTE(:,:)          ! Routing state variables (ensembles, space [reaches]) for the entire river network
  TYPE(STRFLX)    , allocatable  , public :: RCHFLX(:,:)          ! Reach fluxes (ensembles, space [reaches]) for entire river network

  ! lakes data structures
  TYPE(LAKPRP)    , allocatable  , public :: LPARAM(:)            ! Lake parameters
  TYPE(LAKTOPO)   , allocatable  , public :: LKTOPO(:)            ! Lake topology
  TYPE(LKFLX)     , allocatable  , public :: LAKFLX(:,:)          ! Lake fluxes

  ! mapping structures
  type(remap)                    , public :: remap_data           ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
  type(runoff)                   , public :: runoff_data          ! runoff data for one time step for LSM HRUs and River network HRUs

  ! domain data
  type(subbasin_omp), allocatable, public :: river_basin(:)       ! openMP domain decomposition

  ! miscellaneous
  integer(i4b)                   , public :: ixPrint=integerMissing   ! index of desired reach to be on-screen print

end module globalData
