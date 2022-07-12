MODULE globalData
  ! This module includes global data structures

  USE public_var

  ! data types
  USE nrtype

  ! metadata types
  USE dataTypes,  ONLY : struct_info   ! metadata type
  USE dataTypes,  ONLY : dim_info      ! metadata type
  USE dataTypes,  ONLY : var_info      ! metadata type
  USE objTypes,   ONLY : meta_var      ! metadata type

  ! parameter structures
  USE dataTypes,  ONLY : RCHPRP        ! Reach parameters (properties)
  USE dataTypes,  ONLY : RCHTOPO       ! Network topology

  ! routing structures
  USE dataTypes,  ONLY : STRSTA        ! restart state in each reach
  USE dataTypes,  ONLY : STRFLX        ! fluxes in each reach

  ! lake structures
  USE dataTypes,  ONLY : LAKPRP        ! lake properties
  USE dataTypes,  ONLY : LAKTOPO       ! lake topology
  USE dataTypes,  ONLY : LKFLX         ! lake fluxes

  ! remapping structures
  USE dataTypes,  ONLY : remap         ! remapping data type
  USE dataTypes,  ONLY : runoff        ! runoff data type

  ! basin data structure
  USE dataTypes,  ONLY : subbasin_omp  ! mainstem+tributary data structures

  ! time data structure
  USE date_time,  ONLY : datetime      ! time data

  ! time data structure
  USE dataTypes,  ONLY : nc            ! netCDF data

  USE dataTypes,  ONLY : cMolecule     !

  ! data size
  USE var_lookup, ONLY : nStructures   ! number of variables for data structure
  USE var_lookup, ONLY : nDimensions   ! number of variables for data structure
  USE var_lookup, ONLY : nStateDims    ! number of variables for data structure
  USE var_lookup, ONLY : nQdims        ! number of variables for data structure
  USE var_lookup, ONLY : nVarsHRU      ! number of variables for data structure
  USE var_lookup, ONLY : nVarsHRU2SEG  ! number of variables for data structure
  USE var_lookup, ONLY : nVarsSEG      ! number of variables for data structure
  USE var_lookup, ONLY : nVarsNTOPO    ! number of variables for data structure
  USE var_lookup, ONLY : nVarsPFAF     ! number of variables for data structure
  USE var_lookup, ONLY : nVarsRFLX     ! number of variables for data structure
  USE var_lookup, ONLY : nVarsIRFbas   ! number of variables for data structure
  USE var_lookup, ONLY : nVarsBasinQ   ! number of variables for data structure
  USE var_lookup, ONLY : nVarsIRF      ! number of variables for data structure
  USE var_lookup, ONLY : nVarsKWT      ! number of variables for data structure
  USE var_lookup, ONLY : nVarsKW       ! number of variables for data structure
  USE var_lookup, ONLY : nVarsDW       ! number of variables for data structure
  USE var_lookup, ONLY : nVarsMC       ! number of variables for data structure

  implicit none

  save

  ! ---------- MPI/OMP variables -------------------------------------------------------------------
  integer(i4b)                  :: nThreads            ! number of threads

  ! ---------- conversion factors -------------------------------------------------------------------

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
  type(meta_var)                 , public :: meta_rflx   (nVarsRFLX   ) ! reach flux variables
  type(meta_var)                 , public :: meta_irf_bas(nVarsIRFbas ) ! basin IRF routing fluxes/states
  type(meta_var)                 , public :: meta_basinQ (nVarsBasinQ ) ! reach inflow from basin
  type(meta_var)                 , public :: meta_kwt    (nVarsKWT    ) ! KWT routing restart fluxes/states
  type(meta_var)                 , public :: meta_kw     (nVarsKW     ) ! KW routing restart fluxes/states
  type(meta_var)                 , public :: meta_dw     (nVarsDW     ) ! DW routing restart fluxes/states
  type(meta_var)                 , public :: meta_mc     (nVarsMC     ) ! MC routing restart fluxes/states
  type(meta_var)                 , public :: meta_irf    (nVarsIRF    ) ! IRF routing fluxes/states

  ! ---------- data structures ----------------------------------------------------------------------
  integer(i4b)                   , public :: nEns=1               ! number of ensemble
  ! number of spatial elements
  integer(i4b)                   , public :: nHRU                 ! number of HRUs in the whole river network
  integer(i4b)                   , public :: nRch                 ! number of reaches in the whole river network

  ! routing methods
  integer(i4b)                   , public :: nRoutes              ! number of active routing methods
  integer(i4b)    , allocatable  , public :: routeMethods(:)      ! active routing method id
  logical(lgt)                   , public :: onRoute(0:nRouteMethods-1) ! logical to indicate active routing method(s)
  integer(i4b)                   , public :: idxSUM                     ! index of SUM method
  integer(i4b)                   , public :: idxIRF               ! index of IRF method
  integer(i4b)                   , public :: idxKWT               ! index of KWT method
  integer(i4b)                   , public :: idxKW                ! index of KW method
  integer(i4b)                   , public :: idxMC                ! index of MC method
  integer(i4b)                   , public :: idxDW                ! index of DW method

  ! basin and reach IDs (to be removed)
  integer(i8b)    , allocatable  , public :: basinID(:)           ! HRU id
  integer(i4b)    , allocatable  , public :: reachID(:)           ! reach id

  ! DataTime data/variables
  integer(i4b)                   , public :: iTime                ! time index at simulation time step
  real(dp)        , allocatable  , public :: timeVar(:)           ! time variables (unit given by time variable)
  real(dp)                       , public :: TSEC(0:1)            ! begning and end of time step (sec)
  type(datetime)                 , public :: modTime(0:1)         ! previous and current model time (yyyy:mm:dd:hh:mm:ss)
  type(datetime)                 , public :: endCal               ! simulation end date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime)                 , public :: restCal              ! desired restart date/time (yyyy:mm:dd:hh:mm:ss)
  type(datetime)                 , public :: dropCal              ! restart dropoff date/time (yyyy:mm:dd:hh:mm:ss)
  logical(lgt)                   , public :: restartAlarm         ! alarm to triger restart file writing

  ! simulation output netcdf
  type(nc)                       , public :: simout_nc

  ! river topology and parameter structures
  type(RCHPRP)    , allocatable  , public :: RPARAM(:)            ! Reach Parameters for whole domain
  type(RCHTOPO)   , allocatable  , public :: NETOPO(:)            ! River Network topology for whole domain

  ! time delay histogram
  REAL(DP)        , ALLOCATABLE  , public :: FRAC_FUTURE(:)       ! fraction of runoff in future time steps

  ! routing data structures
  TYPE(STRSTA)    , allocatable  , public :: RCHSTA(:,:)          ! Routing state variables (ensembles, space [reaches]) for the entire river network
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
  type(cMolecule)                , public :: nMolecule                ! number of computational molecule (used for KW, MC, DW)

end module globalData
