module globalData
  ! This module includes global data structures

  ! data types
  use nrtype

  ! metadata types
  use dataTypes,  only : struct_info   ! metadata type
  use dataTypes,  only : dim_info      ! metadata type
  use dataTypes,  only : var_info      ! metadata type

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
  use dataTypes,  only : subbasin_mpi  ! reach category (store mainstem code or pfaf code)

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
  USE var_lookup, only : nVarsIRFbas   ! number of variables for data structure
  USE var_lookup, only : nVarsIRF      ! number of variables for data structure
  USE var_lookup, only : nVarsKWT      ! number of variables for data structure

  implicit none

  save

  ! ---------- physical constants -------------------------------------------------------------------

  real(dp)          , parameter  , public :: pi=3.14159265359_dp        ! pi

  ! ---------- constants ----------------------------------------------------------------------------

  ! true/false
  integer(i4b)      , parameter  , public :: true=1001                  ! true
  integer(i4b)      , parameter  , public :: false=1002                 ! false

  ! variable types
  integer(i4b)      , parameter  , public :: varType_integer   = 1001   ! named variable for an integer
  integer(i4b)      , parameter  , public :: varType_double    = 1002   ! named variable for a double precision
  integer(i4b)      , parameter  , public :: varType_character = 1003   ! named variable for a double precision

  ! ---------- conversion factors -------------------------------------------------------------------

  real(dp)                       , public :: time_conv                  ! time conversion factor -- used to convert to mm/s
  real(dp)                       , public :: length_conv                ! length conversion factor -- used to convert to mm/s

  ! ---------- routing parameter names -------------------------------------------------------------------

  real(dp)                      :: fshape                               ! shape parameter in time delay histogram (=gamma distribution) [-]
  real(dp)                      :: tscale                               ! scaling factor for the time delay histogram [sec]
  real(dp)                      :: velo                                 ! velocity [m/s] for Saint-Venant equation   added by NM
  real(dp)                      :: diff                                 ! diffusivity [m2/s] for Saint-Venant equation   added by NM
  real(dp)                      :: mann_n                               ! manning's roughness coefficient [unitless]  added by NM
  real(dp)                      :: wscale                               ! scaling factor for river width [-] added by NM

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
  type(var_info)                 , public :: meta_irf_bas(nVarsIRFbas ) ! basin IRF routing fluxes/states
  type(var_info)                 , public :: meta_kwt    (nVarsKWT    ) ! KWT routing fluxes/states
  type(var_info)                 , public :: meta_irf    (nVarsIRF    ) ! IRF routing fluxes/states
  ! ---------- data structures ----------------------------------------------------------------------

  ! routing parameter structures
  type(RCHPRP)    , allocatable  , public :: RPARAM(:)       ! Reach Parameters
  type(RCHTOPO)   , allocatable  , public :: NETOPO(:)       ! River Network topology

  ! time delay histogram
  REAL(DP)        , ALLOCATABLE  , public :: FRAC_FUTURE(:)  ! fraction of runoff in future time steps

  ! routing data structures
  TYPE(KREACH)    , allocatable  , public :: KROUTE(:,:)     ! Routing state variables (ensembles, space [reaches])
  TYPE(STRFLX)    , allocatable  , public :: RCHFLX(:,:)     ! Reach fluxes (ensembles, space [reaches])

  ! lakes data structures
  TYPE(LAKPRP)    , allocatable  , public :: LPARAM(:)       ! Lake parameters
  TYPE(LAKTOPO)   , allocatable  , public :: LKTOPO(:)       ! Lake topology
  TYPE(LKFLX)     , allocatable  , public :: LAKFLX(:,:)     ! Lake fluxes

  ! mapping structures
  type(remap)                    , public :: remap_data      ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
  type(runoff)                   , public :: runoff_data     ! runoff for one time step for all HRUs

  ! domain data
  type(subbasin_mpi)             , public :: domains(5000)   ! domain decomposition based on subbasin
  integer(i4b)                   , public :: nDomain         ! domain counter

end module globalData
