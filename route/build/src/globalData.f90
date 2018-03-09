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

  ! data size
  USE var_lookup, only : nStructures   ! number of variables for data structure
  USE var_lookup, only : nDimensions   ! number of variables for data structure
  USE var_lookup, only : nVarsHRU      ! number of variables for data structure
  USE var_lookup, only : nVarsHRU2SEG  ! number of variables for data structure
  USE var_lookup, only : nVarsSEG      ! number of variables for data structure
  USE var_lookup, only : nVarsNTOPO    ! number of variables for data structure

  implicit none

  save

  ! ---------- physical constants -------------------------------------------------------------------

  real(dp)          , parameter  , public :: pi=3.14159265359_dp        ! pi

  ! ---------- constants ----------------------------------------------------------------------------

  ! true/false
  integer(i4b)      , parameter  , public :: true=1001                  ! true
  integer(i4b)      , parameter  , public :: false=1002                 ! false

  ! variable types
  integer(i4b)      , parameter  , public :: varType_integer = 1001     ! named variable for an integer
  integer(i4b)      , parameter  , public :: varType_double  = 1002     ! named variable for a double precision

  ! ---------- conversion factors -------------------------------------------------------------------

  real(dp)                       , public :: time_conv                  ! time conversion factor -- used to convert to mm/s
  real(dp)                       , public :: length_conv                ! length conversion factor -- used to convert to mm/s

  ! ---------- general structure information --------------------------------------------------------

  type(struct_info)              , public :: meta_struct(nStructures)   ! metadata on the data structures
  type(dim_info)                 , public :: meta_dims(nDimensions)     ! metadata on the dimensions

  ! ---------- metadata structures ------------------------------------------------------------------

  ! define vectors of metadata
  type(var_info)                 , public :: meta_HRU    (nVarsHRU    ) ! HRU properties
  type(var_info)                 , public :: meta_HRU2SEG(nVarsHRU2SEG) ! HRU-to-segment mapping
  type(var_info)                 , public :: meta_SEG    (nVarsSEG    ) ! stream segment properties
  type(var_info)                 , public :: meta_NTOPO  (nVarsNTOPO  ) ! network topology

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

end module globalData
