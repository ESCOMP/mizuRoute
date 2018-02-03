module globalData
  ! This module includes global data structures

  ! data types
  use nrtype
  use dataTypes, only:struct_info  ! metadata type
  use dataTypes, only:var_info     ! metadata type
  use dataTypes, only:remap        ! remapping data type
  use dataTypes, only:runoff       ! runoff data type

  ! data size
  USE var_lookup, only : nStructures   ! number of variables for data structure
  USE var_lookup, only : nVarsHRU      ! number of variables for data structure
  USE var_lookup, only : nVarsHRU2SEG  ! number of variables for data structure
  USE var_lookup, only : nVarsSEG      ! number of variables for data structure
  USE var_lookup, only : nVarsNTOPO    ! number of variables for data structure


  implicit none

  ! ---------- constants ----------------------------------------------------------------------------

  integer(i4b)      , parameter  , public :: scalar=1001  ! named variable for the scalar
  integer(i4b)      , parameter  , public :: vector=1002  ! named variable for the vector

  integer(i4b)      , parameter  , public :: varType_integer=1001  ! named variable for an integer
  integer(i4b)      , parameter  , public :: varType_double=1002   ! named variable for a double precision

  ! ---------- conversion factors -------------------------------------------------------------------

  real(dp)          , save       , public :: time_conv                  ! time conversion factor -- used to convert to mm/s
  real(dp)          , save       , public :: length_conv                ! length conversion factor -- used to convert to mm/s

  ! ---------- metadata structures ------------------------------------------------------------------

 ! define vectors of metadata
 type(var_info)     , save       , public :: meta_HRU    (nVarsHRU    ) ! HRU properties
 type(var_info)     , save       , public :: meta_HRU2SEG(nVarsHRU2SEG) ! HRU-to-segment mapping
 type(var_info)     , save       , public :: meta_SEG    (nVarsSEG    ) ! stream segment properties
 type(var_info)     , save       , public :: meta_NTOPO  (nVarsNTOPO  ) ! network topology

  ! ---------- structure information ------------------------------------------------------------------
  type(struct_info) , save       , public :: meta_struct(nStructures)   ! metadata on the data structures

  ! ---------- data structures ----------------------------------------------------------------------

  type(remap)       , save       , public :: remap_data                 ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
  type(runoff)      , save       , public :: runoff_data                ! runoff for one time step for all HRUs

end module globalData
