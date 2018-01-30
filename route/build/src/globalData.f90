module globalData
  ! This module includes global data structures
  use nrtype
  use dataTypes, only:remap  ! remapping data type
  use dataTypes, only:runoff ! runoff data type

  implicit none

  ! ---------- conversion factors -------------------------------------------------------------------

  real(dp)      ,save    ,public  :: time_conv              ! time conversion factor -- used to convert to mm/s
  real(dp)      ,save    ,public  :: length_conv            ! length conversion factor -- used to convert to mm/s

  ! ---------- data structures ----------------------------------------------------------------------

  type(remap)   ,save    ,public  :: remap_data             ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
  type(runoff)  ,save    ,public  :: runoff_data            ! runoff for one time step for all HRUs

end module globalData
