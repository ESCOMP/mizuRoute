MODULE lakes_flux
 USE nrtype
 IMPLICIT NONE
 SAVE
 TYPE LKFLX
  REAL(DP)                             :: LAKE_Qav ! lake discharge (average over time step) (m3 s-1)
  REAL(DP)                             :: LAKE_Q ! lake discharge (instantaneous) (m3 s-1)
  REAL(DP)                             :: LAKE_P ! lake precipitation (m3)
  REAL(DP)                             :: LAKE_E ! lake evaporation (m3)
  REAL(DP)                             :: LAKE_I ! inflow to lake (m3 s-1)
 ENDTYPE LKFLX
 TYPE(LKFLX), DIMENSION(:,:), POINTER :: LAKFLX  ! Reach fluxes
END MODULE lakes_flux

