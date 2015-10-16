MODULE nhru2basin
 USE nrtype
 IMPLICIT NONE
 SAVE
 ! --------------------------------------------------------------------------------------
 ! Holds the nhru-basin correspondence
 ! --------------------------------------------------------------------------------------
 ! describe individual HRUs within a given basin
 TYPE IPOINT
  INTEGER(I4B)                           :: hru_ix   ! index of HRU in the nHRU vector
  INTEGER(I4B)                           :: hru_id   ! HRU id
  REAL(DP)                               :: hru_lon  ! longitude
  REAL(DP)                               :: hru_lat  ! latitude
  REAL(DP)                               :: hru_elev ! elevation
  REAL(DP)                               :: hru_area ! area
  REAL(DP)                               :: wght     ! weight assigned to the grid point
 END TYPE IPOINT
 ! Collection of HRUs assigned to a given basin
 TYPE ALL_POINTS
  TYPE(IPOINT), DIMENSION(:), POINTER    :: cHRU
 END TYPE ALL_POINTS
 ! the nhru2basin correspondence
 TYPE(ALL_POINTS), DIMENSION(:), POINTER :: h2b
 ! --------------------------------------------------------------------------------------
END MODULE nhru2basin
