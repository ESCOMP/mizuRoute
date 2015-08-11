MODULE reachstate
 USE nrtype
 IMPLICIT NONE
 SAVE
 ! --------------------------------------------------------------------------------------
 ! Collection of flow particles in each reach (Kinematic wave routing) 
 ! --------------------------------------------------------------------------------------
 ! Individual flow particles
 TYPE FPOINT
  REAL(DP)                             :: QF       ! Flow
  REAL(DP)                             :: QM       ! Modified flow
  REAL(DP)                             :: TI       ! initial time of point in reach
  REAL(DP)                             :: TR       ! time point expected to exit reach
  LOGICAL(LGT)                         :: RF       ! routing flag (T if point has exited)
 END TYPE FPOINT
 ! Collection of flow points within a given reach
 TYPE KREACH
  TYPE(FPOINT), DIMENSION(:), POINTER  :: KWAVE => null()
 END TYPE KREACH
 ! type kreach for all stream segments
 TYPE(KREACH), DIMENSION(:,:), POINTER :: KROUTE => null()
 ! --------------------------------------------------------------------------------------
END MODULE reachstate
