MODULE lake_param
 USE nrtype
 IMPLICIT NONE
 SAVE
 ! Lake Parameters
 TYPE LAKPRP
  REAL(DP)                             :: AREAREF  ! lake area
  REAL(DP)                             :: LAKREFLEV  ! lake elevation
  REAL(DP)                             :: LAKAVGLEV  ! lake average level (for initialization)
  REAL(DP)                             :: HE2AR_C  ! water height-surface area parameter
  REAL(DP)                             :: HE2AR_D  ! water height-surface area parameter
  REAL(DP)                             :: HGHTLOW  ! minimum water height for discharge
  REAL(DP)                             :: HGHTECO  ! minimum height for ecological concerns
  REAL(DP)                             :: HGHTSPL  ! spillway height
  REAL(DP)                             :: DSCHECO  ! discharge at "ecological" height
  REAL(DP)                             :: DSCHSPL  ! discharge at spillway height
  REAL(DP)                             :: RATECVA  ! discharge rating curve parameter
  REAL(DP)                             :: RATECVB  ! discharge rating curve parameter 
 ENDTYPE LAKPRP
 TYPE(LAKPRP), DIMENSION(:), POINTER   :: LPARAM  ! Reach Parameters
 ! Lake topology
 TYPE LAKTOPO
  INTEGER(I4B)                         :: LAKE_IX   ! Lake index (0,1,2,...,nlak-1)
  INTEGER(I4B)                         :: LAKE_ID   ! Lake ID (REC code?)
  REAL(DP)                             :: LAKLAT1   ! Centroid latitude
  REAL(DP)                             :: LAKLAT2   ! Outlet latitude
  REAL(DP)                             :: LAKLON1   ! Centroid longitude
  REAL(DP)                             :: LAKLON2   ! Outlet longitude
  INTEGER(I4B)                         :: DREACHI   ! Downstream reach index
  INTEGER(I4B)                         :: DREACHK   ! Downstream reach ID
  INTEGER(I4B)                         :: DLAKE_I   ! Downstream lake index
  INTEGER(I4B)                         :: DLAKE_K   ! Downstream lake ID
 ENDTYPE LAKTOPO
 TYPE(LAKTOPO), DIMENSION(:), POINTER  :: LKTOPO  ! Lake topology
 ! Reach Parameter multipliers
 TYPE MLAKPAR
  REAL(DP)                             :: MAREAREF
  REAL(DP)                             :: MHE2AR_C
  REAL(DP)                             :: MHE2AR_D
  REAL(DP)                             :: MHGHTLOW
  REAL(DP)                             :: MHGHTECO
  REAL(DP)                             :: MHGHTSPL
  REAL(DP)                             :: MDSCHECO
  REAL(DP)                             :: MDSCHSPL
  REAL(DP)                             :: MRATECVA
  REAL(DP)                             :: MRATECVB
 ENDTYPE MLAKPAR
 TYPE(MLAKPAR), DIMENSION(:), POINTER  :: MULPARM  ! Reach Parameters
END MODULE lake_param
