MODULE reachparam
 USE nrtype
 IMPLICIT NONE
 SAVE
 ! Data for Unit Hydrograph 
 type TDH 
   real(DP), dimension(:), pointer     :: UH_DATA  ! Data for Unit hydrograph (added by NM)
 end type TDH 

 ! Reach Parameters
 type RCHPRP
  real(DP)                             :: R_SLOPE
  real(DP)                             :: R_MAN_N
  real(DP)                             :: R_WIDTH
  real(DP)                             :: RLENGTH
  real(DP)                             :: UPSAREA  ! upstream area (zero if headwater basin)
  real(DP)                             :: BASAREA  ! local basin area
  real(DP)                             :: TOTAREA  ! UPSAREA + BASAREA
  real(DP)                             :: MINFLOW  ! minimum environmental flow
 end type RCHPRP
 type(RCHPRP), dimension(:), pointer   :: RPARAM   ! Reach Parameters

 ! River Network topology
 type RCHTOPO
  integer(I4B)                         :: REACHIX  ! Reach index (0,1,2,...,nrch-1)
  integer(I4B)                         :: REACHID  ! Reach ID (REC code)
  real(DP)                             :: RCHLAT1  ! Start latitude
  real(DP)                             :: RCHLAT2  ! End latitude
  real(DP)                             :: RCHLON1  ! Start longitude
  real(DP)                             :: RCHLON2  ! End longitude
  real(DP),    dimension(:),pointer    :: UPSLENG  ! total upstream length  (added by NM)
  integer(I4B)                         :: DREACHI  ! Immediate Downstream reach index
  integer(I4B)                         :: DREACHK  ! Immediate Downstream reach ID
  integer(I4B),dimension(:),pointer    :: UREACHI  ! Immediate Upstream reach indices
  integer(I4B),dimension(:),pointer    :: UREACHK  ! Immediate Upstream reach IDs
  logical(lgt),dimension(:),pointer    :: goodBas  ! Flag to denote a good basin
  integer(I4B)                         :: RHORDER  ! Processing sequence 
  integer(I4B),dimension(:),pointer    :: RCHLIST  ! List of reaches upstream inices
  type(TDH),   dimension(:),pointer    :: UH       ! Unit hydrograph for upstream (added by NM)
  integer(I4B)                         :: LAKE_IX  ! Lake index (0,1,2,...,nlak-1)
  integer(I4B)                         :: LAKE_ID  ! Lake ID (REC code?)
  real(DP)                             :: BASULAK  ! Area of basin under lake
  real(DP)                             :: RCHULAK  ! Length of reach under lake
  LOGICAL(LGT)                         :: LAKINLT  ! .TRUE. if reach is lake inlet, .FALSE. otherwise
  LOGICAL(LGT)                         :: USRTAKE  ! .TRUE. if user takes from reach, .FALSE. otherwise
 end type RCHTOPO
 type(RCHTOPO), dimension(:), pointer  :: NETOPO   ! River Network topology

 ! Reach Parameter multipliers
 type(RCHPRP), dimension(:,:), pointer :: MULTRCH  ! Multiplier for reach parameters

END MODULE reachparam
