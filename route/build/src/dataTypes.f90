module dataTypes

! used to create/save specific data types

USE nrtype,     only: i4b,i8b,dp,lgt
USE nrtype,     only: strLen   ! string length
USE public_var, only: realMissing
USE public_var, only: integerMissing

implicit none

 ! everything private unless specified otherwise
 private

 ! ---------- information on the data structures -----------------------------------------------------------

 ! data structure information
 type,public :: struct_info
  character(len=strLen)  :: structName = 'empty'          ! name of the data structure
  integer(i4b)           :: varType    = integerMissing   ! type of variable in each data structure
  integer(i4b)           :: nSpace     = integerMissing   ! length of the spatial dimension in each data structure
  integer(i4b)           :: nVars      = integerMissing   ! number of variables in each data structure
 end type struct_info

 ! ---------- information on the dimensions -----------------------------------------------------------------

 ! dimension information
 type,public :: dim_info
  character(len=strLen)  :: dimName   = 'empty'          ! name of the data structure
  integer(i4b)           :: dimId     = integerMissing   ! dimension id
  integer(i4b)           :: dimLength = integerMissing   ! dim length
 end type dim_info

 ! ---------- metadata structures --------------------------------------------------------------------------

 type, public :: var_info
   character(len=strLen)    :: varName  = 'empty'         ! variable name
   character(len=strLen)    :: varDesc  = 'empty'         ! variable description
   character(len=strLen)    :: varUnit  = 'empty'         ! variable units
   integer(i4b)             :: varType  = integerMissing  ! variable type (vectors of different size)
   logical(lgt)             :: varFile  = .true.          ! .true. if the variable should be read from a file
 end type var_info

 ! ---------- time structures ------------------------------------------------------------------------------
 type,public :: time
  integer(i4b)           :: iy       = integerMissing  ! year
  integer(i4b)           :: im       = integerMissing  ! month
  integer(i4b)           :: id       = integerMissing  ! day
  integer(i4b)           :: ih       = integerMissing  ! hour
  integer(i4b)           :: imin     = integerMissing  ! minute
  real(dp)               :: dsec     = realMissing     ! second
 endtype time

 ! ---------- states structure --------------------------------------------------------------------------
 !
 type,public :: var
  integer(i4b),  allocatable  :: array_2d_int(:,:)
  integer(i4b),  allocatable  :: array_3d_int(:,:,:)
  real(dp),      allocatable  :: array_2d_dp(:,:)
  real(dp),      allocatable  :: array_3d_dp(:,:,:)
  logical(lgt),  allocatable  :: array_2d_lgt(:,:)
  logical(lgt),  allocatable  :: array_3d_lgt(:,:,:)
 end type var

 type,public :: states
  type(var),     allocatable :: var(:)
 end type states

 ! ---------- output netcdf structure --------------------------------------------------------------------------
 !
 type,public :: nc
   character(len=strLen)   :: ncname = 'empty'         ! netcdf name
   integer(i4b)            :: ncid   = integerMissing  ! netcdf id
   integer(i4b)            :: status = integerMissing  ! status: 1=defined, 2=open, 3=closed
 end type nc

 ! ---------- basin data structures ----------------------------------------------------------------------
 ! segIndex points to the segment in the entire river network data
 type,public :: subdomain
  character(32)              :: pfaf                  ! subbasin pfaf code - mainstem starting "-"
  integer(i4b)               :: basinType             ! basin type identifier: tributary->1, mainstem->2
  integer(i4b)               :: outletIndex           ! reach index of a domain outlet
  integer(i4b),  allocatable :: segIndex(:)           ! reach indices within a subbasin
  integer(i4b),  allocatable :: hruIndex(:)           ! hru indices within a subbasin
end type subdomain

 type,public :: reach
  integer(i4b), allocatable :: segIndex(:)           ! index of segment index
  integer(i4b)              :: nRch                  ! number of reach
 end type reach

 ! Data structures to hold mainstem and independent tributary reaches separately
 ! Used for openMP
 type,public :: subbasin_omp
   type(reach), allocatable :: branch(:)
 end type subbasin_omp

 ! Data structures to hold mainstem and independent tributary reaches separately
 ! Used for openMP
 type,public :: mslevel
  type(reach), allocatable :: mainstem(:)            ! mainstem reaches
 end type mslevel

 type,public :: basin
  integer(i4b)                 :: outIndex             ! index of outlet segment based on segment array
  type(mslevel), allocatable   :: level(:)             ! mainstem reach
  type(reach), allocatable     :: tributary(:)         ! index of tributary outlet segment
 end type basin

 ! ---------- general data structures ----------------------------------------------------------------------

 ! ** double precision type
 ! dat vector
 type, public :: dlength
  real(dp),allocatable                :: dat(:)    ! dat(:)
 endtype dlength
 ! var vector
 type, public :: var_dlength
  type(dlength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_dlength

 ! ** integer type
 ! dat vector
 type, public :: ilength
  integer(i4b),allocatable            :: dat(:)    ! dat(:)
 endtype ilength
 ! var vector
 type, public :: var_ilength
  type(ilength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_ilength

 ! ** character type
 ! dat vector
 type, public :: clength
  character(len=32),allocatable       :: dat(:)    ! dat(:)
 endtype clength
 ! var vector
 type, public :: var_clength
  type(clength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_clength

 ! ---------- mapping data structures ----------------------------------------------------------------------

 ! data to remap runoff hru to river network hrus
 type, public :: remap
   ! information in the mapping file
   integer(i8b)             , allocatable  :: hru_id(:)    ! Id of hrus associated with river network (="hru")
   integer(i8b)             , allocatable  :: qhru_id(:)   ! Id of hrus associated with runoff simulation (="qhru")
   integer(i4b)             , allocatable  :: num_qhru(:)  ! number of "qhru" within "hru"
   integer(i4b)             , allocatable  :: i_index(:)   ! Index in the y dimension of the runoff grid (starting with 1,1 in LL corner)
   integer(i4b)             , allocatable  :: j_index(:)   ! Index in the x dimension of the runoff grid (starting with 1,1 in LL corner)
   real(dp)                 , allocatable  :: weight(:)    ! area weight of "qhru" intersecting "hru"
   ! ancillary index vectors
   integer(i4b)             , allocatable  :: hru_ix(:)    ! Index of hrus associated with river network (="hru")
   integer(i4b)             , allocatable  :: qhru_ix(:)   ! Index of hrus associated with runoff simulation (="qhru")
 end type remap

 ! simulated runoff data
 type, public :: runoff
   integer(i4b)                            :: nTime         ! number of time steps
   integer(i4b)                            :: nSpace(1:2)   ! number of spatial dimension
   real(dp)                                :: time          ! time variable at one time step
   real(dp)                 , allocatable  :: qsim(:)       ! runoff(HM_HRU) at one time step (size: nSpace(1))
   real(dp)                 , allocatable  :: qsim2D(:,:)   ! runoff(x,y) at one time step (size: /nSpace(1),nSpace(2)/)
   integer(i8b)             , allocatable  :: hru_id(:)     ! id of HM_HRUs or RN_HRUs at which runoff is stored (size: nSpace(1))
   integer(i4b)             , allocatable  :: hru_ix(:)     ! Index of RN_HRUs associated with river network (used only if HM_HRUs = RN_HRUs)
   real(dp)                 , allocatable  :: basinRunoff(:)! remapped river network catchment runoff (size: number of nHRU)
 end type runoff

 ! ---------- reach parameters ----------------------------------------------------------------------------

 ! Reach Parameters
 type, public ::  RCHPRP
  real(DP)                                :: R_SLOPE
  real(DP)                                :: R_MAN_N
  real(DP)                                :: R_WIDTH
  real(DP)                                :: RLENGTH
  real(DP)                                :: UPSAREA  ! upstream area (zero if headwater basin)
  real(DP)                                :: BASAREA  ! local basin area
  real(DP)                                :: TOTAREA  ! UPSAREA + BASAREA
  real(DP)                                :: MINFLOW  ! minimum environmental flow
 end type RCHPRP

 ! River Network topology
 type, public :: RCHTOPO
  integer(I4B)                               :: REACHIX  ! Reach index (0,1,2,...,nrch-1)
  integer(I4B)                               :: REACHID  ! Reach ID (REC code)
  real(DP)                                   :: RCHLAT1  ! Start latitude
  real(DP)                                   :: RCHLAT2  ! End latitude
  real(DP)                                   :: RCHLON1  ! Start longitude
  real(DP)                                   :: RCHLON2  ! End longitude
  real(DP),    dimension(:),allocatable      :: UPSLENG  ! total upstream length
  integer(I4B)                               :: DREACHI  ! Immediate Downstream reach index
  integer(I4B)                               :: DREACHK  ! Immediate Downstream reach ID
  integer(I4B),dimension(:),allocatable      :: UREACHI  ! Immediate Upstream reach indices
  integer(I4B),dimension(:),allocatable      :: UREACHK  ! Immediate Upstream reach IDs
  integer(I4B),dimension(:),allocatable      :: RCHLIST  ! all upstream reach indices
  integer(I4B),dimension(:),allocatable      :: HRUID    ! all contributing HRU IDs
  integer(I4B),dimension(:),allocatable      :: HRUIX    ! all contributing HRU indices
  real(DP),    dimension(:),allocatable      :: HRUWGT   ! areal weight for contributing HRUs
  logical(lgt),dimension(:),allocatable      :: goodBas  ! Flag to denote a good basin
  character(len=32),dimension(:),allocatable :: pfafCode ! pfafstetter code
  integer(I4B)                               :: RHORDER  ! Processing sequence
  real(dp)    ,dimension(:),allocatable      :: UH       ! Unit hydrograph for upstream
  integer(I4B)                               :: LAKE_IX  ! Lake index (0,1,2,...,nlak-1)
  integer(I4B)                               :: LAKE_ID  ! Lake ID (REC code?)
  real(DP)                                   :: BASULAK  ! Area of basin under lake
  real(DP)                                   :: RCHULAK  ! Length of reach under lake
  LOGICAL(LGT)                               :: LAKINLT  ! .TRUE. if reach is lake inlet, .FALSE. otherwise
  LOGICAL(LGT)                               :: USRTAKE  ! .TRUE. if user takes from reach, .FALSE. otherwise
 end type RCHTOPO

 ! ---------- kinematic wave states (collection of particles) ---------------------------------

 ! Individual flow particles
 ! NOTE: type could possibly be private
 TYPE, public :: FPOINT
  REAL(DP)                             :: QF       ! Flow
  REAL(DP)                             :: QM       ! Modified flow
  REAL(DP)                             :: TI       ! initial time of point in reach
  REAL(DP)                             :: TR       ! time point expected to exit reach
  LOGICAL(LGT)                         :: RF       ! routing flag (T if point has exited)
 END TYPE FPOINT

 ! Collection of flow points within a given reach
 TYPE, public :: KREACH
  TYPE(FPOINT),allocatable             :: KWAVE(:)
 END TYPE KREACH

 ! ---------- irf states (future flow series ) ---------------------------------

 ! Future flow series
 TYPE, public :: IRFREACH
  REAL(DP), allocatable                :: qfuture(:)    ! runoff volume in future time steps for IRF routing (m3/s)
 END TYPE IRFREACH

 ! ---------- reach fluxes --------------------------------------------------------------------

 ! fluxes in each reach
 TYPE, public :: STRFLX
  REAL(DP), allocatable                :: QFUTURE(:)        ! runoff volume in future time steps (m3/s)
  REAL(DP), allocatable                :: QFUTURE_IRF(:)    ! runoff volume in future time steps for IRF routing (m3/s)
  REAL(DP)                             :: BASIN_QI          ! instantaneous runoff volume from the local basin (m3/s)
  REAL(DP)                             :: BASIN_QR(0:1)     ! routed runoff volume from the local basin (m3/s)
  REAL(DP)                             :: UPSBASIN_QR       ! routed runoff depth from the upstream basins (m/s)
  REAL(DP)                             :: BASIN_QR_IRF(0:1) ! routed runoff volume from all the upstream basin (m3/s)
  REAL(DP)                             :: REACH_Q           ! time-step average streamflow (m3/s)
  REAL(DP)                             :: REACH_Q_IRF       ! time-step average streamflow (m3/s) from IRF routing
  REAL(DP)                             :: UPSTREAM_QI       ! sum of upstream streamflow (m3/s)
  REAL(DP)                             :: TAKE              ! average take
  logical(lgt)                         :: CHECK_IRF         ! .true. if the reach is routed
 ENDTYPE STRFLX

 ! ---------- lake data types -----------------------------------------------------------------

 ! Lake Parameters
 TYPE, public :: LAKPRP
  REAL(DP)                             :: AREAREF           ! lake area
  REAL(DP)                             :: LAKREFLEV         ! lake elevation
  REAL(DP)                             :: LAKAVGLEV         ! lake average level (for initialization)
  REAL(DP)                             :: HE2AR_C           ! water height-surface area parameter
  REAL(DP)                             :: HE2AR_D           ! water height-surface area parameter
  REAL(DP)                             :: HGHTLOW           ! minimum water height for discharge
  REAL(DP)                             :: HGHTECO           ! minimum height for ecological concerns
  REAL(DP)                             :: HGHTSPL           ! spillway height
  REAL(DP)                             :: DSCHECO           ! discharge at "ecological" height
  REAL(DP)                             :: DSCHSPL           ! discharge at spillway height
  REAL(DP)                             :: RATECVA           ! discharge rating curve parameter
  REAL(DP)                             :: RATECVB           ! discharge rating curve parameter
 ENDTYPE LAKPRP

 ! Lake topology
 TYPE, public :: LAKTOPO
  INTEGER(I4B)                         :: LAKE_IX           ! Lake index (0,1,2,...,nlak-1)
  INTEGER(I4B)                         :: LAKE_ID           ! Lake ID (REC code?)
  REAL(DP)                             :: LAKLAT1           ! Centroid latitude
  REAL(DP)                             :: LAKLAT2           ! Outlet latitude
  REAL(DP)                             :: LAKLON1           ! Centroid longitude
  REAL(DP)                             :: LAKLON2           ! Outlet longitude
  INTEGER(I4B)                         :: DREACHI           ! Downstream reach index
  INTEGER(I4B)                         :: DREACHK           ! Downstream reach ID
  INTEGER(I4B)                         :: DLAKE_I           ! Downstream lake index
  INTEGER(I4B)                         :: DLAKE_K           ! Downstream lake ID
 ENDTYPE LAKTOPO

 ! Lake fluxes
 TYPE, public :: LKFLX
  REAL(DP)                             :: LAKE_Qav          ! lake discharge (average over time step) (m3 s-1)
  REAL(DP)                             :: LAKE_Q            ! lake discharge (instantaneous) (m3 s-1)
  REAL(DP)                             :: LAKE_P            ! lake precipitation (m3)
  REAL(DP)                             :: LAKE_E            ! lake evaporation (m3)
  REAL(DP)                             :: LAKE_I            ! inflow to lake (m3 s-1)
 ENDTYPE LKFLX

END MODULE dataTypes

MODULE objTypes

 USE nrtype,     only: i4b,dp,lgt
 USE nrtype,     only: strLen   ! string length
 USE public_var, only: realMissing
 USE public_var, only: integerMissing

 ! define derived type for model variables, including name, description, and units
 type, public :: var_info_new
   character(len=strLen)    :: varName  = 'empty'         ! variable name
   character(len=strLen)    :: varDesc  = 'empty'         ! variable description
   character(len=strLen)    :: varUnit  = 'empty'         ! variable units
   integer(i4b)             :: varType  = integerMissing  ! variable type
   integer(i4b),allocatable :: varDim(:)                  ! dimension ID associated with variable
   logical(lgt)             :: varFile  = .true.          ! .true. if the variable should be read from a file
 CONTAINS
   procedure, pass :: init
 end type var_info_new

 CONTAINS

  SUBROUTINE init(this, vName, vDesc, vUnit, vType, vDim, vFile)
    implicit none
    class(var_info_new)                 :: this
    character(*),            intent(in) :: vName    ! variable name
    character(*),            intent(in) :: vDesc    ! variable description
    character(*),            intent(in) :: vUnit    ! variable units
    integer(i4b),            intent(in) :: vType    ! variable type
    integer(i4b),            intent(in) :: vDim(:)  ! dimension ID
    logical(lgt),            intent(in) :: vFile    ! .true. if the variable should be read from a file
    integer(i4b)                        :: n        ! size of dimension

    n = size(vDim)
    allocate(this%varDim(n))
    this%varName      = vName
    this%varDesc      = vDesc
    this%varUnit      = vUnit
    this%varType      = vType
    this%varDim(1:n) = vDim(1:n)
    this%varFile      = vFile
  END SUBROUTINE init

END MODULE objTypes
