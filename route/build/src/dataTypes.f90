MODULE dataTypes

! used to create specific data types

USE nrtype
USE public_var, ONLY: realMissing
USE public_var, ONLY: integerMissing
USE datetime_data, ONLY: datetime

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

 ! define derived type for model variables, including name, description, and units
 type,public :: var_info
  character(len=strLen)  :: varName  = 'empty'         ! variable name
  character(len=strLen)  :: varDesc  = 'empty'         ! variable description
  character(len=strLen)  :: varUnit  = 'empty'         ! variable units
  integer(i4b)           :: varType  = integerMissing  ! variable type (vectors of different size)
  logical(lgt)           :: varFile  = .true.          ! .true. if the variable should be read from a file
 endtype var_info

 ! ---------- gauge metadata structures --------------------------------------------------------------------------

 type, public :: gage
   integer(i4b)                   :: nGage
   character(len=30), allocatable :: gageID(:)
   integer(i8b),      allocatable :: reachID(:)
 end type gage

 ! ---------- basin data structures ----------------------------------------------------------------------

 ! segIndex points to the segment in the entire river network data
 ! segOrder is order within subset of mainstem segments or tributary segments
 type,public :: subbasin_mpi
  integer(i4b)               :: basinType             ! basin type identifier: tributary->1, mainstem->2
  integer(i4b)               :: idNode                ! node ID
  integer(i4b)               :: outletIndex           ! reach index of a domain outlet
  integer(i4b),  allocatable :: segIndex(:)           ! reach indices within a subbasin
  integer(i4b),  allocatable :: hruIndex(:)           ! hru indices within a subbasin
 end type subbasin_mpi

 ! -- Data structures to reach indices within each leg of stream order separately
 !    Used for openMP
 type,public :: reach
  integer(i4b), allocatable :: segIndex(:)           ! index of segment index
  integer(i4b)              :: nRch                  ! number of reach
 end type reach

 type,public :: subbasin_omp
   type(reach), allocatable :: branch(:)
 end type subbasin_omp

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


 ! ---------- forcing input file strcuture -----------------------------------------------------------------

 ! -- input file name and strcuture for nc files
  type, public ::  inFileInfo
   integer(i4b)                            :: nTime            ! number of time steps in a netCDF
   integer(i4b)                            :: iTimebound(1:2)  ! cumulative time indices of 1st and last time steps in a netCDf
   real(dp)                 , allocatable  :: timeVar(:)       ! elapsed seconds from reference datetime of netCDF
   type(datetime)                          :: refdatetime      ! reference datetime of netCDF
   character(len=strLen)                   :: infilename       ! name of a netCDF
   character(len=strLen)                   :: calendar         ! calendar used in a netCDF
   character(len=strLen)                   :: unit             ! time unit used in a netCDF
  end type inFileInfo


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

 ! mapping time step between two time series e.g., simulation time step and runoff time step for one simulation time step
 ! runoff value used at each simulaition time step is weighted average of runoff(s) across input time step(s)
 ! (can be one, but several if simulation step is larger than runoff time step)
 ! This store indices of netCDF files in input file streams and time step weight of input variables from each input timestep
 type, public :: map_time
   integer(i4b), allocatable :: iFile(:) ! index of input netCDF file
   integer(i4b), allocatable :: iTime(:) ! index of time steps within the corresponding netCDF
   real(dp),     allocatable :: frac(:)  ! weight
 end type map_time

 ! ---------- input forcing data  ----------------------------------------------------------------------

 type, public :: inputData
   integer(i4b)                            :: nSpace(1:2)     ! number of spatial dimension
   real(dp)                 , allocatable  :: sim(:)          ! flux simulation (HM_HRU) at one time step (size: nSpace(1))
   real(dp)                 , allocatable  :: sim2d(:,:)      ! flux simulation (x,y) at one time step (size: /nSpace(1),nSpace(2)/)
   real(dp)                                :: fillvalue       ! fillvalue
 end type inputData

 type, public, extends(inputData) :: runoff  ! runoff data
   integer(i8b)             , allocatable  :: hru_id(:)       ! id of HM_HRUs or RN_HRUs at which runoff is stored (size: nSpace(1))
   integer(i4b)             , allocatable  :: hru_ix(:)       ! Index of RN_HRUs associated with river network (used only if HM_HRUs = RN_HRUs)
   real(dp)                 , allocatable  :: basinRunoff(:)  ! remapped river network catchment runoff [depth/time] (size: number of nHRU)
   real(dp)                 , allocatable  :: basinEvapo(:)   ! remapped river network catchment evaporation [depth/time] (size: number of nHRU)
   real(dp)                 , allocatable  :: basinPrecip(:)  ! remapped river network catchment precipitation [depth/time] (size: number of nHRU)
 end type runoff

 type, public, extends(inputData) :: wm  ! water-management
   integer(i8b)             , allocatable  :: seg_id(:)       ! id of reach in data (size: nSpace)
   integer(i4b)             , allocatable  :: seg_ix(:)       ! Index of river network reach IDs corresponding reach ID in data
   real(dp)                 , allocatable  :: flux_wm(:)      ! allocated flux to existing river network using sort_flux [m3/s] (size: number of nRCH)
   real(dp)                 , allocatable  :: vol_wm(:)       ! allocated target vol to existing river network using sort_flux [m3/s] (size: number of nRCH)
 end type wm

 ! ---------- reach parameters ----------------------------------------------------------------------------

 ! -- Reach physical parameters
 type, public ::  RCHPRP
  real(dp)                                   :: R_SLOPE        ! channel slope [-]
  real(dp)                                   :: R_MAN_N        ! channel bed manning coefficient [-]
  real(dp)                                   :: R_WIDTH        ! channel width [m]
  real(dp)                                   :: R_DEPTH        ! channel bankfull depth [m]
  real(dp)                                   :: RLENGTH        ! channel length [m]
  real(dp)                                   :: FLDP_SLOPE     ! floodplain slope [-]
  real(dp)                                   :: UPSAREA        ! upstream area (zero if headwater basin)
  real(dp)                                   :: BASAREA        ! local basin area
  real(dp)                                   :: TOTAREA        ! UPSAREA + BASAREA
  real(dp)                                   :: MINFLOW        ! minimum environmental flow
  real(dp)                                   :: D03_MaxStorage ! Doll 2003; maximume storage [m3]
  real(dp)                                   :: D03_Coefficient! Doll 2003; Coefficient [?]
  real(dp)                                   :: D03_Power      ! Doll 2003; Power [-]
  real(dp)                                   :: D03_S0         ! Doll 2003; Additional parameter to represent inactive storage

  real(dp)                                   :: HYP_E_emr      ! HYPE; elevation of emergency spillway [m]
  real(dp)                                   :: HYP_E_lim      ! HYPE; elevation below which primary spillway flow is restrcited [m]
  real(dp)                                   :: HYP_E_min      ! HYPE; elevation below which outflow is zero [m]
  real(dp)                                   :: HYP_E_zero     ! HYPE; elevation at which the reservoir storage is set to 0, lake/reservoir bottom elevation
  real(dp)                                   :: HYP_Qrate_emr  ! HYPE; emergency rate of flow for each unit of elevation above HYP_E_emr [m3/s]
  real(dp)                                   :: HYP_Erate_emr  ! HYPE; power for the rate of flow for each unit of elevation above HYP_E_emr [-]
  real(dp)                                   :: HYP_Qrate_prim ! HYPE; the average yearly or long term output from primary spillway [m3/s]
  real(dp)                                   :: HYP_Qrate_amp  ! HYPE; amplitude of the Qrate_main [-]
  integer(i4b)                               :: HYP_Qrate_phs  ! HYPE; phase of the Qrate_main based on the day of the year [-]; default 100
  logical(lgt)                               :: HYP_prim_F     ! HYPE; if the reservoir has a primary spillway then set to 1 otherwise 0
  real(dp)                                   :: HYP_A_avg      ! HYPE; average area for the lake; this might not be used if bathymetry is provided [m]

  real(dp)                                   :: H06_Smax       ! Hanasaki 2006; maximume reservoir storage [m3]
  real(dp)                                   :: H06_alpha      ! Hanasaki 2006; fraction of active storage compared to total storage [-]
  real(dp)                                   :: H06_envfact    ! Hanasaki 2006; fraction of inflow that can be used to meet demand [-]
  real(dp)                                   :: H06_S_ini      ! Hanasaki 2006; initial storage used for initial estimation of release coefficient [m3]
  real(dp)                                   :: H06_c1         ! Hanasaki 2006; coefficient 1 for target release for irrigation reseroir [-]
  real(dp)                                   :: H06_c2         ! Hanasaki 2006; coefficient 2 for target release for irrigation reseroir [-]
  real(dp)                                   :: H06_exponent   ! Hanasaki 2006; Exponenet of actual release for "within-a-year" reservoir [-]
  real(dp)                                   :: H06_denominator! Hanasaki 2006; Denominator of actual release for "within-a-year" reservoir [-]
  real(dp)                                   :: H06_c_compare  ! Hanasaki 2006; Criterion for distinguish of "within-a-year" or "multi-year" reservoir [-]
  real(dp)                                   :: H06_frac_Sdead ! Hanasaki 2006; Fraction of dead storage to maximume storage [-]
  real(dp)                                   :: H06_E_rel_ini  ! Hanasaki 2006; Initial release coefficient [-]
  real(dp)                                   :: H06_I_Jan      ! Hanasaki 2006; Average January   inflow [m3/s]
  real(dp)                                   :: H06_I_Feb      ! Hanasaki 2006; Average Februrary inflow [m3/s]
  real(dp)                                   :: H06_I_Mar      ! Hanasaki 2006; Average March     inflow [m3/s]
  real(dp)                                   :: H06_I_Apr      ! Hanasaki 2006; Average April     inflow [m3/s]
  real(dp)                                   :: H06_I_May      ! Hanasaki 2006; Average May       inflow [m3/s]
  real(dp)                                   :: H06_I_Jun      ! Hanasaki 2006; Average June      inflow [m3/s]
  real(dp)                                   :: H06_I_Jul      ! Hanasaki 2006; Average July      inflow [m3/s]
  real(dp)                                   :: H06_I_Aug      ! Hanasaki 2006; Average August    inflow [m3/s]
  real(dp)                                   :: H06_I_Sep      ! Hanasaki 2006; Average September inflow [m3/s]
  real(dp)                                   :: H06_I_Oct      ! Hanasaki 2006; Average October   inflow [m3/s]
  real(dp)                                   :: H06_I_Nov      ! Hanasaki 2006; Average November  inflow [m3/s]
  real(dp)                                   :: H06_I_Dec      ! Hanasaki 2006; Average December  inflow [m3/s]
  real(dp)                                   :: H06_D_Jan      ! Hanasaki 2006; Average January   demand [m3/s]
  real(dp)                                   :: H06_D_Feb      ! Hanasaki 2006; Average Februrary demand [m3/s]
  real(dp)                                   :: H06_D_Mar      ! Hanasaki 2006; Average March     demand [m3/s]
  real(dp)                                   :: H06_D_Apr      ! Hanasaki 2006; Average April     demand [m3/s]
  real(dp)                                   :: H06_D_May      ! Hanasaki 2006; Average May       demand [m3/s]
  real(dp)                                   :: H06_D_Jun      ! Hanasaki 2006; Average June      demand [m3/s]
  real(dp)                                   :: H06_D_Jul      ! Hanasaki 2006; Average July      demand [m3/s]
  real(dp)                                   :: H06_D_Aug      ! Hanasaki 2006; Average Agust     demand [m3/s]
  real(dp)                                   :: H06_D_Sep      ! Hanasaki 2006; Average September demand [m3/s]
  real(dp)                                   :: H06_D_Oct      ! Hanasaki 2006; Average October   demand [m3/s]
  real(dp)                                   :: H06_D_Nov      ! Hanasaki 2006; Average November  demand [m3/s]
  real(dp)                                   :: H06_D_Dec      ! Hanasaki 2006; Average December  demand [m3/s]
  integer(i4b)                               :: H06_purpose    ! Hanasaki 2006; reservoir purpose; (0= non-irrigation, 1=irrigation) [-]
  logical(lgt)                               :: H06_I_mem_F    ! Hanasaki 2006; Flag to transition to modelled inflow [-]
  logical(lgt)                               :: H06_D_mem_F    ! Hanasaki 2006; Flag to transition to modelled/provided demand [-]
  integer(i4b)                               :: H06_I_mem_L    ! Hanasaki 2006; Memory length in years for inflow [year]
  integer(i4b)                               :: H06_D_mem_L    ! Hanasaki 2006; Memory length in years for demand [year]
 end type RCHPRP

 ! -- River Network topology
 type, public :: RCHTOPO
  integer(i4b)                               :: REACHIX      ! Reach index (1,2,...,nrch)
  integer(i4b)                               :: REACHID      ! Reach ID (REC code)
  integer(i4b)                               :: DREACHI      ! Immediate Downstream reach index
  integer(i4b)                               :: DREACHK      ! Immediate Downstream reach ID
  integer(i4b),dimension(:),allocatable      :: UREACHI      ! Immediate Upstream reach indices
  integer(i4b),dimension(:),allocatable      :: UREACHK      ! Immediate Upstream reach IDs
  integer(i4b),dimension(:),allocatable      :: RCHLIST      ! all upstream reach indices
  integer(i4b),dimension(:),allocatable      :: HRUID        ! all contributing HRU IDs
  integer(i4b),dimension(:),allocatable      :: HRUIX        ! all contributing HRU indices
  real(dp),    dimension(:),allocatable      :: HRUWGT       ! areal weight for contributing HRUs
  logical(lgt),dimension(:),allocatable      :: goodBas      ! Flag to denote a good basin
  character(len=32),dimension(:),allocatable :: pfafCode     ! pfafstetter code
  integer(i4b)                               :: RHORDER      ! Processing sequence
  real(dp)    ,dimension(:),allocatable      :: UH           ! Unit hydrograph for upstream
  integer(i4b)                               :: LAKE_IX      ! Lake index (1,2,...,nlak)
  integer(i4b)                               :: LAKE_ID      ! Lake ID (REC code)
  real(dp)                                   :: BASULAK      ! Area of basin under lake
  real(dp)                                   :: RCHULAK      ! Length of reach under lake
  logical(lgt)                               :: LAKINLT      ! .TRUE. if reach is lake inlet, .FALSE. otherwise
  logical(lgt)                               :: USRTAKE      ! .TRUE. if user takes from reach, .FALSE. otherwise
  logical(lgt)                               :: ISLAKE       ! .TRUE. if the object is a lake
  logical(lgt)                               :: LAKETARGVOL  ! .TRUE. if the lake follow a given target volume
  integer(i4b)                               :: LAKEMODELTYPE! 1=Doll, 2=Hanasaki, 3=HYPE else=non-parameteric
 end type RCHTOPO

 ! arbitrary element-to-element connection
 type, public :: commLink
   integer(i4b) :: srcTask   ! source task
   integer(i4b) :: destTask  ! destination task
   integer(i4b) :: srcIndex  ! source array index
   integer(i4b) :: destIndex ! destination array index
 end type commLink

 ! ---------- Reach restart variables--------------------------------------------------------------------

 ! - Lagrangian kinematic wave states (collection of particles)
 ! Individual flow particles
 type, public :: FPOINT
  real(dp)                                   :: QF           ! Flow
  real(dp)                                   :: QM           ! Modified flow
  real(dp)                                   :: TI           ! initial time of point in reach
  real(dp)                                   :: TR           ! time point expected to exit reach
  logical(lgt)                               :: RF           ! routing flag (T if point has exited)
 end type FPOINT

 ! Collection of flow points within a given reach
 type, public :: kwtRCH
  type(FPOINT),allocatable             :: KWAVE(:)
 end type kwtRCH

 ! ---------- irf states (future flow series ) ---------------------------------
 ! Future flow series
 type, public :: irfRCH
  real(dp), allocatable                :: qfuture(:)    ! runoff volume in future time steps for IRF routing (m3/s)
 end type irfRCH

 ! ---------- computational node for kw, dw, and mc -----------------------------
 type, public :: cMolecule
   integer(i4b)           :: KW_ROUTE
   integer(i4b)           :: MC_ROUTE
   integer(i4b)           :: DW_ROUTE
 end type cMolecule

 type, public :: SUBRCH
   real(dp), allocatable  :: Q(:)        ! Discharge at sub-reaches at current step (m3/s)
   real(dp), allocatable  :: A(:)        ! Flow area at sub-reach at current step (m2)
   real(dp), allocatable  :: H(:)        ! Flow height at sub-reach at current step (m)
 end type SUBRCH

 type, public :: kwRch
   type(SUBRCH)    :: molecule
 end type kwRCH

 type, public :: mcRch
   type(SUBRCH)    :: molecule
 end type mcRCH

 type, public :: dwRch
   type(SUBRCH)    :: molecule
 end type dwRCH

 type, public :: STRSTA
   type(irfRCH)       :: IRF_ROUTE
   type(kwtRCH)       :: LKW_ROUTE
   type(kwRCH)        :: KW_ROUTE
   type(mcRCH)        :: MC_ROUTE
   type(dwRCH)        :: DW_ROUTE
 end type STRSTA


 ! ---------- reach fluxes --------------------------------------------------------------------
 type, public :: hydraulic
   real(dp)        :: REACH_ELE              ! water height at current time step [m]
   real(dp)        :: REACH_Q                ! discharge at current time step [m3/s]
   real(dp)        :: REACH_VOL(0:1)         ! water volume at previous and current time steps [m3]
   real(dp)        :: REACH_WM_FLUX_actual   ! water management fluxes to and from each reach [m3/s]
   real(dp)        :: WB                     ! reach water balance error [m3]
 end type hydraulic

 ! fluxes and states in each reach
 type, public :: strflx
  real(dp), allocatable                :: QFUTURE(:)             ! runoff volume in future time steps [m3/s]
  real(dp), allocatable                :: QFUTURE_IRF(:)         ! runoff volume in future time steps for IRF routing [m3/s]
  real(dp), allocatable                :: QPASTUP_IRF(:,:)       ! runoff volume in the past time steps for lake upstream [m3/s]
  real(dp), allocatable                :: DEMANDPAST_IRF(:,:)    ! demand volume for lake [m3/s]
  real(dp)                             :: BASIN_QI               ! instantaneous runoff volume from the local basin [m3/s]
  real(dp)                             :: BASIN_QR(0:1)          ! routed runoff volume from the local basin [m3/s]
  type(hydraulic), allocatable         :: ROUTE(:)               ! reach fluxes and states for each routing method
  real(dp)                             :: REACH_WM_FLUX          ! water management fluxes to and from each reach [m3/s]
  real(dp)                             :: REACH_WM_VOL           ! target volume from the second water management file [m3]
  real(dp)                             :: Qobs                   ! observed discharge [m3/s]
  real(dp)                             :: basinEvapo             ! remapped river network catchment Evaporation [unit] (size: number of nHRU)
  real(dp)                             :: basinPrecip            ! remapped river network catchment Precipitation [unit] (size: number of nHRU)
 end type strflx

 ! ---------- lake data types -----------------------------------------------------------------

 ! Lake Parameters
 TYPE, public :: LAKPRP
  real(dp)                             :: AREAREF            ! lake area
  real(dp)                             :: LAKREFLEV          ! lake elevation
  real(dp)                             :: LAKAVGLEV          ! lake average level (for initialization)
  real(dp)                             :: HE2AR_C            ! water height-surface area parameter
  real(dp)                             :: HE2AR_D            ! water height-surface area parameter
  real(dp)                             :: HGHTLOW            ! minimum water height for discharge
  real(dp)                             :: HGHTECO            ! minimum height for ecological concerns
  real(dp)                             :: HGHTSPL            ! spillway height
  real(dp)                             :: DSCHECO            ! discharge at "ecological" height
  real(dp)                             :: DSCHSPL            ! discharge at spillway height
  real(dp)                             :: RATECVA            ! discharge rating curve parameter
  real(dp)                             :: RATECVB            ! discharge rating curve parameter
 END TYPE LAKPRP

 ! Lake topology
 TYPE, public :: LAKTOPO
  integer(i4b)                         :: LAKE_IX           ! Lake index (1,2,...,nlak)
  integer(i4b)                         :: LAKE_ID           ! Lake ID (REC code)
  real(dp)                             :: LAKLAT1           ! Centroid latitude
  real(dp)                             :: LAKLAT2           ! Outlet latitude
  real(dp)                             :: LAKLON1           ! Centroid longitude
  real(dp)                             :: LAKLON2           ! Outlet longitude
  integer(i4b)                         :: DREACHI           ! Downstream reach index
  integer(i4b)                         :: DREACHK           ! Downstream reach ID
  integer(i4b)                         :: DLAKE_I           ! Downstream lake index
  integer(i4b)                         :: DLAKE_K           ! Downstream lake ID
 END TYPE LAKTOPO

 ! Lake fluxes
 TYPE, public :: LKFLX
  real(dp)                             :: LAKE_Qav          ! lake discharge (average over time step) (m3 s-1)
  real(dp)                             :: LAKE_Q            ! lake discharge (instantaneous) (m3 s-1)
  real(dp)                             :: LAKE_P            ! lake precipitation (m3)
  real(dp)                             :: LAKE_E            ! lake evaporation (m3)
  real(dp)                             :: LAKE_I            ! inflow to lake (m3 s-1)
 END TYPE LKFLX

END MODULE dataTypes

MODULE objTypes

 USE nrtype
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
    this%varDim(1:n)  = vDim(1:n)
    this%varFile      = vFile
  END SUBROUTINE init

END MODULE objTypes
