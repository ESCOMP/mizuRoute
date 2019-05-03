MODULE var_lookup
 ! defines named variables used to index array elements
 USE nrtype
 implicit none
 private
 ! define missing value
 integer(i4b),parameter     :: imiss = -999      ! used to initialize named variables
 ! ***********************************************************************************************************
 ! (1) define variables desired for each HRU
 ! ***********************************************************************************************************
 type, public  ::  iLook_HRU
  integer(i4b)     :: xLon = 1          ! longitude
  integer(i4b)     :: yLat = 2          ! latitude
  integer(i4b)     :: elev = 3          ! elevation
  integer(i4b)     :: area = 4          ! basin area
 endtype iLook_HRU
 ! ***********************************************************************************************************
 ! (2) define variables desired for each HRU
 ! ***********************************************************************************************************
 type, public  ::  iLook_SEG
  integer(i4b)     :: bLon   = 1        ! start longitude of segment
  integer(i4b)     :: bLat   = 2        ! start latitude of segment
  integer(i4b)     :: eLon   = 3        ! end longitude of segment
  integer(i4b)     :: eLat   = 4        ! end latitude of segment
  integer(i4b)     :: length = 5        ! length of segment
  integer(i4b)     :: slope  = 6        ! slope of segment
  integer(i4b)     :: upArea = 7        ! upstream area of segment
 endtype iLook_SEG
 ! ***********************************************************************************************************
 ! (3) define variables for the hru2basin mapping
 ! ***********************************************************************************************************
 type, public  ::  iLook_MAP
  integer(i4b)     :: HRUid     = 1     ! HRU id
  integer(i4b)     :: segHRUid  = 2     ! id of the stream segment below each HRU
 endtype iLook_MAP
 ! ***********************************************************************************************************
 ! (4) define variables for the network topology
 ! ***********************************************************************************************************
 type, public  ::  iLook_TOP
  integer(i4b)     :: segid     = 1     ! id of each stream segment
  integer(i4b)     :: toSegment = 2     ! id of the next downstream segment
 endtype iLook_TOP
 ! ***********************************************************************************************************
 ! (X) define size of data vectors
 ! ***********************************************************************************************************
 integer(i4b),parameter,public    :: nVarsHRU=4
 integer(i4b),parameter,public    :: nVarsSEG=7
 integer(i4b),parameter,public    :: nVarsMAP=2
 integer(i4b),parameter,public    :: nVarsTOP=2
 ! ***********************************************************************************************************
 ! (X) define data vectors
 ! ***********************************************************************************************************
 type(iLook_HRU),public,parameter :: ixHRU=iLook_HRU(1,2,3,4)
 type(iLook_SEG),public,parameter :: ixSEG=iLook_SEG(1,2,3,4,5,6,7)
 type(iLook_MAP),public,parameter :: ixMAP=iLook_MAP(1,2)
 type(iLook_TOP),public,parameter :: ixTOP=iLook_TOP(1,2)
 ! ***********************************************************************************************************
END MODULE var_lookup
