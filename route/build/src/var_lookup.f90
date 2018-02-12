MODULE var_lookup
 ! defines named variables used to index array elements
 USE nrtype
 USE nrtype, only: integerMissing  ! missing value for integers
 implicit none
 private
 ! local variables
 integer(i4b),parameter            :: ixVal=1                      ! an example integer
 integer(i4b),parameter            :: iLength=storage_size(ixVal)  ! size of the example integer

 ! ***********************************************************************************************************
 ! ** define data structures
 ! ***********************************************************************************************************
 type, public  ::  iLook_struct
  integer(i4b)     :: HRU          = integerMissing   ! HRU structure
  integer(i4b)     :: HRU2SEG      = integerMissing   ! HRU-to-segment mapping structure
  integer(i4b)     :: SEG          = integerMissing   ! stream segment structure
  integer(i4b)     :: NTOPO        = integerMissing   ! network topology structure
 endtype iLook_struct
 ! ***********************************************************************************************************
 ! ** define dimensions
 ! ***********************************************************************************************************
 type, public  ::  iLook_dims
  integer(i4b)     :: hru          = integerMissing   ! hru vector
  integer(i4b)     :: seg          = integerMissing   ! stream segment vector
  integer(i4b)     :: upHRU        = integerMissing   ! upstream HRUs
  integer(i4b)     :: upSeg        = integerMissing   ! immediate upstream segments
  integer(i4b)     :: upAll        = integerMissing   ! all upstream segments
 endtype iLook_dims
 ! ***********************************************************************************************************
 ! ** define variables desired for each HRU
 ! ***********************************************************************************************************
 type, public  ::  iLook_HRU
  integer(i4b)     :: area         = integerMissing   ! basin area
 endtype iLook_HRU
 ! ***********************************************************************************************************
 ! ** define variables for the hru2segment mapping
 ! ***********************************************************************************************************
 type, public  ::  iLook_HRU2SEG
  integer(i4b)     :: HRUid         = integerMissing  ! unique HRU id
  integer(i4b)     :: HRUindex      = integerMissing  ! HRU index
  integer(i4b)     :: hruSegId      = integerMissing  ! id of the stream segment below each HRU
  integer(i4b)     :: hruSegIndex   = integerMissing  ! index of the stream segment below each HRU
 endtype iLook_HRU2SEG
 ! ***********************************************************************************************************
 ! ** define variables desired for each HRU
 ! ***********************************************************************************************************
 type, public  ::  iLook_SEG
  ! reach properties
  integer(i4b)     :: length        = integerMissing  ! length of segment  (m)
  integer(i4b)     :: slope         = integerMissing  ! slope of segment   (-)
  integer(i4b)     :: width         = integerMissing  ! width of segment   (m)
  integer(i4b)     :: man_n         = integerMissing  ! Manning's n        (weird units)
  ! contributing HRUs
  integer(i4b)     :: hruArea       = integerMissing  ! contributing area for each HRU (m2)
  integer(i4b)     :: weight        = integerMissing  ! weight assigned to each HRU (-)
  ! unit hydrograph routing
  integer(i4b)     :: timeDelayHist = integerMissing  ! time delay histogram for each reach (s)
  integer(i4b)     :: upsLength     = integerMissing  ! length of the vector of reaches above each reach (m)
  integer(i4b)     :: basArea       = integerMissing  ! area of the local HRUs contributing to each reach (m2)
  integer(i4b)     :: upsArea       = integerMissing  ! area above the top of the reach -- zero if headwater (m2)
  integer(i4b)     :: totalArea     = integerMissing  ! basArea + upsArea -- area at the bottom of the reach (m2)
  ! lakes
  integer(i4b)     :: basUnderLake  = integerMissing  ! Area of basin under lake  (m2)
  integer(i4b)     :: rchUnderLake  = integerMissing  ! Length of reach under lake (m)
  ! constraints
  integer(i4b)     :: minFlow       = integerMissing  ! minimum environmental flow
 endtype iLook_SEG
 ! ***********************************************************************************************************
 ! ** define variables for the network topology
 ! ***********************************************************************************************************
 type, public  ::  iLook_NTOPO
  ! upstream HRUs
  integer(i4b)     :: nHRU            = integerMissing  ! number of HRUs that contribute flow to each segment
  integer(i4b)     :: hruContribIx    = integerMissing  ! indices of HRUs that contribute flow to each segment
  integer(i4b)     :: hruContribId    = integerMissing  ! ids of the vector of HRUs that contribute flow to each segment
  ! individual segments
  integer(i4b)     :: segId           = integerMissing  ! unique id of each stream segment
  integer(i4b)     :: segIndex        = integerMissing  ! index of each stream segment (1, 2, 3, ..., n)
  ! downstream segments
  integer(i4b)     :: downSegId       = integerMissing  ! unique id of the next downstream segment
  integer(i4b)     :: downSegIndex    = integerMissing  ! index of downstream reach index
  ! upstream segments
  integer(i4b)     :: upSegIds        = integerMissing  ! ids for the vector of immediate upstream stream segments
  integer(i4b)     :: upSegIndices    = integerMissing  ! indices for the vector of immediate upstream stream segments
  integer(i4b)     :: allUpSegIndices = integerMissing  ! indices of all upstream stream segments
  ! processing sequence
  integer(i4b)     :: rchOrder        = integerMissing  ! order that stream segments are processed
  ! lakes
  integer(i4b)     :: lakeId          = integerMissing  ! unique id of each lake in the river network
  integer(i4b)     :: lakeIndex       = integerMissing  ! index of each lake in the river network
  integer(i4b)     :: isLakeInlet     = integerMissing  ! flag to define if reach is a lake inlet (1=inlet, 0 otherwise)
  ! irrigation
  integer(i4b)     :: userTake        = integerMissing  ! flag to define if user takes water from reach (1=extract, 0 otherwise)
  ! testing
  integer(i4b)     :: goodBasin       = integerMissing  ! flag to define a good basin (1=good, 0=bad)
 endtype iLook_NTOPO
 ! ***********************************************************************************************************
 ! ** define data vectors
 ! ***********************************************************************************************************
 type(iLook_struct),  public,parameter :: ixStruct  = iLook_struct  (1,2,3,4)
 type(iLook_dims),    public,parameter :: ixDims    = iLook_dims    (1,2,3,4,5)
 type(iLook_HRU),     public,parameter :: ixHRU     = iLook_HRU     (1)
 type(iLook_HRU2SEG), public,parameter :: ixHRU2SEG = iLook_HRU2SEG (1,2,3,4)
 type(iLook_SEG),     public,parameter :: ixSEG     = iLook_SEG     (1,2,3,4,5,6,7,8,9,10,11,12,13,14)
 type(iLook_NTOPO)   ,public,parameter :: ixNTOPO   = iLook_NTOPO   (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
 ! ***********************************************************************************************************
 ! ** define size of data vectors
 ! ***********************************************************************************************************
 integer(i4b),parameter,public    :: nStructures  = storage_size(ixStruct )/iLength
 integer(i4b),parameter,public    :: nDimensions  = storage_size(ixDims   )/iLength
 integer(i4b),parameter,public    :: nVarsHRU     = storage_size(ixHRU    )/iLength
 integer(i4b),parameter,public    :: nVarsHRU2SEG = storage_size(ixHRU2SEG)/iLength
 integer(i4b),parameter,public    :: nVarsSEG     = storage_size(ixSEG    )/iLength
 integer(i4b),parameter,public    :: nVarsNTOPO   = storage_size(ixNTOPO  )/iLength
 ! ***********************************************************************************************************

END MODULE var_lookup
