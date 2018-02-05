module popMetadat_module

! data types
USE nrtype
USE dataTypes,  only : struct_info  ! data type for metadata structure
USE dataTypes,  only : var_info     ! data type for metadata structure
USE dataTypes,  only : dim_info     ! data type for metadata structure

! variable type
USE globalData, only : varType_integer  ! named variable for an integer
USE globalData, only : varType_double   ! named variable for a double precision

! metadata on data structures
USE globalData, only : meta_struct  ! structure information
USE globalData, only : meta_dims    ! dimensions
USE globalData, only : meta_HRU     ! HRU properties
USE globalData, only : meta_HRU2SEG ! HRU-to-segment mapping
USE globalData, only : meta_SEG     ! stream segment properties
USE globalData, only : meta_NTOPO   ! network topology

! indices of named variables
USE var_lookup, only : ixStruct  , nStructures   ! index of variables for data structure
USE var_lookup, only : ixDims    , nDimensions   ! index of variables for data structure
USE var_lookup, only : ixHRU     , nVarsHRU      ! index of variables for data structure
USE var_lookup, only : ixHRU2SEG , nVarsHRU2SEG  ! index of variables for data structure
USE var_lookup, only : ixSEG     , nVarsSEG      ! index of variables for data structure
USE var_lookup, only : ixNTOPO   , nVarsNTOPO    ! index of variables for data structure

implicit none

! privacy
private
public::popMetadat
contains

 subroutine popMetadat(err,message)
 implicit none
 ! dummy variables
 integer(i4b),intent(out)       :: err           ! error code
 character(*),intent(out)       :: message       ! error message
 ! initialize error control
 err=0; message='popMetadat/'

 ! ---------- define data structures -----------------------------------------------------------------------------------------------------------

 ! structure index                            name      variable type    length of spatial dimension, number of variables
 meta_struct(ixStruct%HRU    ) = struct_info('HRU',     varType_double,  integerMissing,              nVarsHRU)
 meta_struct(ixStruct%HRU2SEG) = struct_info('HRU2SEG', varType_integer, integerMissing,              nVarsHRU2SEG)
 meta_struct(ixStruct%SEG    ) = struct_info('SEG',     varType_double,  integerMissing,              nVarsSEG)
 meta_struct(ixStruct%NTOPO  ) = struct_info('NTOPO',   varType_integer, integerMissing,              nVarsNTOPO)

 ! ---------- define variable dimensions -------------------------------------------------------------------------------------------------------

 ! structure index                         name     dimension id   dimension length
 meta_dims  (ixDims%hru      ) = dim_info('hru',   integerMissing, integerMissing)  ! hru vector
 meta_dims  (ixDims%seg      ) = dim_info('seg',   integerMissing, integerMissing)  ! stream segment vector
 meta_dims  (ixDims%upHRU    ) = dim_info('upHRU', integerMissing, integerMissing)  ! upstream HRUs
 meta_dims  (ixDims%upSeg    ) = dim_info('upSeg', integerMissing, integerMissing)  ! immediate upstream segments
 meta_dims  (ixDims%upAll    ) = dim_info('upAll', integerMissing, integerMissing)  ! all upstream segments

 ! ---------- populate metadata structures -----------------------------------------------------------------------------------------------------

 ! HRU                                          varName         varDesc                                                varUnit, varType, varFile
 meta_HRU    (ixHRU%area)              = var_info('area'           ,'basin area'                                         ,'m2'    ,ixDims%hru   , .true.)

 ! HRU2SEG                                      varName         varDesc                                                varUnit, varType, varFile
 meta_HRU2SEG(ixHRU2SEG%HRUid        ) = var_info('HRUid'          ,'unique HRU id'                                      ,'-'     ,ixDims%hru   , .true.)
 meta_HRU2SEG(ixHRU2SEG%HRUindex     ) = var_info('HRUindex'       ,'HRU index'                                          ,'-'     ,ixDims%hru   , .false.)
 meta_HRU2SEG(ixHRU2SEG%hruSegId     ) = var_info('hruSegId'       ,'id of the stream segment below each HRU'            ,'-'     ,ixDims%hru   , .true.)
 meta_HRU2SEG(ixHRU2SEG%hruSegIndex  ) = var_info('hruSegIndex'    ,'index of the stream segment below each HRU'         ,'-'     ,ixDims%hru   , .false.)

 ! SEG                                           varName        varDesc                                                varUnit, varType, varFile
 meta_SEG    (ixSEG%length           ) = var_info('length'         , 'length of segment'                                 ,'m'     ,ixDims%seg   , .true.)
 meta_SEG    (ixSEG%slope            ) = var_info('slope'          , 'slope of segment'                                  ,'-'     ,ixDims%seg   , .true.)
 meta_SEG    (ixSEG%width            ) = var_info('width'          , 'width of segment'                                  ,'m'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%man_n            ) = var_info('man_n'          , 'Mannings n'                                        ,'weird' ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%hruArea          ) = var_info('hruArea'        , 'area of each contributing HRU'                     ,'m2'    ,ixDims%upHRU , .false.)
 meta_SEG    (ixSEG%weight           ) = var_info('weight'         , 'weight assigned to each HRU'                       ,'-'     ,ixDims%upHRU , .false.)
 meta_SEG    (ixSEG%timeDelayHist    ) = var_info('timeDelayHist'  , 'time delay histogram for each reach'               ,'s'     ,ixDims%upAll , .false.)
 meta_SEG    (ixSEG%upsLength        ) = var_info('upsLength'      , 'length of the vector of reaches above each reach'  ,'m'     ,ixDims%upAll , .false.)
 meta_SEG    (ixSEG%upsArea          ) = var_info('upsArea'        , 'area above the top of the reach -- 0 if headwater' ,'m2'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%basUnderLake     ) = var_info('basUnderLake'   , 'Area of basin under lake'                          ,'m2'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%rchUnderLake     ) = var_info('rchUnderLake'   , 'Length of reach under lake'                        ,'m'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%minFlow          ) = var_info('minFlow'        , 'minimum environmental flow'                        ,'m s-1' ,ixDims%seg   , .false.)

 ! NTOPO                                         varName        varDesc                                                varUnit, varType, varFile
 meta_NTOPO  (ixNTOPO%nHRU           ) = var_info('nHRU'           , 'number of HRUs contributing flow to each segment'   ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%hruContribIx   ) = var_info('hruContribIx'   , 'indices of HRUs contributing flow to each segment'  ,'-'    ,ixDims%upHRU , .false.)
 meta_NTOPO  (ixNTOPO%hruContribId   ) = var_info('hruContribId'   , 'ids of HRUs that contribute flow to each segment'   ,'-'    ,ixDims%upHRU , .false.)
 meta_NTOPO  (ixNTOPO%segId          ) = var_info('segId'          , 'unique id of each stream segment'                   ,'-'    ,ixDims%seg   , .true.)
 meta_NTOPO  (ixNTOPO%segIndex       ) = var_info('segIndex'       , 'index of each stream segment (1, 2, 3, ..., n)'     ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%downSegId      ) = var_info('downSegId'      , 'unique id of the next downstream segment'           ,'-'    ,ixDims%seg   , .true.)
 meta_NTOPO  (ixNTOPO%downSegIndex   ) = var_info('downSegIndex'   , 'index of downstream reach index'                    ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%upSegIds       ) = var_info('upSegIds'       , 'ids for the immediate upstream stream segments'     ,'-'    ,ixDims%upSeg , .false.)
 meta_NTOPO  (ixNTOPO%upSegIndices   ) = var_info('upSegIndices'   , 'indices for the immediate upstream stream segments' ,'-'    ,ixDims%upSeg , .false.)
 meta_NTOPO  (ixNTOPO%allUpSegIndices) = var_info('allUpSegIndices', 'indices of all upstream stream segments'            ,'-'    ,ixDims%upAll , .false.)
 meta_NTOPO  (ixNTOPO%rchOrder       ) = var_info('rchOrder'       , 'order that stream segments are processed'           ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%lakeId         ) = var_info('lakeId'         , 'unique id of each lake in the river network'        ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%lakeIndex      ) = var_info('lakeIndex'      , 'index of each lake in the river network'            ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%isLakeInlet    ) = var_info('isLakeInlet'    , 'flag to define if a lake inlet (1=true)'            ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%userTake       ) = var_info('userTake'       , 'flag to define if user takes water (1=true)'        ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%goodBasin      ) = var_info('goodBasin'      , 'flag to define a good basin (1=true)'               ,'-'    ,ixDims%seg   , .false.)

 end subroutine popMetadat

end module popMetadat_module
