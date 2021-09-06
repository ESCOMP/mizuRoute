module popMetadat_module

USE pio, ONLY: pio_real, pio_double, pio_int

! common variables
USE public_var, ONLY: MAXQPAR
USE public_var, ONLY: integerMissing
USE public_var, ONLY: varType_integer   ! named variable for an integer
USE public_var, ONLY: varType_double    ! named variable for a double precision
USE public_var, ONLY: varType_character ! named variable for a double precision

! data types
USE nrtype
USE dataTypes,  ONLY: struct_info  ! data type for metadata structure
USE dataTypes,  ONLY: var_info     ! data type for metadata structure
USE dataTypes,  ONLY: dim_info     ! data type for metadata structure

! metadata on data structures
USE globalData, ONLY: meta_struct    ! structure information
USE globalData, ONLY: meta_dims      ! dimensions
USE globalData, ONLY: meta_stateDims ! dimensions for routing states output
USE globalData, ONLY: meta_qDims     ! dimensions for river discharge output
USE globalData, ONLY: meta_HRU       ! HRU properties
USE globalData, ONLY: meta_HRU2SEG   ! HRU-to-segment mapping
USE globalData, ONLY: meta_SEG       ! stream segment properties
USE globalData, ONLY: meta_NTOPO     ! network topology
USE globalData, ONLY: meta_PFAF      ! pfafstetter code

USE globalData, ONLY: meta_rflx      ! reach flux variables
USE globalData, ONLY: meta_basinQ    ! reach inflow from basin
USE globalData, ONLY: meta_irf_bas   ! within-basin irf routing fluxes and states
USE globalData, ONLY: meta_irf       ! irf routing fluxes and states in a segment
USE globalData, ONLY: meta_kwt       ! kinematic wave routing fluxes and states in a segment
USE globalData, ONLY: meta_kwe       ! kinematic wave routing fluxes and states in a segment

! indices of named variables
USE var_lookup, ONLY: ixStruct   , nStructures   ! index of variables for data structure
USE var_lookup, ONLY: ixDims     , nDimensions   ! index of variables for data structure
USE var_lookup, ONLY: ixStateDims, nStateDims    ! index of variables for data structure
USE var_lookup, ONLY: ixQdims    , nQdims        ! index of variables for data structure
USE var_lookup, ONLY: ixHRU      , nVarsHRU      ! index of variables for data structure
USE var_lookup, ONLY: ixHRU2SEG  , nVarsHRU2SEG  ! index of variables for data structure
USE var_lookup, ONLY: ixSEG      , nVarsSEG      ! index of variables for data structure
USE var_lookup, ONLY: ixNTOPO    , nVarsNTOPO    ! index of variables for data structure
USE var_lookup, ONLY: ixPFAF     , nVarsPFAF     ! index of variables for data structure

USE var_lookup, ONLY: ixRFLX                     ! index of variables for data structure
USE var_lookup, ONLY: ixKWT                      ! index of variables for data structure
USE var_lookup, ONLY: ixKWE                      ! index of variables for data structure
USE var_lookup, ONLY: ixIRF                      ! index of variables for data structure
USE var_lookup, ONLY: ixIRFbas                   ! index of variables for data structure
USE var_lookup, ONLY: ixBasinQ                   ! index of variables for data structure

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

 err=0; message='popMetadat/'

 ! ---------- define data structures -----------------------------------------------------------------------------------------------------------

 ! structure index                            name      variable type    length of spatial dimension, number of variables
 meta_struct(ixStruct%HRU    ) = struct_info('HRU',     varType_double,    integerMissing,              nVarsHRU)
 meta_struct(ixStruct%HRU2SEG) = struct_info('HRU2SEG', varType_integer,   integerMissing,              nVarsHRU2SEG)
 meta_struct(ixStruct%SEG    ) = struct_info('SEG',     varType_double,    integerMissing,              nVarsSEG)
 meta_struct(ixStruct%NTOPO  ) = struct_info('NTOPO',   varType_integer,   integerMissing,              nVarsNTOPO)
 meta_struct(ixStruct%PFAF  )  = struct_info('PFAF',    varType_character, integerMissing,              nVarsPFAF)

 ! ---------- define variable dimensions -------------------------------------------------------------------------------------------------------

 ! structure index                         name     dimension id   dimension length
 meta_dims  (ixDims%hru      ) = dim_info('hru',   integerMissing, integerMissing)  ! hru vector
 meta_dims  (ixDims%seg      ) = dim_info('seg',   integerMissing, integerMissing)  ! stream segment vector
 meta_dims  (ixDims%upHRU    ) = dim_info('upHRU', integerMissing, integerMissing)  ! upstream HRUs
 meta_dims  (ixDims%upSeg    ) = dim_info('upSeg', integerMissing, integerMissing)  ! immediate upstream segments
 meta_dims  (ixDims%upAll    ) = dim_info('upAll', integerMissing, integerMissing)  ! all upstream segments
 meta_dims  (ixDims%uh       ) = dim_info('uh'   , integerMissing, integerMissing)  ! all unit hydrograph
 meta_dims  (ixDims%pfaf     ) = dim_info('pfaf' , integerMissing, integerMissing)  ! max pfafstetter code length

 meta_stateDims(ixStateDims%seg     ) = dim_info('seg',     integerMissing, integerMissing)  ! stream segment vector
 meta_stateDims(ixStateDims%time    ) = dim_info('time',    integerMissing, integerMissing)  ! time
 meta_stateDims(ixStateDims%tbound  ) = dim_info('tbound',  integerMissing, 2)               ! time bound (alway 2 - start and end)
 meta_stateDims(ixStateDims%ens     ) = dim_info('ens',     integerMissing, integerMissing)  ! runoff ensemble
 meta_stateDims(ixStateDims%wave    ) = dim_info('wave',    integerMissing, MAXQPAR)         ! reach waves vector (max. number is defined as MAXQPAR)
 meta_stateDims(ixStateDims%fdmesh  ) = dim_info('fdmesh',  integerMissing, 4)               ! KWE finite difference mesh points (1->[t0,x0], 2->[t0,x0], 3->[t0,x0],4->[t1,x1]
 meta_stateDims(ixStateDims%tdh_irf ) = dim_info('tdh_irf', integerMissing, integerMissing)  ! future time steps for irf routing
 meta_stateDims(ixStateDims%tdh     ) = dim_info('tdh',     integerMissing, integerMissing)  ! future time steps for bsasin irf routing

 meta_qDims(ixQdims%time    ) = dim_info('time',    integerMissing, integerMissing)   ! time
 meta_qDims(ixQdims%seg     ) = dim_info('seg',     integerMissing, integerMissing)   ! stream segment vector
 meta_qDims(ixQdims%hru     ) = dim_info('hru',     integerMissing, integerMissing)   ! hru vector
 meta_qDims(ixQdims%ens     ) = dim_info('ens',     integerMissing, integerMissing)   ! ensemble
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
 meta_SEG    (ixSEG%timeDelayHist    ) = var_info('timeDelayHist'  , 'time delay histogram for each reach'               ,'s'     ,ixDims%uh    , .false.)
 meta_SEG    (ixSEG%basArea          ) = var_info('basArea'        , 'total area of the contributing HRUs'               ,'m2'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%upsArea          ) = var_info('upsArea'        , 'area above the top of the reach -- 0 if headwater' ,'m2'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%totalArea        ) = var_info('totalArea'      , 'area above the bottom of the reach -- bas + ups'   ,'m2'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%basUnderLake     ) = var_info('basUnderLake'   , 'Area of basin under lake'                          ,'m2'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%rchUnderLake     ) = var_info('rchUnderLake'   , 'Length of reach under lake'                        ,'m'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%minFlow          ) = var_info('minFlow'        , 'minimum environmental flow'                        ,'m s-1' ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%D03_MaxStorage   ) = var_info('D03_MaxStorage' , 'Doll 2003; maximume storage Doll 2003'             ,'m3'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%D03_Coefficient  ) = var_info('D03_Coefficient', 'Doll 2003; coefficient Doll 2003'                  ,'day-1' ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%D03_Power        ) = var_info('D03_Power'      , 'Doll 2003; power Doll 2003'                        ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_Smax         ) = var_info('H06_Smax'       , 'Hanasaki 2006; maximume reservoir storage'                                               ,'m3'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_alpha        ) = var_info('H06_alpha'      , 'Hanasaki 2006; fraction of active storage compared to total storage'                     ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_envfact      ) = var_info('H06_envfact'    , 'Hanasaki 2006; fraction of inflow that can be used to meet demand'                       ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_S_ini        ) = var_info('H06_S_ini'      , 'Hanasaki 2006; initial storage used for initial estimation of release coefficient'       ,'m3'    ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_c1           ) = var_info('H06_c1'         , 'Hanasaki 2006; coefficient 1 for target release for irrigation reseroir'                 ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_c2           ) = var_info('H06_c2'         , 'Hanasaki 2006; coefficient 2 for target release for irrigation reseroir'                 ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_exponent     ) = var_info('H06_exponent'   , 'Hanasaki 2006; Exponenet of actual release for "within-a-year" reservoir'                ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_denominator  ) = var_info('H06_denominator', 'Hanasaki 2006; Denominator of actual release for "within-a-year" reservoir'              ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_c_compare    ) = var_info('H06_c_compare'  , 'Hanasaki 2006; Criterion for distinguish of "within-a-year" or "multi-year" reservoir'   ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_frac_Sdead   ) = var_info('H06_frac_Sdead' , 'Hanasaki 2006; Fraction of dead storage to maximume storage'                             ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_E_rel_ini    ) = var_info('H06_E_rel_ini'  , 'Hanasaki 2006; Initial release coefficient'                                              ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Jan        ) = var_info('H06_I_Jan'      , 'Hanasaki 2006; Average January   inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Feb        ) = var_info('H06_I_Feb'      , 'Hanasaki 2006; Average Februrary inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Mar        ) = var_info('H06_I_Mar'      , 'Hanasaki 2006; Average March     inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Apr        ) = var_info('H06_I_Apr'      , 'Hanasaki 2006; Average April     inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_May        ) = var_info('H06_I_May'      , 'Hanasaki 2006; Average May       inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Jun        ) = var_info('H06_I_Jun'      , 'Hanasaki 2006; Average June      inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Jul        ) = var_info('H06_I_Jul'      , 'Hanasaki 2006; Average July      inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Aug        ) = var_info('H06_I_Aug'      , 'Hanasaki 2006; Average August    inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Sep        ) = var_info('H06_I_Sep'      , 'Hanasaki 2006; Average September inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Oct        ) = var_info('H06_I_Oct'      , 'Hanasaki 2006; Average October   inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Nov        ) = var_info('H06_I_Nov'      , 'Hanasaki 2006; Average November  inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_Dec        ) = var_info('H06_I_Dec'      , 'Hanasaki 2006; Average December  inflow'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Jan        ) = var_info('H06_D_Jan'      , 'Hanasaki 2006; Average January   demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Feb        ) = var_info('H06_D_Feb'      , 'Hanasaki 2006; Average Februrary demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Mar        ) = var_info('H06_D_Mar'      , 'Hanasaki 2006; Average March     demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Apr        ) = var_info('H06_D_Apr'      , 'Hanasaki 2006; Average April     demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_May        ) = var_info('H06_D_May'      , 'Hanasaki 2006; Average April     demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Jun        ) = var_info('H06_D_Jun'      , 'Hanasaki 2006; Average May       demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Jul        ) = var_info('H06_D_Jul'      , 'Hanasaki 2006; Average June      demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Aug        ) = var_info('H06_D_Aug'      , 'Hanasaki 2006; Average July      demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Sep        ) = var_info('H06_D_Sep'      , 'Hanasaki 2006; Average Agust     demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Oct        ) = var_info('H06_D_Oct'      , 'Hanasaki 2006; Average September demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Nov        ) = var_info('H06_D_Nov'      , 'Hanasaki 2006; Average November  demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_Dec        ) = var_info('H06_D_Dec'      , 'Hanasaki 2006; Average December  demand'                                                 ,'m3 s-1',ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_purpose      ) = var_info('H06intP1'       , 'Hanasaki 2006; reservoir purpose; (0= non-irrigation, 1=irrigation)'                     ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_mem_F      ) = var_info('H06intP2'       , 'Hanasaki 2006; Flag to transition to modelled inflow'                                    ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_mem_F      ) = var_info('H06intP3'       , 'Hanasaki 2006; Flag to transition to modelled/provided demand'                           ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_I_mem_L      ) = var_info('H06intP4'       , 'Hanasaki 2006; Memory length in years for inflow'                                        ,'-'     ,ixDims%seg   , .false.)
 meta_SEG    (ixSEG%H06_D_mem_L      ) = var_info('H06intP5'       , 'Hanasaki 2006; Memory length in years for demand'                                        ,'-'     ,ixDims%seg   , .false.)

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
 meta_NTOPO  (ixNTOPO%streamOrder    ) = var_info('streamOrder'    , 'Shreve Stream Order'                                ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%lakeId         ) = var_info('lakeId'         , 'unique id of each lake in the river network'        ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%lakeIndex      ) = var_info('lakeIndex'      , 'index of each lake in the river network'            ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%isLakeInlet    ) = var_info('isLakeInlet'    , 'flag to define if a lake inlet (1=true)'            ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%userTake       ) = var_info('userTake'       , 'flag to define if user takes water (1=true)'        ,'-'    ,ixDims%seg   , .false.)
 meta_NTOPO  (ixNTOPO%goodBasin      ) = var_info('goodBasin'      , 'flag to define a good basin (1=true)'               ,'-'    ,ixDims%upSeg , .false.)
 meta_NTOPO  (ixNTOPO%islake         ) = var_info('islake'         , 'flag to define if the object is lake (1=true)'      ,'-'    ,ixDims%Seg   , .false.)
 meta_NTOPO  (ixNTOPO%LakeTargVol    ) = var_info('LakeTargVol'    , 'flag to define if lake follow target Vol (1=true)'  ,'-'    ,ixDims%Seg   , .false.)
 meta_NTOPO  (ixNTOPO%LakeModelType  ) = var_info('LakeModelType'  , 'lake model type (1=Doll, 2=Hanasaki, etc=none-para)','-'    ,ixDims%Seg   , .false.)

 ! PFAF CODE                                     varName        varDesc                                                varUnit, varType, varFile
 meta_PFAF  (ixPFAF%code             ) = var_info('code'           , 'pfafstetter code'                                   ,'-'    ,ixDims%seg   , .false.)

 ! ---------- populate segment fluxes/states metadata structures -----------------------------------------------------------------------------------------------------
! Reach Flux                                   varName              varDesc                                 unit,   varType,  varDim,                     writeOut
 call meta_rflx(ixRFLX%basRunoff        )%init('basRunoff'        , 'basin runoff'                        , 'm/s' , pio_real, [ixQdims%hru,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%instRunoff       )%init('instRunoff'       , 'instantaneous runoff in each reach'  , 'm3/s', pio_real, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%dlayRunoff       )%init('dlayRunoff'       , 'delayed runoff in each reach'        , 'm3/s', pio_real, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%sumUpstreamRunoff)%init('sumUpstreamRunoff', 'sum of upstream runoff in each reach', 'm3/s', pio_real, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%KWTroutedRunoff  )%init('KWTroutedRunoff'  , 'KWT routed runoff in each reach'     , 'm3/s', pio_real, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%KWEroutedRunoff  )%init('KWEroutedRunoff'  , 'KWE routed runoff in each reach'     , 'm3/s', pio_real, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%IRFroutedRunoff  )%init('IRFroutedRunoff'  , 'IRF routed runoff in each reach'     , 'm3/s', pio_real, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%IRFlakeVol       )%init('IRFlakeVol'       , 'lake and stream volume for IRF'      , 'm3'  , pio_real, [ixQdims%seg,ixQdims%time], .true.)

 ! Lagrangian Kinematic Wave         varName      varDesc                                           unit,   varType,    varDim,                                             writeOut
 call meta_kwt(ixKWT%tentry   )%init('tentry'   , 'time when a wave enters a segment'             , 's'   , pio_double, [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)
 call meta_kwt(ixKWT%texit    )%init('texit'    , 'time when a wave is expected to exit a segment', 's'   , pio_double, [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)
 call meta_kwt(ixKWT%qwave    )%init('qwave'    , 'flow of a wave'                                , 'm2/s', pio_double, [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)
 call meta_kwt(ixKWT%qwave_mod)%init('qwave_mod', 'modified flow of a wave'                       , 'm2/s', pio_double, [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)
 call meta_kwt(ixKWT%routed   )%init('routed'   , 'routing flag'                                  , '-'   , pio_int,    [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)

 ! Eulerian Kinematic Wave       varName      varDesc                       unit,   varType,    varDim,                                               writeOut
 call meta_kwe(ixKWE%a   )%init('flow_area', 'flow area'                  , 'm2'  , pio_double, [ixStateDims%seg,ixStateDims%fdmesh,ixStateDims%ens], .true.)
 call meta_kwe(ixKWE%q   )%init('kwe_q'    , 'Kinematic wave routed flow' , 'm2/s', pio_double, [ixStateDims%seg,ixStateDims%fdmesh,ixStateDims%ens], .true.)

 ! Impulse Response Function       varName         varDesc              unit,   varType,    varDim,                                                   writeOut
 call meta_irf(ixIRF%qfuture)%init('irf_qfuture', 'future flow series',   'm3/s' ,pio_double, [ixStateDims%seg,ixStateDims%tdh_irf,ixStateDims%ens] , .true.)
 call meta_irf(ixIRF%irfVol)%init ('irfVol'     , 'volume in reach/lake', 'm3'   ,pio_double, [ixStateDims%seg,ixStateDims%tbound, ixStateDims%ens] , .true.)

 ! Basin Impulse Response Function        varName    varDesc               unit,   varType,    varDim,                                           writeOut
 call meta_irf_bas(ixIRFbas%qfuture)%init('qfuture', 'future flow series', 'm3/s' ,pio_double, [ixStateDims%seg,ixStateDims%tdh,ixStateDims%ens], .true.)

! reach inflow from basin                 varName     varDesc              unit,  varType,     varDim,                                           writeOut
 call meta_basinQ(ixBasinQ%q      )%init('basin_q', 'basin routed flow' , 'm3/s' ,pio_double, [ixStateDims%seg,ixStateDims%ens]                , .true.)

 end subroutine popMetadat

end module popMetadat_module
