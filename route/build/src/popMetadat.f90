module popMetadat_module

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
USE globalData, ONLY: meta_irf_bas   ! within-basin irf routing fluxes and states
USE globalData, ONLY: meta_irf       ! irf routing fluxes and states in a segment
USE globalData, ONLY: meta_kwt       ! kinematic wave routing fluxes and states in a segment

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

USE var_lookup, ONLY: ixRFLX     , nVarsRFLX     ! index of variables for data structure
USE var_lookup, ONLY: ixKWT      , nVarsKWT      ! index of variables for data structure
USE var_lookup, ONLY: ixIRF      , nVarsIRF      ! index of variables for data structure
USE var_lookup, ONLY: ixIRFbas   , nVarsIRFbas   ! index of variables for data structure

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
 meta_NTOPO  (ixNTOPO%islake         ) = var_info('islake'         , 'flag to define if the object is lake (1=true)'      ,'-'    ,ixDims%Seg   , .true.) ! has added the lake flag

 ! PFAF CODE                                     varName        varDesc                                                varUnit, varType, varFile
 meta_PFAF  (ixPFAF%code             ) = var_info('code'           , 'pfafstetter code'                                   ,'-'    ,ixDims%seg   , .true.)

 ! ---------- populate segment fluxes/states metadata structures -----------------------------------------------------------------------------------------------------
 ! Reach Flux                                      varName             varDesc                                         unit,     varDim,              writeOut
 meta_rflx   (ixRFLX%basRunoff        ) = var_info('basRunoff'        , 'basin runoff'                                ,'m/s'    ,ixQdims%hru         ,.true.)
 meta_rflx   (ixRFLX%instRunoff       ) = var_info('instRunoff'       , 'instantaneous runoff in each reach'          ,'m3/s'   ,ixQdims%seg         ,.true.)
 meta_rflx   (ixRFLX%dlayRunoff       ) = var_info('dlayRunoff'       , 'delayed runoff in each reach'                ,'m3/s'   ,ixQdims%seg         ,.true.)
 meta_rflx   (ixRFLX%sumUpstreamRunoff) = var_info('sumUpstreamRunoff', 'sum of upstream runoff in each reach'        ,'m3/s'   ,ixQdims%seg         ,.true.)
 meta_rflx   (ixRFLX%KWTroutedRunoff  ) = var_info('KWTroutedRunoff'  , 'KWT routed runoff in each reach'             ,'m3/s'   ,ixQdims%seg         ,.true.)
 meta_rflx   (ixRFLX%IRFroutedRunoff  ) = var_info('IRFroutedRunoff'  , 'IRF routed runoff in each reach'             ,'m3/s'   ,ixQdims%seg         ,.true.)
 meta_rflx   (ixRFLX%IRFlakeVol  )      = var_info('IRFlakeVol'       , 'lake volume for IRF routing'                 ,'m3'     ,ixQdims%seg         ,.true.)

 ! Kinematic Wave                                 varName             varDesc                                          unit,     varDim,              writeOut
 meta_kwt    (ixKWT%tentry          ) = var_info('tentry'          , 'time when a wave enters a segment'              ,'sec'    ,ixStateDims%wave    ,.true.)
 meta_kwt    (ixKWT%texit           ) = var_info('texit'           , 'time when a wave is expected to exit a segment' ,'sec'    ,ixStateDims%wave    ,.true.)
 meta_kwt    (ixKWT%qwave           ) = var_info('qwave'           , 'flow of a wave'                                 ,'m2/sec' ,ixStateDims%wave    ,.true.)
 meta_kwt    (ixKWT%qwave_mod       ) = var_info('qwave_mod'       , 'modified flow of a wave'                        ,'m2/sec' ,ixStateDims%wave    ,.true.)
 meta_kwt    (ixKWT%routed          ) = var_info('routed'          , 'routing flag'                                   ,'-'      ,ixStateDims%wave    ,.true.)
 meta_kwt    (ixKWT%q               ) = var_info('kwt_q'           , 'Kinematic wave routed flow'                     ,'m3/sec' ,ixStateDims%time    ,.true.)

 ! Impulse Response Function                     varName             varDesc                                           unit,     varDim,              writeOut
 meta_irf    (ixIRF%qfuture         ) = var_info('irf_qfuture'     , 'future flow series'                             ,'m3/sec' ,ixStateDims%tdh_irf ,.true.)
 meta_irf    (ixIRF%q               ) = var_info('irf_q'           , 'irf routed flow'                                ,'m3/sec' ,ixStateDims%time    ,.true.)

 ! Basin Impulse Response Function               varName             varDesc                                           unit,     varDim,              writeOut
 meta_irf_bas(ixIRFbas%qfuture      ) = var_info('qfuture'         , 'future flow series'                             ,'m3/sec' ,ixStateDims%tdh     ,.true.)
 meta_irf_bas(ixIRFbas%q            ) = var_info('basin_q'         , 'basin routed flow'                              ,'m3/sec' ,ixStateDims%time    ,.true.)

 end subroutine popMetadat

end module popMetadat_module
