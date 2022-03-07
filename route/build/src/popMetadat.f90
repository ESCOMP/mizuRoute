module popMetadat_module

USE netcdf

! common variables
USE public_var, only : MAXQPAR
USE public_var, only : integerMissing

! data types
USE nrtype
USE dataTypes,  only : struct_info  ! data type for metadata structure
USE dataTypes,  only : var_info     ! data type for metadata structure
USE dataTypes,  only : dim_info     ! data type for metadata structure

! variable type
USE public_var, only : varType_integer   ! named variable for an integer
USE public_var, only : varType_double    ! named variable for a double precision
USE public_var, only : varType_character ! named variable for a double precision

! metadata on data structures
USE globalData, only : meta_struct    ! structure information
USE globalData, only : meta_dims      ! dimensions
USE globalData, only : meta_stateDims ! dimensions for routing states output
USE globalData, only : meta_qDims     ! dimensions for river discharge output
USE globalData, only : meta_HRU       ! HRU properties
USE globalData, only : meta_HRU2SEG   ! HRU-to-segment mapping
USE globalData, only : meta_SEG       ! stream segment properties
USE globalData, only : meta_NTOPO     ! network topology
USE globalData, only : meta_PFAF      ! pfafstetter code

USE globalData, only : meta_rflx      ! reach flux variables
USE globalData, only : meta_irf_bas   ! within-basin irf routing future flow
USE globalData, only : meta_basinQ    ! reach inflow from basin
USE globalData, only : meta_irf       ! irf routing restart states and fluxes in a segment
USE globalData, only : meta_kwt       ! lagrangian kinematic wave routing restart states and fluxes in a segment
USE globalData, only : meta_kw        ! kinematic wave routing restart fluxes and states in a segment
USE globalData, only : meta_dw        ! diffusive wave routing restart fluxes and states in a segment
USE globalData, only : meta_mc        ! muskingum-cunge routing restart fluxes and states in a segment

! indices of named variables
USE var_lookup, only : ixStruct   , nStructures   ! index of variables for data structure
USE var_lookup, only : ixDims     , nDimensions   ! index of variables for data structure
USE var_lookup, only : ixStateDims, nStateDims    ! index of variables for data structure
USE var_lookup, only : ixQdims    , nQdims        ! index of variables for data structure
USE var_lookup, only : ixHRU      , nVarsHRU      ! index of variables for data structure
USE var_lookup, only : ixHRU2SEG  , nVarsHRU2SEG  ! index of variables for data structure
USE var_lookup, only : ixSEG      , nVarsSEG      ! index of variables for data structure
USE var_lookup, only : ixNTOPO    , nVarsNTOPO    ! index of variables for data structure
USE var_lookup, only : ixPFAF     , nVarsPFAF     ! index of variables for data structure

USE var_lookup, only : ixRFLX                     ! index of variables for data structure
USE var_lookup, only : ixKWT                      ! index of variables for data structure
USE var_lookup, only : ixKW                       ! index of variables for data structure
USE var_lookup, only : ixDW                       ! index of variables for data structure
USE var_lookup, only : ixMC                       ! index of variables for data structure
USE var_lookup, only : ixIRF                      ! index of variables for data structure
USE var_lookup, only : ixIRFbas                   ! index of variables for data structure
USE var_lookup, only : ixBasinQ                   ! index of variables for data structure

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
 meta_stateDims(ixStateDims%mol_kw  ) = dim_info('mol_kw',  integerMissing, integerMissing)  ! kw finite difference mesh points
 meta_stateDims(ixStateDims%mol_mc  ) = dim_info('mol_mc',  integerMissing, integerMissing)  ! mc finite difference mesh points
 meta_stateDims(ixStateDims%mol_dw  ) = dim_info('mol_dw',  integerMissing, integerMissing)  ! dw finite difference mesh points
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
 meta_SEG    (ixSEG%Qtake            ) = var_info('Qtake'          , 'target abstraction(-)/injection(+)'                ,'m3 s-1',ixDims%seg   , .false.)
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

 ! PFAF CODE                                     varName        varDesc                                                varUnit, varType, varFile
 meta_PFAF  (ixPFAF%code             ) = var_info('code'           , 'pfafstetter code'                                   ,'-'    ,ixDims%seg   , .false.)

 ! ---------- populate segment fluxes/restart states metadata structures -----------------------------------------------------------------------------------------------------
 ! meta variable = (variable name, description, unit, numeric type, dimemsion, flag to write)
 ! Reach Flux
 call meta_rflx(ixRFLX%basRunoff        )%init('basRunoff'        , 'basin runoff'                                      , 'm/s' , nf90_float, [ixQdims%hru,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%instRunoff       )%init('instRunoff'       , 'instantaneous runoff in each reach'                , 'm3/s', nf90_float, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%dlayRunoff       )%init('dlayRunoff'       , 'delayed runoff in each reach'                      , 'm3/s', nf90_float, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%sumUpstreamRunoff)%init('sumUpstreamRunoff', 'sum of upstream runoff in each reach'              , 'm3/s', nf90_float, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%KWTroutedRunoff  )%init('KWTroutedRunoff'  , 'routed runoff in reach - lagrangian kinematic wave', 'm3/s', nf90_float, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%MCroutedRunoff   )%init('MCroutedRunoff'   , 'routed runoff in reach - muskingum-cunge'          , 'm3/s', nf90_float, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%DWroutedRunoff   )%init('DWroutedRunoff'   , 'routed runoff in reach - diffusive wave'           , 'm3/s', nf90_float, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%KWroutedRunoff   )%init('KWroutedRunoff'   , 'routed runoff in reach - lagrangian kinematic wave', 'm3/s', nf90_float, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%IRFroutedRunoff  )%init('IRFroutedRunoff'  , 'routed runoff in reach - Impulse Response Function', 'm3/s', nf90_float, [ixQdims%seg,ixQdims%time], .true.)
 call meta_rflx(ixRFLX%volume           )%init('volume'           , 'lake and stream volume'                            , 'm3'  , nf90_float, [ixQdims%seg,ixQdims%time], .false.)

 ! Lagrangian kinematic Wave restart state
 call meta_kwt(ixKWT%tentry   )%init('tentry'   , 'time when a wave enters a segment'             , 'sec'   , nf90_double, [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)
 call meta_kwt(ixKWT%texit    )%init('texit'    , 'time when a wave is expected to exit a segment', 'sec'   , nf90_double, [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)
 call meta_kwt(ixKWT%qwave    )%init('qwave'    , 'flow of a wave'                                , 'm2/sec', nf90_double, [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)
 call meta_kwt(ixKWT%qwave_mod)%init('qwave_mod', 'modified flow of a wave'                       , 'm2/sec', nf90_double, [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)
 call meta_kwt(ixKWT%routed   )%init('routed'   , 'routing flag'                                  , '-'     , nf90_int,    [ixStateDims%seg,ixStateDims%wave,ixStateDims%ens], .true.)

 ! Kinematic Wave
 call meta_kw(ixKW%qsub)%init('q_sub_kw', 'flow at computational moelcule', 'm3/s', nf90_double, [ixStateDims%seg,ixStateDims%mol_kw,ixStateDims%ens], .false.)

 ! Diffusive Wave
 call meta_dw(ixDW%qsub)%init('q_sub_dw', 'flow at computational moelcule', 'm3/s', nf90_double, [ixStateDims%seg,ixStateDims%mol_dw,ixStateDims%ens], .false.)

 ! Muskingum-cunge
 call meta_mc(ixMC%qsub)%init('q_sub_mc', 'flow at computational molecule', 'm3/s', nf90_double, [ixStateDims%seg,ixStateDims%mol_mc,ixStateDims%ens], .false.)

 ! Impulse Response Function
 call meta_irf(ixIRF%qfuture)%init('irf_qfuture', 'future flow series', 'm3/sec' ,nf90_double, [ixStateDims%seg,ixStateDims%tdh_irf,ixStateDims%ens] , .true.)
 call meta_irf(ixIRF%irfVol )%init('irf_volume' , 'IRF reach volume'  , 'm3'     ,nf90_double, [ixStateDims%seg,ixStateDims%ens]                     , .true.)

 ! Basin Impulse Response Function
 call meta_irf_bas(ixIRFbas%qfuture)%init('qfuture', 'future flow series', 'm3/sec' ,nf90_double, [ixStateDims%seg,ixStateDims%tdh,ixStateDims%ens], .true.)

 ! reach inflow from basin
 call meta_basinQ(ixBasinQ%q)%init('basin_q', 'inflow into reach from hru', 'm3/sec' ,nf90_double, [ixStateDims%seg,ixStateDims%ens], .true.)

 end subroutine popMetadat

end module popMetadat_module
