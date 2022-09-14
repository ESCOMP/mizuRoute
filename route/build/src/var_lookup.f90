MODULE var_lookup

 ! defines named variables used to index array elements
 USE nrtype,     ONLY: i4b
 USE public_var, ONLY: integerMissing  ! missing value for integers
 implicit none
 private
 ! local variables
 integer(i4b),parameter            :: ixVal=1                      ! an example integer
 integer(i4b),parameter            :: iLength=storage_size(ixVal)  ! size of the example integer

 ! ***********************************************************************************************************
 ! ** define data structures
 ! ***********************************************************************************************************
 type, public  ::  iLook_struct
  integer(i4b)     :: HRU          = integerMissing   ! 1. HRU structure
  integer(i4b)     :: HRU2SEG      = integerMissing   ! 2. HRU-to-segment mapping structure
  integer(i4b)     :: SEG          = integerMissing   ! 3. stream segment structure
  integer(i4b)     :: NTOPO        = integerMissing   ! 4. network topology structure
  integer(i4b)     :: PFAF         = integerMissing   ! 5. Pfafstetter code structure
 endtype iLook_struct
 ! ***********************************************************************************************************
 ! ** define dimensions
 ! ***********************************************************************************************************
 type, public  ::  iLook_dims
  integer(i4b)     :: hru          = integerMissing   ! 1. hru vector
  integer(i4b)     :: seg          = integerMissing   ! 2. stream segment vector
  integer(i4b)     :: upHRU        = integerMissing   ! 3. upstream HRUs
  integer(i4b)     :: upSeg        = integerMissing   ! 4. immediate upstream segments
  integer(i4b)     :: upAll        = integerMissing   ! 5. all upstream segments
  integer(i4b)     :: uh           = integerMissing   ! 6. all segment unit hydrographs
  integer(i4b)     :: pfaf         = integerMissing   ! 7. max pfaf code length
 endtype iLook_dims
 ! For routing state variables
 type, public  ::  iLook_stateDims
  integer(i4b)     :: seg          = integerMissing   !  1. stream segment vector
  integer(i4b)     :: time         = integerMissing   !  2. time
  integer(i4b)     :: tbound       = integerMissing   !  3. 2 elelment time bound vector
  integer(i4b)     :: ens          = integerMissing   !  4. runoff ensemble
  integer(i4b)     :: wave         = integerMissing   !  5. waves in a channel
  integer(i4b)     :: mol_kw       = integerMissing   !  6. kw finite difference computational molecule
  integer(i4b)     :: mol_mc       = integerMissing   !  7. mc finite difference computational molecule
  integer(i4b)     :: mol_dw       = integerMissing   !  8. kw finite difference computational molecule
  integer(i4b)     :: tdh_irf      = integerMissing   !  9. irf routed future channel flow in a segment
  integer(i4b)     :: tdh          = integerMissing   ! 10. uh routed future overland flow
 endtype iLook_stateDims
 ! For river discharge variables
 type, public  ::  iLook_qDims
  integer(i4b)     :: time         = integerMissing   ! 1. time
  integer(i4b)     :: seg          = integerMissing   ! 2. stream segment vector
  integer(i4b)     :: hru          = integerMissing   ! 3. hru vector
  integer(i4b)     :: ens          = integerMissing   ! 4. runoff ensemble
 endtype iLook_qDims
 ! ***********************************************************************************************************
 ! ** define variables desired for each HRU
 ! ***********************************************************************************************************
 type, public  ::  iLook_HRU
  integer(i4b)     :: area         = integerMissing   ! 1. basin area (m2)
 endtype iLook_HRU
 ! ***********************************************************************************************************
 ! ** define variables for the hru2segment mapping
 ! ***********************************************************************************************************
 type, public  ::  iLook_HRU2SEG
  integer(i4b)     :: HRUid         = integerMissing  ! 1. unique HRU id
  integer(i4b)     :: HRUindex      = integerMissing  ! 2. HRU index
  integer(i4b)     :: hruSegId      = integerMissing  ! 3. id of the stream segment below each HRU
  integer(i4b)     :: hruSegIndex   = integerMissing  ! 4. index of the stream segment below each HRU
 endtype iLook_HRU2SEG
 ! ***********************************************************************************************************
 ! ** define variables desired for each HRU
 ! ***********************************************************************************************************
 type, public  ::  iLook_SEG
  ! reach properties
  integer(i4b)     :: length           = integerMissing  !  1. length of segment (m)
  integer(i4b)     :: slope            = integerMissing  !  2. slope of segment (-)
  integer(i4b)     :: width            = integerMissing  !  3. width of segment (m)
  integer(i4b)     :: man_n            = integerMissing  !  4. Manning's n (weird units)
  ! contributing HRUs
  integer(i4b)     :: hruArea          = integerMissing  !  5. contributing area for each HRU (m2)
  integer(i4b)     :: weight           = integerMissing  !  6. weight assigned to each HRU (-)
  ! unit hydrograph routing
  integer(i4b)     :: timeDelayHist    = integerMissing  !  7. time delay histogram for each reach (-)
  integer(i4b)     :: basArea          = integerMissing  !  8. area of the local HRUs contributing to each reach (m2)
  integer(i4b)     :: upsArea          = integerMissing  !  9. area above the top of the reach -- zero if headwater (m2)
  integer(i4b)     :: totalArea        = integerMissing  ! 10. basArea + upsArea -- area at the bottom of the reach (m2)
  ! lakes
  integer(i4b)     :: basUnderLake     = integerMissing  ! 11. Area of basin under lake (m2)
  integer(i4b)     :: rchUnderLake     = integerMissing  ! 12. Length of reach under lake (m)
  ! Doll 2003 parameter (Natural lake outflow)
  integer(i4b)     :: D03_MaxStorage   = integerMissing  ! 13. Lake maximum volume (m3)
  integer(i4b)     :: D03_Coefficient  = integerMissing  ! 14.
  integer(i4b)     :: D03_Power        = integerMissing  ! 15.
  ! HYPE parameter (reservoir outflow)
  integer(i4b)     :: HYP_E_emr        = integerMissing  ! 16.
  integer(i4b)     :: HYP_E_lim        = integerMissing  ! 17.
  integer(i4b)     :: HYP_E_min        = integerMissing  ! 18.
  integer(i4b)     :: HYP_E_zero       = integerMissing  ! 19.
  integer(i4b)     :: HYP_Qrate_emr    = integerMissing  ! 20.
  integer(i4b)     :: HYP_Erate_emr    = integerMissing  ! 21.
  integer(i4b)     :: HYP_Qrate_prim   = integerMissing  ! 22.
  integer(i4b)     :: HYP_Qrate_amp    = integerMissing  ! 23.
  integer(i4b)     :: HYP_Qrate_phs    = integerMissing  ! 24.
  integer(i4b)     :: HYP_prim_F       = integerMissing  ! 25.
  integer(i4b)     :: HYP_A_avg        = integerMissing  ! 26.
  ! Hanasaki 2006 parameter (reservoir outlfow)
  integer(i4b)     :: H06_Smax         = integerMissing  ! 27.
  integer(i4b)     :: H06_alpha        = integerMissing  ! 28.
  integer(i4b)     :: H06_envfact      = integerMissing  ! 29.
  integer(i4b)     :: H06_S_ini        = integerMissing  ! 30.
  integer(i4b)     :: H06_c1           = integerMissing  ! 31.
  integer(i4b)     :: H06_c2           = integerMissing  ! 32.
  integer(i4b)     :: H06_exponent     = integerMissing  ! 33.
  integer(i4b)     :: H06_denominator  = integerMissing  ! 34.
  integer(i4b)     :: H06_c_compare    = integerMissing  ! 35.
  integer(i4b)     :: H06_frac_Sdead   = integerMissing  ! 36.
  integer(i4b)     :: H06_E_rel_ini    = integerMissing  ! 37.
  integer(i4b)     :: H06_I_Jan        = integerMissing  ! 38.
  integer(i4b)     :: H06_I_Feb        = integerMissing  ! 39.
  integer(i4b)     :: H06_I_Mar        = integerMissing  ! 40.
  integer(i4b)     :: H06_I_Apr        = integerMissing  ! 41.
  integer(i4b)     :: H06_I_May        = integerMissing  ! 42.
  integer(i4b)     :: H06_I_Jun        = integerMissing  ! 43.
  integer(i4b)     :: H06_I_Jul        = integerMissing  ! 44.
  integer(i4b)     :: H06_I_Aug        = integerMissing  ! 45.
  integer(i4b)     :: H06_I_Sep        = integerMissing  ! 46.
  integer(i4b)     :: H06_I_Oct        = integerMissing  ! 47.
  integer(i4b)     :: H06_I_Nov        = integerMissing  ! 48.
  integer(i4b)     :: H06_I_Dec        = integerMissing  ! 49.
  integer(i4b)     :: H06_D_Jan        = integerMissing  ! 50.
  integer(i4b)     :: H06_D_Feb        = integerMissing  ! 51.
  integer(i4b)     :: H06_D_Mar        = integerMissing  ! 52.
  integer(i4b)     :: H06_D_Apr        = integerMissing  ! 53.
  integer(i4b)     :: H06_D_May        = integerMissing  ! 54.
  integer(i4b)     :: H06_D_Jun        = integerMissing  ! 55.
  integer(i4b)     :: H06_D_Jul        = integerMissing  ! 56.
  integer(i4b)     :: H06_D_Aug        = integerMissing  ! 57.
  integer(i4b)     :: H06_D_Sep        = integerMissing  ! 58.
  integer(i4b)     :: H06_D_Oct        = integerMissing  ! 59.
  integer(i4b)     :: H06_D_Nov        = integerMissing  ! 60.
  integer(i4b)     :: H06_D_Dec        = integerMissing  ! 61.
  integer(i4b)     :: H06_purpose      = integerMissing  ! 62.
  integer(i4b)     :: H06_I_mem_F      = integerMissing  ! 63.
  integer(i4b)     :: H06_D_mem_F      = integerMissing  ! 64.
  integer(i4b)     :: H06_I_mem_L      = integerMissing  ! 65.
  integer(i4b)     :: H06_D_mem_L      = integerMissing  ! 66.
  ! constraints
  integer(i4b)     :: minFlow       = integerMissing     ! 67. minimum environmental flow (m3/s)
 endtype iLook_SEG
 ! ***********************************************************************************************************
 ! ** define variables for the network topology (all are unitless)
 ! ***********************************************************************************************************
 type, public  ::  iLook_NTOPO
  ! upstream HRUs
  integer(i4b)     :: nHRU            = integerMissing  !  1. number of HRUs that contribute flow to each segment
  integer(i4b)     :: hruContribIx    = integerMissing  !  2. indices of HRUs that contribute flow to each segment
  integer(i4b)     :: hruContribId    = integerMissing  !  3. ids of the vector of HRUs that contribute flow to each segment
  ! individual segments
  integer(i4b)     :: segId           = integerMissing  !  4. unique id of each stream segment
  integer(i4b)     :: segIndex        = integerMissing  !  5. index of each stream segment (1, 2, 3, ..., n)
  ! downstream segments
  integer(i4b)     :: downSegId       = integerMissing  !  6. unique id of the next downstream segment
  integer(i4b)     :: downSegIndex    = integerMissing  !  7. index of downstream reach index
  ! upstream segments
  integer(i4b)     :: upSegIds        = integerMissing  !  8. ids for the vector of immediate upstream stream segments
  integer(i4b)     :: upSegIndices    = integerMissing  !  9. indices for the vector of immediate upstream stream segments
  integer(i4b)     :: allUpSegIndices = integerMissing  ! 10. indices of all upstream stream segments
  ! processing sequence
  integer(i4b)     :: rchOrder        = integerMissing  ! 11. order that stream segments are processed
  ! stream order
  integer(i4b)     :: streamOrder     = integerMissing  ! 12. Shreve Stream Order
  ! lakes
  integer(i4b)     :: lakeId          = integerMissing  ! 13. unique id of each lake in the river network
  integer(i4b)     :: lakeIndex       = integerMissing  ! 14. index of each lake in the river network
  integer(i4b)     :: isLakeInlet     = integerMissing  ! 15. flag to define if reach is a lake inlet (1=inlet, 0 otherwise)
  integer(i4b)     :: islake          = integerMissing  ! 16. flag to define a lake (1=lake, 0=reach)
  integer(i4b)     :: LakeTargVol     = integerMissing  ! 17. flag to define if a lake should follow target volume (1=follow, 0=parameteric)
  integer(i4b)     :: LakeModelType   = integerMissing  ! 18. identifies the lake model (1=Doll, 2=Hanasaki, else=non-parameteric)
  ! irrigation
  integer(i4b)     :: userTake        = integerMissing  ! 19. flag to define if user takes water from reach (1=extract, 0 otherwise)
  integer(i4b)     :: destSegId       = integerMissing  ! 20. id of destination stream segment
  integer(i4b)     :: destSegIndex    = integerMissing  ! 21. index of destination stream segment
  ! testing
  integer(i4b)     :: goodBasin       = integerMissing  ! 22. flag to define a good basin (1=good, 0=bad)
 endtype iLook_NTOPO
 ! ***********************************************************************************************************
 ! ** define variables for Pfafstetter code  (temporary:   separated from NTOPO because type of code is character
 ! ***********************************************************************************************************
 type, public  ::  iLook_PFAF
  ! pfafstetter code
  integer(i4b)     :: code            = integerMissing  ! 1. pfafstetter code
 endtype iLook_PFAF
 ! ***********************************************************************************************************
 ! ** define variables for segment fluxes/states variables
 ! ***********************************************************************************************************
 ! Reach fluxes
 type, public  ::  iLook_RFLX
  integer(i4b)     :: instRunoff        = integerMissing  ! 1. instantaneous runoff at reach outlet (m3/s)
  integer(i4b)     :: dlayRunoff        = integerMissing  ! 2. delayed runoff at reach outlet (m3/s)
  integer(i4b)     :: sumUpstreamRunoff = integerMissing  ! 3. sum of upstream runoff at reach outlet (m3/s)
  integer(i4b)     :: KWTroutedRunoff   = integerMissing  ! 4. Lagrangian KWT routed runoff at reach outlet (m3/s)
  integer(i4b)     :: KWroutedRunoff    = integerMissing  ! 5. KW routed runoff at reach outlet(m3/s)
  integer(i4b)     :: MCroutedRunoff    = integerMissing  ! 6. muskingum-cunge routed runoff at reach outlet (m3/s)
  integer(i4b)     :: DWroutedRunoff    = integerMissing  ! 7. diffusive wave routed runoff at reach outlet (m3/s)
  integer(i4b)     :: IRFroutedRunoff   = integerMissing  ! 8. IRF routed runoff at reach outlet (m3/s)
  integer(i4b)     :: volume            = integerMissing  ! 9. water volume in reach (m3)
 endtype iLook_RFLX
 ! HRU fluxes
 type, public  ::  iLook_HFLX
  integer(i4b)     :: basRunoff         = integerMissing  ! 1. basin runoff (m/s)
 endtype iLook_HFLX
 ! Reach inflow from basin
 type, public  ::  iLook_basinQ
  integer(i4b)     :: q              = integerMissing  ! 1. final discharge (m3/s)
 endtype iLook_basinQ
 ! Basin IRF state/fluxes
 type, public  ::  iLook_IRFbas
  integer(i4b)     :: qfuture        = integerMissing  ! 1. future routed flow (m3/s)
 endtype iLook_IRFbas
 !IRF state/fluxes
 type, public  ::  iLook_IRF
  integer(i4b)     :: qfuture        = integerMissing  ! 1. future routed flow (m3/s)
  integer(i4b)     :: vol            = integerMissing  ! 2. reach volume (m3)
 endtype iLook_IRF
 ! KWT state/fluxes
 type, public  ::  iLook_KWT
  integer(i4b)     :: tentry         = integerMissing  ! 1. wave entry time at a segment (s)
  integer(i4b)     :: texit          = integerMissing  ! 2. wave exit time at a segment (s)
  integer(i4b)     :: qwave          = integerMissing  ! 3. wave flow (m3/s)
  integer(i4b)     :: qwave_mod      = integerMissing  ! 4. wave flow after merged (m3/s)
  integer(i4b)     :: routed         = integerMissing  ! 5. Routed out of a segment or not (-)
 endtype iLook_KWT
 ! KW state/fluxes
 type, public  ::  iLook_KW
  integer(i4b)     :: qsub           = integerMissing  ! 1. discharge (m3/s)
  integer(i4b)     :: vol            = integerMissing  ! 2. reach volume (m3)
 endtype iLook_KW
 ! DW state/fluxes
 type, public  ::  iLook_DW
  integer(i4b)     :: qsub           = integerMissing  ! 1. discharge (m3/s)
  integer(i4b)     :: vol            = integerMissing  ! 2. reach volume (m3)
 endtype iLook_DW
 ! MC state/fluxes
 type, public  ::  iLook_MC
  integer(i4b)     :: qsub           = integerMissing  ! 1. discharge (m3/s)
  integer(i4b)     :: vol            = integerMissing  ! 2. reach volume (m3)
 endtype iLook_MC
 ! ***********************************************************************************************************
 ! ** define data vectors
 ! ***********************************************************************************************************
 type(iLook_struct)   ,public,parameter :: ixStruct    = iLook_struct   ( 1, 2, 3, 4, 5)
 type(iLook_dims)     ,public,parameter :: ixDims      = iLook_dims     ( 1, 2, 3, 4, 5, 6, 7)
 type(iLook_stateDims),public,parameter :: ixStateDims = iLook_stateDims( 1, 2, 3, 4, 5, 6, 7, 8, 9,10)
 type(iLook_qDims)    ,public,parameter :: ixqDims     = iLook_qDims    ( 1, 2, 3, 4)
 type(iLook_HRU)      ,public,parameter :: ixHRU       = iLook_HRU      ( 1)
 type(iLook_HRU2SEG)  ,public,parameter :: ixHRU2SEG   = iLook_HRU2SEG  ( 1, 2, 3, 4)
 type(iLook_SEG)      ,public,parameter :: ixSEG       = iLook_SEG      ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10, &
                                                                         11,12,13,14,15,16,17,18,19,20, &
                                                                         21,22,23,24,25,26,27,28,29,30, &
                                                                         31,32,33,34,35,36,37,38,39,40, &
                                                                         41,42,43,44,45,46,47,48,49,50, &
                                                                         51,52,53,54,55,56,57,58,59,60, &
                                                                         61,62,63,64,65,66,67)
 type(iLook_NTOPO)    ,public,parameter :: ixNTOPO     = iLook_NTOPO    ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10, &
                                                                         11,12,13,14,15,16,17,18,19,20, &
                                                                         21,22)
 type(iLook_PFAF)     ,public,parameter :: ixPFAF      = iLook_PFAF     ( 1)
 type(iLook_RFLX)     ,public,parameter :: ixRFLX      = iLook_RFLX     ( 1, 2, 3, 4, 5, 6, 7, 8, 9)
 type(iLook_HFLX)     ,public,parameter :: ixHFLX      = iLook_HFLX     ( 1)
 type(iLook_basinQ)   ,public,parameter :: ixBasinQ    = iLook_basinQ   ( 1)
 type(iLook_IRFbas)   ,public,parameter :: ixIRFbas    = iLook_IRFbas   ( 1)
 type(iLook_IRF)      ,public,parameter :: ixIRF       = iLook_IRF      ( 1, 2)
 type(iLook_KWT)      ,public,parameter :: ixKWT       = iLook_KWT      ( 1, 2, 3, 4, 5)
 type(iLook_KW)       ,public,parameter :: ixKW        = iLook_KW       ( 1, 2)
 type(iLook_DW)       ,public,parameter :: ixDW        = iLook_DW       ( 1, 2)
 type(iLook_MC)       ,public,parameter :: ixMC        = iLook_MC       ( 1, 2)
 ! ***********************************************************************************************************
 ! ** define size of data vectors
 ! ***********************************************************************************************************
 integer(i4b),parameter,public    :: nStructures  = storage_size(ixStruct   )/iLength
 integer(i4b),parameter,public    :: nDimensions  = storage_size(ixDims     )/iLength
 integer(i4b),parameter,public    :: nStateDims   = storage_size(ixStateDims)/iLength
 integer(i4b),parameter,public    :: nqDims       = storage_size(ixQdims    )/iLength
 integer(i4b),parameter,public    :: nVarsHRU     = storage_size(ixHRU      )/iLength
 integer(i4b),parameter,public    :: nVarsHRU2SEG = storage_size(ixHRU2SEG  )/iLength
 integer(i4b),parameter,public    :: nVarsSEG     = storage_size(ixSEG      )/iLength
 integer(i4b),parameter,public    :: nVarsNTOPO   = storage_size(ixNTOPO    )/iLength
 integer(i4b),parameter,public    :: nVarsPFAF     = storage_size(ixPFAF    )/iLength
 integer(i4b),parameter,public    :: nVarsRFLX     = storage_size(ixRFLX    )/iLength
 integer(i4b),parameter,public    :: nVarsHFLX     = storage_size(ixHFLX    )/iLength
 integer(i4b),parameter,public    :: nVarsKWT      = storage_size(ixKWT     )/iLength
 integer(i4b),parameter,public    :: nVarsKW       = storage_size(ixKW      )/iLength
 integer(i4b),parameter,public    :: nVarsDW       = storage_size(ixDW      )/iLength
 integer(i4b),parameter,public    :: nVarsMC       = storage_size(ixMC      )/iLength
 integer(i4b),parameter,public    :: nVarsIRF      = storage_size(ixIRF     )/iLength
 integer(i4b),parameter,public    :: nVarsIRFbas   = storage_size(ixIRFbas  )/iLength
 integer(i4b),parameter,public    :: nVarsBasinQ   = storage_size(ixBasinQ  )/iLength
 ! ***********************************************************************************************************

END MODULE var_lookup
