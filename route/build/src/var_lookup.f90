MODULE var_lookup

 ! defines named variables used to index array elements
 USE nrtype,     ONLY: i4b
 USE public_var, ONLY: integerMissing  ! missing value for integers
 implicit none
 private
 !
 INTRINSIC :: storage_size
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
  integer(i4b)     :: hru          = integerMissing   !  2. catchment hru vector
  integer(i4b)     :: time         = integerMissing   !  3. time
  integer(i4b)     :: tbound       = integerMissing   !  4. 2 elelment time bound vector
  integer(i4b)     :: ens          = integerMissing   !  5. runoff ensemble
  integer(i4b)     :: wave         = integerMissing   !  6. waves in a channel
  integer(i4b)     :: mol_kw       = integerMissing   !  7. kw finite difference computational molecule
  integer(i4b)     :: mol_mc       = integerMissing   !  8. mc finite difference computational molecule
  integer(i4b)     :: mol_dw       = integerMissing   !  9. kw finite difference computational molecule
  integer(i4b)     :: tdh_irf      = integerMissing   ! 10. irf routed future channel flow in a segment
  integer(i4b)     :: tdh          = integerMissing   ! 11. uh routed future overland flow
  integer(i4b)     :: nchars       = integerMissing   ! 12. number of characters
  integer(i4b)     :: hist_fil     = integerMissing   ! 13. history filenames
 endtype iLook_stateDims
 ! For river discharge variables
 type, public  ::  iLook_qDims
  integer(i4b)     :: time         = integerMissing   ! 1. time stamp
  integer(i4b)     :: tbound       = integerMissing   ! 2. time bound
  integer(i4b)     :: seg          = integerMissing   ! 3. stream segment vector
  integer(i4b)     :: hru          = integerMissing   ! 4. hru vector
  integer(i4b)     :: ens          = integerMissing   ! 5. runoff ensemble
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
  integer(i4b)     :: depth            = integerMissing  !  4. bankfull depth (m)
  integer(i4b)     :: sideSlope        = integerMissing  !  5. bankfull side slope, v:h=1:sideSlope (-)
  integer(i4b)     :: storage          = integerMissing  !  6. channel storage (m3)
  integer(i4b)     :: man_n            = integerMissing  !  7. Manning's n (weird units)
  integer(i4b)     :: floodplainSlope  = integerMissing  !  8. floodplain slope,v:h=1:floodplainSlope (-)
  ! contributing HRUs
  integer(i4b)     :: hruArea          = integerMissing  !  9. contributing area for each HRU (m2)
  integer(i4b)     :: weight           = integerMissing  ! 10. weight assigned to each HRU (-)
  ! unit hydrograph routing
  integer(i4b)     :: timeDelayHist    = integerMissing  ! 11. time delay histogram for each reach (-)
  integer(i4b)     :: basArea          = integerMissing  ! 12. area of the local HRUs contributing to each reach (m2)
  integer(i4b)     :: upsArea          = integerMissing  ! 13. area above the top of the reach -- zero if headwater (m2)
  integer(i4b)     :: totalArea        = integerMissing  ! 14. basArea + upsArea -- area at the bottom of the reach (m2)
  ! lakes
  integer(i4b)     :: basUnderLake     = integerMissing  ! 15. Area of basin under lake (m2)
  integer(i4b)     :: rchUnderLake     = integerMissing  ! 16. Length of reach under lake (m)
  ! Doll 2003 parameter (Natural lake outflow)
  integer(i4b)     :: D03_MaxStorage   = integerMissing  ! 17. Lake maximum volume (m3)
  integer(i4b)     :: D03_Coefficient  = integerMissing  ! 18.
  integer(i4b)     :: D03_Power        = integerMissing  ! 19.
  integer(i4b)     :: D03_S0           = integerMissing  ! 20.
  ! HYPE parameter (reservoir outflow)
  integer(i4b)     :: HYP_E_emr        = integerMissing  ! 21.
  integer(i4b)     :: HYP_E_lim        = integerMissing  ! 22.
  integer(i4b)     :: HYP_E_min        = integerMissing  ! 23.
  integer(i4b)     :: HYP_E_zero       = integerMissing  ! 24.
  integer(i4b)     :: HYP_Qrate_emr    = integerMissing  ! 25.
  integer(i4b)     :: HYP_Erate_emr    = integerMissing  ! 26.
  integer(i4b)     :: HYP_Qrate_prim   = integerMissing  ! 27.
  integer(i4b)     :: HYP_Qrate_amp    = integerMissing  ! 28.
  integer(i4b)     :: HYP_Qrate_phs    = integerMissing  ! 29.
  integer(i4b)     :: HYP_prim_F       = integerMissing  ! 30.
  integer(i4b)     :: HYP_A_avg        = integerMissing  ! 31.
  integer(i4b)     :: HYP_Qsim_mode    = integerMissing  ! 32.
  ! Hanasaki 2006 parameter (reservoir outlfow)
  integer(i4b)     :: H06_Smax         = integerMissing  ! 33.
  integer(i4b)     :: H06_alpha        = integerMissing  ! 34.
  integer(i4b)     :: H06_envfact      = integerMissing  ! 35.
  integer(i4b)     :: H06_S_ini        = integerMissing  ! 36.
  integer(i4b)     :: H06_c1           = integerMissing  ! 37.
  integer(i4b)     :: H06_c2           = integerMissing  ! 38.
  integer(i4b)     :: H06_exponent     = integerMissing  ! 39.
  integer(i4b)     :: H06_denominator  = integerMissing  ! 40.
  integer(i4b)     :: H06_c_compare    = integerMissing  ! 41.
  integer(i4b)     :: H06_frac_Sdead   = integerMissing  ! 42.
  integer(i4b)     :: H06_E_rel_ini    = integerMissing  ! 43.
  integer(i4b)     :: H06_I_Jan        = integerMissing  ! 44.
  integer(i4b)     :: H06_I_Feb        = integerMissing  ! 45.
  integer(i4b)     :: H06_I_Mar        = integerMissing  ! 46.
  integer(i4b)     :: H06_I_Apr        = integerMissing  ! 47.
  integer(i4b)     :: H06_I_May        = integerMissing  ! 48.
  integer(i4b)     :: H06_I_Jun        = integerMissing  ! 49.
  integer(i4b)     :: H06_I_Jul        = integerMissing  ! 50.
  integer(i4b)     :: H06_I_Aug        = integerMissing  ! 51.
  integer(i4b)     :: H06_I_Sep        = integerMissing  ! 52.
  integer(i4b)     :: H06_I_Oct        = integerMissing  ! 53.
  integer(i4b)     :: H06_I_Nov        = integerMissing  ! 54.
  integer(i4b)     :: H06_I_Dec        = integerMissing  ! 55.
  integer(i4b)     :: H06_D_Jan        = integerMissing  ! 56.
  integer(i4b)     :: H06_D_Feb        = integerMissing  ! 57.
  integer(i4b)     :: H06_D_Mar        = integerMissing  ! 58.
  integer(i4b)     :: H06_D_Apr        = integerMissing  ! 59.
  integer(i4b)     :: H06_D_May        = integerMissing  ! 60.
  integer(i4b)     :: H06_D_Jun        = integerMissing  ! 61.
  integer(i4b)     :: H06_D_Jul        = integerMissing  ! 62.
  integer(i4b)     :: H06_D_Aug        = integerMissing  ! 63.
  integer(i4b)     :: H06_D_Sep        = integerMissing  ! 64.
  integer(i4b)     :: H06_D_Oct        = integerMissing  ! 65.
  integer(i4b)     :: H06_D_Nov        = integerMissing  ! 66.
  integer(i4b)     :: H06_D_Dec        = integerMissing  ! 67.
  integer(i4b)     :: H06_purpose      = integerMissing  ! 68.
  integer(i4b)     :: H06_I_mem_F      = integerMissing  ! 69.
  integer(i4b)     :: H06_D_mem_F      = integerMissing  ! 70.
  integer(i4b)     :: H06_I_mem_L      = integerMissing  ! 71.
  integer(i4b)     :: H06_D_mem_L      = integerMissing  ! 72.
  ! constraints
  integer(i4b)     :: minFlow       = integerMissing     ! 73. minimum environmental flow (m3/s)
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
 ! ** define history output fluxes/states variables
 ! ***********************************************************************************************************
 ! Reach/lake fluxes
 type, public  ::  iLook_RFLX
  integer(i4b)     :: instRunoff        = integerMissing  !  1. instantaneous runoff at reach outlet (m3/s)
  integer(i4b)     :: dlayRunoff        = integerMissing  !  2. delayed runoff at reach outlet (m3/s)
  integer(i4b)     :: sumUpstreamRunoff = integerMissing  !  3. sum of upstream runoff at reach outlet (m3/s)
  integer(i4b)     :: IRFroutedRunoff   = integerMissing  !  4. IRF routed runoff at reach outlet (m3/s)
  integer(i4b)     :: KWTroutedRunoff   = integerMissing  !  5. Lagrangian KWT routed runoff at reach outlet (m3/s)
  integer(i4b)     :: KWroutedRunoff    = integerMissing  !  6. KW routed runoff at reach outlet(m3/s)
  integer(i4b)     :: MCroutedRunoff    = integerMissing  !  7. muskingum-cunge routed runoff at reach outlet (m3/s)
  integer(i4b)     :: DWroutedRunoff    = integerMissing  !  8. diffusive wave routed runoff at reach outlet (m3/s)
  integer(i4b)     :: IRFvolume         = integerMissing  !  9. IRF water volume in reach/lake (m3)
  integer(i4b)     :: KWTvolume         = integerMissing  ! 10. KWT water volume in reach/lake (m3)
  integer(i4b)     :: KWvolume          = integerMissing  ! 11. KW water volume in reach/lake (m3)
  integer(i4b)     :: MCvolume          = integerMissing  ! 12. MC water volume in reach/lake (m3)
  integer(i4b)     :: DWvolume          = integerMissing  ! 13. DW water volume in reach/lake (m3)
  integer(i4b)     :: IRFheight         = integerMissing  ! 14. KW water volume in floodplain (m3)
  integer(i4b)     :: KWTheight         = integerMissing  ! 15. KW water volume in floodplain (m3)
  integer(i4b)     :: KWheight          = integerMissing  ! 16. KW water height from bottom in reach/lake (m)
  integer(i4b)     :: MCheight          = integerMissing  ! 17. MC water height from bottom in reach/lake (m)
  integer(i4b)     :: DWheight          = integerMissing  ! 18. DW water height from bottom in reach/lake (m)
  integer(i4b)     :: IRFfloodVolume    = integerMissing  ! 19. KW water volume in floodplain (m3)
  integer(i4b)     :: KWTfloodVolume    = integerMissing  ! 20. KW water volume in floodplain (m3)
  integer(i4b)     :: KWfloodVolume     = integerMissing  ! 21. KW water volume in floodplain (m3)
  integer(i4b)     :: MCfloodVolume     = integerMissing  ! 22. MC water volume in floodplain (m3)
  integer(i4b)     :: DWfloodVolume     = integerMissing  ! 23. DW water volume in floodplain (m3)
  integer(i4b)     :: IRFinflow         = integerMissing  ! 24. IRF inflow from upstreams into reach/lake (m3/s)
  integer(i4b)     :: KWTinflow         = integerMissing  ! 25. KWT inflow flow upstreams into reach/lake (m3/s)
  integer(i4b)     :: KWinflow          = integerMissing  ! 26. KW inflow from upstreams into reach/lake (m3/s)
  integer(i4b)     :: MCinflow          = integerMissing  ! 27. MC inflow from upstreams into reach/lake (m3/s)
  integer(i4b)     :: DWinflow          = integerMissing  ! 28. DW inflow from upstreams into reach/lake (m3/s)
  integer(i4b)     :: localCC           = integerMissing  ! 29. concentration from local basin (g/m3)
  integer(i4b)     :: DWsoluteFlux      = integerMissing  ! 30. DW routed solute flux from upstreams into reach/lake (mg/s)
  integer(i4b)     :: DWsoluteMass      = integerMissing  ! 31. DW routed solute mass in reach/lake (mg)
 endtype iLook_RFLX
 ! HRU fluxes
 type, public  ::  iLook_HFLX
  integer(i4b)     :: basRunoff         = integerMissing  ! 1. basin runoff (m/s)
 endtype iLook_HFLX
 ! ***********************************************************************************************************
 ! ** define restart variables
 ! ***********************************************************************************************************
 ! Reach inflow from basin
 type, public  ::  iLook_basinQ
  integer(i4b)     :: q              = integerMissing  ! 1. final discharge (m3/s)
 endtype iLook_basinQ
 ! Reach constituent from basin
 type, public  ::  iLook_basinC
  integer(i4b)     :: c              = integerMissing  ! 1. constituent flux (g/m3/s)
 endtype iLook_basinC
 ! Basin IRF state/fluxes
 type, public  ::  iLook_IRFbas
  integer(i4b)     :: qfuture        = integerMissing  ! 1. future routed flow (m3/s)
 endtype iLook_IRFbas
 !IRF state/fluxes
 type, public  ::  iLook_IRF
  integer(i4b)     :: qfuture        = integerMissing  ! 1. future routed flow (m3/s)
  integer(i4b)     :: vol            = integerMissing  ! 2. reach volume (m3)
  integer(i4b)     :: qerror         = integerMissing  ! 3. discharge error (m3/s): optional
 endtype iLook_IRF
 ! KWT state/fluxes
 type, public  ::  iLook_KWT
  integer(i4b)     :: tentry         = integerMissing  ! 1. wave entry time at a segment (s)
  integer(i4b)     :: texit          = integerMissing  ! 2. wave exit time at a segment (s)
  integer(i4b)     :: qwave          = integerMissing  ! 3. wave flow (m3/s)
  integer(i4b)     :: qwave_mod      = integerMissing  ! 4. wave flow after merged (m3/s)
  integer(i4b)     :: routed         = integerMissing  ! 5. Routed out of a segment or not (-)
  integer(i4b)     :: vol            = integerMissing  ! 6. reach volume (m3)
 endtype iLook_KWT
 ! KW state/fluxes
 type, public  ::  iLook_KW
  integer(i4b)     :: qsub           = integerMissing  ! 1. discharge (m3/s)
  integer(i4b)     :: vol            = integerMissing  ! 2. reach volume (m3)
  integer(i4b)     :: qerror         = integerMissing  ! 3. discharge error (m3/s): optional
 endtype iLook_KW
 ! DW state/fluxes
 type, public  ::  iLook_DW
  integer(i4b)     :: qsub           = integerMissing  ! 1. discharge (m3/s)
  integer(i4b)     :: vol            = integerMissing  ! 2. reach volume (m3)
  integer(i4b)     :: qerror         = integerMissing  ! 3. discharge error (m3/s): optional
 endtype iLook_DW
 ! MC state/fluxes
 type, public  ::  iLook_MC
  integer(i4b)     :: qsub           = integerMissing  ! 1. discharge (m3/s)
  integer(i4b)     :: vol            = integerMissing  ! 2. reach volume (m3)
  integer(i4b)     :: qerror         = integerMissing  ! 3. discharge error (m3/s): optional
 endtype iLook_MC
 ! ***********************************************************************************************************
 ! ** define data vectors
 ! ***********************************************************************************************************
 type(iLook_struct)   ,public,parameter :: ixStruct    = iLook_struct   ( 1, 2, 3, 4, 5)
 type(iLook_dims)     ,public,parameter :: ixDims      = iLook_dims     ( 1, 2, 3, 4, 5, 6, 7)
 type(iLook_stateDims),public,parameter :: ixStateDims = iLook_stateDims( 1, 2, 3, 4, 5, 6, 7, 8, 9,10, &
                                                                         11, 12, 13)
 type(iLook_qDims)    ,public,parameter :: ixQdims     = iLook_qDims    ( 1, 2, 3, 4, 5)
 type(iLook_HRU)      ,public,parameter :: ixHRU       = iLook_HRU      ( 1)
 type(iLook_HRU2SEG)  ,public,parameter :: ixHRU2SEG   = iLook_HRU2SEG  ( 1, 2, 3, 4)
 type(iLook_SEG)      ,public,parameter :: ixSEG       = iLook_SEG      ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10, &
                                                                         11,12,13,14,15,16,17,18,19,20, &
                                                                         21,22,23,24,25,26,27,28,29,30, &
                                                                         31,32,33,34,35,36,37,38,39,40, &
                                                                         41,42,43,44,45,46,47,48,49,50, &
                                                                         51,52,53,54,55,56,57,58,59,60, &
                                                                         61,62,63,64,65,66,67,68,69,70, &
                                                                         71,72,73)
 type(iLook_NTOPO)    ,public,parameter :: ixNTOPO     = iLook_NTOPO    ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10, &
                                                                         11,12,13,14,15,16,17,18,19,20, &
                                                                         21,22)
 type(iLook_PFAF)     ,public,parameter :: ixPFAF      = iLook_PFAF     ( 1)
 type(iLook_RFLX)     ,public,parameter :: ixRFLX      = iLook_RFLX     ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10, &
                                                                         11,12,13,14,15,16,17,18,19,20, &
                                                                         21,22,23,24,25,26,27,28,29,30, &
                                                                         31)
 type(iLook_HFLX)     ,public,parameter :: ixHFLX      = iLook_HFLX     ( 1)
 type(iLook_basinQ)   ,public,parameter :: ixBasinQ    = iLook_basinQ   ( 1)
 type(iLook_basinC)   ,public,parameter :: ixBasinC    = iLook_basinC   ( 1)
 type(iLook_IRFbas)   ,public,parameter :: ixIRFbas    = iLook_IRFbas   ( 1)
 type(iLook_IRF)      ,public,parameter :: ixIRF       = iLook_IRF      ( 1, 2, 3)
 type(iLook_KWT)      ,public,parameter :: ixKWT       = iLook_KWT      ( 1, 2, 3, 4, 5, 6)
 type(iLook_KW)       ,public,parameter :: ixKW        = iLook_KW       ( 1, 2, 3)
 type(iLook_DW)       ,public,parameter :: ixDW        = iLook_DW       ( 1, 2, 3)
 type(iLook_MC)       ,public,parameter :: ixMC        = iLook_MC       ( 1, 2, 3)
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
 integer(i4b),parameter,public    :: nVarsBasinC   = storage_size(ixBasinC  )/iLength
 ! ***********************************************************************************************************

END MODULE var_lookup
