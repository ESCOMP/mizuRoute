MODULE read_streamSeg

! data types
USE nrtype,    ONLY: i4b,dp,lgt
USE nrtype,    ONLY: strLen               ! string length
USE dataTypes, ONLY: var_ilength          ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_dlength          ! double precision type: var(:)%dat
USE dataTypes, ONLY: var_clength          ! character type:        var(:)%dat
USE dataTypes, ONLY: var_info             ! metadata

! global data
USE public_var

! metadata on data structures
USE globalData, ONLY: meta_struct         ! structure information
USE globalData, ONLY: meta_HRU            ! HRU properties
USE globalData, ONLY: meta_HRU2SEG        ! HRU-to-segment mapping
USE globalData, ONLY: meta_SEG            ! stream segment properties
USE globalData, ONLY: meta_NTOPO          ! network topology
USE globalData, ONLY: meta_PFAF           ! network topology

! named variables
USE var_lookup, ONLY: ixStruct, nStructures  ! index of data structures
USE var_lookup, ONLY: ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup, ONLY: ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup, ONLY: ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup, ONLY: ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology
USE var_lookup, ONLY: ixPFAF,   nVarsPFAF    ! index of variables for the pfafstetter code

! netcdf modules
USE netcdf

! external utilities
USE nr_utility_module, ONLY: arth    ! Num. Recipies utilities
USE allocation,        ONLY: alloc_struct

implicit none

! privacy
private
public::getData
contains

! *********************************************************************
! new subroutine: get ancillary data for HRUs and stream segments
! *********************************************************************
subroutine getData(&
                  ! input
                  fname,        & ! input: file name
                  dname_nhru,   & ! input: dimension name of the HRUs
                  dname_sseg,   & ! input: dimension name of the stream segments
                  ! output: model control
                  nHRU_in,      & ! output: number of HRUs
                  nRch_in,      & ! output: number of stream segments
                  ! output: populate data structures
                  structHRU,    & ! ancillary data for HRUs
                  structSeg,    & ! ancillary data for stream segments
                  structHRU2seg,& ! ancillary data for mapping hru2basin
                  structNTOPO,  & ! ancillary data for network toopology
                  structPFAF,   & ! ancillary data for pfafstetter code
                  ! output: error control
                  ierr,message)   ! output: error control

  implicit none
  ! input variables
  character(*)      , intent(in)               :: fname             ! filename
  character(*)      , intent(in)               :: dname_nhru        ! dimension name for HRUs
  character(*)      , intent(in)               :: dname_sseg        ! dimension name for stream segments
  ! output: model control
  integer(i4b)      , intent(out)              :: nHRU_in           ! number of HRUs
  integer(i4b)      , intent(out)              :: nRch_in           ! number of stream segments
  ! output: data structures
  type(var_dlength) , intent(out), allocatable :: structHRU(:)      ! HRU properties
  type(var_dlength) , intent(out), allocatable :: structSeg(:)      ! stream segment properties
  type(var_ilength) , intent(out), allocatable :: structHRU2seg(:)  ! HRU-to-segment mapping
  type(var_ilength) , intent(out), allocatable :: structNTOPO(:)    ! network topology
  type(var_clength) , intent(out), allocatable :: structPFAF(:)     ! network topology
  ! output: error control
  integer(i4b)      , intent(out)              :: ierr              ! error code
  character(*)      , intent(out)              :: message           ! error message
  ! ==========================================================================================================
  ! local variables
  integer(i4b)                                 :: iStruct                 ! structure index
  integer(i4b)                                 :: iSpace                  ! spatial index
  integer(i4b)                                 :: iVar                    ! variable index
  integer(i4b)                                 :: ncid                    ! NetCDF file ID
  integer(i4b)                                 :: idimID_nHRU             ! dimension ID for HRUs
  integer(i4b)                                 :: idimID_sseg             ! dimension ID for stream segments
  integer(i4b)                                 :: iVarID                  ! variable ID
  integer(i4b)                                 :: jxStart                 ! Start index for a given reach
  integer(i4b)                                 :: jxCount                 ! Number of elements for a given reach
  integer(i4b),              allocatable       :: ixStart(:)              ! Start index for each reach
  integer(i4b),              allocatable       :: ixCount(:)              ! Number of elements in each reach
  integer(i4b),              allocatable       :: iTemp(:)                ! temporary integer vector
  real(dp),                  allocatable       :: dTemp(:)                ! temporary double precision vector
  character(Len=maxPfafLen), allocatable       :: cTemp(:)                ! temporary charactervector
  integer(i4b)                                 :: dimLength               ! dimension length
  logical(lgt)                                 :: isVarDesired            ! .true. if the variable is desired
  character(len=strLen)                        :: varName                 ! variable name
  character(len=strLen)                        :: cmessage                ! error message of downwind routine
  integer(i4b),              allocatable       :: islake_local(:)         ! local array to save islake flag
  integer(i4b),              allocatable       :: LakeTargVol_local(:)    ! local array to save LakeTargetVol flag
  integer(i4b),              allocatable       :: LakeModelType_local(:)  ! local array to save LakeModelType flag
  ! logical(lgt)                                 :: Doll_is_called          ! if Doll model is called in lake type varibale
  ! logical(lgt)                                 :: Hanasaki_is_called      ! if Hanasaki model is called in lake type varibale
  ! logical(lgt)                                 :: HYPE_is_called          ! if HYPE model is called in lake type varibale
  integer(i4b)                                 :: i                       ! counter
  !logical(lgt)                                 :: lake_model_conflict     ! if both parameteric model and non parameteric models are on for a lake
  integer(i4b)                                 :: number_lakes            ! number of lakes in network topology
  integer(i4b)                                 :: number_Endorheic        ! number of Endorheic lakes
  integer(i4b)                                 :: number_Doll             ! number of lakes with parameteric Doll 2003 formulation
  integer(i4b)                                 :: number_Hanasaki         ! number of lakes with parameteric Hanasaki 2006 formulation
  integer(i4b)                                 :: number_HYPE             ! number of lakes with parameteric Hanasaki 2006 formulation
  integer(i4b)                                 :: number_TargVol          ! number of lakes with target volume


  ierr=0; message='getData, read_segment/'

  ! ---------- initial reading of dimensions ------------------------------------------------------------------------

  ! open file for reading
  ierr = nf90_open(fname, nf90_nowrite, ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; file='//trim(fname); return; endif

  ! get the ID of the HRU dimension
  ierr = nf90_inq_dimid(ncid, dname_nhru, idimID_nHRU)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_nhru); return; endif

  ! get the length of the HRU dimension
  ierr = nf90_inquire_dimension(ncid, idimID_nHRU, len=nHRU_in)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get the ID of the stream segment dimension
  ierr = nf90_inq_dimid(ncid, dname_sseg, idimID_sseg)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_sseg); return; endif

  ! get the length of the stream segment dimension
  ierr = nf90_inquire_dimension(ncid, idimID_sseg, len=nRch_in)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! ---------- allocate space for higher-level structure components -------------------------------------------------
  call alloc_struct(&
                    nHRU_in,      & ! output: number of HRUs
                    nRch_in,      & ! output: number of stream segments
                    structHRU,    & ! inout: ancillary data for HRUs
                    structSeg,    & ! inout: ancillary data for stream segments
                    structHRU2seg,& ! inout: ancillary data for mapping hru2basin
                    structNTOPO,  & ! inout: ancillary data for network toopology
                    structPFAF,   & ! inout: ancillary data for pfafstetter code
                    ierr,cmessage)  ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initial allocation of the temporary vectors
  allocate(iTemp(nHRU_in), dTemp(nHRU_in), cTemp(nHRU_in), stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating temporary vectors'; return; endif

  ! -----------------------------------------------------------------------------------------------------------------
  ! ---------- read in data -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------------------------

  ! set flags if we want to read hdraulic geometry from file
  if(hydGeometryOption==readFromFile)then
    meta_SEG(ixSEG%width)%varFile = .true.
    meta_SEG(ixSEG%man_n)%varFile = .true.
  endif

  ! set flags if we simulate lake and need lake parameters; was set to false in pop_metadata.f90
  if(is_lake_sim)then

    meta_NTOPO(ixNTOPO%islake)%varFile          = .true.       ! if the object is lake should be provided
    meta_NTOPO(ixNTOPO%lakeModelType)%varFile   = .true.       ! if the object is lake, lake type should be provided

    ! read the lake type variable
    varName = trim(meta_NTOPO(ixNTOPO%lakeModelType)%varName)  ! get the varibale name of lakeModelType
    ierr = nf90_inq_varid(ncid, varName, ivarID)               ! get the id of the lakeModelType variable
    if(ierr/=0)then; ierr=20; message=trim(message)//'problem reading variable ID for lakeModelType'; return; endif
    allocate (LakeModelType_local(nRch_in), stat=ierr)         ! allocate the lakeModelType_local
    if(ierr/=0)then; ierr=20; message=trim(message)//'not able to allocate LakeModelType_local'; return; endif
    ierr = nf90_get_var(ncid, ivarID, LakeModelType_local)     ! read the variable
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; varname='//trim(varName); return; endif

    ! read the is_lake flag
    varName = trim(meta_NTOPO(ixNTOPO%islake)%varName)         ! get the varibale name of islake
    ierr = nf90_inq_varid(ncid, varName, ivarID)               ! get the id of the islake variable
    if(ierr/=0)then; ierr=20; message=trim(message)//'problem reading variable ID for islake'; return; endif
    allocate (islake_local(nRch_in), stat=ierr)                ! allocate the islake_local
    if(ierr/=0)then; ierr=20; message=trim(message)//'not able to allocate islake_local'; return; endif
    ierr = nf90_get_var(ncid, ivarID, islake_local)            ! read the variable
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; varname='//trim(varName); return; endif

    ! read lake target volume
    if (is_vol_wm) then
      meta_NTOPO(ixNTOPO%laketargvol)%varFile   = .true.         ! if target volume flag is on, varibale should be provided
      ! read the is_lake flag
      varName = trim(meta_NTOPO(ixNTOPO%laketargvol)%varName)    ! get the varibale name of target volume
      ierr = nf90_inq_varid(ncid, varName, ivarID)               ! get the id of the target volume variable
      if(ierr/=0)then; ierr=20; message=trim(message)//'problem reading variable ID for lake target volume'; return; endif
      allocate (LakeTargVol_local(nRch_in), stat=ierr)           ! allocate the LakeTargVol_local
      if(ierr/=0)then; ierr=20; message=trim(message)//'not able to allocate LakeTargVol_local'; return; endif
      ierr = nf90_get_var(ncid, ivarID, LakeTargVol_local)       ! read the variable
      if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; varname='//trim(varName); return; endif
      ! check if the size of islake_local and LakeTargVol_local
      if (size(islake_local)/=size(LakeTargVol_local))then
        ierr=20; message=trim(message)//'lake flag and lake target volume do not have similar length'; return;
      endif
    endif

    ! check if the size of islake_local and LakeModelType_local
    if (size(islake_local)/=size(LakeModelType_local))then
      ierr=20; message=trim(message)//'lake flag and lake types do not have similar length'; return;
    endif

    ! assign the initial values
    number_lakes      =  0
    number_Endorheic  =  0
    number_Doll       =  0
    number_Hanasaki   =  0
    number_HYPE       =  0
    number_TargVol    =  0

    !is_flux_wm
    !is_vol_wm

    ! specifying which lake models are called and if there is conflict between lake model type and data driven flag
    do i = 1, size(islake_local)
      if (islake_local(i) == 1) then ! if the segement is flagged as lake
        number_lakes = number_lakes + 1 ! total number of lakes
        if (is_vol_wm) then
          if (LakeTargVol_local(i) == 1) then ! if lake is not target volume then it should be parametertic
            number_TargVol = number_TargVol + 1 ! add number of target volume case
            select case(LakeModelType_local(i))
              case(0); ierr=20; message=trim(message)//'both data driven (follow target volume) and Endorheic lake are activated for a lake'; return
              case(1); ierr=20; message=trim(message)//'both data driven (follow target volume) and Doll lake formulation are activated for a lake'; return
              case(2); ierr=20; message=trim(message)//'both data driven (follow target volume) and Hanasaki lake formulation are activated for a lake'; return
              case(3); ierr=20; message=trim(message)//'both data driven (follow target volume) and HYPY lake formulation are activated for a lake'; return
            end select
          else
            select case(LakeModelType_local(i))
              case(0);                              number_Endorheic = number_Endorheic + 1; ! add number of Endorheic lakes
              case(1); lake_model_D03     = .true.; number_Doll      = number_Doll      + 1; ! add number of Doll lakes
              case(2); lake_model_H06     = .true.; number_Hanasaki  = number_Hanasaki  + 1; ! add number of Hanasaki lakes
              case(3); lake_model_HYPE    = .true.; number_HYPE      = number_HYPE      + 1; ! add number of HYPE lakes
              case default; ierr=20; message=trim(message)//'unable to identify the lake model type'; return
            end select
          endif
        else
          select case(LakeModelType_local(i))
            case(0);                              number_Endorheic = number_Endorheic + 1; ! add number of Endorheic lakes
            case(1); lake_model_D03     = .true.; number_Doll      = number_Doll      + 1; ! add number of Doll lakes
            case(2); lake_model_H06     = .true.; number_Hanasaki  = number_Hanasaki  + 1; ! add number of Hanasaki lakes
            case(3); lake_model_HYPE    = .true.; number_HYPE      = number_HYPE      + 1; ! add number of HYPE lakes
            case default; ierr=20; message=trim(message)//'unable to identify the lake model type'; return
          end select
        endif
      endif
    enddo

    ! print the numbers
    print*, "total number of lakes             = ", number_lakes
    print*, "total number of Endorheic lakes   = ", number_Endorheic
    print*, "lakes with Doll formulation       = ", number_Doll
    print*, "lakes with Hanasaki formulation   = ", number_Hanasaki
    print*, "lakes with HYPE formulation       = ", number_HYPE
    print*, "lakes with target volume          = ", number_TargVol

    ! check is the number of parameteric lakes and target volume sums up to the total number of lakes
    if (number_Endorheic+number_Doll+number_Hanasaki+number_HYPE+number_TargVol == number_lakes) then
      print*, "number of lake models and target volume models matches the total number of lakes; should be good to go!"
    else
      ! print*, "number of lake models and target volume models do not match the total number of lakes"
      ierr=20; message=trim(message)//'number of lake models and target volume models do not match the total number of lakes'; return
    endif

    ! warning about the Endorheic lakes and absence of evaporation
    if ((number_Endorheic>0).and.(suppress_P_Ep)) then
      print*, "The river network topology includes Endorheic lakes while no precipitation and evaporation is provided"
      print*, "In absence of extraction from lake the Endorheic lake volumns will be always increasing"
    endif

    ! deallocate the variables
    deallocate(islake_local, LakeModelType_local, stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem deallocating islake_local and LakeModelType_local'; return; endif

    ! print flags
    print*, "Doll is activated                 = ", lake_model_D03
    print*, "HYPE is activated                 = ", lake_model_HYPE
    print*, "Hanasaki is activated             = ", lake_model_H06


    if (lake_model_D03) then
      meta_SEG(ixSEG%D03_MaxStorage)%varFile    = .true.    ! Doll parameter
      meta_SEG(ixSEG%D03_coefficient)%varFile   = .true.    ! Doll parameter
      meta_SEG(ixSEG%D03_power)%varFile         = .true.    ! Doll parameter
    endif

    if (lake_model_HYPE) then
      meta_SEG(ixSEG%HYP_E_emr)%varFile         = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_E_lim)%varFile         = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_E_min)%varFile         = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_E_zero)%varFile        = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_Qrate_emr)%varFile     = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_Erate_emr)%varFile     = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_Qrate_prim)%varFile    = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_Qrate_amp)%varFile     = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_Qrate_phs)%varFile     = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_prim_F)%varFile        = .true.    ! HYPE parameter
      meta_SEG(ixSEG%HYP_A_avg)%varFile         = .true.    ! HYPE parameter
    endif

    if (lake_model_H06) then
      meta_SEG(ixSEG%H06_Smax)%varFile          = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_alpha)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_envfact)%varFile       = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_S_ini)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_c1)%varFile            = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_c2)%varFile            = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_exponent)%varFile      = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_denominator)%varFile   = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_c_compare)%varFile     = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_frac_Sdead)%varFile    = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_E_rel_ini)%varFile     = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Jan)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Feb)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Mar)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Apr)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_May)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Jun)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Jul)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Aug)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Sep)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Oct)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Nov)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_Dec)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Jan)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Feb)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Mar)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Apr)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_May)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Jun)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Jul)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Aug)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Sep)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Oct)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Nov)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_Dec)%varFile         = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_purpose)%varFile       = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_mem_F)%varFile       = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_mem_F)%varFile       = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_I_mem_L)%varFile       = .true.    ! Hanasaki parameter
      meta_SEG(ixSEG%H06_D_mem_L)%varFile       = .true.    ! Hanasaki parameter
    endif

  endif

  ! loop through data structures
  write(iulog,'(2a)') new_line('a'), '---- Read river network data --- '
  do iStruct=1,nStructures

    ! loop through the variables
    do iVar=1,meta_struct(iStruct)%nVars

      ! ---------- get information on variable ------------------------------------------------------------------------

      ! get the variable name
      select case(iStruct)
        case(ixStruct%HRU    ); varName=trim(meta_HRU(    ivar)%varName) ; isVarDesired=(meta_HRU(    ivar)%varFile)
        case(ixStruct%HRU2SEG); varName=trim(meta_HRU2SEG(ivar)%varName) ; isVarDesired=(meta_HRU2SEG(ivar)%varFile)
        case(ixStruct%SEG    ); varName=trim(meta_SEG(    ivar)%varName) ; isVarDesired=(meta_SEG(    ivar)%varFile)
        case(ixStruct%NTOPO  ); varName=trim(meta_NTOPO(  ivar)%varName) ; isVarDesired=(meta_NTOPO(  ivar)%varFile)
        case(ixStruct%PFAF   ); varName=trim(meta_PFAF(   ivar)%varName) ; isVarDesired=(meta_PFAF(   ivar)%varFile)
        case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
      end select

      ! check need to read the variable
      if(topoNetworkOption==readFromFile) isVarDesired=.true.
      if(.not.isVarDesired) cycle

      ! get the netCDF variable ID
      ierr = nf90_inq_varid(ncid, trim(varName), ivarID)
      if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(varName); return; endif

      ! get subset of indices for each reach
      call getSubetIndices(&
                          ! input
                          ncid,                       &  ! netCDF file id
                          ivarID,                     &  ! netCDF variable id
                          meta_struct(iStruct)%nSpace,&  ! length of the spatial dimension
                          ! output
                          ixStart,                    &  ! vector of start indices
                          ixCount,                    &  ! vector defining number of elements in each reach
                          dimLength,                  &  ! dimension length
                          ierr,cmessage)                 ! error control
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! skip for cases where te dimension length is zero
      if(dimLength==0) cycle

      ! ---------- read data into temporary structures ----------------------------------------------------------------

      ! print progress
      write(iulog,'(2x,a)') 'Reading '//trim(varName)//' into structure '//trim(meta_struct(iStruct)%structName)

      ! read data from NetCDF files
      select case(iStruct)

        ! integer vector
        case(ixStruct%HRU2SEG,ixStruct%NTOPO)

          ! allocate space
          if(size(iTemp)/=dimLength)then
            deallocate(iTemp,stat=ierr);          if(ierr/=0)then; message=trim(message)//'problem deallocating iTemp'; return; endif
            allocate(iTemp(dimLength),stat=ierr); if(ierr/=0)then; message=trim(message)//'problem allocating iTemp'; return; endif
          endif

          ! read data
          ierr = nf90_get_var(ncid, ivarID, iTemp)
          if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; varname='//trim(varName); return; endif

        ! double precision vector
        case(ixStruct%HRU,ixStruct%SEG)

          ! allocate space
          if(size(dTemp)/=dimLength)then
            deallocate(dTemp,stat=ierr);          if(ierr/=0)then; message=trim(message)//'problem deallocating dTemp'; return; endif
            allocate(dTemp(dimLength),stat=ierr); if(ierr/=0)then; message=trim(message)//'problem allocating dTemp'; return; endif
          endif

          ! read data
          ierr = nf90_get_var(ncid, ivarID, dTemp)
          if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; varname='//trim(varName); return; endif

        ! character vector 2D character array array (character-dimension x data dimension)
        case(ixStruct%PFAF)

          ! allocate space
          if(size(cTemp)/=dimLength)then
            deallocate(cTemp,stat=ierr);          if(ierr/=0)then; message=trim(message)//'problem deallocating cTemp'; return; endif
            allocate(cTemp(dimLength),stat=ierr); if(ierr/=0)then; message=trim(message)//'problem allocating cTemp'; return; endif
          endif

         ! read data
         ierr = nf90_get_var(ncid, ivarID, cTemp, start=(/1,1/), count=(/maxPfafLen, dimLength/))
         if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; varname='//trim(varName); return; endif

        ! check errors
        case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
      end select   ! selection of the data structure

      ! ---------- put the temporary data into data structures  -------------------------------------------------------

      ! loop through the spatial elements
      do iSpace=1,meta_struct(iStruct)%nSpace

        ! define the start and count index for a given spatial element
        jxStart = ixStart(iSpace)
        jxCount = ixCount(iSpace)

        ! allocate space for the data
        ! NOTE: if dimension length is equal to the size of the structure, then the data is a scalar and already allocated
        select case(iStruct)
          case(ixStruct%HRU    ); if(size(structHRU)    /=dimLength) allocate(structHRU(    iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
          case(ixStruct%HRU2SEG); if(size(structHRU2seg)/=dimLength) allocate(structHRU2seg(iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
          case(ixStruct%SEG    ); if(size(structSeg)    /=dimLength) allocate(structSeg(    iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
          case(ixStruct%NTOPO  ); if(size(structNTOPO)  /=dimLength) allocate(structNTOPO(  iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
          case(ixStruct%PFAF   ); if(size(structPFAF)   /=dimLength) allocate(structPFAF(   iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
          case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
        end select
        if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for the data vectors'; return; endif

        ! get the data
        select case(iStruct)
          case(ixStruct%HRU    ); structHRU(    iSpace)%var(iVar)%dat(1:jxCount) = dTemp(jxStart:jxStart+jxCount-1)  ! dp
          case(ixStruct%SEG    ); structSeg(    iSpace)%var(iVar)%dat(1:jxCount) = dTemp(jxStart:jxStart+jxCount-1)  ! dp
          case(ixStruct%HRU2SEG); structHRU2seg(iSpace)%var(iVar)%dat(1:jxCount) = iTemp(jxStart:jxStart+jxCount-1)  ! i4b
          case(ixStruct%NTOPO  ); structNTOPO(  iSpace)%var(iVar)%dat(1:jxCount) = iTemp(jxStart:jxStart+jxCount-1)  ! i4b
          case(ixStruct%PFAF   ); structPFAF(   iSpace)%var(iVar)%dat(1:jxCount) = cTemp(jxStart:jxStart+jxCount-1)  ! character
          case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
        end select

      end do  ! looping through space

      ! deallocate space
      deallocate(ixStart, ixCount, stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem deallocating space for array indices'; return; endif

    end do  ! looping through variables
  end do  ! looping through the structures

  ! deallocate space
  deallocate(iTemp, dTemp, cTemp, stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem deallocating space for temporary vectors'; return; endif

  ! close the NetCDF file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

!  ! ---------- check the consistency of the lake models  ----------------------------------------------------------
!
!  if (is_lake_sim) then
!    ! check which lakes model are specified in segment
!    do iStruct=1,nStructures
!      do iVar=1,meta_struct(iStruct)%nVars
!        do iSpace=1,meta_struct(iStruct)%nSpace
!          select case(iStruct)
!            case(ixStruct%HRU    )
!            case(ixStruct%HRU2SEG)
!            case(ixStruct%SEG    )
!            case(ixStruct%PFAF   )
!            case(ixStruct%NTOPO) ! just the NTOPO case
!              !print*, meta_struct(iStruct)%nVars, iVar, trim(meta_NTOPO(ivar)%varName), meta_struct(iStruct)%nSpace
!              !print*, trim(meta_NTOPO(ivar)%varName)
!              select case (iVar) ! get the index of the varibale
!                case(ixNTOPO%islake) !iVar is euqal to the location of islake
!                  !print*, trim(meta_NTOPO(ivar)%varName), ivar
!                  if (allocated(islake_local)) then
!                    islake_local(iSpace) = structNTOPO(iSpace)%var(iVar)%dat(1)
!                  else
!                    allocate(islake_local(meta_struct(iStruct)%nSpace),stat=ierr)
!                    if(ierr/=0)then; message=trim(message)//'problem allocating islake_local'; return; endif
!                    islake_local(iSpace) = structNTOPO(iSpace)%var(iVar)%dat(1)
!                  endif
!                case(ixNTOPO%laketargvol) !iVar is euqal to the location of laketargvol
!                  if (is_vol_wm) then
!                    !print*, trim(meta_NTOPO(ivar)%varName), ivar
!                    if (allocated(LakeTargVol_local)) then
!                      LakeTargVol_local(iSpace) = structNTOPO(iSpace)%var(iVar)%dat(1)
!                    else
!                      allocate(LakeTargVol_local(meta_struct(iStruct)%nSpace),stat=ierr)
!                      if(ierr/=0)then; message=trim(message)//'problem allocating LakeTargVol_local'; return; endif
!                      LakeTargVol_local(iSpace) = structNTOPO(iSpace)%var(iVar)%dat(1)
!                    endif
!                  else
!                    if (allocated(LakeTargVol_local)) then
!                      LakeTargVol_local(iSpace) = 0
!                    else
!                      allocate(LakeTargVol_local(meta_struct(iStruct)%nSpace),stat=ierr)
!                      if(ierr/=0)then; message=trim(message)//'problem allocating LakeTargVol_local'; return; endif
!                      LakeTargVol_local(iSpace) = 0
!                    endif
!                  endif
!                case(ixNTOPO%LakeModelType)
!                  if ((lake_model_D03).or.(lake_model_H06).or.(lake_model_HYPE)) then
!                    !print*, trim(meta_NTOPO(ivar)%varName), ivar
!                    if (allocated(LakeModelType_local)) then
!                      LakeModelType_local(iSpace) = structNTOPO(iSpace)%var(iVar)%dat(1)
!                    else
!                      allocate(LakeModelType_local(meta_struct(iStruct)%nSpace),stat=ierr)
!                      if(ierr/=0)then; message=trim(message)//'problem allocating LakeModelType_local'; return; endif
!                      LakeModelType_local(iSpace) = structNTOPO(iSpace)%var(iVar)%dat(1)
!                    endif
!                  else
!                    if (allocated(LakeModelType_local)) then
!                      LakeModelType_local(iSpace) = 0
!                    else
!                      allocate(LakeModelType_local(meta_struct(iStruct)%nSpace),stat=ierr)
!                      if(ierr/=0)then; message=trim(message)//'problem allocating LakeModelType_local'; return; endif
!                      LakeModelType_local(iSpace) = 0
!                    endif
!                  endif
!              end select
!            case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
!          end select
!        end do
!      end do
!    end do
!
!    write(iulog,'(2a)') new_line('a'), '---- Check the lake model types and flags for each lake --- '
!
!    ! initializing the flags
!    Doll_is_called        =  .false. ! initialize the flag to false
!    Hanasaki_is_called    =  .false. ! initialize the flag to false
!    HYPE_is_called        =  .false. ! initialize the flag to false
!    lake_model_conflict   =  .false. ! initialize the flag to false
!
!    ! check if the length of the local read varibales are equal
!    if (.not.(size (islake_local)==size (LakeTargVol_local)).and.(size (LakeTargVol_local)==size (LakeModelType_local)) ) then
!      ierr=20; message=trim(message)//'the lake flags, lake target volume flags and lake model type flags are with different dimensions'; return
!    endif
!
!    ! assign the initial values
!    number_lakes      =  0
!    number_Doll       =  0
!    number_Hanasaki   =  0
!    number_HYPE       =  0
!    number_TargVol    =  0
!
!    ! specifying which lake models are called and if there is conflict between lake model type and data driven flag
!    do i = 1, size(islake_local)
!      if (islake_local(i) == 1) then ! if the segement is flagged as lake
!        number_lakes = number_lakes + 1 ! total number of lakes
!        ! check if the lakes are data driven or parameteric
!        if (LakeTargVol_local(i) /= 1) then ! if lake is not target volume then it should be parametertic
!          select case(LakeModelType_local(i))
!            case(1); Doll_is_called     = .true.; number_Doll     = number_Doll     + 1; ! add number of Doll lakes
!            case(2); Hanasaki_is_called = .true.; number_Hanasaki = number_Hanasaki + 1; ! add number of Hanasaki lakes
!            case(3); HYPE_is_called     = .true.; number_HYPE     = number_HYPE     + 1; ! add number of HYPE lakes
!            case default; ierr=20; message=trim(message)//'unable to identify the lake model type'; return
!          end select
!        else ! if for a lake both parameteric and non parameteric are set to correct
!          number_TargVol = number_TargVol + 1 ! add number of target volume case
!          select case(LakeModelType_local(i))
!            case(1); ierr=20; message=trim(message)//'both data driven (follow target volume) and Doll lake formulation are activated for a lake'; return
!            case(2); ierr=20; message=trim(message)//'both data driven (follow target volume) and Hanasaki lake formulation are activated for a lake'; return
!            case(3); ierr=20; message=trim(message)//'both data driven (follow target volume) and HYPY lake formulation are activated for a lake'; return
!          end select
!        endif
!      endif
!    enddo
!
!    ! check the local lake model flags with the global flags
!    if ((Doll_is_called).and.(.not.(lake_model_D03))) then
!      ierr=20; message=trim(message)//'Doll 2003 is called in the netwrok topology but the flag is not set to true in control file; set the flag to true in control file'; return
!    endif
!
!    ! check the local lake model type with the global flags
!    if ((Hanasaki_is_called).and.(.not.(lake_model_H06))) then
!      ierr=20; message=trim(message)//'Hanasaki 2006 is called in the netwrok topology but the flag is not set to true in control file; set the flag to true in control file'; return
!    endif
!
!    ! check the local lake model type with the global flags
!    if ((HYPE_is_called).and.(.not.(lake_model_HYPE))) then
!      ierr=20; message=trim(message)//'HYPE is called in the netwrok topology but the flag is not set to true in control file; set the flag to true in control file'; return
!    endif
!
!    ! print the numbers
!    print*, "total number of lakes             = ", number_lakes
!    print*, "lakes with Doll formulation       = ", number_Doll
!    print*, "lakes with Hanasaki formulation   = ", number_Hanasaki
!    print*, "lakes with HYPE formulation       = ", number_HYPE
!    print*, "lakes with target volume          = ", number_TargVol
!
!    ! check is the number of parameteric lakes and target volume sums up to the total number of lakes
!    if (number_Doll+number_Hanasaki+number_HYPE+number_TargVol == number_lakes) then
!      print*, "number of lake models and target volume models matches the total number of lakes; should be good to go!"
!    else
!      ! print*, "number of lake models and target volume models do not match the total number of lakes"
!      ierr=20; message=trim(message)//'number of lake models and target volume models do not match the total number of lakes'; return
!    endif
!  endif

end subroutine getData

! ==========================================================================================================
! ==========================================================================================================
! ==========================================================================================================
! ==========================================================================================================
! ==========================================================================================================

! *********************************************************************
! private subroutine: get start and count vectors
! *********************************************************************
subroutine getSubetIndices(&
                          ! input
                          ncid,           &  ! netCDF file id
                          ivarid,         &  ! netCDF variable id
                          nSpace,         &  ! length of the spatial dimension
                          ! output
                          ixStart,        &  ! vector of start indices
                          ixCount,        &  ! vector defining number of elements in each reach
                          dimLength,      &  ! dimension length
                          ierr,message)      ! error control
  implicit none
  ! input variables
  integer(i4b)  , intent(in)                :: ncid           ! netCDF file id
  integer(i4b)  , intent(in)                :: ivarid         ! netCDF variable id
  integer(i4b)  , intent(in)                :: nSpace         ! length of the spatial dimension
  ! output variables
  integer(i4b)  , intent(out) , allocatable :: ixStart(:)     ! vector of start indices
  integer(i4b)  , intent(out) , allocatable :: ixCount(:)     ! vector defining number of elements in each reach
  integer(i4b)  , intent(out)               :: dimLength      ! dimension length
  integer(i4b)  , intent(out)               :: ierr           ! error code
  character(*)  , intent(out)               :: message        ! error message
  ! local variables
  integer(i4b), allocatable                 :: ncDimIDs(:)    ! dimension IDs for a given variable
  character(len=strLen)                     :: dimName        ! dimension name
  logical(lgt)                              :: isRaggedArray  ! logical flag to denote a ragged array
  integer(i4b)                              :: nDims          ! number of dimensions in a variable
  integer(i4b)                              :: iStartID       ! ID for start of ragged array
  integer(i4b)                              :: iCountID       ! ID for count of ragged array

  ierr=0; message='getSubetIndices/'

  ! get the variable type
  ierr = nf90_inquire_variable(ncid, ivarID, ndims=nDims)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  allocate(ncDimIDs(nDims))

  ! get the dimension ID -- vector of length=1
  ierr = nf90_inquire_variable(ncid, ivarID, dimids=ncDimIDs)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! get the dimension name and length
  ierr = nf90_inquire_dimension(ncid, ncDimIDs(nDims), dimName, dimLength)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dimName); return; endif

  ! allocate space for the ragged arrays
  allocate(ixStart(nSpace), ixCount(nSpace), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for ragged array indices'; return; endif

  ! check if it is a ragged array
  isRaggedArray=(dimLength/=nSpace)
  if(isRaggedArray)then

    ! get the start index in the ragged arrays
    ! -- variable ID
    ierr = nf90_inq_varid(ncid, trim(dimName)//'_start', iStartID)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dimName)//'_start'; return; endif
    ! -- variable data
    ierr = nf90_get_var(ncid, iStartID, ixStart)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dimName)//'_start'; return; endif

    ! get the count of the ragged arrays
    ! -- variable ID
    ierr = nf90_inq_varid(ncid, trim(dimName)//'_count', iCountID)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dimName)//'_count'; return; endif
    ! -- variable data
    ierr = nf90_get_var(ncid, iCountID, ixCount)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dimName)//'_count'; return; endif

    ! not a ragged array
  else

    ! define array indices
    ixStart(:) = arth(1,1,nSpace)
    ixCount(:) = 1

  endif

end subroutine getSubetIndices

end module read_streamSeg
