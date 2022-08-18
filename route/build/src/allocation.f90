MODULE alloc_data

USE nrtype
! data types
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
USE var_lookup, ONLY: ixDims,   nDimensions  ! index of dimensions
USE var_lookup, ONLY: ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup, ONLY: ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup, ONLY: ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup, ONLY: ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology
USE var_lookup, ONLY: ixPFAF,   nVarsPFAF    ! index of variables for the pfafstetter code

implicit none

private
public::alloc_struct

CONTAINS

 ! *********************************************************************
 ! new subroutine: get ancillary data for HRUs and stream segments
 ! *********************************************************************
 SUBROUTINE alloc_struct(nHRU,         & ! input: number of HRUs
                         nSeg,         & ! input: number of stream segments
                         structHRU,    & ! input/output: ancillary data for HRUs
                         structSeg,    & ! input/output: ancillary data for stream segments
                         structHRU2seg,& ! input/output: ancillary data for mapping hru2basin
                         structNTOPO,  & ! input/output: ancillary data for network toopology
                         structPFAF,   & ! input/output: ancillary data for pfafstetter code
                         ierr,message)   ! output: error control
 implicit none
 ! Argument variables
 integer(i4b)      , intent(in)                 :: nHRU             ! number of HRUs
 integer(i4b)      , intent(in)                 :: nSeg             ! number of stream segments
 type(var_dlength) , intent(inout), allocatable :: structHRU(:)     ! HRU properties
 type(var_dlength) , intent(inout), allocatable :: structSeg(:)     ! stream segment properties
 type(var_ilength) , intent(inout), allocatable :: structHRU2seg(:) ! HRU-to-segment mapping
 type(var_ilength) , intent(inout), allocatable :: structNTOPO(:)   ! network topology
 type(var_clength) , intent(inout), allocatable :: structPFAF(:)    ! network topology
 integer(i4b)      , intent(out)                :: ierr             ! error code
 character(*)      , intent(out)                :: message          ! error message
 ! local variables
 integer(i4b)                                   :: iStruct      ! structure index
 integer(i4b)                                   :: iSpace       ! spatial index
 integer(i4b)                                   :: iHRU         ! HRU index
 integer(i4b)                                   :: iSeg         ! segment index
 integer(i4b)                                   :: iVar         ! variable index
 logical(lgt)                                   :: isDimScalar  ! .true. if the dimension is a scalar

 ierr=0; message='alloct_struc/'

 ! ---------- allocate space for higher-level structure components -------------------------------------------------

 ! allocate the spatial dimension in all data structures
 allocate(structHRU(nHRU), structHRU2seg(nHRU), structSeg(nSeg), structNTOPO(nSeg), structPFAF(nSeg), stat=ierr)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating [structHRU,structHRU2seg,structNTOPO,structPFAF]'; return; endif

 ! allocate the variable dimension in the data structures with length nHRU
 do iHRU=1,nHRU
  allocate(structHRU(iHRU)%var(nVarsHRU), structHRU2seg(iHRU)%var(nVarsHRU2SEG), stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating variables for HRUs'; return; endif
 end do

 ! allocate the variable dimension in the data structures with length nSeg
 do iSeg=1,nSeg
  allocate(structSeg(iSeg)%var(nVarsSEG), structNTOPO(iSeg)%var(nVarsNTOPO), structPFAF(iSeg)%var(nVarsPFAF),stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating variables for stream segments'; return; endif
 end do

 ! ---------- allocate space for the scalar variables --------------------------------------------------------------

 ! loop through data structures
 do iStruct=1,nStructures

  ! populate the spatial dimension
  select case(iStruct)
   case(ixStruct%HRU, ixStruct%HRU2SEG);              meta_struct(iStruct)%nSpace=nHRU
   case(ixStruct%SEG, ixStruct%NTOPO, ixStruct%PFAF); meta_struct(iStruct)%nSpace=nSeg
   case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
  end select

  ! loop through the spatial elements
  do iSpace=1,meta_struct(iStruct)%nSpace

   ! loop through the variables
   do iVar=1,meta_struct(iStruct)%nVars

    ! define the need to allocate
    select case(iStruct)
     case(ixStruct%HRU    ); isDimScalar = ( meta_HRU(    ivar)%varType==ixDims%hru .or. meta_HRU(    ivar)%varType==ixDims%seg )
     case(ixStruct%HRU2SEG); isDimScalar = ( meta_HRU2SEG(ivar)%varType==ixDims%hru .or. meta_HRU2SEG(ivar)%varType==ixDims%seg )
     case(ixStruct%SEG    ); isDimScalar = ( meta_SEG(    ivar)%varType==ixDims%hru .or. meta_SEG(    ivar)%varType==ixDims%seg )
     case(ixStruct%NTOPO  ); isDimScalar = ( meta_NTOPO(  ivar)%varType==ixDims%hru .or. meta_NTOPO(  ivar)%varType==ixDims%seg )
     case(ixStruct%PFAF   ); isDimScalar = ( meta_PFAF(   ivar)%varType==ixDims%hru .or. meta_PFAF(   ivar)%varType==ixDims%seg )
     case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
    end select

    ! allocate space for the data
    select case(iStruct)
     case(ixStruct%HRU    ); if(isDimScalar) allocate(structHRU(    iSpace)%var(iVar)%dat(1), stat=ierr)
     case(ixStruct%HRU2SEG); if(isDimScalar) allocate(structHRU2seg(iSpace)%var(iVar)%dat(1), stat=ierr)
     case(ixStruct%SEG    ); if(isDimScalar) allocate(structSeg(    iSpace)%var(iVar)%dat(1), stat=ierr)
     case(ixStruct%NTOPO  ); if(isDimScalar) allocate(structNTOPO(  iSpace)%var(iVar)%dat(1), stat=ierr)
     case(ixStruct%PFAF   ); if(isDimScalar) allocate(structPFAF(   iSpace)%var(iVar)%dat(1), stat=ierr)
     case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
    end select
    if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for the data vectors'; return; endif

   end do  ! loop through variab;es
  end do  ! loop through space
 end do  ! loop through structures
 ! ---------- initialize variables ---------------------------------------------------------------------------

 ! loop through stream segments
 do iSeg=1,nSeg

  ! initialize variables not yet computed / assigned
  structSEG(iSeg)%var(ixSEG%width          )%dat(1) = realMissing
  structSEG(iSeg)%var(ixSEG%man_n          )%dat(1) = realMissing
  structSEG(iSeg)%var(ixSEG%upsArea        )%dat(1) = realMissing
  structSEG(iSeg)%var(ixSEG%basUnderLake   )%dat(1) = realMissing
  structSEG(iSeg)%var(ixSEG%rchUnderLake   )%dat(1) = realMissing
  structSEG(iSeg)%var(ixSEG%minFlow        )%dat(1) = realMissing

  ! initialize variables not yet computed / assigned
  structNTOPO(iSeg)%var(ixNTOPO%rchOrder   )%dat(1) = integerMissing
  structNTOPO(iSeg)%var(ixNTOPO%lakeId     )%dat(1) = integerMissing
  structNTOPO(iSeg)%var(ixNTOPO%lakeIndex  )%dat(1) = integerMissing
  structNTOPO(iSeg)%var(ixNTOPO%isLakeInlet)%dat(1) = integerMissing
  structNTOPO(iSeg)%var(ixNTOPO%userTake   )%dat(1) = integerMissing

 end do  ! looping through stream segments

 END SUBROUTINE alloc_struct

END MODULE alloc_data
