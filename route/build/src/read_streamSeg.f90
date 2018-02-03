module read_streamSeg

! data types
USE nrtype
USE nrtype,    only : strLen            ! string length
USE nrtype,    only : integerMissing    ! missing value for integers
USE dataTypes, only : var_ilength       ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength       ! double precision type: var(:)%dat

! old data structures
USE nhru2basin, only : all_points        ! derived data types for hru2segment mapping
USE reachparam, only : RCHTOPO           ! Network topology structure
USE reachparam, only : RCHPRP            ! Reach Parameter structure

! metadata on data structures
USE globalData, only : meta_struct  ! structure information
USE globalData, only : meta_HRU     ! HRU properties
USE globalData, only : meta_HRU2SEG ! HRU-to-segment mapping
USE globalData, only : meta_SEG     ! stream segment properties
USE globalData, only : meta_NTOPO   ! network topology

! named variables
USE globalData, only : scalar       ! scalar variable
USE globalData, only : vector       ! vector variable

! named variables
USE var_lookup,only:ixStruct, nStructures  ! index of data structures
USE var_lookup,only:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup,only:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup,only:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup,only:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology

! old named variables
USE var_lookup,only:ixMAP,nVarsMAP   ! index of variables for the hru2segment mapping
USE var_lookup,only:ixTOP,nVarsTOP   ! index of variables for the network topology

! old data structures
USE dataTypes,  only : namepvar,nameivar      ! provide access to data types

! global data
USE public_var

! netcdf modules
USE netcdf

! external utilities
USE nr_utility_module, ONLY: indexx  ! Num. Recipies utilities
USE nr_utility_module, ONLY: arth    ! Num. Recipies utilities

implicit none

! define module-level constants
integer(i4b),parameter  :: down2noSegment=0     ! index in the input file if the HRU does not drain to a segment

! privacy
private
public::getData
public::hru2segment
public::up2downSegment
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
                    nHRU,         & ! output: number of HRUs
                    nSeg,         & ! output: number of stream segments
                    ! output: populate data structures
                    structHRU,    & ! ancillary data for HRUs
                    structSeg,    & ! ancillary data for stream segments
                    structHRU2seg,& ! ancillary data for mapping hru2basin
                    structNTOPO,  & ! ancillary data for network toopology
                    ! output: error control
                    ierr,message)   ! output: error control
 implicit none
 ! input variables
 character(*)      , intent(in)               :: fname            ! filename
 character(*)      , intent(in)               :: dname_nhru       ! dimension name for HRUs
 character(*)      , intent(in)               :: dname_sseg       ! dimension name for stream segments
 ! output: model control
 integer(i4b)      , intent(out)              :: nHRU             ! number of HRUs
 integer(i4b)      , intent(out)              :: nSeg             ! number of stream segments
 ! output: data structures
 type(var_dlength) , intent(out), allocatable :: structHRU(:)     ! HRU properties
 type(var_dlength) , intent(out), allocatable :: structSeg(:)     ! stream segment properties
 type(var_ilength) , intent(out), allocatable :: structHRU2seg(:) ! HRU-to-segment mapping
 type(var_ilength) , intent(out), allocatable :: structNTOPO(:)   ! network topology
 ! output: error control
 integer(i4b)      , intent(out)              :: ierr             ! error code
 character(*)      , intent(out)              :: message          ! error message
 ! ==========================================================================================================
 ! local variables
 integer(i4b)                        :: iStruct      ! structure index
 integer(i4b)                        :: iSpace       ! spatial index
 integer(i4b)                        :: iHRU         ! HRU index
 integer(i4b)                        :: iSeg         ! segment index
 integer(i4b)                        :: iVar         ! variable index
 integer(i4b)                        :: ncid         ! NetCDF file ID
 integer(i4b)                        :: idimID_nHRU  ! dimension ID for HRUs
 integer(i4b)                        :: idimID_sseg  ! dimension ID for stream segments
 integer(i4b)                        :: iVarID       ! variable ID
 integer(i4b)                        :: jxStart      ! Start index for a given reach
 integer(i4b)                        :: jxCount      ! Number of elements for a given reach
 integer(i4b), allocatable           :: ixStart(:)   ! Start index for each reach
 integer(i4b), allocatable           :: ixCount(:)   ! Number of elements in each reach
 integer(i4b), allocatable           :: iTemp(:)     ! temporary integer vector
 real(dp),     allocatable           :: dTemp(:)     ! temporary double precision vector
 integer(i4b)                        :: dimLength    ! dimension length
 logical(lgt)                        :: isVarDesired ! .true. if the variable is desired
 character(len=strLen)               :: varName      ! variable name
 character(len=strLen)               :: cmessage     ! error message of downwind routine
 ! initialize error control
 ierr=0; message='getData/'

 ! ---------- initial reading of dimensions ------------------------------------------------------------------------

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; file='//trim(fname); return; endif

 ! get the ID of the HRU dimension
 ierr = nf90_inq_dimid(ncid, dname_nhru, idimID_nHRU)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_nhru); return; endif

 ! get the length of the HRU dimension
 ierr = nf90_inquire_dimension(ncid, idimID_nHRU, len=nHRU)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the ID of the stream segment dimension
 ierr = nf90_inq_dimid(ncid, dname_sseg, idimID_sseg)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_sseg); return; endif

 ! get the length of the stream segment dimension
 ierr = nf90_inquire_dimension(ncid, idimID_sseg, len=nSeg)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! ---------- allocate space for higher-level structure components -------------------------------------------------

 ! print progress
 print*, 'Allocating space for the higher-level structure components'; call flush(6)

 ! allocate the spatial dimension in all data structures
 allocate(structHRU(nHRU), structHRU2seg(nHRU), structSeg(nSeg), structNTOPO(nSeg), stat=ierr)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating spatial dimension for data structures'; return; endif

 ! allocate the variable dimension in the data structures with length nHRU
 do iHRU=1,nHRU
  allocate(structHRU(iHRU)%var(nVarsHRU), structHRU2seg(iHRU)%var(nVarsHRU2SEG), stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating variables for HRUs'; return; endif
 end do

 ! allocate the variable dimension in the data structures with length nSeg
 do iSeg=1,nSeg
  allocate(structSeg(iSeg)%var(nVarsSEG), structNTOPO(iSeg)%var(nVarsNTOPO), stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating variables for stream segments'; return; endif
 end do

 ! initial allocation of the temporary vectors
 allocate(iTemp(nHRU), dTemp(nHRU), stat=ierr)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating temporary vectors'; return; endif

 ! ---------- allocate space for the scalar variables --------------------------------------------------------------

 ! print progress
 print*, 'Allocating space for the scalar variables'; call flush(6)

 ! loop through data structures
 do iStruct=1,nStructures

  ! populate the spatial dimension
  select case(iStruct)
   case(ixStruct%HRU, ixStruct%HRU2SEG); meta_struct(iStruct)%nSpace=nHRU
   case(ixStruct%SEG, ixStruct%NTOPO  ); meta_struct(iStruct)%nSpace=nSeg
   case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
  end select

  ! loop through the spatial elements
  do iSpace=1,meta_struct(iStruct)%nSpace

   ! loop through the variables
   do iVar=1,meta_struct(iStruct)%nVars

    ! allocate space for the data
    select case(iStruct)
     case(ixStruct%HRU    ); if(meta_HRU(    ivar)%varType==scalar) allocate(structHRU(    iSpace)%var(iVar)%dat(1), stat=ierr)
     case(ixStruct%HRU2SEG); if(meta_HRU2SEG(ivar)%varType==scalar) allocate(structHRU2seg(iSpace)%var(iVar)%dat(1), stat=ierr)
     case(ixStruct%SEG    ); if(meta_SEG(    ivar)%varType==scalar) allocate(structSeg(    iSpace)%var(iVar)%dat(1), stat=ierr)
     case(ixStruct%NTOPO  ); if(meta_NTOPO(  ivar)%varType==scalar) allocate(structNTOPO(  iSpace)%var(iVar)%dat(1), stat=ierr)
     case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
    end select
    if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for the data vectors'; return; endif

   end do  ! loop through variab;es
  end do  ! loop through space
 end do  ! loop through structures

 ! -----------------------------------------------------------------------------------------------------------------
 ! ---------- read in data -----------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------

 ! loop through data structures
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
    case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
   end select

   ! check need to read the variable
   if(topoNetworkOption/=0 .and. .not.isVarDesired) cycle
   print*, 'Reading '//trim(varName)//' into structure '//trim(meta_struct(iStruct)%structName)
   call flush(6)

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

   ! ---------- read data into temporary structures ----------------------------------------------------------------

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
    select case(iStruct)
     case(ixStruct%HRU    ); if(meta_HRU(    ivar)%varType/=scalar) allocate(structHRU(    iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
     case(ixStruct%HRU2SEG); if(meta_HRU2SEG(ivar)%varType/=scalar) allocate(structHRU2seg(iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
     case(ixStruct%SEG    ); if(meta_SEG(    ivar)%varType/=scalar) allocate(structSeg(    iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
     case(ixStruct%NTOPO  ); if(meta_NTOPO(  ivar)%varType/=scalar) allocate(structNTOPO(  iSpace)%var(iVar)%dat( ixCount(iSpace) ), stat=ierr)
     case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
    end select
    if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for the data vectors'; return; endif

    ! get the data
    select case(iStruct)
     case(ixStruct%HRU    ); structHRU(    iSpace)%var(iVar)%dat(1:jxCount) = dTemp(jxStart:jxStart+jxCount-1)  ! dp
     case(ixStruct%SEG    ); structSeg(    iSpace)%var(iVar)%dat(1:jxCount) = dTemp(jxStart:jxStart+jxCount-1)  ! dp
     case(ixStruct%HRU2SEG); structHRU2seg(iSpace)%var(iVar)%dat(1:jxCount) = iTemp(jxStart:jxStart+jxCount-1)  ! i4b
     case(ixStruct%NTOPO  ); structNTOPO(  iSpace)%var(iVar)%dat(1:jxCount) = iTemp(jxStart:jxStart+jxCount-1)  ! i4b
     case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
    end select

   end do  ! looping through space

   ! deallocate space
   deallocate(ixStart, ixCount, stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem deallocating space for array indices'; return; endif

  end do  ! looping through variables
 end do  ! looping through the structures

 ! deallocate space
 deallocate(iTemp, dTemp, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for temporary vectors'; return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine getData

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
 integer(i4b), dimension(1)                :: ncDimIDs       ! dimension IDs for a given variable
 character(len=strLen)                     :: dimName        ! dimension name
 logical(lgt)                              :: isRaggedArray  ! logical flag to denote a ragged array
 integer(i4b)                              :: iStartID       ! ID for start of ragged array
 integer(i4b)                              :: iCountID       ! ID for count of ragged array
 ! initialize error control
 ierr=0; message='getSubetIndices/'

 ! get the dimension ID -- vector of length=1
 ierr = nf90_inquire_variable(ncid, ivarID, dimids=ncDimIDs)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the dimension name and length
 ierr = nf90_inquire_dimension(ncid, ncDimIDs(1), dimName, dimLength)
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


 ! *********************************************************************
 ! new subroutine: compute correspondence between HRUs and segments
 ! *********************************************************************
 subroutine hru2segment(&
                        ! input
                        nHRU,       &   ! input: number of HRUs
                        nSeg,       &   ! input: number of stream segments
                        ! input-output: data structures
                        structHRU,     & ! ancillary data for HRUs
                        structSeg,     & ! ancillary data for stream segments
                        structHRU2seg, & ! ancillary data for mapping hru2basin
                        structNTOPO,   & ! ancillary data for network toopology
                        ! output
                        total_hru,  &   ! output: total number of HRUs that drain into any segments
                        ierr, message)  ! output: error control
 implicit none
 ! input variables
 integer(i4b), intent(in)                      :: nHRU              ! number of HRUs
 integer(i4b), intent(in)                      :: nSeg              ! number of stream segments
 ! input-output: data structures
 type(var_dlength), intent(inout), allocatable :: structHRU(:)      ! HRU properties
 type(var_dlength), intent(inout), allocatable :: structSeg(:)      ! stream segment properties
 type(var_ilength), intent(inout), allocatable :: structHRU2seg(:)  ! HRU-to-segment mapping
 type(var_ilength), intent(inout), allocatable :: structNTOPO(:)    ! network topology
 ! output variables
 integer(i4b), intent(out)                     :: total_hru         ! total number of HRUs that drain into any segments
 integer(i4b), intent(out)                     :: ierr              ! error code
 character(*), intent(out)                     :: message           ! error message
 ! local variables
 logical(lgt),parameter          :: checkMap=.true.   ! flag to check the mapping
 character(len=strLen)           :: cmessage          ! error message of downwind routine
 integer(i4b)                    :: hruIndex          ! index of HRU (from another data structure)
 integer(i4b)                    :: iHRU              ! index of HRU
 integer(i4b)                    :: iSeg              ! index of stream segment
 integer(i4b)                    :: segId(nSeg)       ! unique identifier of the HRU
 integer(i4b)                    :: hruSegId(nHRU)    ! unique identifier of the segment where HRU drains
 integer(i4b)                    :: segHRUix(nHRU)    ! index of segment where HRU drains
 integer(i4b)                    :: nHRU2seg(nSeg)    ! number of HRUs that drain into a given segment
 real(dp)                        :: totarea           ! total area of all HRUs feeding into a given stream segment (m2)
 !integer*8                       :: time0,time1       ! times

 ! initialize error control
 ierr=0; message='hru2segment/'

 !print*, 'PAUSE: start of '//trim(message); read(*,*)

 ! initialize timing
 !call system_clock(time0)

 ! ---------- get the index of the stream segment that a given HRU drains into ------------------------------

 ! get input vectors
 do iSeg=1,nSeg; segId(iSeg)    = structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1); end do
 do iHRU=1,nHRU; hruSegId(iHRU) = structHRU2seg(iHRU)%var(ixHRU2seg%hruSegId)%dat(1); end do

 call downReachIndex(&
                     ! input
                     nHRU,          & ! number of upstream elements
                     nSeg,          & ! number of stream segments
                     segId,         & ! unique identifier of the stream segments
                     hruSegId,      & ! unique identifier of the segment where water drains
                     ! output
                     segHRUix,      & ! index of downstream stream segment
                     nHRU2seg,      & ! number of elements that drain into each segment
                     ierr,cmessage)   ! error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the total number of HRUs that drain into any segments
 total_hru = sum(nHRU2seg)

 ! ---------- allocate space for the mapping structures -----------------------------------------------------

 ! loop through stream segments
 do iSeg=1,nSeg
  ! allocate space (number of elements that drain into each segment)
  allocate(structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat( nHRU2seg(iSeg) ), &
           structNTOPO(iSeg)%var(ixNTOPO%hruContribId)%dat( nHRU2seg(iSeg) ), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for hru2seg structure component'; return; endif
  ! initialize the number of HRUs
  structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1) = 0
 end do

 ! get timing
 !call system_clock(time1)
 !print*, 'timing: allocate space = ', time1-time0

 ! ---------- populate structure components for HRU-2-Segment mapping ---------------------------------------

 ! loop through HRUs
 do iHRU=1,nHRU

  ! identify the index of the stream segment that the HRU drains into
  iSeg = segHRUix(iHRU)

  ! if there is no stream segment associated with current hru
  if (iSeg == integerMissing) cycle

  ! associate variables in data structure
  associate(nContrib       => structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1),      & ! contributing HRUs
            hruContribIx   => structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat, & ! index of contributing HRU
            hruContribId   => structNTOPO(iSeg)%var(ixNTOPO%hruContribId)%dat  ) ! unique ids of contributing HRU

  ! increment the HRU counter
  nContrib = nContrib + 1

  ! populate structure components
  hruContribIx(nContrib)   = iHRU
  hruContribId(nContrib)   = structHRU2seg(iHRU)%var(ixHRU2SEG%HRUid)%dat(1)

  ! end associations
  end associate

 end do ! looping through HRUs

 ! check
 if(checkMap)then
  do iSeg=1,nSeg
   if(nHRU2seg(iSeg)/=structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1))then
    message=trim(message)//'problems identifying the HRUs draining into stream segment'
    ierr=20; return
   endif
  end do
 endif

 ! get timing
 !call system_clock(time1)
 !print*, 'timing: populate structure components = ', time1-time0

 ! ---------- compute HRU weights ---------------------------------------------------------------------------

 ! loop through segments
 do iSeg=1,nSeg

   ! compute total area of the HRUs draining to the stream segment
   totarea = 0._dp
   do iHRU=1,structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1)
    hruIndex = structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat(iHRU)
    totarea  = totarea + structHRU(hruIndex)%var(ixHRU%area)%dat(1)
   end do

   ! compute the weights
   structHRU(iHRU)%var(ixHRU%weight)%dat(1) = structHRU(iHRU)%var(ixHRU%area)%dat(1) / totarea

 end do  ! (looping thru stream segments)

 ! get timing
 !call system_clock(time1)
 !print*, 'timing: compute HRU weights = ', time1-time0
 !print*, 'PAUSE: end of '//trim(message); read(*,*)

 end subroutine hru2segment


 ! *********************************************************************
 ! new subroutine: mapping between upstream and downstream segments
 ! *********************************************************************
 subroutine up2downSegment(&
                           ! input
                           nRch,         &    ! input: number of stream segments
                           sseg_acil,    &    ! input: stream segment parameters
                           ntop_acil,    &    ! input: network topology
                           hru2seg,      &    ! input: hru-segment mapping structure
                           ! output
                           NETOPO,       &    ! output: River Network topology
                           RPARAM,       &    ! output: Reach Parameters
                           total_upseg,  &    ! output: sum of immediate upstream segments
                           ierr, message)     ! output (error control)
 implicit none
 ! input variables
 integer(i4b),   intent(in)                 :: nRch         ! number of reaches
 type(namepvar), intent(in)                 :: sseg_acil(:) ! ancillary data for stream segments
 type(nameivar), intent(in)                 :: ntop_acil(:) ! ancillary data for the network topology
 type(all_points), intent(in)               :: hru2seg(:)   ! hru-segment mapping structure
 ! output variables
 type(RCHTOPO),  intent(out) , allocatable  :: NETOPO(:)    ! River Network topology
 type(RCHPRP),   intent(out) , allocatable  :: RPARAM(:)    ! Reach Parameters
 integer(i4b),   intent(out)                :: total_upseg  ! sum of immediate upstream segments
 integer(i4b),   intent(out)                :: ierr         ! error code
 character(*),   intent(out)                :: message      ! error message
 ! local variables
 logical(lgt),parameter          :: checkMap=.true.     ! flag to check the mapping
 character(len=strLen)           :: cmessage            ! error message of downwind routine
 integer(i4b)                    :: iRch                ! reach index
 integer(i4b)                    :: ixDownRch           ! index of the downstream reach
 integer(i4b)                    :: downIndex(nRch)     ! index of downstream stream segment
 integer(i4b)                    :: nUpstream(nRch)     ! number of elements that drain into each segment
 integer(i4b)                    :: mUpstream(nRch)     ! number of elements that drain into each segment
 real(dp),parameter              :: min_slope=1.e-6_dp  ! minimum slope
 ! initialize error control
 ierr=0; message='up2downSegment/'

 ! ---------- initialization ---------------------------------------------------------------------------------

 ! allocate space for the reach parameter structure
 allocate(RPARAM(nRch),NETOPO(nRch),stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for reach parameter structures'; return; endif

 ! transfer information to the network topology structures
 NETOPO(:)%REACHIX = arth(1,1,nRch)
 NETOPO(:)%REACHID = ntop_acil(ixTOP%segid    )%varData(:)
 NETOPO(:)%DREACHK = ntop_acil(ixTOP%toSegment)%varData(:)

 ! transfer information to the reach structures
 RPARAM(:)%RLENGTH = sseg_acil(ixSEG%length   )%varData(:)
 RPARAM(:)%R_SLOPE = sseg_acil(ixSEG%slope    )%varData(:)

 ! loop through reaches
 do iRch=1,nRch

  ! compute area draining to each stream segment
  RPARAM(iRch)%BASAREA = sum(hru2seg(iRch)%cHRU(:)%hru_area)

  ! ensure that slope exceeds minimum slope
  if(RPARAM(iRch)%R_SLOPE < min_slope) RPARAM(iRch)%R_SLOPE = min_slope

 end do  ! looping through reaches

 ! just to be safe, specify some things that we should not need
 NETOPO(:)%RCHLAT1 =  huge(kind(dp))    ! Start latitude
 NETOPO(:)%RCHLAT2 =  huge(kind(dp))    ! End latitude
 NETOPO(:)%RCHLON1 =  huge(kind(dp))    ! Start longitude
 NETOPO(:)%RCHLON2 =  huge(kind(dp))    ! End longitude
 NETOPO(:)%LAKE_IX =  -1                ! Lake index (0,1,2,...,nlak-1)
 NETOPO(:)%LAKE_ID =  -1                ! Lake ID (REC code?)
 NETOPO(:)%BASULAK =   0._dp            ! Area of basin under lake
 NETOPO(:)%RCHULAK =   0._dp            ! Length of reach under lake
 NETOPO(:)%LAKINLT = .false.            ! .TRUE. if reach is lake inlet, .FALSE. otherwise
 NETOPO(:)%USRTAKE = .false.            ! .TRUE. if user takes from reach, .FALSE. otherwise

 ! ---------- define the index of the downstream reach ID ----------------------------------------------------

 call downReachIndex(&
                     ! input
                     nRch,                              & ! number of upstream elements
                     nRch,                              & ! number of stream segments
                     ntop_acil(ixTOP%segid)%varData,    & ! unique identifier of the stream segments
                     ntop_acil(ixTOP%toSegment)%varData,& ! unique identifier of the segment where water drains
                     ! output
                     downIndex,                         & ! index of downstream stream segment
                     nUpstream,                         & ! number of elements that drain into each segment
                     ierr,cmessage)                       ! error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the total number of HRUs that drain into any segments
 total_upseg = sum(nUpstream)

 ! ---------- allocate space for the number of upstream reaches ---------------------------------------------

 ! loop through the reaches
 do iRch=1,nRch
  allocate(NETOPO(iRch)%UREACHI(nUpstream(iRch)),NETOPO(iRch)%UREACHK(nUpstream(iRch)), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for upstream reaches'; return; endif
 end do

 ! ---------- populate data structures for the upstream reaches ----------------------------------------------

 ! initialize the number of upstream elements in each reach
 mUpstream(:)=0

 ! loop through the reaches
 do iRch=1,nRch

  ! identify the index of the downstream segment
  ixDownRch = downIndex(iRch)
  if(ixDownRch == integerMissing) cycle

  ! increment the number of upstream elements in the downstream segment
  mUpstream(ixDownRch) = mUpstream(ixDownRch)+1
  if(mUpstream(ixDownRch)>nUpstream(ixDownRch))then
   message=trim(message)//'upstream index exceeds dimension'
   ierr=20; return
  endif

  ! populate the structure components
  NETOPO(ixDownRch)%UREACHI( mUpstream(ixDownRch) ) = NETOPO(iRch)%REACHIX
  NETOPO(ixDownRch)%UREACHK( mUpstream(ixDownRch) ) = NETOPO(iRch)%REACHID

 end do  ! looping through reaches

 ! set missing values to -1
 ! NOTE: check if the -1 is special and if the replacement with integer missing is necessary
 where(downIndex==integerMissing) downIndex=-1

 ! populate data structures
 NETOPO(:)%DREACHI = downIndex(:)

 end subroutine up2downSegment

 ! *********************************************************************
 ! new subroutine: define index of downstream reach
 ! *********************************************************************
 subroutine downReachIndex(&
                           ! input
                           nUp,          & ! number of upstream elements
                           nSeg,         & ! number of stream segments
                           segId,        & ! unique identifier of the stream segments
                           downId,       & ! unique identifier of the segment where water drains
                           ! output
                           downSegIndex, & ! index of downstream stream segment
                           nElement2Seg, & ! number of elements that drain into each segment
                           ierr,message)
 ! external modules
 USE nr_utility_module, ONLY: indexx  ! Num. Recipies utilities
 implicit none
 ! input variables
 integer(i4b), intent(in)        :: nUp             ! number of upstream elements
 integer(i4b), intent(in)        :: nSeg            ! number of stream segments
 integer(i4b), intent(in)        :: segId(:)        ! unique identifier of the stream segments
 integer(i4b), intent(in)        :: downId(:)       ! unique identifier of the segment where water drains
 ! output variables
 integer(i4b), intent(out)       :: downSegIndex(:) ! index of downstream stream segment
 integer(i4b), intent(out)       :: nElement2Seg(:) ! number of elements that drain into each segment
 integer(i4b), intent(out)       :: ierr            ! error code
 character(*), intent(out)       :: message         ! error message
 ! local variables
 integer(i4b)                    :: iUp                 ! index of upstream element
 integer(i4b)                    :: iSeg,jSeg           ! index of stream segment
 integer(i4b)                    :: rankSegId           ! ranked Id of the stream segment
 integer(i4b)                    :: rankDownId          ! ranked Id of the downstream stream segment
 integer(i4b)                    :: rankSeg(nSeg)       ! rank index of each segment in the nRch vector
 integer(i4b)                    :: rankDownSeg(nUp)    ! rank index of each downstream stream segment
 logical(lgt),parameter          :: checkMap=.true.     ! flag to check the mapping
 ! initialize error control
 ierr=0; message='downReachIndex/'

 ! initialize output
 nElement2Seg(:) = 0
 downSegIndex(:) = integerMissing

 ! rank the stream segments
 call indexx(segId, rankSeg)

 ! rank the drainage segment
 call indexx(downId, rankDownSeg)

 iSeg=1  ! second counter
 ! loop through the upstream elements
 do iUp=1,nUp

  ! get Ids for the up2seg mapping vector
  rankDownId = downId( rankDownSeg(iUp) )
  if (rankDownId == down2noSegment) cycle ! upstream element does not drain into any stream segment (closed basin or coastal HRU)

  ! keep going until found the index
  do jSeg=iSeg,nSeg ! normally a short loop

   ! get Id of the stream segment
   rankSegId = segId( rankSeg(jSeg) )
   !print*, 'iUp, iSeg, jSeg, rankDownId, rankSegId, rankDownSeg(iUp), rankSeg(jSeg) = ', &
   !         iUp, iSeg, jSeg, rankDownId, rankSegId, rankDownSeg(iUp), rankSeg(jSeg)

   ! define the index where we have a match
   if(rankDownId==rankSegId)then

    ! identify the index of the segment that the HRU drains into
    downSegIndex( rankDownSeg(iUp) ) = rankSeg(jSeg)
    nElement2Seg( rankSeg(jSeg)    ) = nElement2Seg( rankSeg(jSeg) ) + 1

    ! check if we should increment the stream segment
    ! NOTE: we can have multiple upstream elements draining into the same segment
    !        --> in this case, we want to keep the segment the same
    if(iUp<nUp .and. jSeg<nSeg)then
     if(downId( rankDownSeg(iUp+1) ) >= segId( rankSeg(jSeg+1) ) ) iSeg=jSeg+1
    endif

    ! identified the segment so exit the segment loop and evaluate the next upstream element
    exit

   endif  ! match between the upstream drainage segment and the stream segment
  end do  ! skipping segments that have no input

 end do  ! looping through upstream elements

 ! check
 if(checkMap)then
  do iUp=1,nUp
   if(downId(iUp) == down2noSegment .or. downSegIndex(iUp)==integerMissing) cycle
   if(downId(iUp) /= segId( downSegIndex(iUp) ) )then
    message=trim(message)//'problems identifying the index of the stream segment that a given HRU drains into'
    ierr=20; return
   endif
  end do
 endif

 end subroutine downReachIndex

end module read_streamSeg
