module read_streamSeg

! data types
USE nrtype,    only : i4b,dp,lgt
USE nrtype,    only : strLen               ! string length
USE dataTypes, only : var_ilength          ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength          ! double precision type: var(:)%dat
USE dataTypes, only : var_clength          ! character type:        var(:)%dat
USE dataTypes, only : var_info             ! metadata

! global data
USE public_var

! metadata on data structures
USE globalData, only : meta_struct         ! structure information
USE globalData, only : meta_HRU            ! HRU properties
USE globalData, only : meta_HRU2SEG        ! HRU-to-segment mapping
USE globalData, only : meta_SEG            ! stream segment properties
USE globalData, only : meta_NTOPO          ! network topology
USE globalData, only : meta_PFAF           ! network topology

! named variables
USE var_lookup,only:ixStruct, nStructures  ! index of data structures
USE var_lookup,only:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup,only:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup,only:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup,only:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology
USE var_lookup,only:ixPFAF,   nVarsPFAF    ! index of variables for the pfafstetter code

! netcdf modules
USE netcdf

! external utilities
USE nr_utility_module, ONLY: arth    ! Num. Recipies utilities

USE alloc_data,        ONLY: alloc_struct

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
 character(*)      , intent(in)               :: fname            ! filename
 character(*)      , intent(in)               :: dname_nhru       ! dimension name for HRUs
 character(*)      , intent(in)               :: dname_sseg       ! dimension name for stream segments
 ! output: model control
 integer(i4b)      , intent(out)              :: nHRU_in          ! number of HRUs
 integer(i4b)      , intent(out)              :: nRch_in          ! number of stream segments
 ! output: data structures
 type(var_dlength) , intent(out), allocatable :: structHRU(:)     ! HRU properties
 type(var_dlength) , intent(out), allocatable :: structSeg(:)     ! stream segment properties
 type(var_ilength) , intent(out), allocatable :: structHRU2seg(:) ! HRU-to-segment mapping
 type(var_ilength) , intent(out), allocatable :: structNTOPO(:)   ! network topology
 type(var_clength) , intent(out), allocatable :: structPFAF(:)    ! network topology
 ! output: error control
 integer(i4b)      , intent(out)              :: ierr             ! error code
 character(*)      , intent(out)              :: message          ! error message
 ! ==========================================================================================================
 ! local variables
 integer(i4b)                           :: iStruct      ! structure index
 integer(i4b)                           :: iSpace       ! spatial index
 integer(i4b)                           :: iVar         ! variable index
 integer(i4b)                           :: ncid         ! NetCDF file ID
 integer(i4b)                           :: idimID_nHRU  ! dimension ID for HRUs
 integer(i4b)                           :: idimID_sseg  ! dimension ID for stream segments
 integer(i4b)                           :: iVarID       ! variable ID
 integer(i4b)                           :: jxStart      ! Start index for a given reach
 integer(i4b)                           :: jxCount      ! Number of elements for a given reach
 integer(i4b),              allocatable :: ixStart(:)   ! Start index for each reach
 integer(i4b),              allocatable :: ixCount(:)   ! Number of elements in each reach
 integer(i4b),              allocatable :: iTemp(:)     ! temporary integer vector
 real(dp),                  allocatable :: dTemp(:)     ! temporary double precision vector
 character(Len=maxPfafLen), allocatable :: cTemp(:)     ! temporary charactervector
 integer(i4b)                           :: dimLength    ! dimension length
 logical(lgt)                           :: isVarDesired ! .true. if the variable is desired
 character(len=strLen)                  :: varName      ! variable name
 character(len=strLen)                  :: cmessage     ! error message of downwind routine
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
 ! set flags if we want to turn on abstraction/injection option (require Qtake in network data)
 if(qtakeOption)then
  meta_SEG(ixSEG%Qtake)%varFile = .true.
 endif

 ! set flags if we want to read hdraulic geometry from file
 if(hydGeometryOption==readFromFile)then
  meta_SEG(ixSEG%width)%varFile = .true.
  meta_SEG(ixSEG%man_n)%varFile = .true.
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
 ! initialize error control
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
