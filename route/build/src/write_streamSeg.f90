module write_streamSeg

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
USE globalData, only : meta_dims           ! dimension information
USE globalData, only : meta_HRU            ! HRU properties
USE globalData, only : meta_HRU2SEG        ! HRU-to-segment mapping
USE globalData, only : meta_SEG            ! stream segment properties
USE globalData, only : meta_NTOPO          ! network topology
USE globalData, only : meta_PFAF           ! network topology

! named variables
USE var_lookup,only:ixStruct, nStructures  ! index of data structures
USE var_lookup,only:ixDims,   nDimensions  ! index of dimensions
USE var_lookup,only:ixHRU,    nVarsHRU     ! index of variables for the HRUs
USE var_lookup,only:ixSEG,    nVarsSEG     ! index of variables for the stream segments
USE var_lookup,only:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
USE var_lookup,only:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology
USE var_lookup,only:ixPFAF,   nVarsPFAF    ! index of variables for the pfafstetter code

! netcdf modules
USE netcdf
USE io_netcdf, only: open_nc               ! open netcdf
USE io_netcdf, only: close_nc              ! close netcdf
USE io_netcdf, only: write_nc

! external utilities
USE nr_utility_module, ONLY: indexx  ! Num. Recipies utilities
USE nr_utility_module, ONLY: arth    ! Num. Recipies utilities

implicit none

! privacy
private
public::writeData
contains

 ! *********************************************************************
 ! new subroutine: write ancillary data for HRUs and stream segments
 ! *********************************************************************
 subroutine writeData(&
                      ! input
                      fname,         & ! input: file name
                      ! input: model control
                      tot_hru,       & ! input: total number of all the upstream hrus for all stream segments
                      tot_upseg,     & ! input: total number of immediate upstream segments for all  stream segments
                      tot_upstream,  & ! input: total number of all of the upstream stream segments for all stream segments
                      tot_uh,        & ! input: total number of all of unit hydrographs for all the segments
                      ! input: reach masks
                      ixHRU_desired, & ! input: indices of desired hrus
                      ixSeg_desired, & ! input: indices of desired reaches
                      ! input: data structures
                      structHRU,     & ! input: ancillary data for HRUs
                      structSeg,     & ! input: ancillary data for stream segments
                      structHRU2seg, & ! input: ancillary data for mapping hru2basin
                      structNTOPO,   & ! input: ancillary data for network toopology
                      structPFAF,    & ! input: ancillary data for pfafstetter code
                      ! output: error control
                      ierr,message)    ! output: error control
 implicit none
 ! input variables
 character(*)      , intent(in)      :: fname            ! filename
 ! input: model control
 integer(i4b)      , intent(in)      :: tot_hru          ! total number of all the upstream hrus for all stream segments
 integer(i4b)      , intent(in)      :: tot_upseg        ! total number of immediate upstream segments for all  stream segments
 integer(i4b)      , intent(in)      :: tot_upstream     ! total number of all of the upstream stream segments for all stream segments
 integer(i4b)      , intent(in)      :: tot_uh           ! total number of all of the unit hydrograph for all stream segments
 ! input: reach masks
 integer(i4b)      , intent(in)      :: ixHRU_desired(:) ! indices of desired hrus
 integer(i4b)      , intent(in)      :: ixSeg_desired(:) ! indices of desired reaches
 ! input: data structures
 type(var_dlength) , intent(in)      :: structHRU(:)     ! HRU properties
 type(var_dlength) , intent(in)      :: structSeg(:)     ! stream segment properties
 type(var_ilength) , intent(in)      :: structHRU2seg(:) ! HRU-to-segment mapping
 type(var_ilength) , intent(in)      :: structNTOPO(:)   ! network topology
 type(var_clength) , intent(in)      :: structPFAF(:)    ! network topology
 ! output: error control
 integer(i4b)      , intent(out)     :: ierr             ! error code
 character(*)      , intent(out)     :: message          ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                        :: nSpace           ! number of spatial elements
 integer(i4b)                        :: iStruct          ! structure index
 integer(i4b)                        :: iVar             ! variable index
 logical(lgt)                        :: dimCheck(nDimensions)  ! dimension used for output
 character(len=strLen)               :: cmessage         ! error message of downwind routine
 ! initialize error control
 ierr=0; message='writeData/'

 ! ---------- Find dimension needed for output --------------------------------------------------------------------
 dimCheck(:) = .true.
 ! check dimensions need for output (only subset mode)
 if (idSegOut>0) then
   dimCheck(:) = .false.
   do iVar=1,size(meta_HRU)
     if (meta_HRU(iVar)%varFile) dimCheck(meta_HRU(iVar)%vartype) = .true.
   end do
   do iVar=1,size(meta_HRU2SEG)
     if (meta_HRU2SEG(iVar)%varFile) dimCheck(meta_HRU2SEG(iVar)%vartype) = .true.
   end do
   do iVar=1,size(meta_SEG)
     if (meta_SEG(iVar)%varFile) dimCheck(meta_SEG(iVar)%vartype) = .true.
   end do
   do iVar=1,size(meta_NTOPO)
     if (meta_NTOPO(iVar)%varFile) dimCheck(meta_NTOPO(iVar)%vartype) = .true.
   end do
   do iVar=1,size(meta_PFAF)
     if (meta_PFAF(iVar)%varFile) &
        dimCheck(meta_PFAF(iVar)%vartype) = .true.
        dimCheck(ixDims%pfaf) = .true.
   end do
 endif

 ! define dimension lengths
 meta_dims(ixDims%hru  )%dimLength =  size(ixHRU_desired)   ! hru vector
 meta_dims(ixDims%seg  )%dimLength =  size(ixSeg_desired)   ! stream segment vector
 meta_dims(ixDims%upHRU)%dimLength =  tot_hru               ! upstream HRUs
 meta_dims(ixDims%upSeg)%dimLength =  tot_upseg             ! immediate upstream segments
 meta_dims(ixDims%upAll)%dimLength =  tot_upstream          ! all upstream segments
 meta_dims(ixDims%uh   )%dimLength =  tot_uh                ! all unit hydrograph
 meta_dims(ixDims%pfaf )%dimLength =  maxPfafLen            ! maximum pfaf code length

 ! create NetCDF file
 call createFile(trim(fname), dimCheck, ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ---------- write data to the NetCDF file  ---------------------------------------------------------------------

 ! loop through data structures
 do iStruct=1,nStructures

  ! get the length of the spatial dimension
  nSpace = meta_struct(iStruct)%nSpace

  ! write data from each structure
  select case(iStruct)
   case(ixStruct%HRU    ); call writeVar_dp  (trim(fname), nSpace, meta_HRU    , structHRU    , ixHRU_desired, ierr, cmessage)       ! HRU properties
   case(ixStruct%SEG    ); call writeVar_dp  (trim(fname), nSpace, meta_SEG    , structSeg    , ixSeg_desired, ierr, cmessage)       ! stream segment properties
   case(ixStruct%HRU2SEG); call writeVar_i4b (trim(fname), nSpace, meta_HRU2SEG, structHRU2seg, ixHRU_desired, ierr, cmessage)       ! HRU-to-segment mapping
   case(ixStruct%NTOPO  ); call writeVar_i4b (trim(fname), nSpace, meta_NTOPO  , structNTOPO  , ixSeg_desired, ierr, cmessage)       ! network topology
   case(ixStruct%PFAF   ); call writeVar_char(trim(fname), nSpace, meta_PFAF   , structPFAF   , ixSeg_desired, ierr, cmessage)       ! pfafstetter code
   case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end do  ! looping through data structures

 end subroutine writeData


 ! *********************************************************************
 ! new subroutine: create NetCDF output file
 ! *********************************************************************
 subroutine createFile(fname, dimCheck, ierr,message)
 implicit none
 ! dummy variables
 character(*)      , intent(in)      :: fname            ! filename
 logical(lgt)                        :: dimCheck(:)      ! dimension used for output
 integer(i4b)      , intent(out)     :: ierr             ! error code
 character(*)      , intent(out)     :: message          ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                        :: ncid             ! NetCDF file ID
 integer(i4b)                        :: jDim             ! dimension index
 integer(i4b)                        :: iStruct          ! structure index
 integer(i4b),parameter              :: nVars=30         ! number of variables
 character(len=strLen)               :: cmessage         ! error message of downwind routine
 ! initialize error control
 ierr=0; message='createFile/'

 ! ---------- create file ----------------------------------------------------------------------------------------

 ! create file
 ierr = nf90_create(trim(fname), NF90_64BIT_OFFSET, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! ---------- define dimensions ----------------------------------------------------------------------------------

 ! define dimensions
 do jDim=1,size(meta_dims)
  if (.not.dimCheck(jDim)) cycle
  ierr = nf90_def_dim(ncid, trim(meta_dims(jDim)%dimName), meta_dims(jDim)%dimLength, meta_dims(jDim)%dimId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name = '//trim(meta_dims(jDim)%dimName); return; endif
 end do  ! looping through dimensions

 ! ---------- define start index and count for the ragged arrays -------------------------------------------------

 ! loop through dimensions
 do jDim=1,size(meta_dims)

  ! HRUs and stream segments are not ragged, so cycle
  if(jDim==ixDims%hru .or. jDim==ixDims%seg .or. jDim==ixDims%pfaf) cycle
  ! if dimension is not used
  if (.not.dimCheck(jDim)) cycle

  ! define the start index
  call varDefine(ncid, trim(meta_dims(jDim)%dimName)//'_start', 'start index in ragged array', '-', nf90_int, ixDims%seg, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! define the count
  call varDefine(ncid, trim(meta_dims(jDim)%dimName)//'_count', 'count of spatial element in ragged array', '-', nf90_int, ixDims%seg, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end do  ! looping through dimensions

 ! ---------- define variables -----------------------------------------------------------------------------------

 ! loop through data structures
 do iStruct=1,nStructures

  ! process metadata from each structure
  select case(iStruct)
   case(ixStruct%HRU    ); call defineVar(ncid, meta_HRU    , nf90_double, ierr, cmessage)       ! HRU properties
   case(ixStruct%SEG    ); call defineVar(ncid, meta_SEG    , nf90_double, ierr, cmessage)       ! stream segment properties
   case(ixStruct%HRU2SEG); call defineVar(ncid, meta_HRU2SEG, nf90_int,    ierr, cmessage)       ! HRU-to-segment mapping
   case(ixStruct%NTOPO  ); call defineVar(ncid, meta_NTOPO  , nf90_int,    ierr, cmessage)       ! network topology
   case(ixStruct%PFAF   ); call defineVar(ncid, meta_PFAF   , nf90_char,   ierr, cmessage)       ! pfafstetter code
   case default; ierr=20; message=trim(message)//'unable to identify data structure'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end do  ! looping through data structures

 ! ---------- complete file definition ---------------------------------------------------------------------------

 ! exit define mode
 ierr = nf90_enddef(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close netCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine createFile

 ! *********************************************************************
 ! new subroutine: define variable in NetCDF file
 ! *********************************************************************
 subroutine defineVar(&
                      ! input
                      ncid,         & ! input: NetCDF id
                      meta,         & ! input: metadata structure
                      ivtype,       & ! input: variable type
                      ! output: error control
                      ierr,message)   ! output: error control
 implicit none
 ! input variables
 integer(i4b)      , intent(in)      :: ncid      ! netcdf id
 type(var_info)    , intent(in)      :: meta(:)   ! metadata structure
 integer(i4b)      , intent(in)      :: ivtype    ! variable type
 ! output: error control
 integer(i4b)      , intent(out)     :: ierr      ! error code
 character(*)      , intent(out)     :: message   ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 integer(i4b)                        :: iVar      ! variable index
 character(len=strLen)               :: cmessage  ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------------------------
 ierr=0; message='defineVar/'

 ! loop through variables
 do iVar=1,size(meta)
   ! if this is a subset mode (idSegOut>0) then only write variables where meta(iVar)%varFile = .true.
   if(idSegOut>0 .and. .not.meta(iVar)%varFile) cycle
   call varDefine(ncid, & ! NetCDF ID
                  trim(meta(iVar)%varName), trim(meta(iVar)%varDesc), trim(meta(iVar)%varUnit),&  ! name, description, units
                  ivtype, meta(iVar)%varType, ierr, cmessage)                                     ! variable type, dimension index
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do

 end subroutine defineVar

 ! *********************************************************************
 ! new subroutine: define variable in NetCDF file
 ! *********************************************************************
 subroutine varDefine(&
                      ! input
                      ncid,         & ! input: NetCDF id
                      varName,      & ! input: variable name
                      varDesc,      & ! input: variable description
                      varUnit,      & ! input: variable units
                      ivtype,       & ! input: variable type
                      ivdim,        & ! input: variable dimension
                      ! output: error control
                      ierr,message)   ! output: error control
 implicit none
 ! input variables
 integer(i4b)      , intent(in)      :: ncid      ! netcdf id
 character(*)      , intent(in)      :: varName   ! variable name
 character(*)      , intent(in)      :: varDesc   ! variable description
 character(*)      , intent(in)      :: varUnit   ! variable units
 integer(i4b)      , intent(in)      :: ivtype    ! variable type
 integer(i4b)      , intent(in)      :: ivdim     ! variable dimension
 ! output: error control
 integer(i4b)      , intent(out)     :: ierr      ! error code
 character(*)      , intent(out)     :: message   ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 integer(i4b)                        :: iVarID    ! variable ID
 ! ---------------------------------------------------------------------------------------------------------------
 ierr=0; message='varDefine/'

 ! define variable
 if (ivtype /= 2) then ! character array
  ierr = nf90_def_var(ncid,trim(varName), ivtype, meta_dims(ivdim)%dimId, iVarId)
 else
  ierr = nf90_def_var(ncid,trim(varName), ivtype, (/meta_dims(ixDims%pfaf)%dimId, meta_dims(ivdim)%dimId/), iVarId)
 end if
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'name = '//trim(varName); return; endif

 ! add variable description
 ierr = nf90_put_att(ncid,iVarId,'long_name',trim(varDesc))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'name = '//trim(varName); return; endif

 ! add variable units
 ierr = nf90_put_att(ncid,iVarId,'units',trim(varUnit))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'name = '//trim(varName); return; endif

 ! add missing values
 select case(ivtype)
  case(nf90_int);    ierr = nf90_put_att(ncid,iVarId,'_FillValue',integerMissing)
  case(nf90_double); ierr = nf90_put_att(ncid,iVarId,'_FillValue',realMissing)
  !case(nf90_char);   ierr = nf90_put_att(ncid,iVarId,'_FillValue',characterMissing)
 endselect
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'name = '//trim(varName); return; endif

 end subroutine varDefine

 ! *********************************************************************
 ! new subroutine: write variables to NetCDF file (double precision)
 ! *********************************************************************
 subroutine writeVar_dp(&
                        ! input
                        fname,        & ! input: filename
                        nSpace,       & ! input: number of spatial elements
                        meta,         & ! input: metadata structure
                        struct,       & ! input: variable type
                        ixDesired,    & ! input: vector of desired spatial elements
                        ! output: error control
                        ierr,message)   ! output: error control
 implicit none
 ! input variables
 character(*)      , intent(in)      :: fname            ! filename
 integer(i4b)      , intent(in)      :: nSpace           ! number of spatial elements
 type(var_info)    , intent(in)      :: meta(:)          ! metadata structure
 type(var_dlength) , intent(in)      :: struct(:)        ! data structure
 integer(i4b)      , intent(in)      :: ixDesired(:)     ! vector of desired spatial elements
 ! output: error control
 integer(i4b)      , intent(out)     :: ierr             ! error code
 character(*)      , intent(out)     :: message          ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 integer(i4b)                        :: ncid             ! netCDF ID
 integer(i4b)                        :: ix,jx,nx         ! write indices
 integer(i4b)                        :: iVar             ! variable index
 integer(i4b)                        :: jDim             ! dimension index
 integer(i4b)                        :: iSpace           ! space index
 integer(i4b)                        :: jSpace           ! space index
 integer(i4b)                        :: nSubset          ! spatial subset
 integer(i4b)                        :: dimLength        ! dimension length
 integer(i4b)                        :: ixStart(nSpace)  ! start index in the ragged arrays
 integer(i4b)                        :: ixCount(nSpace)  ! count index in the ragged arrays
 logical(lgt)                        :: isRaggedArray    ! logical flag to define a ragged array
 character(len=strLen)               :: cmessage         ! error message of downwind routine
 real(dp)          , allocatable     :: tempVec(:)       ! temporary vector
 ! ---------------------------------------------------------------------------------------------------------------
 ierr=0; message='writeVar_dp/'

 ! get the size of the spatial subset
 nSubset = size(ixDesired)

 ! initial allocation of the temporary vectors
 allocate(tempVec(nSubset), stat=ierr)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating temporary vector'; return; endif

 call open_nc(fname, 'w', ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! loop through variables
 do iVar=1,size(meta)

  ! ---------- get the data vector for a given variable -----------------------------------------------------------

  ! if we created a subset (idSegOut>0) then only write variables where meta(iVar)%varFile = .true.
  if(idSegOut>0 .and. .not.meta(iVar)%varFile) cycle

  ! print progress
  print*, 'Writing '//trim(meta(iVar)%varName)//' to file '//trim(fname)

  ! save the dimension length
  jDim      = meta(iVar)%varType        ! dimension index
  dimLength = meta_dims(jDim)%dimLength ! dimension length
  if(dimLength==0) cycle

  ! reallocate the temporary vector
  if(size(tempVec)/=dimLength)then
   deallocate(tempVec,         stat=ierr); if(ierr/=0)then; message=trim(message)//'problem deallocating tempVec'; return; endif
   allocate(tempVec(dimLength),stat=ierr); if(ierr/=0)then; message=trim(message)//'problem allocating tempVec'; return; endif
  endif

  ! check if the variable should be output in a ragged array
  isRaggedArray = (jDim/=ixDims%hru .and. jdim/=ixDims%seg)

  ! ragged array
  ix = 1
  do jSpace=1,nSubset

   ! get the spatial index
   iSpace = ixDesired(jSpace)

   ! ***** case 1: ragged array
   if(isRaggedArray)then  ! check the need to create the ragged array

    ! get count and end index
    nx = size(struct(iSpace)%var(iVar)%dat)
    jx = ix + nx -1

    ! write vector
    tempVec(ix:jx) = struct(iSpace)%var(iVar)%dat(1:nx)

    ! update indices
    ixStart(jSpace) = ix
    ixCount(jSpace) = nx
    ix = ix + nx

   ! ***** case 2: regular array
   else
    tempVec(jSpace) = struct(iSpace)%var(iVar)%dat(1)
   endif

  end do  ! looping through space

  ! ---------- write the data vector to the NetCDF file -----------------------------------------------------------

  ! write data
  call write_nc(ncid, trim(meta(iVar)%varName), tempVec, (/1/), (/dimLength/), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write ragged array
  ! NOTE: only need to do once
  if(isRaggedArray)then

   ! write start index
   call write_nc(ncid, trim(meta_dims(jDim)%dimName)//'_start', ixStart, (/1/), (/nSubset/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! write count
   call write_nc(ncid, trim(meta_dims(jDim)%dimName)//'_count', ixCount, (/1/), (/nSubset/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif   ! if a ragged array

 end do  ! looping through variables

 call close_nc(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ---------- clean up ------------------------------------------------------------------------------------------

 ! deallocate space for the ragged arrays
 deallocate(tempVec, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for the temporary vector'; return; endif

 end subroutine writeVar_dp

 ! *********************************************************************
 ! new subroutine: write variables to NetCDF file (integer)
 ! *********************************************************************
 subroutine writeVar_i4b(&
                         ! input
                         fname,        & ! input: filename
                         nSpace,       & ! input: number of spatial elements
                         meta,         & ! input: metadata structure
                         struct,       & ! input: variable type
                         ixDesired,    & ! input: vector of desired spatial elements
                         ! output: error control
                         ierr,message)   ! output: error control
 implicit none
 ! input variables
 character(*)      , intent(in)      :: fname            ! filename
 integer(i4b)      , intent(in)      :: nSpace           ! number of spatial elements
 type(var_info)    , intent(in)      :: meta(:)          ! metadata structure
 type(var_ilength) , intent(in)      :: struct(:)        ! data structure
 integer(i4b)      , intent(in)      :: ixDesired(:)     ! vector of desired spatial elements
 ! output: error control
 integer(i4b)      , intent(out)     :: ierr             ! error code
 character(*)      , intent(out)     :: message          ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 integer(i4b)                        :: ncid             ! netCDF ID
 integer(i4b)                        :: ix,jx,nx         ! write indices
 integer(i4b)                        :: iVar             ! variable index
 integer(i4b)                        :: jDim             ! dimension index
 integer(i4b)                        :: iSpace           ! space index
 integer(i4b)                        :: jSpace           ! space index
 integer(i4b)                        :: nSubset          ! spatial subset
 integer(i4b)                        :: dimLength        ! dimension length
 integer(i4b)                        :: ixStart(nSpace)  ! start index in the ragged arrays
 integer(i4b)                        :: ixCount(nSpace)  ! count index in the ragged arrays
 logical(lgt)                        :: isRaggedArray    ! logical flag to define a ragged array
 character(len=strLen)               :: cmessage         ! error message of downwind routine
 integer(i4b)      , allocatable     :: tempVec(:)       ! temporary vector
 ! ---------------------------------------------------------------------------------------------------------------
 ierr=0; message='writeVar_i4b/'

 ! get the size of the spatial subset
 nSubset = size(ixDesired)

 ! initial allocation of the temporary vectors
 allocate(tempVec(nSubset), stat=ierr)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating temporary vector'; return; endif

 call open_nc(fname, 'w', ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! loop through variables
 do iVar=1,size(meta)

  ! ---------- get the data vector for a given variable -----------------------------------------------------------

  ! if we created a subset (idSegOut>0) then only write variables where meta(iVar)%varFile = .true.
  if(idSegOut>0 .and. .not.meta(iVar)%varFile) cycle

  ! print progress
  print*, 'Writing '//trim(meta(iVar)%varName)//' to file '//trim(fname)

  ! save the dimension length
  jDim      = meta(iVar)%varType        ! dimension index
  dimLength = meta_dims(jDim)%dimLength ! dimension length
  if(dimLength==0) cycle

  ! reallocate the temporary vector
  if(size(tempVec)/=dimLength)then
   deallocate(tempVec,         stat=ierr); if(ierr/=0)then; message=trim(message)//'problem deallocating tempVec'; return; endif
   allocate(tempVec(dimLength),stat=ierr); if(ierr/=0)then; message=trim(message)//'problem allocating tempVec'; return; endif
  endif

  ! check if the variable should be output in a ragged array
  isRaggedArray = (jDim/=ixDims%hru .and. jdim/=ixDims%seg)

  ! ragged array
  ix = 1
  do jSpace=1,nSubset

   ! get the spatial index
   iSpace = ixDesired(jSpace)

   ! ***** case 1: ragged array
   if(isRaggedArray)then  ! check the need to create the ragged array

    ! get count and end index
    nx = size(struct(iSpace)%var(iVar)%dat)
    jx = ix + nx -1

    ! write vector
    tempVec(ix:jx) = struct(iSpace)%var(iVar)%dat(1:nx)

    ! update indices
    ixStart(jSpace) = ix
    ixCount(jSpace) = nx
    ix = ix + nx

   ! ***** case 2: regular array
   else
    tempVec(jSpace) = struct(iSpace)%var(iVar)%dat(1)
   endif

  end do  ! looping through space

  ! ---------- write the data vector to the NetCDF file -----------------------------------------------------------

  ! write data
  call write_nc(ncid, trim(meta(iVar)%varName), tempVec, (/1/), (/dimLength/), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write ragged array
  ! NOTE: only need to do once
  if(isRaggedArray)then

   ! write start index
   call write_nc(ncid, trim(meta_dims(jDim)%dimName)//'_start', ixStart, (/1/), (/nSubset/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! write count
   call write_nc(ncid, trim(meta_dims(jDim)%dimName)//'_count', ixCount, (/1/), (/nSubset/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif   ! if a ragged array

 end do  ! looping through variables

 call close_nc(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ---------- clean up ------------------------------------------------------------------------------------------

 ! deallocate space for the ragged arrays
 deallocate(tempVec, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for the temporary vector'; return; endif

 end subroutine writeVar_i4b


 ! *********************************************************************
 ! new subroutine: write variables to NetCDF file (character array)
 ! *********************************************************************
 subroutine writeVar_char(&
                         ! input
                         fname,        & ! input: filename
                         nSpace,       & ! input: number of spatial elements
                         meta,         & ! input: metadata structure
                         struct,       & ! input: variable type
                         ixDesired,    & ! input: vector of desired spatial elements
                         ! output: error control
                         ierr,message)   ! output: error control
 implicit none
 ! input variables
 character(*)      , intent(in)         :: fname            ! filename
 integer(i4b)      , intent(in)         :: nSpace           ! number of spatial elements
 type(var_info)    , intent(in)         :: meta(:)          ! metadata structure
 type(var_clength) , intent(in)         :: struct(:)        ! data structure
 integer(i4b)      , intent(in)         :: ixDesired(:)     ! vector of desired spatial elements
 ! output: error control
 integer(i4b)      , intent(out)        :: ierr             ! error code
 character(*)      , intent(out)        :: message          ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 integer(i4b)                           :: ncid             ! netCDF ID
 integer(i4b)                           :: ix,jx,nx         ! write indices
 integer(i4b)                           :: iVar             ! variable index
 integer(i4b)                           :: jDim             ! dimension index
 integer(i4b)                           :: iSpace           ! space index
 integer(i4b)                           :: jSpace           ! space index
 integer(i4b)                           :: nSubset          ! spatial subset
 integer(i4b)                           :: dimLength        ! dimension length
 integer(i4b)                           :: ixStart(nSpace)  ! start index in the ragged arrays
 integer(i4b)                           :: ixCount(nSpace)  ! count index in the ragged arrays
 logical(lgt)                           :: isRaggedArray    ! logical flag to define a ragged array
 character(len=strLen)                  :: cmessage         ! error message of downwind routine
 character(len=maxPfafLen), allocatable :: tempVec(:)       ! temporary vector
 ! ---------------------------------------------------------------------------------------------------------------
 ierr=0; message='writeVar_char/'

 ! get the size of the spatial subset
 nSubset = size(ixDesired)

 ! initial allocation of the temporary vectors
 allocate(tempVec(nSubset), stat=ierr)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating temporary vector'; return; endif

 call open_nc(fname, 'w', ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! loop through variables
 do ivar=1,size(meta)

  ! ---------- get the data vector for a given variable -----------------------------------------------------------

  ! if we created a subset (idSegOut>0) then only write variables where meta(ivar)%varFile = .true.
  if(idSegOut>0 .and. .not.meta(ivar)%varFile) cycle

  ! print progress
  print*, 'Writing '//trim(meta(ivar)%varName)//' to file '//trim(fname)

  ! save the dimension length
  jDim      = meta(ivar)%varType        ! dimension index
  dimLength = meta_dims(jDim)%dimLength ! dimension length
  if(dimLength==0) cycle

  ! reallocate the temporary vector
  if(size(tempVec)/=dimLength)then
   deallocate(tempVec,         stat=ierr); if(ierr/=0)then; message=trim(message)//'problem deallocating tempVec'; return; endif
   allocate(tempVec(dimLength),stat=ierr); if(ierr/=0)then; message=trim(message)//'problem allocating tempVec'; return; endif
  endif

  ! check if the variable should be output in a ragged array
  isRaggedArray = (jDim/=ixDims%hru .and. jdim/=ixDims%seg)

  ! ragged array
  ix = 1
  do jSpace=1,nSubset

   ! get the spatial index
   iSpace = ixDesired(jSpace)

   ! ***** case 1: ragged array
   if(isRaggedArray)then  ! check the need to create the ragged array

    ! get count and end index
    nx = size(struct(iSpace)%var(ivar)%dat)
    jx = ix + nx -1

    ! write vector
    tempVec(ix:jx) = struct(iSpace)%var(ivar)%dat(1:nx)

    ! update indices
    ixStart(jSpace) = ix
    ixCount(jSpace) = nx
    ix = ix + nx

   ! ***** case 2: regular array
   else
    tempVec(jSpace) = struct(iSpace)%var(ivar)%dat(1)
   endif

  end do  ! looping through space

  ! ---------- write the data vector to the NetCDF file -----------------------------------------------------------

  ! write data
  call write_nc(ncid, trim(meta(ivar)%varName), tempVec, (/1,1/), (/maxPfafLen, dimLength/), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write ragged array
  ! NOTE: only need to do once
  if(isRaggedArray)then

   ! write start index
   call write_nc(ncid, trim(meta_dims(jDim)%dimName)//'_start', ixStart, (/1,1/), (/maxPfafLen, nSubset/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! write count
   call write_nc(ncid, trim(meta_dims(jDim)%dimName)//'_count', ixCount, (/1,1/), (/maxPfafLen, nSubset/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif   ! if a ragged array

 end do  ! looping through variables

 call close_nc(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ---------- clean up ------------------------------------------------------------------------------------------

 ! deallocate space for the ragged arrays
 deallocate(tempVec, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for the temporary vector'; return; endif

 end subroutine writeVar_char

end module write_streamSeg
