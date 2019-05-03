module read_streamSeg
USE nrtype
USE netcdf
implicit none
private
public::getData
public::hru2basin
public::assign_reachparam
contains

 ! *********************************************************************
 ! new subroutine: get ancillary data for HRUs and stream segments
 ! *********************************************************************
 subroutine getData(fname,      &   ! input: file name
                    dname_nhru, &   ! input: dimension name for HRUs
                    dname_sseg, &   ! input: dimension name for stream segments
                    nhru_acil,  &   ! input-output: ancillary data for HRUs
                    sseg_acil,  &   ! input-output: ancillary data for stream segments
                    imap_acil,  &   ! input-output: ancillary data for mapping hru2basin
                    ntop_acil,  &   ! input-output: ancillary data for network topology
                    nHRU,       &   ! output: number of HRUs
                    nSeg,       &   ! output: number of stream segments
                    ierr, message)  ! output: error control
 USE dataTypes,only:namepvar,nameivar    ! provide access to data types
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: dname_nhru   ! dimension name for HRUs
 character(*), intent(in)        :: dname_sseg   ! dimension name for stream segments
 ! input-output
 type(namepvar), intent(inout)   :: nhru_acil(:) ! ancillary data for the HRUs
 type(namepvar), intent(inout)   :: sseg_acil(:) ! ancillary data for the stream segments
 type(nameivar), intent(inout)   :: imap_acil(:) ! ancillary data for the hru2basin mapping
 type(nameivar), intent(inout)   :: ntop_acil(:) ! ancillary data for the network topology
 ! output variables
 integer(i4b), intent(out)       :: nHRU         ! number of HRUs
 integer(i4b), intent(out)       :: nSeg         ! number of stream segments
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: iVar         ! variable index
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: idimID_nHRU  ! dimension ID for HRUs
 integer(i4b)                    :: idimID_sseg  ! dimension ID for stream segments
 integer(i4b)                    :: iVarID       ! variable ID
 ! initialize error control
 ierr=0; message='getData/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

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

 ! ** read in HRU variables
 do iVar=1,size(nhru_acil)

  ! allocate space for the structure component
  allocate(nhru_acil(ivar)%varData(nHRU),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for nhru data structure'; return; endif

  ! get the variable ID
  ierr = nf90_inq_varid(ncid, trim(nhru_acil(ivar)%varName), ivarID)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(nhru_acil(ivar)%varName); return; endif

  ! get the data
  ierr = nf90_get_var(ncid, ivarID, nhru_acil(ivar)%varData)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end do  ! (looping through variables)

 ! ** read in stream segment variables
 do iVar=1,size(sseg_acil)

  ! allocate space for the structure component
  allocate(sseg_acil(ivar)%varData(nSEG),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for sseg data structure'; return; endif

  ! get the variable ID
  ierr = nf90_inq_varid(ncid, trim(sseg_acil(ivar)%varName), ivarID)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(sseg_acil(ivar)%varName); return; endif

  ! get the data
  ierr = nf90_get_var(ncid, ivarID, sseg_acil(ivar)%varData)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end do  ! (looping through variables)

 ! ** read in hru2basin mapping
 do iVar=1,size(imap_acil)

  ! allocate space for the structure component
  allocate(imap_acil(ivar)%varData(nHRU),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for imap data structure'; return; endif

  ! get the variable ID
  ierr = nf90_inq_varid(ncid, trim(imap_acil(ivar)%varName), ivarID)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(imap_acil(ivar)%varName); return; endif

  ! get the data
  ierr = nf90_get_var(ncid, ivarID, imap_acil(ivar)%varData)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end do  ! (looping through variables)

 ! ** read in network topology
 do iVar=1,size(ntop_acil)

  ! allocate space for the structure component
  allocate(ntop_acil(ivar)%varData(nSEG),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for ntop data structure'; return; endif

  ! get the variable ID
  ierr = nf90_inq_varid(ncid, trim(ntop_acil(ivar)%varName), ivarID)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(ntop_acil(ivar)%varName); return; endif

  ! get the data
  ierr = nf90_get_var(ncid, ivarID, ntop_acil(ivar)%varData)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end do  ! (looping through variables)

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine getData


 ! *********************************************************************
 ! new subroutine: compute correspondence between HRUs and basins
 ! *********************************************************************
 subroutine hru2basin(nHRU,       &   ! input: number of HRUs
                      nSeg,       &   ! input: number of stream segments
                      nhru_acil,  &   ! input: ancillary data for HRUs
                      imap_acil,  &   ! input: ancillary data for mapping hru2basin
                      ntop_acil,  &   ! input: ancillary data for network topology
                      ierr, message)  ! output: error control
 USE nrutil,only:arth                 ! use to build vectors with regular increments
 USE dataTypes,only:namepvar,nameivar ! provide access to data types
 USE var_lookup,only:ixHRU,nVarsHRU   ! index of variables for the HRUs
 USE var_lookup,only:ixSEG,nVarsSEG   ! index of variables for the stream segments
 USE var_lookup,only:ixMAP,nVarsMAP   ! index of variables for the hru2basin mapping
 USE var_lookup,only:ixTOP,nVarsTOP   ! index of variables for the network topology
 USE nhru2basin  ! data structures holding the nhru2basin correspondence
 implicit none
 ! input variables
 integer(i4b), intent(in)        :: nHRU         ! number of HRUs
 integer(i4b), intent(in)        :: nSeg         ! number of stream segments
 type(namepvar), intent(in)      :: nhru_acil(:) ! ancillary data for the HRUs
 type(nameivar), intent(in)      :: imap_acil(:) ! ancillary data for the hru2basin mapping
 type(nameivar), intent(in)      :: ntop_acil(:) ! ancillary data for the network topology
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: iSeg         ! index of stream segment
 integer(i4b)                    :: jSeg         ! id of stream segment
 integer(i4b)                    :: HRUix(nHRU)  ! index of HRU in the nHRU vector
 integer(i4b)                    :: iDrain(nHRU) ! id of segment where HRU drains
 logical(lgt)                    :: hFlag(nHRU)  ! flag to denote if HRU drains into the current segment
 integer(i4b)                    :: nDrain       ! number of HRUs that drain into the current segment
 real(dp)                        :: totarea      ! total area of all HRUs feeding into a given stream segment (m2)

 ! initialize error control
 ierr=0; message='hru2basin/'

 ! allocate space for the hru2basin mapping
 allocate(h2b(nSeg),stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for h2b structure'; return; endif

 ! define indices for each HRU
 HRUix = arth(1,1,nHRU)

 ! define id of segment where HRU drains
 iDrain(:) = imap_acil(ixMAP%segHRUid)%varData(:)

 ! loop through stream segments
 do iSeg=1,nSeg

  ! define segment id
  jSeg = ntop_acil(ixTOP%segid)%varData(iSeg)

  ! identify the HRUs that drain into the current segment
  hFlag=.false. ! initialize flag to denote if HRU drains into the current segment
  where(iDrain == jSeg) hFlag=.true.

  ! get number of HRUs that drain into the current segment
  nDrain = count(hFlag)

  ! allocate space for the structure component
  allocate(h2b(iSeg)%cHRU(nDrain), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for h2b structure component'; return; endif

  ! populate structure components
  h2b(iSeg)%cHRU(:)%hru_ix   = pack(HRUix,hFlag)
  h2b(iSeg)%cHRU(:)%hru_id   = pack(imap_acil(ixMAP%HRUid)%varData,hFlag)
  h2b(iSeg)%cHRU(:)%hru_lon  = pack(nhru_acil(ixHRU%xLon )%varData,hFlag)
  h2b(iSeg)%cHRU(:)%hru_lat  = pack(nhru_acil(ixHRU%yLat )%varData,hFlag)
  h2b(iSeg)%cHRU(:)%hru_elev = pack(nhru_acil(ixHRU%elev )%varData,hFlag)
  h2b(iSeg)%cHRU(:)%hru_area = pack(nhru_acil(ixHRU%area )%varData,hFlag)

  ! compute total area of the HRUs draining to the stream segment
  totarea = sum(h2b(iSeg)%cHRU(:)%hru_area)

  ! compute the weights
  h2b(iSeg)%cHRU(:)%wght = h2b(iSeg)%cHRU(:)%hru_area / totarea

  ! check
  !print*, 'jSeg, nDrain, totarea = ', jSeg, nDrain, totarea

 end do  ! (looping thru stream segments)

 end subroutine hru2basin


 ! *********************************************************************
 ! new subroutine: assign data to the reach parameter structure
 ! *********************************************************************
 subroutine assign_reachparam(nRch,         &    ! input: number of stream segments
                              sseg_acil,    &    ! input: stream segment parameters
                              ntop_acil,    &    ! input: network topology
                              ierr, message)     ! output (error control)
 USE dataTypes,only:namepvar,nameivar ! provide access to data types
 USE var_lookup,only:ixHRU,nVarsHRU   ! index of variables for the HRUs
 USE var_lookup,only:ixSEG,nVarsSEG   ! index of variables for the stream segments
 USE var_lookup,only:ixMAP,nVarsMAP   ! index of variables for the hru2basin mapping
 USE var_lookup,only:ixTOP,nVarsTOP   ! index of variables for the network topology
 USE nhru2basin                       ! data structures holding the nhru2basin correspondence
 USE reachparam                       ! reach parameter structure
 implicit none
 ! input variables
 integer(i4b), intent(in)        :: nRch         ! number of reaches
 type(namepvar), intent(in)      :: sseg_acil(:) ! ancillary data for stream segments
 type(nameivar), intent(in)      :: ntop_acil(:) ! ancillary data for the network topology
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: iRch         ! reach index
 real(dp),parameter              :: min_slope=1.e-6_dp  ! minimum slope
 integer(i4b),parameter          :: strLen=256   ! length of character string
 character(len=strLen)           :: cmessage     ! error message of downwind routine

 ! initialize error control
 ierr=0; message='assign_reachparam/'

 ! allocate space for the reach parameter structure
 allocate(RPARAM(nRch),NETOPO(nRch),stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for reach parameter structures'; return; endif

 ! transfer information to the network topology structures
 NETOPO(:)%REACHID = ntop_acil(ixTOP%segid    )%varData(:)
 NETOPO(:)%DREACHK = ntop_acil(ixTOP%toSegment)%varData(:)

 ! transfer information to the reach structures
 RPARAM(:)%UPSAREA = sseg_acil(ixSEG%upArea   )%varData(:)
 RPARAM(:)%RLENGTH = sseg_acil(ixSEG%length   )%varData(:)
 RPARAM(:)%R_SLOPE = sseg_acil(ixSEG%slope    )%varData(:)

 ! compute area draining to each stream segment
 do iRch=1,nRch
  RPARAM(iRch)%BASAREA = sum(h2b(iRch)%cHRU(:)%hru_area)
 end do

 ! compute total area above the bottom of each reach (m2)
 RPARAM(:)%TOTAREA = RPARAM(:)%UPSAREA + RPARAM(:)%BASAREA

 ! define additional aspects of the network topology
 ! (populates reachparam)
 call define_topology(nrch, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ensure that slope exceeds minimum slope
 do iRch=1,nRch
  if(RPARAM(iRch)%R_SLOPE < min_slope) RPARAM(iRch)%R_SLOPE = min_slope
 end do

 !! specify some additional parameters (temporary "fix")
 !RPARAM(:)%R_WIDTH =  widthScale * sqrt(RPARAM(:)%TOTAREA)    ! channel width (m)
 !!RPARAM(:)%R_MAN_N =  0.05_dp  ! Manning's "n" paramater (unitless)
 !RPARAM(:)%R_MAN_N =  0.01_dp  ! Manning's "n" paramater (unitless)

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

 end subroutine assign_reachparam


 ! *********************************************************************
 ! new subroutine: define additional network topology
 ! *********************************************************************
 subroutine define_topology(nRch, &            ! input
                            ierr, message)     ! output (error control)
 USE nrutil, ONLY: arth       ! Num. Recipies utilities
 USE reachparam               ! reach parameter structure
 implicit none
 ! input variables
 integer(i4b), intent(in)        :: nRch         ! number of reaches
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: iRch,jRch    ! reach indices
 integer(i4b)                    :: iUps         ! index of upstream reach
 integer(i4b)                    :: nUps         ! number of upstream reaches
 logical(lgt)                    :: uFlag(nRch)  ! flag for upstream reaches
 logical(lgt)                    :: hFlag(nRch)  ! flag for headwater reaches
 real(dp),parameter              :: smallArea=1.e-2_dp  ! very small area (m2)
 ! initialize error control
 ierr=0; message='define_topology/'

 ! define reach indices
 NETOPO(:)%REACHIX = arth(1,1,nRch)
 NETOPO(:)%DREACHI = -1

 hFlag(:) = .false.
 ! first identify non-headwater basins with no upstream area
 do iRch=1,nRch
  uFlag(:) = .false.; where(NETOPO(:)%DREACHK == NETOPO(iRch)%REACHID) uFlag(:) = .true.
  ! remove links where there is no upstream area
  if(RPARAM(iRch)%UPSAREA < smallArea .and. count(uFlag) > 0)then
   where(uFlag) NETOPO(:)%DREACHK = 0
  endif
 end do

 ! define additional indices
 do iRch=1,nRch
  ! get the number of upstream reaches
  nUps=0
  do jRch=1,nRch
   ! check if another reach flows into the current reach
   if(NETOPO(jRch)%DREACHK == NETOPO(iRch)%REACHID)then
    NETOPO(jRch)%DREACHI=NETOPO(iRch)%REACHIX  ! define indices of the downstream reaches
    nUps = nUps + 1                            ! define the number of upstream reaches
   endif
  end do
  ! allocate space for the number of upstream reaches
  allocate(NETOPO(iRch)%UREACHI(nUps),NETOPO(iRch)%UREACHK(nUps), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for upstream reaches'; return; endif
  ! define the upstream reaches
  iUps=0
  do jRch=1,nRch
   ! check if another reach flows into the current reach
   if(NETOPO(jRch)%DREACHK == NETOPO(iRch)%REACHID)then
    iUps = iUps + 1
    if(iUps > size(NETOPO(iRch)%UREACHI))then; message=trim(message)//'upstream index exceeds dimension'; return; endif
    NETOPO(iRch)%UREACHI(iUps) = NETOPO(jRch)%REACHIX
    NETOPO(iRch)%UREACHK(iUps) = NETOPO(jRch)%REACHID
   endif   ! (if another reach flows into the current reach)
  end do  ! (looping through all other reaches)
  ! nullify the reach list (just to be safe)
  NETOPO(iRch)%RCHLIST => null()            ! List of reaches upstream
  !print*, irch, NETOPO(iRch)%UREACHK
 end do  ! (looping through reaches)

 ! check there is only one reach with the index of -1
 !if(count(NETOPO(:)%DREACHI < 0) /= 1)then
 ! print*, 'number of outlets = ', count(NETOPO(:)%DREACHI < 0)
 ! ierr=20; message=trim(message)//'more than one outlet'; return
 !endif

 end subroutine define_topology


end module read_streamSeg
