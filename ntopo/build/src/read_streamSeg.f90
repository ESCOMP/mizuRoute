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
                      total_hru,  &   ! output: total number of HRUs that drain into any segments
                      ierr, message)  ! output: error control
 USE nr_utility_module,ONLY: arth     ! Num. Recipies utilities
 USE nr_utility_module,ONLY: indexx   ! Num. Recipies utilities
 USE dataTypes,only:namepvar,nameivar ! provide access to data types
 USE var_lookup,only:ixHRU,nVarsHRU   ! index of variables for the HRUs
 USE var_lookup,only:ixSEG,nVarsSEG   ! index of variables for the stream segments
 USE var_lookup,only:ixMAP,nVarsMAP   ! index of variables for the hru2basin mapping
 USE var_lookup,only:ixTOP,nVarsTOP   ! index of variables for the network topology
 USE nhru2basin  ! data structures holding the nhru2basin correspondence
 implicit none
 ! input variables
 integer(i4b), intent(in)        :: nHRU              ! number of HRUs
 integer(i4b), intent(in)        :: nSeg              ! number of stream segments
 type(namepvar), intent(in)      :: nhru_acil(:)      ! ancillary data for the HRUs
 type(nameivar), intent(in)      :: imap_acil(:)      ! ancillary data for the hru2basin mapping
 type(nameivar), intent(in)      :: ntop_acil(:)      ! ancillary data for the network topology
 ! output variables
 integer(i4b), intent(out)       :: total_hru         ! total number of HRUs that drain into any segments (nHRU minus number of HRU not draining int any segments)
 integer(i4b), intent(out)       :: ierr              ! error code
 character(*), intent(out)       :: message           ! error message
 ! local variables
 integer(i4b)                    :: iHRU              ! index of HRU
 integer(i4b)                    :: iSeg              ! index of stream segment
 integer(i4b)                    :: jSeg              ! id of stream segment
 integer(i4b)                    :: mapCOMid          ! COM Id for the hru2seg mapping vector
 integer(i4b)                    :: segCOMid          ! COM Id for the stream segment
 integer(i4b)                    :: rankSeg(nSeg)     ! index of segment in the nSeg vector
 integer(i4b)                    :: rankSegHRU(nHRU)  ! index of segment in the hru2seg mapping vector
 integer(i4b)                    :: segHRUix(nHRU)    ! index of segment where HRU drains
 integer(i4b)                    :: nHRU2seg(nSeg)    ! number of HRUs that drain into a given segment
 logical(lgt),parameter          :: checkMap=.true.   ! flag to check the mapping
 real(dp)                        :: totarea           ! total area of all HRUs feeding into a given stream segment (m2)
 integer*8                       :: time0,time1       ! times

 ! initialize error control
 ierr=0; message='hru2basin/'

 !print*, 'PAUSE: start of '//trim(message); read(*,*)

 ! initialize timing
 call system_clock(time0)

 ! allocate space for the hru2basin mapping
 allocate(h2b(nSeg),stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for h2b structure'; return; endif

 ! ---------- get the index of the stream segment that a given HRU drains into ------------------------------

 ! initialize the HRUs that drain into a given segment
 nHRU2seg(:) = 0
 segHRUix(:) = -999

 ! initialize total number of HRUs that drain into any segments
 total_hru = 0

 ! rank the stream segments
 call indexx(ntop_acil(ixTOP%segid)%varData, rankSeg)

 ! rank the drainage segment
 call indexx(imap_acil(ixMAP%segHRUid)%varData, rankSegHRU)

 iSeg=1  ! second counter
 ! loop through the HRUs
 do iHRU=1,nHRU

  mapCOMid = imap_acil(ixMAP%segHRUid)%varData( rankSegHRU(iHRU) )
  if (mapCOMid == 0) cycle

  ! keep going until found the index
  do jSeg=iSeg,nSeg ! normally a short loop

   ! get the COM Ids for the hru2seg mapping vector and the segments
   segCOMid = ntop_acil(ixTOP%segid)%varData( rankSeg(jSeg) )
   !print*, 'iHRU, iSeg, jSeg, mapCOMid, segCOMid, rankSegHRU(iHRU), rankSeg(jSeg) = ', &
   !         iHRU, iSeg, jSeg, mapCOMid, segCOMid, rankSegHRU(iHRU), rankSeg(jSeg)

   ! define the index where we have a match
   if(mapCOMid==segCOMid)then
    segHRUix( rankSegHRU(iHRU) ) = rankSeg(jSeg)
    nHRU2seg( rankSeg(jSeg) )    = nHRU2seg( rankSeg(jSeg) ) + 1
    total_hru = total_hru + 1
    !iSeg=jSeg+1
    exit
   endif

  end do  ! skipping segments that have no HRU input
 end do  ! looping through HRUs

 ! check
 if(checkMap)then
  do iHRU=1,nHRU
   if (imap_acil(ixMAP%segHRUid)%varData(iHRU) == 0) cycle
   if( imap_acil(ixMAP%segHRUid)%varData(iHRU) /= ntop_acil(ixTOP%segid)%varData( segHRUix(iHRU) ) )then
    message=trim(message)//'problems identifying the index of the stream segment that a given HRU drains into'
    ierr=20; return
   endif
  end do
 endif

 ! ---------- allocate space for the mapping structures -----------------------------------------------------

 ! loop through stream segments
 do iSeg=1,nSeg
  ! allocate space
  allocate(h2b(iSeg)%cHRU( nHRU2seg(iSeg) ), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for h2b structure component'; return; endif
  ! initialize the number of HRUs
  h2b(iSeg)%nHRU = 0
 end do

 ! get timing
 call system_clock(time1)
 print*, 'timing: allocate space = ', time1-time0

 ! ---------- populate structure components for HRU-2-Segment mapping ---------------------------------------

 ! loop through HRUs
 do iHRU=1,nHRU

  ! identify the index of the stream segment that the HRU drains into
  iSeg = segHRUix(iHRU)

  ! if there is no stream segment associated with current hru
  if (iSeg ==-999) cycle

  ! increment the HRU counter
  h2b(iSeg)%nHRU = h2b(iSeg)%nHRU+1

  ! populate structure components
  h2b(iSeg)%cHRU( h2b(iSeg)%nHRU )%hru_ix   = iHRU
  h2b(iSeg)%cHRU( h2b(iSeg)%nHRU )%hru_id   = imap_acil(ixMAP%HRUid)%varData(iHRU)
  h2b(iSeg)%cHRU( h2b(iSeg)%nHRU )%hru_lon  = nhru_acil(ixHRU%xLon )%varData(iHRU)
  h2b(iSeg)%cHRU( h2b(iSeg)%nHRU )%hru_lat  = nhru_acil(ixHRU%yLat )%varData(iHRU)
  h2b(iSeg)%cHRU( h2b(iSeg)%nHRU )%hru_elev = nhru_acil(ixHRU%elev )%varData(iHRU)
  h2b(iSeg)%cHRU( h2b(iSeg)%nHRU )%hru_area = nhru_acil(ixHRU%area )%varData(iHRU)

 end do ! looping through HRUs

 ! check
 if(checkMap)then
  do iSeg=1,nSeg
   if(nHRU2seg(iSeg)/=h2b(iSeg)%nHRU)then
    message=trim(message)//'problems identifying the HRUs draining into stream segment'
    ierr=20; return
   endif
  end do
 endif

 ! get timing
 call system_clock(time1)
 print*, 'timing: populate structure components = ', time1-time0

 ! ---------- compute HRU weights ---------------------------------------------------------------------------

 ! loop through segments
 do iSeg=1,nSeg

   ! compute total area of the HRUs draining to the stream segment
   totarea = sum(h2b(iSeg)%cHRU(:)%hru_area)

   ! compute the weights
   h2b(iSeg)%cHRU(:)%wght = h2b(iSeg)%cHRU(:)%hru_area / totarea

 end do  ! (looping thru stream segments)

 ! get timing
 call system_clock(time1)
 print*, 'timing: compute HRU weights = ', time1-time0
 print*, 'PAUSE: end of '//trim(message); read(*,*)

 end subroutine hru2basin


 ! *********************************************************************
 ! new subroutine: assign data to the reach parameter structure
 ! *********************************************************************
 subroutine assign_reachparam(nRch,         &    ! input: number of stream segments
                              sseg_acil,    &    ! input: stream segment parameters
                              ntop_acil,    &    ! input: network topology
                              total_upseg,  &    ! output: sum of immediate upstream segments
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
 integer(i4b), intent(out)       :: total_upseg  ! sum of immediate upstream segments
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
 RPARAM(:)%RLENGTH = sseg_acil(ixSEG%length   )%varData(:)
 RPARAM(:)%R_SLOPE = sseg_acil(ixSEG%slope    )%varData(:)

 ! compute area draining to each stream segment
 do iRch=1,nRch
  RPARAM(iRch)%BASAREA = sum(h2b(iRch)%cHRU(:)%hru_area)
 end do

 ! define additional aspects of the network topology
 ! (populates reachparam)
 call define_topology(nrch, total_upseg, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ensure that slope exceeds minimum slope
 do iRch=1,nRch
  if(RPARAM(iRch)%R_SLOPE < min_slope) RPARAM(iRch)%R_SLOPE = min_slope
 end do

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
                            total_upseg, &     ! output: sum of immediate upstream segments
                            ierr, message)     ! output (error control)
 USE nr_utility_module, ONLY: arth             ! Num. Recipies utilities
 USE reachparam                                ! reach parameter structure
 implicit none
 ! input variables
 integer(i4b), intent(in)        :: nRch         ! number of reaches
 ! output variables
 integer(i4b), intent(out)       :: total_upseg  ! total number of immediate upstream reaches
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: iRch,jRch    ! reach indices
 integer(i4b)                    :: iUps         ! index of upstream reach
 integer(i4b)                    :: nUps         ! number of upstream reaches
 ! initialize error control
 ierr=0; message='define_topology/'

 ! define reach indices
 NETOPO(:)%REACHIX = arth(1,1,nRch)
 NETOPO(:)%DREACHI = -1

 ! define additional indices
 total_upseg = 0
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
  total_upseg = total_upseg+nUps
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
