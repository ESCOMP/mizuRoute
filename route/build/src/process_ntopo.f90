module process_ntopo

! data types
USE nrtype,    only : i4b,dp,lgt          ! variable types, etc.
USE nrtype,    only : strLen              ! length of characters
USE dataTypes, only : var_ilength         ! integer type:          var(:)%dat
USE dataTypes, only : var_clength         ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength,dlength ! double precision type: var(:)%dat, or dat

! global vars
USE public_var, only : idSegOut           ! ID for stream segment at the bottom of the subset

! options
USE public_var, only : topoNetworkOption  ! option to compute network topology
USE public_var, only : computeReachList   ! option to compute reach list
USE public_var, only : hydGeometryOption  ! option to obtain routing parameters
USE public_var, only : routOpt            ! option for desired routing method
USE public_var, only : allRoutingMethods  ! option for routing methods - all the methods
USE public_var, only : kinematicWave      ! option for routing methods - kinematic wave only
USE public_var, only : impulseResponseFunc! option for routing methods - IRF only

! named variables
USE globalData, only : true,false         ! named integers for true/false

! named variables
USE var_lookup,only:ixSEG                 ! index of variables for the stream segments
USE var_lookup,only:ixNTOPO               ! index of variables for the network topology
USE var_lookup,only:ixPFAF                ! index of variables for the pfafstetter code

! common variables
USE public_var, only : compute            ! compute given variable
USE public_var, only : doNotCompute       ! do not compute given variable
USE public_var, only : readFromFile       ! read given variable from a file
USE public_var, only : realMissing         ! missing value for real
USE public_var, only : integerMissing      ! missing value for integers

implicit none

! privacy -- everything private unless declared explicitly
private
public::augment_ntopo
public::put_data_struct

contains

 ! *********************************************************************
 ! public subroutine: augment river network data
 ! *********************************************************************
 subroutine augment_ntopo(&
                  ! input: model control
                  nHRU,             & ! number of HRUs
                  nSeg,             & ! number of stream segments
                  ! inout: populate data structures
                  structHRU,        & ! ancillary data for HRUs
                  structSEG,        & ! ancillary data for stream segments
                  structHRU2seg,    & ! ancillary data for mapping hru2basin
                  structNTOPO,      & ! ancillary data for network toopology
                  ! output:
                  tot_hru,          & ! total number of all the upstream hrus for all stream segments
                  tot_upseg,        & ! total number of immediate upstream segments for all  stream segments
                  tot_upstream,     & ! total number of all the upstream segments for all stream segments
                  tot_uh,           & ! total number of unit hydrograph from all the stream segments
                  ixHRU_desired,    & ! indices of desired hrus
                  ixSeg_desired,    & ! indices of desired reaches
                  ! output: error control
                  ierr, message)

 ! external subroutines/data
 ! network topology routine
 use network_topo,     only:hru2segment           ! get the mapping between HRUs and segments
 use network_topo,     only:up2downSegment        ! get the mapping between upstream and downstream segments
 use network_topo,     only:reachOrder            ! define the processing order
 use network_topo,     only:reach_list            ! reach list
 use network_topo,     only:reach_mask            ! identify all reaches upstream of a given reach
 ! Routing parameter estimation routine
 use routing_param,    only:basinUH               ! construct basin unit hydrograph
 use routing_param,    only:make_uh               ! construct reach unit hydrograph
 ! routing spatial constant parameters
 use globalData,       only:fshape, tscale        ! basin IRF routing parameters (Transfer function parameters)
 use globalData,       only:mann_n, wscale        ! KWT routing parameters (Transfer function parameters)
 use globalData,       only:velo, diff            ! IRF routing parameters (Transfer function parameters)

 USE public_var, only : dt                        ! simulation time step [sec]

 ! This subroutine populate river network topology data strucutres
 implicit none
 ! output: model control
 integer(i4b)      , intent(in)                 :: nHRU             ! number of HRUs
 integer(i4b)      , intent(in)                 :: nSeg             ! number of stream segments
 ! inout: populate data structures
 type(var_dlength) , intent(inout), allocatable :: structHRU(:)     ! HRU properties
 type(var_dlength) , intent(inout), allocatable :: structSEG(:)     ! stream segment properties
 type(var_ilength) , intent(inout), allocatable :: structHRU2seg(:) ! HRU-to-segment mapping
 type(var_ilength) , intent(inout), allocatable :: structNTOPO(:)   ! network topology
 ! output:
 integer(i4b),       intent(out)                :: tot_upstream     ! total number of all of the upstream stream segments for all stream segments
 integer(i4b),       intent(out)                :: tot_upseg        ! total number of immediate upstream segments for all  stream segments
 integer(i4b),       intent(out)                :: tot_hru          ! total number of all the upstream hrus for all stream segments
 integer(i4b),       intent(out)                :: tot_uh           ! total number of unit hydrograph from all the stream segments
 integer(i4b),       intent(out),   allocatable :: ixHRU_desired(:) ! indices of desired hrus
 integer(i4b),       intent(out),   allocatable :: ixSeg_desired(:) ! indices of desired reaches
 ! output: error control
 integer(i4b)      , intent(out)                :: ierr             ! error code
 character(*)      , intent(out)                :: message          ! error message
 ! --------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=strLen)                          :: cmessage           ! error message of downwind routine
 integer(i4b)                                   :: iSeg               ! indices for stream segment
 integer(i4b)   , parameter                     :: maxUpstreamFile=10000000 ! 10 million: maximum number of upstream reaches to enable writing
 integer*8                                      :: time0,time1        ! times
 real(dp)     , allocatable                     :: seg_length(:)      ! temporal array for segment length
 type(dlength), allocatable                     :: temp_dat(:)        ! temporal data storage

 ! initialize error control
 ierr=0; message='augment_ntopo/'

 ! initialize times
 call system_clock(time0)

 ! ---------- get the mapping between HRUs and segments ------------------------------------------------------

 ! check the need to compute network topology
 if(topoNetworkOption==compute)then

  ! get the mapping between HRUs and basins
  call hru2segment(&
                   ! input
                   nHRU,          & ! input: number of HRUs
                   nSeg,          & ! input: number of stream segments
                   ! input-output: data structures
                   structHRU,     & ! ancillary data for HRUs
                   structSEG,     & ! ancillary data for stream segments
                   structHRU2seg, & ! ancillary data for mapping hru2basin
                   structNTOPO,   & ! ancillary data for network toopology
                   ! output
                   tot_hru,    &   ! output: total number of all the upstream hrus for all stream segments
                   ierr, cmessage) ! output: error control

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,i20)') 'after hru2segment: time = ', time1-time0
  !print*, trim(message)//'PAUSE : '; read(*,*)

 endif  ! if need to compute network topology

 ! ---------- get the mapping between upstream and downstream segments ---------------------------------------

 ! check the need to compute network topology
 if(topoNetworkOption==compute)then

  ! get the mapping between upstream and downstream segments
  call up2downSegment(&
                      ! input
                      nSeg,          & ! input: number of stream segments
                      ! input-output: data structures
                      structNTOPO,   & ! ancillary data for network toopology
                      ! output
                      tot_upseg,     & ! output: sum of immediate upstream segments
                      ierr, cmessage)  ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,i20)') 'after up2downSegment: time = ', time1-time0
  !print*, trim(message)//'PAUSE : '; read(*,*)

 endif  ! if need to compute network topology

 ! ---------- get the processing order -----------------------------------------------------------------------

 ! check the need to compute network topology
 if(topoNetworkOption==compute)then

  ! defines the processing order for the individual stream segments in the river network
  call REACHORDER(nSeg,         &   ! input:        number of reaches
                  structNTOPO,  &   ! input:output: network topology
                  ierr, cmessage)   ! output:       error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,i20)') 'after reachOrder: time = ', time1-time0
  !print*, trim(message)//'PAUSE : '; read(*,*)

 endif  ! if need to compute network topology

 ! ---------- get the list of all upstream reaches above a given reach ---------------------------------------

 ! get the list of all upstream reaches above a given reach
 call reach_list(&
                 ! input
                 nSeg,                        & ! Number of reaches
                 (computeReachList==compute), & ! flag to compute the reach list
                 structNTOPO,                 & ! Network topology
                 ! output
                 structSEG,                   & ! input: ancillary data for stream segments
                 tot_upstream,                & ! Total number of upstream reaches for all reaches
                 ierr, cmessage)                ! Error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 !write(*,'(a,1x,i20)') 'after reach_list: time = ', time1-time0
 !print*, trim(message)//'PAUSE : '; read(*,*)

 ! ---------- Compute routing parameters  --------------------------------------------------------------------

 ! get lag times in the basin unit hydrograph
 call basinUH(dt, fshape, tscale, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! compute hydraulic geometry (width and Manning's "n")
 if(hydGeometryOption==compute)then

  ! (hydraulic geometry only needed for the kinematic wave method)
  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
   do iSeg=1,nSeg
    structSEG(iSeg)%var(ixSEG%width)%dat(1) = wscale * sqrt(structSEG(iSeg)%var(ixSEG%totalArea)%dat(1))  ! channel width (m)
    structSEG(iSeg)%var(ixSEG%man_n)%dat(1) = mann_n                                                      ! Manning's "n" paramater (unitless)
   end do
  end if

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,i20)') 'after river geometry : time = ', time1-time0
  !print*, trim(message)//'PAUSE : '; read(*,*)

 endif  ! computing hydraulic geometry

 ! get the channel unit hydrograph
 if(topoNetworkOption==compute)then

  ! (channel unit hydrograph is only needed for the impulse response function)
  if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then

   ! extract the length information from the structure and place it in a vector
   allocate(seg_length(nSeg), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': seg_length'; return; endif
   forall(iSeg=1:nSeg) seg_length(iSeg) = structSEG(iSeg)%var(ixSEG%length)%dat(1)

   ! compute lag times in the channel unit hydrograph
   call make_uh(seg_length, dt, velo, diff, temp_dat, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! put the lag times in the data structures
   tot_uh = 0
   do iSeg=1,nSeg
    allocate(structSEG(iSeg)%var(ixSEG%timeDelayHist)%dat(size(temp_dat(iSeg)%dat)), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//': structSEG%var(ixSEG%uh)%dat'; return; endif
    structSEG(iSeg)%var(ixSEG%timeDelayHist)%dat(:) = temp_dat(iSeg)%dat(:)
    tot_uh = tot_uh+size(temp_dat(iSeg)%dat)
   enddo

  endif ! if using the impulse response function

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,i20)') 'after topoNetwork : time = ', time1-time0
  !print*, trim(message)//'PAUSE : '; read(*,*)

 endif ! if there is a need to compute the channel unit hydrograph

 ! ---------- get the mask of all upstream reaches above a given reach ---------------------------------------

 ! get the mask of all upstream reaches above a given reach
 call reach_mask(&
                 ! input
                 idSegOut,      &  ! input: reach index
                 structNTOPO,   &  ! input: network topology structures
                 structSeg,     &  ! input: river reach properties
                 nHRU,          &  ! input: number of HRUs
                 nSeg,          &  ! input: number of reaches
                 ! output: updated dimensions
                 tot_hru,       &  ! input+output: total number of all the upstream hrus for all stream segments
                 tot_upseg,     &  ! input+output: sum of immediate upstream segments
                 tot_upstream,  &  ! input+output: total number of upstream reaches for all reaches
                 tot_uh,        &  ! input+output: total number of unit hydrograph dimensions
                 ! output: dimension masks
                 ixHRU_desired, &  ! output: indices of desired hrus
                 ixSeg_desired, &  ! output: indices of desired reaches
                 ! output: error control
                 ierr, cmessage )  ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 !write(*,'(a,1x,i20)') 'after reach_mask: time = ', time1-time0
 !print*, trim(message)//'PAUSE : '; read(*,*)

end subroutine augment_ntopo

 ! *********************************************************************
 ! public subroutine: populate old data strucutures
 ! *********************************************************************
 ! ---------- temporary code: populate old data structures --------------------------------------------------
 subroutine put_data_struct(nSeg, structSEG, structNTOPO, ierr, message)
  ! Updating global parameters
  use globalData, only : RPARAM             ! Reach parameters
  use globalData, only : NETOPO             ! Network topology

  USE public_var, only : min_slope          ! minimum slope

  implicit none
  ! input
  integer(i4b)      , intent(in)                 :: nSeg             ! number of stream segments
  ! inout: populate data structures
  type(var_dlength) , intent(in)                 :: structSEG(:)     ! stream segment properties
  type(var_ilength) , intent(in)                 :: structNTOPO(:)   ! network topology
  ! output: error control
  integer(i4b)      , intent(out)                :: ierr             ! error code
  character(*)      , intent(out)                :: message          ! error message
  ! local varialbles
  character(len=strLen)                          :: cmessage         ! error message of downwind routine
  integer(i4b)                                   :: nUps             ! number of upstream segments for a segment
  integer(i4b)                                   :: iSeg,iUps        ! loop indices

  ! initialize error control
  ierr=0; message='put_data_struct/'

  ! allocate space
  allocate(RPARAM(nSeg), NETOPO(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for old data structures'; return; endif

  ! loop through stream segments
  do iSeg=1,nSeg

   ! print progress
   if(mod(iSeg,1000000)==0) print*, 'Copying to the old data structures: iSeg, nSeg = ', iSeg, nSeg

   ! ----- reach parameters -----

   ! copy data into the reach parameter structure
   RPARAM(iSeg)%RLENGTH =     structSEG(iSeg)%var(ixSEG%length)%dat(1)
   RPARAM(iSeg)%R_SLOPE = max(structSEG(iSeg)%var(ixSEG%slope)%dat(1), min_slope)
   RPARAM(iSeg)%R_MAN_N =     structSEG(iSeg)%var(ixSEG%man_n)%dat(1)
   RPARAM(iSeg)%R_WIDTH =     structSEG(iSeg)%var(ixSEG%width)%dat(1)

   ! compute variables
   RPARAM(iSeg)%BASAREA = structSEG(iSeg)%var(ixSEG%basArea)%dat(1)
   RPARAM(iSeg)%UPSAREA = structSEG(iSeg)%var(ixSEG%upsArea)%dat(1)
   RPARAM(iSeg)%TOTAREA = structSEG(iSeg)%var(ixSEG%totalArea)%dat(1)

   ! NOT USED: MINFLOW -- minimum environmental flow
   RPARAM(iSeg)%MINFLOW = structSEG(iSeg)%var(ixSEG%minFlow)%dat(1)

   ! ----- network topology -----

   ! reach indices
   NETOPO(iSeg)%REACHIX = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)     ! reach index (1, 2, 3, ..., nSeg)
   NETOPO(iSeg)%REACHID = structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1)        ! reach ID (unique reach identifier)

   ! downstream reach indices
   NETOPO(iSeg)%DREACHI = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1) ! Immediate Downstream reach index
   NETOPO(iSeg)%DREACHK = structNTOPO(iSeg)%var(ixNTOPO%downSegId)%dat(1)    ! Immediate Downstream reach ID

   ! allocate space for immediate upstream reach indices
   nUps = size(structNTOPO(iSeg)%var(ixNTOPO%upSegIds)%dat)
   allocate(NETOPO(iSeg)%UREACHI(nUps), NETOPO(iSeg)%UREACHK(nUps), NETOPO(iSeg)%goodBas(nUps), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate space for upstream structures'; return; endif

   ! populate immediate upstream data structures
   if(nUps>0)then
    do iUps=1,nUps   ! looping through upstream reaches
     NETOPO(iSeg)%UREACHI(iUps) = structNTOPO(iSeg)%var(ixNTOPO%upSegIndices)%dat(iUps)      ! Immediate Upstream reach indices
     NETOPO(iSeg)%UREACHK(iUps) = structNTOPO(iSeg)%var(ixNTOPO%upSegIds    )%dat(iUps)      ! Immediate Upstream reach Ids
     NETOPO(iSeg)%goodBas(iUps) = (structNTOPO(iSeg)%var(ixNTOPO%goodBasin)%dat(iUps)==true) ! "good" basin
    end do  ! Loop through upstream reaches
   endif

   ! define the reach order
   NETOPO(iSeg)%RHORDER = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)  ! Processing sequence

   ! allocate space for contributing HRUs
   nUps = structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1)
   allocate(NETOPO(iSeg)%HRUID(nUps), NETOPO(iSeg)%HRUIX(nUps), NETOPO(iSeg)%HRUWGT(nUps), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate space for contributing HRUs'; return; endif

   ! HRU2SEG topology
   if(nUps>0)then
     do iUps = 1, nUps
       NETOPO(iSeg)%HRUID(iUps) = structNTOPO(iSeg)%var(ixNTOPO%hruContribId)%dat(iUps)
       NETOPO(iSeg)%HRUIX(iUps) = structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat(iUps)
       NETOPO(iSeg)%HRUWGT(iUps) = structSEG(iSeg)%var(ixSEG%weight)%dat(iUps)
     end do  ! Loop through contributing HRU loop
   end if

   ! NOT USED: lake parameters
   NETOPO(iSeg)%LAKE_IX = integerMissing  ! Lake index (0,1,2,...,nlak-1)
   NETOPO(iSeg)%LAKE_ID = integerMissing  ! Lake ID (REC code?)
   NETOPO(iSeg)%BASULAK = realMissing     ! Area of basin under lake
   NETOPO(iSeg)%RCHULAK = realMissing     ! Length of reach under lake
   NETOPO(iSeg)%LAKINLT = .false.         ! .TRUE. if reach is lake inlet, .FALSE. otherwise
   NETOPO(iSeg)%USRTAKE = .false.         ! .TRUE. if user takes from reach, .FALSE. otherwise

   ! NOT USED: Location (available in the input files)
   NETOPO(iSeg)%RCHLAT1 = realMissing     ! Start latitude
   NETOPO(iSeg)%RCHLAT2 = realMissing     ! End latitude
   NETOPO(iSeg)%RCHLON1 = realMissing     ! Start longitude
   NETOPO(iSeg)%RCHLON2 = realMissing     ! End longitude

   ! reach unit hydrograph
   if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
     allocate(NETOPO(iSeg)%UH(size(structSEG(iSeg)%var(ixSEG%timeDelayHist)%dat)), stat=ierr, errmsg=cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage)//': NETOPO(iSeg)%UH'; return; endif
     NETOPO(iSeg)%UH(:) =  structSEG(iSeg)%var(ixSEG%timeDelayHist)%dat(:)
   end if

   ! upstream reach list
   allocate(NETOPO(iSeg)%RCHLIST(size(structNTOPO(iSeg)%var(ixNTOPO%allUpSegIndices)%dat)), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': NETOPO(iSeg)%RCHLIST'; return; endif
   NETOPO(iSeg)%RCHLIST(:) =  structNTOPO(iSeg)%var(ixNTOPO%allUpSegIndices)%dat(:)

  end do  ! looping through stream segments

 end subroutine put_data_struct

end module process_ntopo
