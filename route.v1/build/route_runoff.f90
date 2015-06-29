program route_runoff

! ******
! provide access to desired modules
! ************************
USE nrtype                                    ! variable types, etc.
USE dataTypes,only:namepvar,nameivar          ! provide access to data types
USE var_lookup,only:ixHRU,nVarsHRU            ! index of variables for the HRUs
USE var_lookup,only:ixSEG,nVarsSEG            ! index of variables for the stream segments
USE var_lookup,only:ixMAP,nVarsMAP            ! index of variables for the hru2basin mapping
USE var_lookup,only:ixTOP,nVarsTOP            ! index of variables for the network topology
USE reachparam                                ! reach parameters 
USE reachstate                                ! reach states
USE reach_flux                                ! fluxes in each reach
USE nhru2basin                                ! data structures holding the nhru2basin correspondence
USE nrutil,only:arth                          ! use to build vectors with regular increments
USE ascii_util_module,only:file_open          ! open file (performs a few checks as well)
USE ascii_util_module,only:get_vlines         ! get a list of character strings from non-comment lines
USE read_simoutput,only:get_qDims             ! get the dimensions from the runoff file
USE read_simoutput,only:get_qMeta             ! get the metadata from the runoff file
USE read_simoutput,only:getVarQsim            ! get the data from the runoff file
USE read_ntopo,only:get_vec_ivar              ! get the data from the ntopo file
USE read_ntopo,only:get_scl_ivar              ! get the data from the ntopo file
USE read_ntopo,only:get_scl_dvar              ! get the data from the ntopo file
USE read_ntopo,only:get_vec_dvar              ! get the data from the ntopo file
USE read_ntopo,only:get_Vec_dim               ! get the data from the ntopo file
USE write_simoutput,only:defineFile           ! define output file
USE write_simoutput,only:write_iVec           ! write an integer vector
USE write_simoutput,only:write_dVec           ! write a double precision vector
USE network_route,only:reachorder             ! define the processing order for the stream segments
USE network_route,only:qroute_rch             ! route kinematic waves through the river network
USE vic_route,only:upstrm_length              ! Compute total upstream length  NM
USE vic_route,only:make_uh                    ! Compute upstream segment UH for segment NM
USE vic_route,only:get_upsbas_qr              ! Compute total upstream basin q NM
USE vic_route,only:conv_upsbas_qr             ! Compute convoluted upstream basin q NM

! ******
! define variables
! ************************
implicit none
! index for printing (set to negative to supress printing
integer(i4b),parameter     :: ixPrint = -9999     ! index for printing
! useful constants
real(dp),parameter         :: secprhour=3600._dp  ! number of seconds in an hour
real(dp),parameter         :: secprday=86400._dp  ! number of seconds in a day
logical(lgt),parameter     :: doVicRoute=.True.  ! switch to turn off/on VIC routing added by NM
logical(lgt),parameter     :: doKWaveRoute=.True.
real(dp),parameter         :: verySmall=tiny(1.0_dp)  ! a very small number 
! general guff
integer(i4b),parameter     :: strLen=256          ! length of character string
integer(i4b)               :: ierr                ! error code
character(len=strLen)      :: cmessage            ! error message of downwind routine
integer(i4b)               :: iTime               ! loop through time
character(len=strLen)      :: str                 ! miscellaneous string
! read control file
character(len=strLen)      :: cfile_name          ! name of the control file
character(len=strLen)      :: param_nml           ! name of the namelist file 
integer(i4b)               :: iunit               ! file unit
character(len=strLen),allocatable :: cLines(:)    ! vector of character strings
integer(i4b)               :: ipos                ! index of character string
integer(i4b)               :: iLine               ! index of line in cLines
integer(i4b)               :: ibeg_name           ! start index of variable name in string cLines(iLine)
integer(i4b)               :: iend_name           ! end index of variable name in string cLines(iLine)
integer(i4b)               :: iend_data           ! end index of data in string cLines(iLine)
character(len=strLen)      :: cName,cData         ! name and data from cLines(iLine)
! define output file
character(len=strLen)      :: fname_output        ! name of output file
integer(i4b),parameter     :: outunit=31          ! unit for output file
! define directories
character(len=strLen)      :: ancil_dir           ! directory containing ancillary data
character(len=strLen)      :: input_dir           ! directory containing input data
character(len=strLen)      :: output_dir          ! directory containing output data
! define stream segment information
integer(i4b),target        :: nSeg                ! number of all the stream segments
integer(i4b),pointer       :: nSegRoute           ! number of stream segments to be routed
integer(i4b)               :: nUpstream           ! number of reaches upstream of each stream segment
integer(i4b)               :: iSeg                ! index of stream segment
integer(i4b)               :: jSeg                ! index of stream segment
integer(i4b)               :: iSegOut             ! index of outlet stream segment
integer(i4b)               :: iSelect(1)          ! index of desired stream segment (iSegOut) from the minloc operation
integer(i4b)               :: iSegDesire          ! index of desired stream segment -- de-vectorized version of iSelect(1)
integer(i4b)               :: iUps                ! index of upstream stream segment added by NM
!integer(i4b)               :: jUps               ! index of upstream stream segment added by NM
integer(i4b)               :: iStart              ! start index of the ragged array
integer(i4b),dimension(1)  :: iDesire             ! index of stream segment with maximum upstream area (vector)
integer(i4b)               :: ixDesire            ! index of stream segment with maximum upstream area (scalar)
!integer(i4b)               :: iTDH               ! index of unit hydrograph data element 
! define stream network information
character(len=strLen)      :: fname_ntop          ! filename containing stream network topology information
integer(i4b),allocatable   :: REACHIDGV(:)
integer(i4b),allocatable   :: RCHIXLIST(:)
integer(i4b)               :: nTotal              ! total number of upstream segments for all stream segments
integer(i4b)               :: iRchStart   
integer(i4b)               :: iRchStart1   
integer(i4b),target        :: nRchCount   
integer(i4b)               :: nRchCount1   
integer(i4b)               :: iUpRchStart 
integer(i4b)               :: nUpRchCount 
integer(i4b)               :: iUpHruStart 
integer(i4b)               :: nUpHruCount 
integer(i4b),allocatable   :: upStrmRchList(:)
! define metadata from model output file
character(len=strLen)      :: fname_qsim          ! filename containing simulated runoff
character(len=strLen)      :: vname_qsim          ! variable name for simulated runoff
character(len=strLen)      :: vname_hruid         ! coordinate name for the HRUid
character(len=strLen)      :: vname_time          ! coordinate name for time
character(len=strLen)      :: units_qsim          ! units of simulated runoff data
character(len=strLen)      :: units_time          ! time units
integer(i4b),allocatable   :: qsimHRUid(:)        ! vector of HRU ids from simulated runoff file
logical(lgt),allocatable   :: qsimHRUid_mask(:)   ! vector of HRU mask that is in ustream of outlet
real(dp),allocatable       :: qsimHRUarea(:)         ! vector of HRU area from simulated runoff file
!integer(i4b),allocatable   :: HRUix(:)           ! vector of HRU indices
integer(i4b)               :: nTime               ! number of time elements
integer(i4b)               :: nHRU_data           ! number of HRUs in the data file
integer(i4b)               :: iRch                ! index in reach structures
!integer(i4b)               :: jRch               ! index in reach structures
! define simulated runoff data at the HRUs
real(dp)                   :: dt                  ! time step (seconds)
real(dp)                   :: dTime               ! time variable (in units units_time)
real(dp), allocatable      :: qsim_hru(:)         ! simulated runoff at the HRUs
character(len=strLen)      :: cLength,cTime       ! length and time units
real(dp)                   :: time_conv           ! time conversion factor -- used to convert to mm/s
real(dp)                   :: length_conv         ! length conversion factor -- used to convert to mm/s
real(dp),parameter         :: runoffMin=1.e-15_dp ! minimum runoff from each basin
! interpolate simulated runoff data to the basins
integer(i4b)               :: ibas                ! index of the basins
integer(i4b)               :: iHRU                ! index of the HRUs associated to the basin
integer(i4b)               :: nDrain              ! number of HRUs that drain into a given stream segment
integer(i4b)               :: ix                  ! index of the HRU assigned to a given basin
real(dp), allocatable      :: qsim_basin(:)       ! simulated runoff at the basins
! route simulated runoff through the local basin
real(dp)                   :: scale_time          ! scaling factor for the time delay histogram (seconds)
integer(i4b)               :: jtim                ! index of the time delay vectors
integer(i4b)               :: ntdh                ! number of elements in the time delay histogram
! route delaied runoff through river network with St.Venant UH
real(dp)                   :: velo                ! velocity [m/s] for Saint-Venant equation   added by NM
real(dp)                   :: diff                ! diffusivity [m2/s] for Saint-Venant equation   added by NM
integer(i4b)               :: nUH_DATA_MAX        ! maximum number of elements in the UH data among all the upstreamfs for a segment 
! compute total instantaneous runoff upstream of each reach
integer(i4b),allocatable   :: iUpstream(:)        ! indices for all reaches upstream
real(dp),allocatable       :: qUpstream(:)        ! streamflow for all reaches upstream
! route kinematic waves through the river network
integer(i4b), parameter    :: nens=1              ! number of ensemble members
integer(i4b)               :: iens                ! index of ensemble member
real(dp)                   :: T0,T1               ! start and end of the time step (seconds)
real(dp)                   :: mann_n              ! manning's roughness coefficient [unitless]  added by NM
real(dp)                   :: width_scale         ! scaling factor for river width [-] added by NM
integer(i4b)               :: LAKEFLAG            ! >0 if processing lakes
! desired variables when printing progress
real(dp)                   :: qDomain_hru         ! domain: total runoff at the HRUs (mm/s)
real(dp)                   :: qDomain_basin       ! domain: total instantaneous basin runoff (mm/s)
real(dp)                   :: qDomain_reach       ! domain: total instantaneous reach runoff (mm/s)
real(dp)                   :: qDesire_instant     ! desire: sum of instantaneous runoff for all basins upstream of the desired reach (mm/s)
real(dp)                   :: qDesire_routed      ! desire: routed runoff for the desired reach (mm/s)
namelist /ST_VENANT_UH/velo,diff
namelist /PRTCL_TRCK/mann_n,width_scale

! *****
! (0) Read control file...
! ************************
! get command-line argument defining the full path to the control file
call getarg(1,cfile_name)
if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')

! *** get a list of character strings from non-comment lines ****
! open file (also returns un-used file unit used to open the file)
call file_open(trim(cfile_name),iunit,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
! get a list of character strings from non-comment lines
call get_vlines(iunit,cLines,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
! close the file unit
close(iunit)

! loop through the non-comment lines in the input file, and extract the name and the information
do iLine=1,size(cLines)
  ! identify start and end of the name and the data
  ibeg_name = index(cLines(iLine),'<'); if(ibeg_name==0) ierr=20
  iend_name = index(cLines(iLine),'>'); if(iend_name==0) ierr=20
  iend_data = index(cLines(iLine),'!'); if(iend_data==0) ierr=20
  if(ierr/=0) call handle_err(60,'problem disentangling cLines(iLine) [string='//trim(cLines(iLine))//']')
  ! extract name of the information, and the information itself
  cName = adjustl(cLines(iLine)(ibeg_name:iend_name))
  cData = adjustl(cLines(iLine)(iend_name+1:iend_data-1))
  print*, trim(cName), ' --> ', trim(cData)
  ! populate variables
  select case(trim(cName))
    ! define directories
    case('<ancil_dir>');    ancil_dir    = trim(cData)           ! directory containing ancillary data
    case('<input_dir>');    input_dir    = trim(cData)           ! directory containing input data
    case('<output_dir>');   output_dir   = trim(cData)           ! directory containing output data
    ! define variables for the network topology
    case('<fname_ntop>');   fname_ntop   = trim(cData)           ! name of file containing stream network topology information
    case('<seg_outlet>'   ); read(cData,*,iostat=ierr) iSegOut                ! seg_id of outlet streamflow segment 
      if(ierr/=0) call handle_err(70,'problem with internal read of iSegOut info, read from control file')
    ! define the file and variable name for runoff netCDF 
    case('<fname_qsim>');   fname_qsim   = trim(cData)    ! name of file containing the runoff
    case('<vname_qsim>');   vname_qsim   = trim(cData)    ! name of runoff variable
    case('<vname_hruid>');  vname_hruid  = trim(cData)    ! name of the HRU id
    case('<vname_time>');   vname_time   = trim(cData)    ! name of time variable in the runoff file
    case('<units_qsim>');   units_qsim   = trim(cData)    ! units of runoff
    case('<dt_qsim>');      read(cData,*,iostat=ierr) dt  ! time interval of the gridded runoff
      if(ierr/=0) call handle_err(70,'problem with internal read of dt info, read from control file')
    ! define the output filename
    case('<fname_output>'); fname_output = trim(cData)    ! filename for the model output
    ! define namelist name for routing parameters
    case('<param_nml>');    param_nml    = trim(cData)     ! name of namelist including routing parameter value 
  end select
end do  ! looping through lines in the control file

! find the position of the "/" character
ipos = index(trim(units_qsim),'/')
if(ipos==0) call handle_err(80,'expect the character "/" exists in the units string [units='//trim(units_qsim)//']')
! get the length and time units
cLength = units_qsim(1:ipos-1)
cTime   = units_qsim(ipos+1:len_trim(units_qsim))
! get the conversion factor for length
select case(trim(cLength))
 case('m');  length_conv = 1._dp
 case('mm'); length_conv = 1._dp/1000._dp
 case default; call handle_err(81,'expect the length units to be "m" or "mm" [units='//trim(cLength)//']')
endselect
! get the conversion factor for time
select case(trim(cTime))
 case('d','day');    time_conv = 1._dp/secprday
 case('h','hour');   time_conv = 1._dp/secprhour
 case('s','second'); time_conv = 1._dp
 case default; call handle_err(81,'cannot identify the time units [time units = '//trim(cTime)//']')
endselect

! *****
! (0.b) Read Namelist...
! ************************
call file_open(trim(param_nml),iunit,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
read(iunit, nml=ST_VENANT_UH)
read(iunit, nml=PRTCL_TRCK)

! get the number of stream segments (needed for allocate statements)
call get_vec_dim(trim(ancil_dir)//trim(fname_ntop), &
                 'sSeg',                            &
                 nSeg,                              &
                 ierr, cmessage)
call handle_err(ierr,cmessage)
print*,'Number of nSeg in Network topology: ', nSeg

! *****
! (1.2) Read in the network topology information...
! *********************************************
! Data type to be populated 
! NETOPO(:)%REACHIX       index
! NETOPO(:)%REACHID       id
! NETOPO(:)%DREACHK       id
! NETOPO(:)%DREACHI       index
! RPARAM(:)%R_SLOPE
! RPARAM(:)%BASAREA
! RPARAM(:)%UPSAREA
! RPARAM(:)%TOTAREA
! NETOPO(:)%RCHLAT1
! NETOPO(:)%RCHLAT2
! NETOPO(:)%RCHLON1
! NETOPO(:)%RCHLON2
! NETOPO(:)%RCHLIST(:)    index
! NETOPO(:)%UPSLENG(:)
! NETOPO(:)%UREACHK(:)    id
! NETOPO(:)%UREACHI(:)    index
! h2b(:)%cHRU(:)%hru_ix   index
! h2b(:)%cHRU(:)%hru_id   id 
! h2b(:)%cHRU(:)%hru_lon 
! h2b(:)%cHRU(:)%hru_lat 
! h2b(:)%cHRU(:)%hru_elev
! h2b(:)%cHRU(:)%hru_area
! h2b(:)%cHRU(:)%wght    

! Read global reach id  (input = filename, variable name, variable vector, start index; output = error control)
allocate(REACHIDGV(nSeg), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for REACHIDGV')
call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop), 'reachID', REACHIDGV, 1, nSeg, ierr, cmessage); call handle_err(ierr,cmessage)

if ( iSegOut /= -9999 ) then
  print*, 'Outlet segment = ', iSegOut
  ! Identify index of the desired stream segment from reachID vector (dimension size: nSeg)
  iSelect = minloc(abs(REACHIDGV - iSegOut))
  iSegDesire = iSelect(1)  ! de-vectorize the desired stream segment
  if(REACHIDGV(iSegDesire) /= iSegOut)&
    call handle_err(20,'unable to find desired stream segment')

  ! Read the start index and the count for lagged array - all the upstream segments, immediate upstream segment, immediate upstream HRUs
  call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachStart', iRchStart, iSegDesire, ierr, cmessage); call handle_err(ierr,cmessage)
  call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachCount', nRchCount, iSegDesire, ierr, cmessage); call handle_err(ierr,cmessage)
  print*,'iRchStart   = ',iRchStart
  print*,'Number of upstream segment from outlet segment (nRchCount): ',nRchCount
  
  ! Read reach list of index from global segments (all the upstream reachs for each segment) 
  allocate(upStrmRchList(nRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for upStrmRchList')
  call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop), 'reachList', upStrmRchList, iRchStart,nRchCount, ierr, cmessage); call handle_err(ierr,cmessage)
  
  ! Reach upstream segment and associated HRU infor from non-ragged vector 
  allocate(NETOPO(nRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO')
  allocate(RPARAM(nRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for RPARAM')
  allocate(h2b(nRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for h2b')
  
  ! Create REACH index for local segments
  NETOPO(:)%REACHIX=arth(1,1,nRchCount)
  
  do iSeg=1,nRchCount
    ! Reach reach topology and parameters (integer)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachID',     NETOPO(iSeg)%REACHID,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'downReachID', NETOPO(iSeg)%DREACHK,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    ! Reach reach topology and parameters (double precision)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachSlope',  RPARAM(iSeg)%R_SLOPE,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLength', RPARAM(iSeg)%RLENGTH,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'basinArea',   RPARAM(iSeg)%BASAREA,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'upstreamArea',RPARAM(iSeg)%UPSAREA,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'totalArea',   RPARAM(iSeg)%TOTAREA,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLat1',   NETOPO(iSeg)%RCHLAT1,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLat2',   NETOPO(iSeg)%RCHLAT2,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLon1',   NETOPO(iSeg)%RCHLON1,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLon2',   NETOPO(iSeg)%RCHLON2,upStrmRchList(iSeg),ierr,cmessage); call handle_err(ierr,cmessage)
  enddo
  
  ! Recompute downstream segment index as local segment list, or NETOPO(:)%REACHID
  do iSeg=1,nRchCount
    ! Assign downstream segment ID = 0 at desired outlet segment
    if (NETOPO(iSeg)%REACHID == iSegOut) then 
      NETOPO(iSeg)%DREACHK = 0
    else
      ! Identify the index of the desired stream segment from reachID vector (dimension size: nSeg)
      iSelect = minloc(abs(NETOPO(:)%REACHID - NETOPO(iSeg)%DREACHK))
      NETOPO(iSeg)%DREACHI = iSelect(1)  ! de-vectorize the desired stream segment
      if (NETOPO(NETOPO(iSeg)%DREACHI)%REACHID /= NETOPO(iSeg)%DREACHK) then 
        print*,'iSeg = ', iSeg 
        print*,'NETOPO(iSeg)%DREACHK = ', NETOPO(iSeg)%DREACHK
        print*,'NETOPO(NETOPO(iSeg)%DREACHI)%REACHID = ', NETOPO(NETOPO(iSeg)%DREACHI)%REACHID
        call handle_err(20,'unable to find desired downstream segment')
      endif
    endif
  enddo
  
  ! Reach upstream segment and associated HRU infor from ragged vector 
  nTotal=0
  do iSeg=1,nRchCount
    ! sAll dimension 
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachStart', iRchStart1, upStrmRchList(iSeg), ierr, cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachCount', nRchCount1, upStrmRchList(iSeg), ierr, cmessage); call handle_err(ierr,cmessage)
  
    allocate(NETOPO(iSeg)%UPSLENG(nRchCount1), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%UPSLENG')
    allocate(NETOPO(iSeg)%RCHLIST(nRchCount1), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%RCHLIST')
    allocate(RCHIXLIST(nRchCount1), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for RCHIXLIST(nRchCount)')
  
    call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop),'reachList'         ,RCHIXLIST,iRchStart1,nRchCount1,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'upReachTotalLength',NETOPO(iSeg)%UPSLENG(:),iRchStart1,nRchCount1,ierr,cmessage); call handle_err(ierr,cmessage)
  
    ! Recompute all the upstream segment indices as local segment list = NETOPO(:)%REACHID
    nTotal = nTotal + nRchCount1
    do jSeg=1,nRchCount1
      ! Identify the index of the desired stream segment from reachID vector (dimension size: nSeg)
      iSelect = minloc( abs( NETOPO(:)%REACHID - REACHIDGV(RCHIXLIST(jSeg)) ) )
      NETOPO(iSeg)%RCHLIST(jSeg) = iSelect(1)  ! de-vectorize the desired stream segment
    enddo
    !print*,'NETOPO(iSeg)%RCHLIST(:) = ',NETOPO(iSeg)%RCHLIST(:)
    deallocate(RCHIXLIST, stat=ierr)
  
    ! sUps dimension 
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'upReachStart', iUpRchStart, upStrmRchList(iSeg), ierr, cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'upReachCount', nUpRchCount, upStrmRchList(iSeg), ierr, cmessage); call handle_err(ierr,cmessage)
  
    allocate(NETOPO(iSeg)%UREACHI(nUpRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%UREACHI')
    allocate(NETOPO(iSeg)%UREACHK(nUpRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%UREACHK')
    allocate(NETOPO(iSeg)%goodBas(nUpRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%goodBas')
    if (nUpRchCount > 0) then
  
      call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop),'upReachID'   ,NETOPO(iSeg)%UREACHK(:),iUpRchStart,nUpRchCount,ierr,cmessage); call handle_err(ierr,cmessage)
      do jSeg=1,nUpRchCount
        ! Identify the index of the desired stream segment from reachID vector (dimension size: nSeg)
        iSelect = minloc(abs(NETOPO(:)%REACHID - NETOPO(iSeg)%UREACHK(jSeg)))
        NETOPO(iSeg)%UREACHI(jSeg) = iSelect(1)  ! de-vectorize the desired stream segment
        ! check that we identify the correct upstream reach
        if (NETOPO(NETOPO(iSeg)%UREACHI(jSeg))%REACHID /= NETOPO(iSeg)%UREACHK(jSeg)) then 
          print*,'iSeg = ', iSeg 
          print*,'NETOPO(iSeg)%UREACHK(jSeg) = ', NETOPO(iSeg)%UREACHK(jSeg)
          print*,'NETOPO(NETOPO(iSeg)%UREACHI(jSeg))%REACHID = ', NETOPO(NETOPO(iSeg)%UREACHI(jSeg))%REACHID
          call handle_err(20,'unable to find desired immediate upstream segment')
        endif
  
        ! check that the upstream reach has a basin area > 0
        if(RPARAM(NETOPO(iSeg)%UREACHI(jSeg))%TOTAREA > verySmall)then
         NETOPO(iSeg)%goodBas(jSeg) = .true.
        else
         NETOPO(iSeg)%goodBas(jSeg) = .false.
        endif
  
      enddo ! looping through the immediate upstream reaches
    endif  ! if not a headwater
    !print*,'NETOPO(iSeg)%UREACHI(:) = ',NETOPO(iSeg)%UREACHI(:)
    !print*,'NETOPO(iSeg)%REACHID = ',NETOPO(iSeg)%REACHID
    !print*,'NETOPO(iSeg)%UREACHK(:) = ',NETOPO(iSeg)%UREACHK(:)
    
    ! sHrus dimension 
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'upHruStart', iUpHruStart, upStrmRchList(iSeg), ierr, cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'upHruCount', nUpHruCount, upStrmRchList(iSeg), ierr, cmessage); call handle_err(ierr,cmessage)
    allocate(h2b(iSeg)%cHRU(nUpHruCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for h2b(iSeg)%cHRU(nUpHruCount)')
    if (nUpHrucount /= 0) then
      call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop),'hru_id',    h2b(iSeg)%cHRU(:)%hru_id,  iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_lon',   h2b(iSeg)%cHRU(:)%hru_lon, iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_lat',   h2b(iSeg)%cHRU(:)%hru_lat, iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_elev',  h2b(iSeg)%cHRU(:)%hru_elev,iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_area',  h2b(iSeg)%cHRU(:)%hru_area,iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_weight',h2b(iSeg)%cHRU(:)%wght,    iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
    endif
    !print*,'h2b(iSeg)%cHRU(:)%hru_id = ',h2b(iSeg)%cHRU(:)%hru_id
  enddo  ! looping through the stream segments within the model domain
  nSegRoute => nRchCount

else ! if the entire river network routing is selected
  print*, 'Route all the segments included in network topology'
  ! Populate sSeg dimensioned variable 
  allocate(NETOPO(nSeg), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO')
  allocate(RPARAM(nSeg), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for RPARAM')
  do iSeg=1,nSeg
    ! Reach reach topology and parameters (integer)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachIndex',    NETOPO(iSeg)%REACHIX,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachID',       NETOPO(iSeg)%REACHID,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLat1',     NETOPO(iSeg)%RCHLAT1,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLat2',     NETOPO(iSeg)%RCHLAT2,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLon1',     NETOPO(iSeg)%RCHLON1,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLon2',     NETOPO(iSeg)%RCHLON2,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'downReachIndex',NETOPO(iSeg)%DREACHI,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'downReachID',   NETOPO(iSeg)%DREACHK,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    ! Reach reach topology and parameters (double precision)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachSlope',    RPARAM(iSeg)%R_SLOPE,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'reachLength',   RPARAM(iSeg)%RLENGTH,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'upstreamArea',  RPARAM(iSeg)%UPSAREA,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'basinArea',     RPARAM(iSeg)%BASAREA,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_scl_dvar(trim(ancil_dir)//trim(fname_ntop),'totalArea',     RPARAM(iSeg)%TOTAREA,iSeg,ierr,cmessage); call handle_err(ierr,cmessage)
  enddo
  ! Populate sAll dimensioned variable 
  ! NETOPO%RCHLIST - upstream reach list
  ! NETOPO%UPSLENG - total upstream reach length 
  nTotal=0
  do iSeg=1,nSeg
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachStart', iRchStart1, iSeg, ierr, cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'reachCount', nRchCount1, iSeg, ierr, cmessage); call handle_err(ierr,cmessage)
    allocate(NETOPO(iSeg)%UPSLENG(nRchCount1), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%UPSLENG')
    allocate(NETOPO(iSeg)%RCHLIST(nRchCount1), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%RCHLIST')
  
    call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'upReachTotalLength',NETOPO(iSeg)%UPSLENG(:),iRchStart1,nRchCount1,ierr,cmessage); call handle_err(ierr,cmessage)
    call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop),'reachList',NETOPO(iSeg)%RCHLIST(:),iRchStart1,nRchCount1,ierr,cmessage); call handle_err(ierr,cmessage)
    nTotal = nTotal + nRchCount1
  enddo 
  ! Populate sUps dimensioned variable 
  ! NETOPO%UREACHI - Immediate upstream reach index list
  ! NETOPO%UREACHK - Immediate upstream reach ID list 
  do iSeg=1,nSeg
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'upReachStart', iUpRchStart, iSeg, ierr, cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'upReachCount', nUpRchCount, iSeg, ierr, cmessage); call handle_err(ierr,cmessage)

    allocate(NETOPO(iSeg)%UREACHI(nUpRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%UREACHI')
    allocate(NETOPO(iSeg)%UREACHK(nUpRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%UREACHK')
    allocate(NETOPO(iSeg)%goodBas(nUpRchCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for NETOPO%goodBas')
    if (nUpRchCount > 0) then
      call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop),'upReachID'   ,NETOPO(iSeg)%UREACHK(:),iUpRchStart,nUpRchCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop),'upReachIndex',NETOPO(iSeg)%UREACHI(:),iUpRchStart,nUpRchCount,ierr,cmessage); call handle_err(ierr,cmessage)
      do jSeg=1,nUpRchCount
        ! check that we identify the correct upstream reach
        if (NETOPO(NETOPO(iSeg)%UREACHI(jSeg))%REACHID /= NETOPO(iSeg)%UREACHK(jSeg)) then 
          print*,'iSeg = ', iSeg 
          print*,'NETOPO(iSeg)%UREACHK(jSeg) = ', NETOPO(iSeg)%UREACHK(jSeg)
          print*,'NETOPO(NETOPO(iSeg)%UREACHI(jSeg))%REACHID = ', NETOPO(NETOPO(iSeg)%UREACHI(jSeg))%REACHID
          call handle_err(20,'unable to find desired immediate upstream segment')
        endif
        ! check that the upstream reach has a basin area > 0
        if(RPARAM(NETOPO(iSeg)%UREACHI(jSeg))%TOTAREA > verySmall)then
         NETOPO(iSeg)%goodBas(jSeg) = .true.
        else
         NETOPO(iSeg)%goodBas(jSeg) = .false.
        endif
      enddo ! looping through the immediate upstream reaches
    endif  ! if not a headwater
  enddo 
  ! Populate sHrus dimensioned variable 
  allocate(h2b(nSeg), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for h2b')
  do iSeg=1,nSeg
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'upHruStart', iUpHruStart, iSeg, ierr, cmessage); call handle_err(ierr,cmessage)
    call get_scl_ivar(trim(ancil_dir)//trim(fname_ntop),'upHruCount', nUpHruCount, iSeg, ierr, cmessage); call handle_err(ierr,cmessage)

    allocate(h2b(iSeg)%cHRU(nUpHruCount), stat=ierr); if(ierr/=0) call handle_err(ierr,'problem allocating space for h2b(iSeg)%cHRU(nUpHruCount)')
    if (nUpHrucount /= 0) then
      call get_vec_ivar(trim(ancil_dir)//trim(fname_ntop),'hru_id',    h2b(iSeg)%cHRU(:)%hru_id,  iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_lon',   h2b(iSeg)%cHRU(:)%hru_lon, iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_lat',   h2b(iSeg)%cHRU(:)%hru_lat, iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_elev',  h2b(iSeg)%cHRU(:)%hru_elev,iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_area',  h2b(iSeg)%cHRU(:)%hru_area,iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
      call get_vec_dvar(trim(ancil_dir)//trim(fname_ntop),'hru_weight',h2b(iSeg)%cHRU(:)%wght,    iUpHruStart,nUpHruCount,ierr,cmessage); call handle_err(ierr,cmessage)
    endif
  enddo
  nSegRoute => nSeg
endif ! outlet segment choice

! compute the time-delay histogram (to route runoff within basins)
scale_time = 1._dp*dt
call qtimedelay(dt,scale_time,ierr,cmessage)
call handle_err(ierr, cmessage)

if (doKWaveRoute) then
  ! define processing order of the reaches
  call reachorder(nSegRoute, ierr, cmessage); call handle_err(ierr, cmessage)
  !! check
  !do iRch=1,nSegRoute
  ! jRch = NETOPO(iRch)%RHORDER
  ! write(*,'(a,1x,2(i4,1x),f20.2)') 'iRch, NETOPO(jRch)%DREACHI, RPARAM(jRch)%TOTAREA/1000000._dp  = ', iRch, NETOPO(jRch)%DREACHI, RPARAM(jRch)%TOTAREA/1000000._dp
  !end do 
  !pause ' check reachorder'
end if

! identify the stream segment with the largest upstream area
iDesire = maxLoc(RPARAM(:)%TOTAREA)
ixDesire= iDesire(1)
print*, 'maximum upstream area = ', RPARAM(ixDesire)%TOTAREA, size(NETOPO(ixDesire)%RCHLIST)

! set the downstream index of the outlet reach to negative (the outlet reach does not flow into anything)
NETOPO(ixDesire)%DREACHI = -9999

if (doVicRoute) then
  ! For VIC routing scheme
  ! Compute total length from the bottom of segment to each upstream segment
  !  call upstrm_length(nSeg, ierr, cmessage);
  ! Compute total length from the bottom of segment to each upstream segment
  call make_uh(nSegRoute, dt, velo, diff, ierr, cmessage); call handle_err(ierr, cmessage)
  !check
 ! !do iSeg=1949,1949
 ! do iSeg=9,9
 !   nUpstream = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment 
 !   do iUps=1,nUpstream
 !     jUps=NETOPO(iSeg)%RCHLIST(iUps)
 !     nTDH= size(NETOPO(iSeg)%UH(iUps)%UH_DATA) ! size of UH data 
 !     write(*,'(a,1x,i4,1x,i4.1,1x,f10.2)') 'upstrm index, upstrmID,length = ', NETOPO(iSeg)%RCHLIST(iUps), NETOPO(jUps)%REACHID, NETOPO(iSeg)%UPSLENG(iUps)
 !     write(*,'(96f6.3)') (NETOPO(iSeg)%UH(iUps)%UH_DATA(iTDH), iTDH=1,nTDH)
 !   enddo
 ! enddo
end if

! *****
! (2) Read in metadata for the runoff file...
! *******************************************
! read in dimensions for the runoff file
call get_qDims(trim(input_dir)//trim(fname_qsim), &  ! input: filename
               vname_hruid,                       &  ! input: name of coordinate dimension HRUid
               vname_time,                        &  ! input: name of coordinate dimension time
               units_time,                        &  ! output: time units
               nTime,                             &  ! output: number of time elements
               nHRU_data,                         &  ! output: number of HRUs in the runoff data file
               ierr, cmessage)                       ! output: error control
call handle_err(ierr, cmessage)

! allocate space for the HRUs
allocate(qsimHRUid(nHRU_data), stat=ierr)
if(ierr/=0) call handle_err(ierr,'problem allocating space for qsimHRUid')
allocate(qsimHRUid_mask(nHRU_data), stat=ierr)
if(ierr/=0) call handle_err(ierr,'problem allocating space for qsimHRUid_mask')
qsimHRUid_mask(:)=.false.
allocate(qsimHRUarea(nHRU_data), stat=ierr)
if(ierr/=0) call handle_err(ierr,'problem allocating space for qsimHRUarea')
qsimHRUarea(:)=-999

! read in metadata for the runoff file
call get_qMeta(trim(input_dir)//trim(fname_qsim), &  ! input: filename
               vname_hruid,                       &  ! input: name of coordinate dimension HRUid
               qsimHRUid,                         &  ! output: HRUid in the simulations
               ierr, cmessage)                       ! output: error control
call handle_err(ierr, cmessage)
print*,'size(qsimHRUid) = ', size(qsimHRUid)

! check all the upstream hrus at the desired outlet exist in runoff file 
! Assign hru_ix based on order of hru in runoff file 
do iSeg=1,nSegRoute ! (loop through stream segments)
  nDrain = size(h2b(iSeg)%cHRU)  ! (number of HRUs that drain into a given stream segment)
  if(nDrain > 0)then
  do iHRU=1,nDrain ! (loop through HRUs that drain into a given stream segment)
    ! check existence of upstream hrus in funoff file
    if( minval(abs(qsimHRUid - h2b(iSeg)%cHRU(iHRU)%hru_id)) /= 0 )then
      write (str,'(I10)') h2b(iSeg)%cHRU(iHRU)%hru_id
      call handle_err(20,'runoff file does not include runoff time series at HRU'//trim(str))
    end if  
    ! Assign hru index  
    iSelect = minloc(abs(qsimHRUid - h2b(iSeg)%cHRU(iHRU)%hru_id))
    h2b(iSeg)%cHRU(iHRU)%hru_ix=iSelect(1)
    if(h2b(iSeg)%cHRU(iHRU)%hru_id /= qsimHRUid(h2b(iSeg)%cHRU(iHRU)%hru_ix))then
      call handle_err(20,'mis-match in HRUs')
    endif
    qsimHRUid_mask(iSelect(1)) = .true. 
    qsimHRUarea(iSelect(1)) = h2b(iSeg)%cHRU(iHRU)%hru_area
  end do  ! (loop through HRUs that drain into a given stream segment) 
 endif  ! (if HRUs drain into the stream segment)
end do ! (loop through stream segments)

! *****
! (3) Define NetCDF output file and write ancillary data...
! *********************************************************
! create NetCDF file
call defineFile(trim(output_dir)//trim(fname_output),  &  ! input: file name
                nSegRoute,                             &  ! input: number of stream segments
                nTotal,                                &  ! input: total number of upstream reaches for all reaches
                units_time,                            &  ! input: time units
                ierr, cmessage)                           ! output: error control
call handle_err(ierr, cmessage)

! write network toplogy (input = filename, variable name, variable vector, start index; output = error control)
call write_iVec(trim(output_dir)//trim(fname_output), 'reachID',    NETOPO(:)%REACHID, 1,  ierr, cmessage); call handle_err(ierr,cmessage)
call write_iVec(trim(output_dir)//trim(fname_output), 'reachOrder', NETOPO(:)%RHORDER, 1,  ierr, cmessage); call handle_err(ierr,cmessage)

iStart=1  ! initialize the start index of the ragged array
! write list of reaches upstream of each reach (ragged array)
do iSeg=1,nSegRoute
 ! get the number of reaches
 nUpstream = size(NETOPO(iSeg)%RCHLIST)
 ! write the vector to the ragged array
 call write_iVec(trim(output_dir)//trim(fname_output), 'reachList', NETOPO(iSeg)%RCHLIST(:), iStart, ierr, cmessage)
 call handle_err(ierr,cmessage)
 ! write the start index and the count (NOTE: pass as a vector)
 call write_iVec(trim(output_dir)//trim(fname_output), 'listStart', (/iStart/),    iSeg, ierr, cmessage); call handle_err(ierr,cmessage)
 call write_iVec(trim(output_dir)//trim(fname_output), 'listCount', (/nUpstream/), iSeg, ierr, cmessage); call handle_err(ierr,cmessage)
 ! update the start index
 iStart = iStart + nUpstream
end do

! write reach parameters
call write_dVec(trim(output_dir)//trim(fname_output), 'basinArea',    RPARAM(:)%BASAREA, (/1/), (/nSegRoute/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_dVec(trim(output_dir)//trim(fname_output), 'upstreamArea', RPARAM(:)%TOTAREA, (/1/), (/nSegRoute/), ierr, cmessage); call handle_err(ierr,cmessage)


! *****
! (4) Prepare for the routing simulations...
! *******************************************
! allocate space for the simulated runoff at the HRUs
allocate(qsim_hru(nHRU_data), stat=ierr)
if(ierr/=0) call handle_err(ierr,'problem allocating space for simulated runoff at the HRUs')

! allocate space for the simulated runoff at basins and reaches
allocate(qsim_basin(nSegRoute),RCHFLX(nens,nSegRoute), KROUTE(nens,nSegRoute), stat=ierr)
if(ierr/=0) call handle_err(ierr,'problem allocating space for simulated runoff at the basins')

! initialize the routed elements
RCHFLX(:,:)%BASIN_QR(1) = 0._dp
RCHFLX(:,:)%BASIN_QR_VIC(1) = 0._dp 

! initialize the time-delay histogram
do iens=1,nens
 do ibas=1,nSegRoute
  ! allocate space for the delayed runoff
  allocate(RCHFLX(iens,ibas)%QFUTURE(size(FRAC_FUTURE)), stat=ierr)
  call handle_err(ierr, 'problem allocating space for QFUTURE element')
  ! initialize to zeroes
  RCHFLX(iens,ibas)%QFUTURE(:) = 0._dp
  
  if (doVicRoute) then
    ! allocate space for the delayed runoff
    nUpstream = size(NETOPO(ibas)%RCHLIST) ! size of upstream segment 
    nUH_DATA_MAX=0
    upstrms_loop: do iUps=1,nUpstream
      nUH_DATA_MAX= max(nUH_DATA_MAX ,size(NETOPO(ibas)%UH(iUps)%UH_DATA))
    enddo upstrms_loop 
    allocate(RCHFLX(iens,ibas)%QFUTURE_VIC(nUH_DATA_MAX), stat=ierr)
    call handle_err(ierr, 'problem allocating space for QFUTURE_VIC element')
    ! initialize to zeroes
    RCHFLX(iens,ibas)%QFUTURE_VIC(:) = 0._dp
  end if

 end do
end do

! define flags
LAKEFLAG=0  ! no lakes in the river network

! define time
T0 = 0._dp
T1 = dt

! *****
! (5) Perform the routing...
! **************************
! loop through time
do iTime=1,nTime

 ! loop through ensemble members
 do iens=1,nens

  ! *****
  ! (5a) Get the simulated runoff for the current time step...
  ! **********************************************************
  ! get the simulated runoff for the current time step
  call getVarQsim(trim(input_dir)//trim(fname_qsim), & ! input: filename
                  vname_time,                        & ! input: name of coordinate dimension time
                  vname_qsim,                        & ! input: name of runoff variable
                  iTime,                             & ! input: time index
                  dTime,                             & ! output: time
                  qsim_hru,                          & ! output: simulated runoff
                  ierr, cmessage)                      ! output: error control
  call handle_err(ierr, cmessage)
  ! ensure that simulated runoff is non-zero
  where(qsim_hru < runoffMin) qsim_hru=runoffMin

  ! write time -- note time is just carried across from the input
  call write_dVec(trim(output_dir)//trim(fname_output), 'time', (/dTime/), (/iTime/), (/1/), ierr, cmessage)
  call handle_err(ierr,cmessage)

  ! *****
  ! (5b) Interpolate simulated runoff to local basins...
  ! ****************************************************
  ! interpolate the data to the basins
  do ibas=1,nSegRoute
   ! intialize the basin runoff
   qsim_basin(ibas) = 0._dp
   ! get the number of HRUs that drain into the basin
   nDrain = size(h2b(ibas)%cHRU)
   ! * case where HRUs drain into the segment
   if(nDrain > 0)then
    ! loop through the HRUs
    do iHRU=1,nDrain
     ! get the HRU index
     ix = h2b(ibas)%cHRU(iHRU)%hru_ix

     !Error check - runoff depth cannot be negative (no missing value)
     if( qsim_hru(ix) < 0._dp )then
       write (str,'(I10)') h2b(ibas)%cHRU(iHRU)%hru_id
       call handle_err(20,'negaive runoff value is invalid')
     endif
     
     ! compute the weighted average
     qsim_basin(ibas) = qsim_basin(ibas) + h2b(ibas)%cHRU(iHRU)%wght*qsim_hru(ix)*time_conv*length_conv  ! ensure m/s
    end do  ! (looping through gridpoints associated with each basin)
   ! * special case where no HRUs drain into the segment
   else
    qsim_basin(ibas) = 0._dp
   endif
   ! convert runoff to m3/s
   RCHFLX(iens,ibas)%BASIN_QI = qsim_basin(ibas)*RPARAM(ibas)%BASAREA
  end do  ! (looping through basins)
  !print*,'RCHFLX(iens,:)%BASIN_QI = ',RCHFLX(iens,:)%BASIN_QI

  ! write instantaneous local runoff in each stream segment (m3/s)
  call write_dVec(trim(output_dir)//trim(fname_output), 'instBasinRunoff', RCHFLX(iens,:)%BASIN_QI, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
  call handle_err(ierr,cmessage)
  
  ! *****
  ! (5c) Delay runoff within local basins...
  ! ****************************************
  ! identify the number of future time steps for a given basin
  ntdh = size(FRAC_FUTURE)

  ! route streamflow through the basin
  do ibas=1,nSegRoute ! place a fraction of runoff in future time steps
    do JTIM=1,NTDH
     RCHFLX(iens,ibas)%QFUTURE(JTIM) = RCHFLX(iens,ibas)%QFUTURE(JTIM) + FRAC_FUTURE(JTIM)*RCHFLX(iens,ibas)%BASIN_QI
    end do
    ! save the routed runoff
    RCHFLX(iens,ibas)%BASIN_QR(0) = RCHFLX(iens,ibas)%BASIN_QR(1)  ! (save the runoff from the previous time step)
    RCHFLX(iens,ibas)%BASIN_QR(1) = RCHFLX(iens,ibas)%QFUTURE(1)
    ! move array back
    do JTIM=2,NTDH
     RCHFLX(iens,ibas)%QFUTURE(JTIM-1) = RCHFLX(iens,ibas)%QFUTURE(JTIM)
    end do 
    RCHFLX(iens,ibas)%QFUTURE(NTDH)    = 0._dp
  end do  ! (looping through basins)

  ! write routed local runoff in each stream segment (m3/s)
  call write_dVec(trim(output_dir)//trim(fname_output), 'dlayBasinRunoff', RCHFLX(iens,:)%BASIN_QR(1), (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
  call handle_err(ierr,cmessage)
 
  ! *****
  ! (5d) Route streamflow for each upstream segment through the river network...
  ! **************************************************
  if (doVicRoute) then
    ! Get upstreme routed runoff depth at top of each segment for VIC routing 
    call get_upsbas_qr(nSegRoute,iens,ierr,cmessage)
    ! write routed runoff (m3/s)
    call write_dVec(trim(output_dir)//trim(fname_output), 'UpBasRoutedRunoff', RCHFLX(iens,:)%UPSBASIN_QR, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    call handle_err(ierr,cmessage)
    ! Get upstreme routed runoff depth at top of each segment for VIC routing 
    call conv_upsbas_qr(nSegRoute,iens,ierr,cmessage)
    ! write routed runoff (m3/s)
    call write_dVec(trim(output_dir)//trim(fname_output), 'VICroutedRunoff', RCHFLX(iens,:)%REACH_Q_VIC, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    call handle_err(ierr,cmessage)
  endif

  ! *****
  ! (5e) Compute total instantaneous runoff from all upstream basins...
  ! *******************************************************************
  ! compute the sum of all upstream runoff at each point in the river network
  do iSeg=1,nSegRoute
    ! identify how many reaches are upstream
    nUpstream = size(NETOPO(iSeg)%RCHLIST)
    ! allocate space for upstream vectors
    allocate(iUpstream(nUpstream), qUpstream(nUpstream), stat=ierr)
    if(ierr/=0) call handle_err(ierr,'problem allocating vectors for all upstream basins')
    ! get indices for all reaches upstream
    iUpstream = NETOPO(iSeg)%RCHLIST(1:nUpstream)
    ! get streamflow for all reaches upstream
    qUpstream = RCHFLX(iens,iUpstream(1:nUpstream))%BASIN_QR(1)
    ! get mean streamflow 
    RCHFLX(IENS,iSeg)%UPSTREAM_QI = sum(qUpstream)
    ! test
    if(NETOPO(iSeg)%REACHID == ixPrint)then
     print*, 'ixUpstream = ', NETOPO(iUpstream(1:nUpstream))%REACHIX
     print*, 'idUpstream = ', NETOPO(iUpstream(1:nUpstream))%REACHID
     print*, 'qUpstream = ', qUpstream
    endif
    ! deallocate space for upstream vectors
    deallocate(iUpstream,qUpstream, stat=ierr)
    if(ierr/=0) call handle_err(ierr,'problem deallocating vectors for all upstream basins')
  end do  ! looping through stream segments

  ! write instantaneous runoff from all upstream basins (m3/s)
  call write_dVec(trim(output_dir)//trim(fname_output), 'sumUpstreamRunoff', RCHFLX(iens,:)%UPSTREAM_QI, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
  call handle_err(ierr,cmessage)

  ! *****
  ! (5f) Route streamflow through the river network...
  ! **************************************************
  if (doKWaveRoute) then
    ! specify some additional parameters (temporary "fix")
    RPARAM(:)%R_WIDTH =  width_scale * sqrt(RPARAM(:)%TOTAREA)    ! channel width (m)
    RPARAM(:)%R_MAN_N =  mann_n  ! Manning's "n" paramater (unitless)

    ! route streamflow through the river network
    do iSeg=1,nSegRoute
      ! identify reach to process
      irch = NETOPO(iSeg)%RHORDER
      !print*, 'irch, ixDesire = ', irch, ixDesire

      ! route kinematic waves through the river network
      CALL QROUTE_RCH(IENS,irch,    & ! input: array indices
                      ixDesire,     & ! input: index of the outlet reach
                      T0,T1,        & ! input: start and end of the time step
                      LAKEFLAG,     & ! input: flag if lakes are to be processed
                      ierr,cmessage)  ! output: error control
      call handle_err(ierr,cmessage)
      !if(iRch==5) pause 'finished stream segment'
    end do  ! (looping through stream segments)
    !pause 'first time step'
    
    ! write routed runoff (m3/s)
    call write_dVec(trim(output_dir)//trim(fname_output), 'routedRunoff', RCHFLX(iens,:)%REACH_Q, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    call handle_err(ierr,cmessage)
  end if

  ! *****
  ! (5g) print progress...
  ! **********************
  ! extract desired variables
  !qDomain_hru     = sum(qsim_hru(:)*nhru_acil(ixHRU%area)%varData(:))/sum(nhru_acil(ixHRU%area)%varData(:)) ! domain: total runoff at the HRUs (mm/s)
  qDomain_hru     = sum(pack(qsim_hru,qsimHRUid_mask)*pack(qsimHRUarea,qsimHRUid_mask))/sum(pack(qsimHRUarea,qsimHRUid_mask)) ! domain: total runoff at the HRUs (mm/s)
  qDomain_basin   = 1000._dp*sum(qsim_basin*RPARAM(:)%BASAREA) / sum(RPARAM(:)%BASAREA)                  ! domain: total instantaneous basin runoff (mm/s)
  qDomain_reach   = 1000._dp*sum(RCHFLX(iens,:)%BASIN_QI)/sum(RPARAM(:)%BASAREA)                            ! domain: total instantaneous reach runoff (mm/s)
  qDesire_instant = 1000._dp*RCHFLX(iens,ixDesire)%UPSTREAM_QI/RPARAM(ixDesire)%TOTAREA                     ! desire: sum of instantaneous runoff for all basins upstream of the desired reach (mm/s)
  qDesire_routed  = 1000._dp*RCHFLX(iens,ixDesire)%REACH_Q/RPARAM(ixDesire)%TOTAREA                         ! desire: routed runoff for the desired reach (mm/s)

  ! print progress
  write(*,'(i5,1x,10(e20.10,1x))') iTime, qDomain_hru, qDomain_basin, qDomain_reach, qDesire_instant, qDesire_routed

 end do  ! (looping through ensemble members)

 ! increment time bounds
 T0 = T0 + dt
 T1 = T0 + dt

 !pause 'finished time step'
end do  ! (looping through time)

stop

contains

 subroutine handle_err(err,message)
 ! handle error codes
 implicit none
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 if(err/=0)then
  print*,'FATAL ERROR: '//trim(message)
  stop
 endif
 end subroutine handle_err


 subroutine findix(scl,vec,iSelect)
  ! Find vec index where the value match up with scl  
  implicit none
      
  integer(i4b),intent(in)              :: scl 
  integer(i4b),intent(in),allocatable  :: vec(:)
  integer(i4b),intent(out)             :: iSelect
 
  integer(i4b)                         :: i(1)
  i = minloc(abs(vec - scl))
  iSelect = i(1)  ! de-vectorize the desired stream segment
  if(vec(iSelect) /= scl)&
   call handle_err(20,'unable to find matched value')
 end subroutine findix 

end

