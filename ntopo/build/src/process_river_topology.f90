program process_river_topology

! provide access to desired modules
USE nrtype                                    ! variable types, etc.
USE dataTypes,only:namepvar,nameivar          ! provide access to data types
USE var_lookup,only:ixHRU,nVarsHRU            ! index of variables for the HRUs
USE var_lookup,only:ixSEG,nVarsSEG            ! index of variables for the stream segments
USE var_lookup,only:ixMAP,nVarsMAP            ! index of variables for the hru2basin mapping
USE var_lookup,only:ixTOP,nVarsTOP            ! index of variables for the network topology
USE reachparam                                ! reach parameters
USE nhru2basin                                ! data structures holding the nhru2basin correspondence
USE nr_utility_module,only:arth               ! use to build vectors with regular increments
USE ascii_util_module,only:file_open          ! open file (performs a few checks as well)
USE ascii_util_module,only:get_vlines         ! get a list of character strings from non-comment lines
USE read_streamSeg,only:getData               ! get the ancillary data
USE read_streamSeg,only:hru2basin             ! get the mapping between HRUs and basins
USE read_streamSeg,only:assign_reachparam     ! assign reach parameters
USE write_ntopo,only:defineFile               ! define output file
USE write_netcdf,only:write_nc                ! write netcdf
USE network_topo,only:reach_list              ! identify all reaches upstream of the each reach
USE network_topo,only:upstrm_length           ! Compute total upstream length  NM

! define variables
implicit none
integer(i4b),parameter     :: missingInt  = -9999   ! missing value (integer)
real(dp),parameter         :: missingReal = -9999.0 ! missing value (real)
! general guff
integer(i4b),parameter     :: strLen=256         ! length of character string
integer(i4b)               :: ierr               ! error code
character(len=strLen)      :: cmessage           ! error message of downwind routine
! read control file
character(len=strLen)      :: cfile_name      ! name of the control file
integer(i4b)               :: iunit           ! file unit
character(len=strLen),allocatable :: cLines(:) ! vector of character strings
integer(i4b)               :: iLine               ! index of line in cLines
integer(i4b)               :: ibeg_name           ! start index of variable name in string cLines(iLine)
integer(i4b)               :: iend_name           ! end index of variable name in string cLines(iLine)
integer(i4b)               :: iend_data           ! end index of data in string cLines(iLine)
character(len=strLen)      :: cName,cData         ! name and data from cLines(iLine)
! define output file
character(len=strLen)      :: fname_output        ! name of output file
integer(i4b),parameter     :: outunit=31          ! unit for output file
! define directories
character(len=strLen)      :: input_dir           ! directory containing input data
character(len=strLen)      :: output_dir          ! directory containing output data
! define stream segment information
character(len=strLen)      :: fname_sseg          ! filename containing stream segment information
character(len=strLen)      :: dname_sseg          ! dimension name of the stream segments
character(len=strLen)      :: dname_nhru          ! dimension name of the HRUs
integer(i4b)               :: nHRU                ! number of HRUs
integer(i4b)               :: nSeg                ! number of stream segments
integer(i4b)               :: tot_all_upseg       ! total number of all the upstream segments for all stream segments
integer(i4b)               :: tot_upseg           ! total number of immediate upstream segments for all  stream segments
integer(i4b)               :: tot_hru             ! total number of all the upstream hrus for all stream segments
integer(i4b)               :: n_all_upseg         ! number of all pstream reaches of each stream segment
integer(i4b)               :: n_upseg             ! number of reaches immediate upstream of each stream segment
integer(i4b)               :: n_uphru             ! number of HRH draining to each stream segment
type(namepvar)             :: nhru_acil(nVarsHRU) ! ancillary data for HRUs
type(namepvar)             :: sseg_acil(nVarsSEG) ! ancillary data for stream segments
type(nameivar)             :: imap_acil(nVarsMAP) ! ancillary data for mapping hru2basin
type(nameivar)             :: ntop_acil(nVarsTOP) ! ancillary data for network toopology
integer(i4b)               :: iSeg                ! index of stream segment
integer(i4b)               :: iUps                ! index of upstream stream segment
!integer(i4b)               :: jUps               ! index of upstream stream segment added by NM
integer(i4b)               :: iStartAll           ! start index of the ragged array for all upstream reach vector
integer(i4b)               :: iStartUps           ! start index of the ragged array for immediate upstream vector
integer(i4b)               :: iStartHru           ! start index of the ragged array for immediate upstream Hru vector
!integer(i4b)               :: iRch!,jRch       ! index in reach structures
integer*8                  :: time0,time1         ! times
integer(i4b),allocatable   :: reachStart(:)
integer(i4b),allocatable   :: reachCount(:)
integer(i4b),allocatable   :: upReachStart(:)
integer(i4b),allocatable   :: upReachCount(:)
integer(i4b),allocatable   :: upHruStart(:)
integer(i4b),allocatable   :: upHruCount(:)
integer(i4b),allocatable   :: reachList(:)
real(dp),allocatable   :: upReachTotalLength(:)
integer(i4b),allocatable   :: upReachIndex(:)
integer(i4b),allocatable   :: upReachID(:)
integer(i4b),allocatable   :: hruIndex(:)
integer(i4b),allocatable   :: hru_id(:)
real(dp),allocatable   :: hru_lat(:)
real(dp),allocatable   :: hru_lon(:)
real(dp),allocatable   :: hru_elev(:)
real(dp),allocatable   :: hru_area(:)
real(dp),allocatable   :: hru_weight(:)

! *****
! (0) Read control file...
! ************************

! initialize times
call system_clock(time0)

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
  case('<input_dir>');    input_dir    = trim(cData)           ! directory containing input data
  case('<output_dir>');   output_dir   = trim(cData)           ! directory containing output data
  ! define the file and the dimensions
  case('<fname_sseg>');   fname_sseg   = trim(cData)           ! name of file containing stream segment information
  case('<dname_nhru>');   dname_nhru   = trim(cData)           ! dimension name of the HRUs
  case('<dname_sseg>');   dname_sseg   = trim(cData)           ! dimension name of the stream segments
  ! define the variables desired for each HRU
  case('<HRUname_xLon>');      nhru_acil(ixHRU%xLon     )%varName = trim(cData) ! name of variable holding HRU longitude
  case('<HRUname_yLat>');      nhru_acil(ixHRU%yLat     )%varName = trim(cData) ! name of variable holding HRU latitude
  case('<HRUname_elev>');      nhru_acil(ixHRU%elev     )%varName = trim(cData) ! name of variable holding HRU elevation
  case('<HRUname_area>');      nhru_acil(ixHRU%area     )%varName = trim(cData) ! name of variable holding HRU area
  ! define the variables desired for each stream segment
  case('<SEGname_bLon>'  );    sseg_acil(ixSEG%bLon     )%varName = trim(cData) ! name of variable holding start longitude
  case('<SEGname_bLat>'  );    sseg_acil(ixSEG%bLat     )%varName = trim(cData) ! name of variable holding start latitude
  case('<SEGname_eLon>'  );    sseg_acil(ixSEG%eLon     )%varName = trim(cData) ! name of variable holding end longitude
  case('<SEGname_eLat>'  );    sseg_acil(ixSEG%eLat     )%varName = trim(cData) ! name of variable holding end latitude
  case('<SEGname_length>');    sseg_acil(ixSEG%length   )%varName = trim(cData) ! name of variable holding segment length
  case('<SEGname_slope>' );    sseg_acil(ixSEG%slope    )%varName = trim(cData) ! name of variable holding segment slope
  ! define the variables desired for mapping hru2basin
  case('<HRUname_HRUid>'    ); imap_acil(ixMAP%HRUid    )%varName = trim(cData) ! name of variable holding HRU id
  case('<HRUname_segHRUid>' ); imap_acil(ixMAP%segHRUid )%varName = trim(cData) ! name of variable holding the stream segment below each HRU
  ! define variables for the network topology
  case('<SEGname_segid>'    ); ntop_acil(ixTOP%segid    )%varName = trim(cData) ! name of variable holding the ID of each stream segment
  case('<SEGname_toSegment>'); ntop_acil(ixTOP%toSegment)%varName = trim(cData) ! name of variable holding the ID of the next downstream segment
  ! define the output filename
  case('<fname_output>'); fname_output = trim(cData)    ! filename for the model output
 end select
end do  ! looping through lines in the control file

! *****
! (1) Read in the stream segment information...
! *********************************************
! get the number of HRUs and stream segments (needed for allocate statements)
call getData(trim(input_dir)//trim(fname_sseg), & ! input: file name
             dname_nhru, &  ! input: dimension name of the HRUs
             dname_sseg, &  ! input: dimension name of the stream segments
             nhru_acil,  &  ! input-output: ancillary data for HRUs
             sseg_acil,  &  ! input-output: ancillary data for stream segments
             imap_acil,  &  ! input-output: ancillary data for mapping hru2basin
             ntop_acil,  &  ! input-output: ancillary data for network topology
             nHRU,       &  ! output: number of HRUs
             nSeg,       &  ! output: number of stream segments
             ierr,cmessage) ! output: error control
call handle_err(ierr, cmessage)

! get timing
call system_clock(time1)
write(*,'(a,1x,i20)') 'after getData: time = ', time1-time0

! get the mapping between HRUs and basins
call hru2basin(nHRU,       &   ! input: number of HRUs
               nSeg,       &   ! input: number of stream segments
               nhru_acil,  &   ! input: ancillary data for HRUs
               imap_acil,  &   ! input: ancillary data for mapping hru2basin
               ntop_acil,  &   ! input: ancillary data for network topology
               tot_hru,    &   ! output: total number of all the upstream hrus for all stream segments
               ierr, cmessage) ! output: error control
call handle_err(ierr, cmessage)

! get timing
call system_clock(time1)
write(*,'(a,1x,i20)') 'after hru2basin: time = ', time1-time0

! put data in structures
call assign_reachparam(nSeg,         & ! input: number of stream segments
                       sseg_acil,    & ! input: ancillary data for stream segments
                       ntop_acil,    & ! input: ancillary data for network topology
                       tot_upseg,    & ! output: sum of number of immediate upstream reaches
                       ierr, cmessage) ! output: error control
call handle_err(ierr, cmessage)

! get timing
call system_clock(time1)
write(*,'(a,1x,i20)') 'after assign_reachparam: time = ', time1-time0

!stop

! check
!do iRch=1,nSeg
! jRch = NETOPO(iRch)%RHORDER
! write(*,'(a,1x,2(i4,1x),f20.2)') 'iRch, NETOPO(jRch)%DREACHI, RPARAM(jRch)%TOTAREA/1000000._dp  = ', iRch, NETOPO(jRch)%DREACHI, RPARAM(jRch)%TOTAREA/1000000._dp
!end do

! identify all reaches upstream of each reach
call reach_list(nSeg, tot_all_upseg, ierr, cmessage); call handle_err(ierr, cmessage)

! Compute UPSAREA and TOTAREA
do iSeg=1,nSeg
  print*,'--------------------------------------------------------------'
  print*,'NETOPO(iSeg)%REACHID(:)',NETOPO(iSeg)%REACHID
  ! Count how many upstream reaches
  n_all_upseg = size(NETOPO(iSeg)%RCHLIST(:))
  ! Initialize UPSAREA for current segment
  RPARAM(iSeg)%UPSAREA = 0._dp
  do iUps = 1,n_all_upseg
    if (NETOPO(NETOPO(iSeg)%RCHLIST(iUps))%REACHID/=NETOPO(iSeg)%REACHID) then
      print*,'RCHLIST(iUps) = ',NETOPO(NETOPO(iSeg)%RCHLIST(iUps))%REACHID
      RPARAM(iSeg)%UPSAREA = RPARAM(iSeg)%UPSAREA +RPARAM(NETOPO(iSeg)%RCHLIST(iUps))%BASAREA
    else
      print*,'local reach,',NETOPO(NETOPO(iSeg)%RCHLIST(iUps))%REACHID
    endif
  enddo
enddo

! compute total area above the bottom of each reach (m2)
RPARAM(:)%TOTAREA = RPARAM(:)%UPSAREA + RPARAM(:)%BASAREA

!do iSeg=1,nSeg
!  print*,'------------------------------'
!  print*,'NETOPO(iSeg)%REACHID = ',NETOPO(iSeg)%REACHID
!  print*,'NETOPO(iSeg)%UREACHK(:) =',NETOPO(iSeg)%UREACHK(:)
!  print*,'RPARAM(iSeg)%UPSAREA, RPARAM(iSeg)%TOTAREA =',RPARAM(iSeg)%UPSAREA,RPARAM(iSeg)%TOTAREA
!end do

! Compute total length from the bottom of segment to each upstream segment
call upstrm_length(nSeg, ierr, cmessage);

!do iSeg=9,9
!  n_upseg = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment
!  do iUps=1,n_upseg
!    jUps=NETOPO(iSeg)%RCHLIST(iUps)
!    write(*,'(a,1x,i4,1x,i10,1x,f10.2)') 'upstrm index, upstrmID,length = ', NETOPO(iSeg)%RCHLIST(iUps), NETOPO(jUps)%REACHID, NETOPO(iSeg)%UPSLENG(iUps)
!  enddo
!enddo

iStartAll=1  ! initialize the start index of the ragged array
iStartUps=1  ! initialize the start index of the ragged array
iStartHru=1  ! initialize the start index of the ragged array

allocate(reachStart(nSeg), stat=ierr)
allocate(reachCount(nSeg), stat=ierr)
allocate(upReachStart(nSeg), stat=ierr)
allocate(upReachCount(nSeg), stat=ierr)
allocate(upHruStart(nSeg), stat=ierr)
allocate(upHruCount(nSeg), stat=ierr)
allocate(reachList(tot_all_upseg), stat=ierr)
allocate(upReachTotalLength(tot_all_upseg), stat=ierr)
allocate(upReachIndex(tot_upseg), stat=ierr)
allocate(upReachID(tot_upseg), stat=ierr)
allocate(hruIndex(tot_hru), stat=ierr)
allocate(hru_id(tot_hru), stat=ierr)
allocate(hru_lat(tot_hru), stat=ierr)
allocate(hru_lon(tot_hru), stat=ierr)
allocate(hru_elev(tot_hru), stat=ierr)
allocate(hru_area(tot_hru), stat=ierr)
allocate(hru_weight(tot_hru), stat=ierr)
! write attributes of upstream reach, immediate upstream reachs and immediate hrus of each reach (ragged array)
do iSeg=1,nSeg
  ! get the number of all upstream reaches
  n_all_upseg = size(NETOPO(iSeg)%RCHLIST)
  ! get the number of immediate upstream reaches
  n_upseg    = size(NETOPO(iSeg)%UREACHI)
  ! get the number of immediate upstream hru
  n_uphru       = size(h2b(iSeg)%cHRU)

  reachList(iStartAll:iStartAll+n_all_upseg-1) = NETOPO(iSeg)%RCHLIST(:)
  upReachTotalLength(iStartAll:iStartAll+n_all_upseg) = NETOPO(iSeg)%UPSLENG(:)
  reachStart(iSeg) = iStartAll
  reachCount(iSeg) = n_all_upseg
  if ( n_upseg/=0 ) then
    upReachIndex(iStartUps:iStartUps+n_upseg-1) = NETOPO(iSeg)%UREACHI(:)
    upReachID(iStartUps:iStartUps+n_upseg-1) = NETOPO(iSeg)%UREACHK(:)
    upReachStart(iSeg) = iStartUps
    upReachCount(iSeg) = n_upseg
    iStartUps = iStartUps + n_upseg ! update the start index for immediate upstream array
  else
    upReachStart(iSeg) = missingInt
    upReachCount(iSeg) = 0
  endif
  if ( n_uphru/=0 ) then
    hruIndex(istartHru:istartHru+n_uphru-1) =  h2b(iSeg)%cHRU(:)%hru_ix
    hru_id(istartHru:istartHru+n_uphru-1) = h2b(iSeg)%cHRU(:)%hru_id
    hru_lon(istartHru:istartHru+n_uphru-1) = h2b(iSeg)%cHRU(:)%hru_lon
    hru_lat(istartHru:istartHru+n_uphru-1) = h2b(iSeg)%cHRU(:)%hru_lat
    hru_elev(istartHru:istartHru+n_uphru-1) = h2b(iSeg)%cHRU(:)%hru_elev
    hru_area(istartHru:istartHru+n_uphru-1) = h2b(iSeg)%cHRU(:)%hru_area
    hru_weight(istartHru:istartHru+n_uphru-1) = h2b(iSeg)%cHRU(:)%wght
    upHruStart(iSeg) = iStartHru
    upHruCount(iSeg) = n_uphru
    iStartHru = iStartHru + n_uphru ! update the start index for immediate upstream array
  else
    upHruStart(iSeg) = missingInt
    upHruCount(iSeg) = 0
  endif
  ! update the start index for all upstream array
  iStartAll = iStartAll + n_all_upseg
end do

! *****
! (3) Define NetCDF output file and write network topology data...
! *********************************************************
! create NetCDF file
call defineFile(trim(output_dir)//trim(fname_output),  & ! input: file name
                nSeg,                                  & ! input: number of stream segments
                tot_all_upseg,                         & ! input: total number of upstream reaches for all reaches
                tot_upseg,                             & ! input: sum of number of immediate upstream reachs for all reaches
                tot_hru,                              & ! input: sum of number of immediate upstream HRUs for all reaches
                ierr, cmessage)                          ! output: error control
call handle_err(ierr, cmessage)
! write network toplogy (input = filename, variable name, variable vector, start index; output = error control)
call write_nc(trim(output_dir)//trim(fname_output), 'reachID',        NETOPO(:)%REACHID,(/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'reachIndex',     NETOPO(:)%REACHIX,(/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'downReachIndex', NETOPO(:)%DREACHI,(/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'downReachID',    NETOPO(:)%DREACHK,(/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
! write reach parameters
call write_nc(trim(output_dir)//trim(fname_output), 'reachSlope',   RPARAM(:)%R_SLOPE, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'reachLength',  RPARAM(:)%RLENGTH, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'basinArea',    RPARAM(:)%BASAREA, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'upstreamArea', RPARAM(:)%UPSAREA, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'totalArea',    RPARAM(:)%TOTAREA, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'reachLat1',    NETOPO(:)%RCHLAT1, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'reachLat2',    NETOPO(:)%RCHLAT2, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'reachLon1',    NETOPO(:)%RCHLON1, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output), 'reachLon2',    NETOPO(:)%RCHLON2, (/1/), (/nSeg/), ierr, cmessage); call handle_err(ierr,cmessage)
! Write all upstream reach attribute variables
call write_nc(trim(output_dir)//trim(fname_output),'reachList',  reachList, (/1/), (/tot_all_upseg/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'upReachTotalLength',  upReachTotalLength, (/1/), (/tot_all_upseg/), ierr,cmessage); call handle_err(ierr,cmessage)
! Write positional information all upstream reach attribute variables
call write_nc(trim(output_dir)//trim(fname_output),'reachStart',  reachStart, (/1/), (/nSeg/), ierr,cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'reachCount',  reachCount, (/1/), (/nSeg/), ierr,cmessage); call handle_err(ierr,cmessage)
! Write attributes of immediate upstream reaches
call write_nc(trim(output_dir)//trim(fname_output),'upReachIndex', upReachIndex, (/1/), (/tot_upseg/),ierr,cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'upReachID',    upReachID,    (/1/), (/tot_upseg/),ierr,cmessage); call handle_err(ierr,cmessage)
! Write positional information immediate upstream reach attribute variables
call write_nc(trim(output_dir)//trim(fname_output),'upReachStart', upReachStart, (/1/), (/nSeg/), ierr,cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'upReachCount', upReachCount, (/1/), (/nSeg/), ierr,cmessage); call handle_err(ierr,cmessage)
! Write attributes of immediate upstream hrus
call write_nc(trim(output_dir)//trim(fname_output),'hruIndex',  hruIndex,(/1/), (/tot_hru/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'hru_id',    hru_id,  (/1/), (/tot_hru/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'hru_lon',   hru_lon, (/1/), (/tot_hru/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'hru_lat',   hru_lat, (/1/), (/tot_hru/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'hru_elev',  hru_elev,(/1/), (/tot_hru/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'hru_area',  hru_area,(/1/), (/tot_hru/), ierr, cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'hru_weight',hru_weight,(/1/), (/tot_hru/), ierr, cmessage); call handle_err(ierr,cmessage)
! Write positional information hru attribute variables
call write_nc(trim(output_dir)//trim(fname_output),'upHruStart', upHruStart, (/1/), (/nSeg/), ierr,cmessage); call handle_err(ierr,cmessage)
call write_nc(trim(output_dir)//trim(fname_output),'upHruCount', upHruCount, (/1/), (/nSeg/), ierr,cmessage); call handle_err(ierr,cmessage)

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

end
