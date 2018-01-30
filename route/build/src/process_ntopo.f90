module process_ntopo

! provide access to desired modules
USE nrtype                                    ! variable types, etc.
USE public_var

private
public::ntopo

contains

 ! *********************************************************************
 ! public subroutine: read and process river network data
 ! *********************************************************************
 subroutine ntopo(ierr, message)

 ! This subroutine 1) read river network data and 2) populate river network topology data strucutres
 ! Data structure populated
 ! 1. h2b
 ! 2. NETOPO
 ! 3. RPARAM

 use read_streamSeg, only:getData               ! get the ancillary data
 use read_streamSeg, only:hru2basin             ! get the mapping between HRUs and basins
 use read_streamSeg, only:assign_reachparam     ! assign reach parameters
 use network_topo,   only:reach_list            ! identify all reaches upstream of the each reach
 use network_topo,   only:upstrm_length         ! Compute total upstream length  NM
 use reachparam                                 ! reach parameters
 use nhru2basin                                 ! data structures holding the nhru2basin correspondence

 implicit none
 ! Output variables
 integer(i4b), intent(out)      :: ierr           ! error code
 character(*), intent(out)      :: message        ! error message
 ! Local variables
 character(len=strLen)          :: cmessage       ! error message of downwind routine
 integer(i4b)                   :: nHRU           ! number of HRUs
 integer(i4b)                   :: nSeg           ! number of stream segments
 integer(i4b)                   :: tot_all_upseg  ! total number of all the upstream segments for all stream segments
 integer(i4b)                   :: tot_upseg      ! total number of immediate upstream segments for all  stream segments
 integer(i4b)                   :: tot_hru        ! total number of all the upstream hrus for all stream segments
 integer(i4b)                   :: n_all_upseg    ! number of all pstream reaches of each stream segment
 integer(i4b)                   :: n_upseg        ! number of reaches immediate upstream of each stream segment
 integer(i4b)                   :: n_uphru        ! number of HRH draining to each stream segment
 integer(i4b)                   :: iSeg           ! index of stream segment
 integer(i4b)                   :: iUps           ! index of upstream stream segment
 !integer(i4b)                   :: jUps           ! index of upstream stream segment added by NM
 !integer(i4b)                   :: iRch!,jRch     ! index in reach structures
 integer*8                      :: time0,time1    ! times

 ! initialize error control
 ierr=0; message='ntopo/'

 ! initialize times
 call system_clock(time0)

 ! *****
 ! (1) Read in the stream segment information...
 ! *********************************************
 ! get the number of HRUs and stream segments (needed for allocate statements)
 call getData(trim(ancil_dir)//trim(fname_ntop), & ! input: file name
              dname_nhru, &  ! input: dimension name of the HRUs
              dname_sseg, &  ! input: dimension name of the stream segments
              nhru_acil,  &  ! input-output: ancillary data for HRUs
              sseg_acil,  &  ! input-output: ancillary data for stream segments
              imap_acil,  &  ! input-output: ancillary data for mapping hru2basin
              ntop_acil,  &  ! input-output: ancillary data for network topology
              nHRU,       &  ! output: number of HRUs
              nSeg,       &  ! output: number of stream segments
              ierr,cmessage) ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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

 ! get timing
 call system_clock(time1)
 write(*,'(a,1x,i20)') 'after hru2basin: time = ', time1-time0

 ! put data in structures
 call assign_reachparam(nSeg,         & ! input: number of stream segments
                        sseg_acil,    & ! input: ancillary data for stream segments
                        ntop_acil,    & ! input: ancillary data for network topology
                        tot_upseg,    & ! output: sum of number of immediate upstream reaches
                        ierr, cmessage) ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 write(*,'(a,1x,i20)') 'after assign_reachparam: time = ', time1-time0
 print*, 'PAUSE: '; read(*,*)

 !stop

 ! identify all reaches upstream of each reach
 call reach_list(nSeg, tot_all_upseg, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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
 call upstrm_length(nSeg, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 !do iSeg=9,9
 !  n_upseg = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment
 !  do iUps=1,n_upseg
 !    jUps=NETOPO(iSeg)%RCHLIST(iUps)
 !    write(*,'(a,1x,i4,1x,i10,1x,f10.2)') 'upstrm index, upstrmID,length = ', NETOPO(iSeg)%RCHLIST(iUps), NETOPO(jUps)%REACHID, NETOPO(iSeg)%UPSLENG(iUps)
 !  enddo
 !enddo
 end subroutine

end module process_ntopo
