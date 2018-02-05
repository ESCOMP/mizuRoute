program route_runoff

! ******
! provide access to desired modules
! ************************

! variable types
USE nrtype                                    ! variable types, etc.
USE dataTypes, only : var_ilength             ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength             ! double precision type: var(:)%dat

! global data
USE public_var
!USE globalData, only:time_conv,length_conv    ! conversion factors

! general data structures
!USE globalData, only:remap_data
!USE globalData, only:runoff_data

! metadata
USE var_lookup, only : ixHRU     , nVarsHRU      ! index of variables for data structure
USE var_lookup, only : ixHRU2SEG , nVarsHRU2SEG  ! index of variables for data structure
USE var_lookup, only : ixSEG     , nVarsSEG      ! index of variables for data structure
USE var_lookup, only : ixNTOPO   , nVarsNTOPO    ! index of variables for data structure

! subroutines: populate metadata
USE popMetadat_module, only : popMetadat      ! populate metadata

! subroutines: model control
USE read_control_module, only : read_control  ! read the control file
USE ascii_util_module, only : file_open       ! open file (performs a few checks as well)

! subroutines: get network topology
USE process_ntopo, only: ntopo                ! process the network topology
!USE reach_mask_module, only:reach_mask        ! identify all reaches upstream of a given reach

! ******
! define variables
! ************************
implicit none

! index for printing (set to negative to supress printing
integer(i4b),parameter        :: ixPrint = -9999     ! index for printing

! error control
integer(i4b)                  :: ierr                ! error code
character(len=strLen)         :: cmessage            ! error message of downwind routine

!  network topology data structures
type(var_dlength),allocatable :: structHRU(:)        ! HRU properties
type(var_dlength),allocatable :: structSeg(:)        ! stream segment properties
type(var_ilength),allocatable :: structHRU2seg(:)    ! HRU-to-segment mapping
type(var_ilength),allocatable :: structNTOPO(:)      ! network topology

! read control file
character(len=strLen)         :: cfile_name          ! name of the control file
integer(i4b)                  :: iunit               ! file unit

! define desired reaches
integer(i4b)                  :: nHRU                ! number of HRUs
integer(i4b)                  :: nRch                ! number of desired reaches

! namelist parameters
real(dp)                   :: fshape              ! shape parameter in time delay histogram (=gamma distribution) [-]
real(dp)                   :: tscale              ! scaling factor for the time delay histogram [sec]
real(dp)                   :: velo                ! velocity [m/s] for Saint-Venant equation   added by NM
real(dp)                   :: diff                ! diffusivity [m2/s] for Saint-Venant equation   added by NM
real(dp)                   :: mann_n              ! manning's roughness coefficient [unitless]  added by NM
real(dp)                   :: wscale              ! scaling factor for river width [-] added by NM
namelist /HSLOPE/fshape,tscale  ! route simulated runoff through the local basin
namelist /IRF_UH/velo,diff      ! route delayed runoff through river network with St.Venant UH
namelist /KWT/mann_n,wscale     ! route kinematic waves through the river network

! *****
! *** Populate metadata...
! ************************

! populate the metadata files
call popMetadat(ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! *****
! *** Read control files...
! *************************

! get command-line argument defining the full path to the control file
call getarg(1,cfile_name)
if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')

! read the control file
call read_control(trim(cfile_name), ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

! read the name list
call file_open(trim(param_nml),iunit,ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)
read(iunit, nml=HSLOPE)
read(iunit, nml=IRF_UH)
read(iunit, nml=KWT)
close(iunit)

! *****
! *** Process the network topology...
! ***********************************

! get the network topology
call ntopo(&
           ! output: model control
           nHRU,             & ! number of HRUs
           nRch,             & ! number of stream segments
           ! output: populate data structures
           structHRU,        & ! ancillary data for HRUs
           structSeg,        & ! ancillary data for stream segments
           structHRU2seg,    & ! ancillary data for mapping hru2basin
           structNTOPO,      & ! ancillary data for network toopology
           ! output: error control
           ierr, cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

print*, 'PAUSE: after network topology'; read(*,*)

! *****
! *** Route runoff...
! *******************

! ---------- temporary code: populate data structures ---------------------------------------

 ! Reach Parameters
    !    R_SLOPE
    !    R_MAN_N
    !    R_WIDTH
    !    RLENGTH
    !    UPSAREA  ! upstream area (zero if headwater basin)
    !    BASAREA  ! local basin area
    !    TOTAREA  ! UPSAREA + BASAREA
    !    MINFLOW  ! minimum environmental flow

















    !
    !          ! define flags
    !          LAKEFLAG=0  ! no lakes in the river network
    !
    !          ! define time
    !          T0 = 0._dp
    !          T1 = dt
    !          ! (5) Perform the routing...
    !          ! **************************
    !          ! loop through time
    !          do iTime=1,nTime
    !
    !           ! loop through ensemble members
    !           do iens=1,nens
    !
    !            ! *****
    !            ! (5a) Get the simulated runoff for the current time step...
    !            ! **********************************************************
    !            ! get the simulated runoff for the current time step
    !            call get_runoff(trim(input_dir)//trim(fname_qsim), & ! input: filename
    !          				  iTime,                             & ! input: time index
    !          				  dTime,                             & ! output: time
    !          				  ierr, cmessage)                      ! output: error control
    !            call handle_err(ierr, cmessage)
    !            ! ensure that simulated runoff is non-zero
    !            where(runoff_data%qsim < runoffMin) runoff_data%qsim=runoffMin
    !
    !            if (is_remap) then
    !          	call remap_runoff(rn_hru_id, runoff_data%qsim, runoff_data%hru_id, qsim_hru, ierr, cmessage)
    !          	if(ierr/=0) call handle_err(ierr,cmessage)
    !            else
    !          	qsim_hru=runoff_data%qsim
    !            end if
    !
    !            ! write time -- note time is just carried across from the input
    !            call write_nc(trim(output_dir)//trim(fname_output), 'time', (/dTime/), (/iTime/), (/1/), ierr, cmessage)
    !            call handle_err(ierr,cmessage)
    !
    !            ! *****
    !            ! (5b) Interpolate simulated runoff to local basins...
    !            ! ****************************************************
    !            ! interpolate the data to the basins
    !            do ibas=1,nSegRoute
    !             ! intialize the basin runoff
    !             qsim_basin(ibas) = 0._dp
    !             ! get the number of HRUs that drain into the basin
    !             nDrain = size(h2b(ibas)%cHRU)
    !             ! * case where HRUs drain into the segment
    !             if(nDrain > 0)then
    !          	! loop through the HRUs
    !          	do iHRU=1,nDrain
    !          	 ! get the HRU index
    !          	 ix = h2b(ibas)%cHRU(iHRU)%hru_ix
    !
    !          	 !Error check - runoff depth cannot be negative (no missing value)
    !          	 if( qsim_hru(ix) < 0._dp )then
    !          	   write (str,'(I10)') h2b(ibas)%cHRU(iHRU)%hru_id
    !          	   call handle_err(20,'negaive runoff value is invalid')
    !          	 endif
    !
    !          	 ! compute the weighted average
    !          	 qsim_basin(ibas) = qsim_basin(ibas) + h2b(ibas)%cHRU(iHRU)%wght*qsim_hru(ix)*time_conv*length_conv  ! ensure m/s
    !          	end do  ! (looping through gridpoints associated with each basin)
    !             ! * special case where no HRUs drain into the segment
    !             else
    !          	qsim_basin(ibas) = 0._dp
    !             endif
    !             ! convert runoff to m3/s
    !             RCHFLX(iens,ibas)%BASIN_QI = qsim_basin(ibas)*RPARAM(ibas)%BASAREA
    !            end do  ! (looping through basins)
    !
    !            ! write instantaneous local runoff in each stream segment (m3/s)
    !            call write_nc(trim(output_dir)//trim(fname_output), 'instBasinRunoff', RCHFLX(iens,:)%BASIN_QI, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    !            call handle_err(ierr,cmessage)
    !
    !            ! *****
    !            ! (5c) Delay runoff within local basins (hill-slope routing) ...
    !            ! ****************************************
    !
    !            ! route streamflow through the basin
    !            do ibas=1,nSegRoute ! place a fraction of runoff in future time steps
    !          	do jtim=1,ntdh
    !          	 RCHFLX(iens,ibas)%QFUTURE(jtim) = RCHFLX(iens,ibas)%QFUTURE(jtim) + FRAC_FUTURE(jtim)*RCHFLX(iens,ibas)%BASIN_QI
    !          	end do
    !          	! save the routed runoff
    !          	RCHFLX(iens,ibas)%BASIN_QR(0) = RCHFLX(iens,ibas)%BASIN_QR(1)  ! (save the runoff from the previous time step)
    !          	RCHFLX(iens,ibas)%BASIN_QR(1) = RCHFLX(iens,ibas)%QFUTURE(1)
    !          	! move array back
    !          	do jtim=2,ntdh
    !          	 RCHFLX(iens,ibas)%QFUTURE(jtim-1) = RCHFLX(iens,ibas)%QFUTURE(jtim)
    !          	end do
    !          	RCHFLX(iens,ibas)%QFUTURE(ntdh)    = 0._dp
    !            end do  ! (looping through basins)
    !
    !            ! write routed local runoff in each stream segment (m3/s)
    !            call write_nc(trim(output_dir)//trim(fname_output), 'dlayBasinRunoff', RCHFLX(iens,:)%BASIN_QR(1), (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    !            call handle_err(ierr,cmessage)
    !
    !            ! *****
    !            ! (5d) Route streamflow for each upstream segment through the river network with "IRF routing scheme"...
    !            ! **************************************************
    !            if (routOpt==0 .or. routOpt==1) then
    !          	! Get upstreme routed runoff depth at top of each segment for IRF routing
    !          	call get_upsbas_qr(nSegRoute,iens,ierr,cmessage)
    !          	! write routed runoff (m3/s)
    !          	call write_nc(trim(output_dir)//trim(fname_output), 'UpBasRoutedRunoff', RCHFLX(iens,:)%UPSBASIN_QR, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    !          	call handle_err(ierr,cmessage)
    !          	! Get upstreme routed runoff depth at top of each segment for IRF routing
    !          	call conv_upsbas_qr(nSegRoute,iens,ierr,cmessage)
    !          	! write routed runoff (m3/s)
    !          	call write_nc(trim(output_dir)//trim(fname_output), 'IRFroutedRunoff', RCHFLX(iens,:)%REACH_Q_IRF, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    !          	call handle_err(ierr,cmessage)
    !            endif
    !
    !            ! *****
    !            ! (5e) Compute total instantaneous runoff from all upstream basins...
    !            ! *******************************************************************
    !            ! compute the sum of all upstream runoff at each point in the river network
    !            do iSeg=1,nSegRoute
    !          	! identify how many reaches are upstream
    !          	nUpstream = size(NETOPO(iSeg)%RCHLIST)
    !          	! allocate space for upstream vectors
    !          	allocate(iUpstream(nUpstream), qUpstream(nUpstream), stat=ierr)
    !          	if(ierr/=0) call handle_err(ierr,'problem allocating vectors for all upstream basins')
    !          	! get indices for all reaches upstream
    !          	iUpstream(1:nUpstream) = NETOPO(iSeg)%RCHLIST(1:nUpstream)
    !          	! get streamflow for all reaches upstream
    !          	qUpstream(1:nUpstream) = RCHFLX(iens,iUpstream(1:nUpstream))%BASIN_QR(1)
    !          	! get mean streamflow
    !          	RCHFLX(iens,iSeg)%UPSTREAM_QI = sum(qUpstream)
    !          	! test
    !          	if(NETOPO(iSeg)%REACHID == ixPrint)then
    !          	 print*, 'ixUpstream = ', NETOPO(iUpstream(1:nUpstream))%REACHIX
    !          	 print*, 'idUpstream = ', NETOPO(iUpstream(1:nUpstream))%REACHID
    !          	 print*, 'qUpstream = ', qUpstream
    !          	endif
    !          	! deallocate space for upstream vectors
    !          	deallocate(iUpstream,qUpstream, stat=ierr)
    !          	if(ierr/=0) call handle_err(ierr,'problem deallocating vectors for all upstream basins')
    !            end do  ! looping through stream segments
    !
    !            ! write instantaneous runoff from all upstream basins (m3/s)
    !            call write_nc(trim(output_dir)//trim(fname_output), 'sumUpstreamRunoff', RCHFLX(iens,:)%UPSTREAM_QI, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    !            call handle_err(ierr,cmessage)
    !
    !            ! *****
    !            ! (5f) Route streamflow through the river network with "Kinematic Wave Routing scheme" ...
    !            ! **************************************************
    !            if (routOpt==0 .or. routOpt==2) then
    !          	! route streamflow through the river network
    !          	do iSeg=1,nSegRoute
    !          	  ! identify reach to process
    !          	  irch = NETOPO(iSeg)%RHORDER
    !          	  ! route kinematic waves through the river network
    !          	  call QROUTE_RCH(iens,irch,           & ! input: array indices
    !          					  ixDesire,            & ! input: index of the outlet reach
    !          					  T0,T1,               & ! input: start and end of the time step
    !          					  MAXQPAR,             & ! input: maximum number of particle in a reach
    !          					  LAKEFLAG,            & ! input: flag if lakes are to be processed
    !          					  ierr,cmessage)         ! output: error control
    !          	  if (ierr/=0) call handle_err(ierr,cmessage)
    !          	end do  ! (looping through stream segments)
    !          	do iSeg=1,nSegRoute
    !          	  wavesize(iens,iSeg) =size(KROUTE(iens,iSeg)%KWAVE)
    !          	  wavestate(iens,iSeg)=KROUTE(iens,iSeg)
    !          	end do  ! (looping through stream segments)
    !
    !          	! write routed runoff (m3/s)
    !          	call write_nc(trim(output_dir)//trim(fname_output), 'KWTroutedRunoff', RCHFLX(iens,:)%REACH_Q, (/1,iTime/), (/nSegRoute,1/), ierr, cmessage)
    !          	call handle_err(ierr,cmessage)
    !
    !            end if
    !
    !            ! *****
    !            ! (5g) print progress...
    !            ! **********************
    !            ! extract desired variables
    !            qDomain_hru     = sum(pack(qsim_hru,rn_hru_id_mask)*pack(rn_hru_area,rn_hru_id_mask))/sum(pack(rn_hru_area,rn_hru_id_mask)) ! domain: total runoff at the HRUs (mm/s)
    !            qDomain_basin   = 1000._dp*sum(qsim_basin*RPARAM(:)%BASAREA) / sum(RPARAM(:)%BASAREA)                  ! domain: total instantaneous basin runoff (mm/s)
    !            qDomain_reach   = 1000._dp*sum(RCHFLX(iens,:)%BASIN_QI)/sum(RPARAM(:)%BASAREA)                            ! domain: total instantaneous reach runoff (mm/s)
    !            qDesire_instant = 1000._dp*RCHFLX(iens,ixDesire)%UPSTREAM_QI/RPARAM(ixDesire)%TOTAREA                     ! desire: sum of instantaneous runoff for all basins upstream of the desired reach (mm/s)
    !            qDesire_routed  = 1000._dp*RCHFLX(iens,ixDesire)%REACH_Q/RPARAM(ixDesire)%TOTAREA                         ! desire: routed runoff for the desired reach (mm/s)
    !
    !            ! print progress
    !            write(*,'(i5,1x,10(e20.10,1x))') iTime, qDomain_hru, qDomain_basin, qDomain_reach, qDesire_instant, qDesire_routed
    !
    !           end do  ! (looping through ensemble members)
    !
    !           ! increment time bounds
    !           T0 = T0 + dt
    !           T1 = T0 + dt
    !
    !          end do  ! (looping through time)
    !
    !          ! *****
    !          ! (6) write state variable
    !          ! **********************
    !          ! create States NetCDF file
    !          call defineStateFile(trim(output_dir)//trim(fname_state_out), &  ! input: file name
    !          					 nEns,                                    &  ! input: dimension size: number of ensembles
    !          					 nSegRoute,                               &  ! input: dimension size: number of stream segments
    !          					 ntdh,                                    &  ! input: dimension size: number of future time steps for hillslope UH routing
    !          					 maxval(irfsize),                         &  ! input: dimension size: number of future time steps for irf UH routing (max. of all the segments)
    !          					 maxval(wavesize),                        &  ! input: dimension size: number of waves (max. of all the segments)
    !          					 routOpt,                                 &  ! input: routing scheme options
    !          					 ierr, cmessage)                             ! output: error control
    !          if(ierr/=0) call handle_err(ierr, cmessage)
    !          call write_nc(trim(output_dir)//trim(fname_state_out),'reachID', NETOPO(:)%REACHID, (/1/), (/size(NETOPO)/), ierr, cmessage);
    !          if(ierr/=0) call handle_err(ierr,cmessage)
    !          call write_nc(trim(output_dir)//trim(fname_state_out),'time_bound', (/T0,T1/), (/1/), (/2/), ierr, cmessage)
    !          if(ierr/=0) call handle_err(ierr,cmessage)
    !
    !          allocate(qfuture_array(nSegRoute,ntdh,nens),stat=ierr)
    !          if(ierr/=0) call handle_err(ierr,'problem allocating qfuture_array')
    !          allocate(basin_qr_array(nSegRoute,nens),stat=ierr)
    !          if(ierr/=0) call handle_err(ierr,'problem allocating basin_qr_array')
    !          allocate(qfuture_irf_array(nSegRoute,maxval(irfsize),nens),stat=ierr)
    !          if(ierr/=0) call handle_err(ierr,'problem allocating qfuture_irf_array')
    !          allocate(qf_array(nSegRoute,maxval(wavesize),nens),stat=ierr)
    !          if(ierr/=0) call handle_err(ierr,'problem allocating qf_array')
    !          allocate(qm_array(nSegRoute,maxval(wavesize),nens),stat=ierr)
    !          if(ierr/=0) call handle_err(ierr,'problem allocating qm_array')
    !          allocate(ti_array(nSegRoute,maxval(wavesize),nens),stat=ierr)
    !          if(ierr/=0) call handle_err(ierr,'problem allocating ti_array')
    !          allocate(tr_array(nSegRoute,maxval(wavesize),nens),stat=ierr)
    !          if(ierr/=0) call handle_err(ierr,'problem allocating tr_array')
    !          allocate(rf_array(nSegRoute,maxval(wavesize),nens),stat=ierr)
    !          if(ierr/=0) call handle_err(ierr,'problem allocating rf_array')
    !
    !          do iens=1,nens
    !            do iSeg=1,nSegRoute
    !          	qfuture_array(iSeg,1:ntdh,iens) = RCHFLX(iens,iSeg)%QFUTURE(:)
    !          	basin_qr_array(iSeg,iens) = RCHFLX(iens,iSeg)%BASIN_QR(1)
    !          	if (routOpt==0 .or. routOpt==1) then
    !          	  qfuture_irf_array(iSeg,1:irfsize(iens,iSeg),iens) = RCHFLX(iens,iSeg)%QFUTURE_IRF(:)
    !          	end if
    !          	if (routOpt==0 .or. routOpt==2) then
    !          	  if (allocated(RFvec)) deallocate(RFvec,stat=ierr)
    !          	  allocate(RFvec(wavesize(iens,iSeg)),stat=ierr)
    !          	  RFvec=0_i4b
    !          	  where (wavestate(iens,iSeg)%KWAVE(:)%RF) RFvec=1_i4b
    !          	  qf_array(iSeg,1:wavesize(iens,iSeg),iens) = wavestate(iens,iSeg)%KWAVE(:)%QF
    !          	  qm_array(iSeg,1:wavesize(iens,iSeg),iens) = wavestate(iens,iSeg)%KWAVE(:)%QM
    !          	  ti_array(iSeg,1:wavesize(iens,iSeg),iens) = wavestate(iens,iSeg)%KWAVE(:)%TI
    !          	  tr_array(iSeg,1:wavesize(iens,iSeg),iens) = wavestate(iens,iSeg)%KWAVE(:)%TR
    !          	  rf_array(iSeg,1:wavesize(iens,iSeg),iens) = RFvec
    !          	end if
    !            enddo
    !          enddo
    !
    !          ! write hill-slope routing state
    !          call write_nc(trim(output_dir)//trim(fname_state_out), 'QFUTURE',qfuture_array, (/1,1,1/), (/nSegRoute,ntdh,nens/), ierr, cmessage)
    !          if(ierr/=0) call handle_err(ierr,cmessage)
    !          call write_nc(trim(output_dir)//trim(fname_state_out),'BASIN_QR',BASIN_QR_array, (/1,1/), (/nSegRoute,nens/), ierr, cmessage)
    !          if(ierr/=0) call handle_err(ierr,cmessage)
    !          ! write IRF routing  state for restart
    !          if (routOpt==0 .or. routOpt==1) then
    !            call write_nc(trim(output_dir)//trim(fname_state_out), 'irfsize', irfsize, (/1,1/), (/nSegRoute,nens/), ierr, cmessage)
    !            if(ierr/=0) call handle_err(ierr,cmessage)
    !            call write_nc(trim(output_dir)//trim(fname_state_out), 'QFUTURE_IRF', qfuture_irf_array, (/1,1,1/), (/nSegRoute,maxval(irfsize),nens/), ierr, cmessage)
    !            if(ierr/=0) call handle_err(ierr,cmessage)
    !          end if
    !          !write KWT routing state for restart
    !          if (routOpt==0 .or. routOpt==2) then
    !            call write_nc(trim(output_dir)//trim(fname_state_out), 'wavesize', wavesize, (/1,1/), (/nSegRoute,nens/), ierr, cmessage)
    !            if(ierr/=0) call handle_err(ierr,cmessage)
    !            call write_nc(trim(output_dir)//trim(fname_state_out), 'QF', QF_array, (/1,1,1/), (/nSegRoute,maxval(wavesize),nens/), ierr, cmessage)
    !            if(ierr/=0) call handle_err(ierr,cmessage)
    !            call write_nc(trim(output_dir)//trim(fname_state_out), 'QM', QM_array, (/1,1,1/), (/nSegRoute,maxval(wavesize),nens/), ierr, cmessage)
    !            if(ierr/=0) call handle_err(ierr,cmessage)
    !            call write_nc(trim(output_dir)//trim(fname_state_out), 'TI', TI_array, (/1,1,1/), (/nSegRoute,maxval(wavesize),nens/), ierr, cmessage)
    !            if(ierr/=0) call handle_err(ierr,cmessage)
    !            call write_nc(trim(output_dir)//trim(fname_state_out), 'TR', TR_array, (/1,1,1/), (/nSegRoute,maxval(wavesize),nens/), ierr, cmessage)
    !            if(ierr/=0) call handle_err(ierr,cmessage)
    !            call write_nc(trim(output_dir)//trim(fname_state_out), 'RF', rf_array, (/1,1,1/), (/nSegRoute,maxval(wavesize),nens/), ierr, cmessage)
    !            if(ierr/=0) call handle_err(ierr,cmessage)
    !          end if

stop

contains

 subroutine handle_err(err,message)
 ! handle error codes
 implicit none
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 if(err/=0)then
  print*,'FATAL ERROR: '//trim(message)
  call flush(6)
  stop
 endif
 end subroutine handle_err

end program route_runoff
