module getAncillary_module

! data types
use nrtype
USE dataTypes,  only : remap                 ! remapping data type
USE dataTypes,  only : runoff                ! runoff data type

! global data
USE public_var

! external subroutines
USE read_remap, only : get_runoff_metadata   ! runoff metadata
USE read_remap, only : get_remap_data        ! get remap data

implicit none

! privacy
private
public::getAncillary

contains

 ! *****
 ! public subroutine: get mapping data between runoff hru and river network hru
 ! *********************************************************************
 subroutine getAncillary(&
                         ierr, message)     ! output: error control
 implicit none
 ! dummy variables
 integer(i4b), intent(out)          :: ierr            ! error code
 character(*), intent(out)          :: message         ! error message
 ! local variables
 integer(i4b)                       :: nTime           ! number of time steps
 integer(i4b)                       :: nSpatial        ! number of spatial elements
 character(len=strLen)              :: units_time      ! time units
 type(remap)                        :: remap_data      ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
 type(runoff)                       :: runoff_data     ! runoff for one time step for all HRUs
 character(len=strLen)              :: cmessage        ! error message from subroutine

 ! get runoff metadata
 call get_runoff_metadata(trim(input_dir)//trim(fname_qsim), & ! input: filename
                          nSpatial, runoff_data%hru_id,      & ! output: number spatial elements and HRU ids
                          nTime,    units_time,              & ! output: number of time steps and time units
                          ierr,     cmessage)                  ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 print*, 'nSpatial, nTime, trim(units_time) = ', nSpatial, nTime, trim(units_time)
 print*, 'runoff_data%hru_id = ', runoff_data%hru_id
 print*, trim(message)//'PAUSE : '; read(*,*)


 ! get runoff mapping file
 call get_remap_data(trim(ancil_dir)//trim(fname_remap),     & ! input: file name
                     remap_data,                             & ! output: data structure to remap data from a polygon
                     ierr, cmessage)                           ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif





 end subroutine getAncillary

end module getAncillary_module
        !          ! *****
        !          ! (3) get river network hru id
        !          ! *******************************************
        !
        !          ! get the number of HRUs
        !          if (is_remap) then
        !            nHRU_rn = size(remap_data%hru_id)
        !          else
        !            nHRU_rn = size(runoff_data%hru_id)
        !          endif
        !
        !          ! allocate space for the HRUids
        !          allocate(rn_hru_id(nHRU_rn), stat=ierr)
        !          if(ierr/=0) call handle_err(ierr,'problem allocating space for rn_hru_id')
        !
        !          ! get the HRU ids
        !          if (is_remap) then
        !            rn_hru_id(:) = remap_data%hru_id(:)
        !          else
        !            rn_hru_id(:) = runoff_data%hru_id(:)
        !          endif
        !
        !          ! allocate space for the HRU mask
        !          allocate(rn_hru_id_mask(nHRU_rn), stat=ierr)
        !          if(ierr/=0) call handle_err(ierr,'problem allocating space for rn_hru_id_mask')
        !
        !          ! allocate space for the HRU area
        !          allocate(rn_hru_area(nHRU_rn), stat=ierr)
        !          if(ierr/=0) call handle_err(ierr,'problem allocating space for rn_hru_area')
        !
        !          ! initialize mask and area
        !          rn_hru_id_mask(:)=.false.
        !          rn_hru_area(:)=-999
        !
        !          ! check all the upstream hrus at the desired outlet exist in runoff file
        !          ! Assign hru_ix based on order of hru in runoff file
        !          do iSeg=1,nSegRoute ! (loop through stream segments)
        !            nDrain = size(h2b(iSeg)%cHRU)  ! (number of HRUs that drain into a given stream segment)
        !            if(nDrain > 0)then
        !            do iHRU=1,nDrain ! (loop through HRUs that drain into a given stream segment)
        !          	! check existence of upstream hrus in funoff file
        !          	if( minval(abs(rn_hru_id - h2b(iSeg)%cHRU(iHRU)%hru_id)) /= 0 )then
        !          	  write (str,'(I10)') h2b(iSeg)%cHRU(iHRU)%hru_id
        !          	  call handle_err(20,'runoff file does not include runoff time series at HRU'//trim(str))
        !          	end if
        !          	! Assign hru index
        !          	iSelect = minloc(abs(rn_hru_id - h2b(iSeg)%cHRU(iHRU)%hru_id))
        !          	h2b(iSeg)%cHRU(iHRU)%hru_ix=iSelect(1)
        !          	if(h2b(iSeg)%cHRU(iHRU)%hru_id /= rn_hru_id(h2b(iSeg)%cHRU(iHRU)%hru_ix)) call handle_err(20,'mismatch in HRUs')
        !          	rn_hru_id_mask(iSelect(1)) = .true.
        !          	rn_hru_area(iSelect(1)) = h2b(iSeg)%cHRU(iHRU)%hru_area
        !            end do  ! (loop through HRUs that drain into a given stream segment)
        !           endif  ! (if HRUs drain into the stream segment)
        !          end do ! (loop through stream segments)
        !
        !          ! *****
        !          ! (3) Define NetCDF output file and write ancillary data...
        !          ! *********************************************************
        !          ! create NetCDF file
        !          call defineFile(trim(output_dir)//trim(fname_output),  &  ! input: file name
        !          				nSegRoute,                             &  ! input: number of stream segments
        !          				nTotal,                                &  ! input: total number of upstream reaches for all reaches
        !          				units_time,                            &  ! input: time units
        !          				ierr, cmessage)                           ! output: error control
        !          call handle_err(ierr, cmessage)
        !
        !          ! write network toplogy (input = filename, variable name, variable vector, start index; output = error control)
        !          call write_nc(trim(output_dir)//trim(fname_output), 'reachID',    NETOPO(:)%REACHID, (/1/), (/size(NETOPO)/), ierr, cmessage); call handle_err(ierr,cmessage)
        !          call write_nc(trim(output_dir)//trim(fname_output), 'reachOrder', NETOPO(:)%RHORDER, (/1/), (/size(NETOPO)/), ierr, cmessage); call handle_err(ierr,cmessage)
        !
        !          iStart=1  ! initialize the start index of the ragged array
        !          ! write list of reaches upstream of each reach (ragged array)
        !          do iSeg=1,nSegRoute
        !           ! get the number of reaches
        !           nUpstream = size(NETOPO(iSeg)%RCHLIST)
        !           ! write the vector to the ragged array
        !           call write_nc(trim(output_dir)//trim(fname_output), 'reachList', NETOPO(iSeg)%RCHLIST(:), (/iStart/), (/size(NETOPO(iSeg)%RCHLIST)/), ierr, cmessage)
        !           call handle_err(ierr,cmessage)
        !           ! write the start index and the count (NOTE: pass as a vector)
        !           call write_nc(trim(output_dir)//trim(fname_output), 'listStart', (/iStart/),    (/iSeg/), (/1/), ierr, cmessage); call handle_err(ierr,cmessage)
        !           call write_nc(trim(output_dir)//trim(fname_output), 'listCount', (/nUpstream/), (/iSeg/), (/1/), ierr, cmessage); call handle_err(ierr,cmessage)
        !           ! update the start index
        !           iStart = iStart + nUpstream
        !          end do
        !          ! write reach parameters
        !          call write_nc(trim(output_dir)//trim(fname_output), 'basinArea',    RPARAM(:)%BASAREA, (/1/), (/nSegRoute/), ierr, cmessage); call handle_err(ierr,cmessage)
        !          call write_nc(trim(output_dir)//trim(fname_output), 'upstreamArea', RPARAM(:)%TOTAREA, (/1/), (/nSegRoute/), ierr, cmessage); call handle_err(ierr,cmessage)
        !
        !          ! *****
        !          ! (4) Prepare for the routing simulations...
        !          ! *******************************************
        !          ! allocate space for the simulated runoff at the HRUs
        !          allocate(qsim_hru(nHRU_rn), stat=ierr)
        !          if(ierr/=0) call handle_err(ierr,'problem allocating space for simulated runoff at the HRUs')
        !
        !          ! allocate space for the simulated runoff at basins and reaches
        !          allocate(qsim_basin(nSegRoute),RCHFLX(nens,nSegRoute), KROUTE(nens,nSegRoute), wavestate(nens,nSegRoute), stat=ierr)
        !          if(ierr/=0) call handle_err(ierr,'problem allocating space for simulated runoff at the basins')
        !
        !          allocate(irfsize(nens,nSegRoute), wavesize(nens,nSegRoute), stat=ierr)
        !          if(ierr/=0) call handle_err(ierr,'problem allocating space for simulated runoff at the basins')
        !
        !          ! initialize the routed elements
        !          RCHFLX(:,:)%BASIN_QR(1) = 0._dp
        !          RCHFLX(:,:)%BASIN_QR_IRF(1) = 0._dp
        !
        !          ! initialize the time-delay histogram
        !          ! identify the number of future time steps for a given basin
        !          ntdh = size(FRAC_FUTURE)
        !          do iens=1,nens
        !            nUH_DATA_MAX=0
        !            do ibas=1,nSegRoute
        !          	! allocate space for the delayed runoff for Hillslope routing
        !          	allocate(RCHFLX(iens,ibas)%QFUTURE(ntdh), stat=ierr)
        !          	call handle_err(ierr, 'problem allocating space for QFUTURE element')
        !          	! initialize to zeroes
        !          	RCHFLX(iens,ibas)%QFUTURE(:) = 0._dp
        !
        !            !  allocate space for the delayed runoff for IRF routing
        !          	if ((routOpt==0 .or. routOpt==1) .and. .not.(isRestart)) then
        !          	  nUpstream = size(NETOPO(ibas)%RCHLIST) ! size of upstream segment
        !          	  nUH_DATA_MAX= max(nUH_DATA_MAX ,size(NETOPO(ibas)%UH))
        !          	  allocate(RCHFLX(iens,ibas)%QFUTURE_IRF(nUH_DATA_MAX), stat=ierr)
        !          	  call handle_err(ierr, 'problem allocating space for QFUTURE_IRF element')
        !          	  ! initialize to zeroes
        !          	  RCHFLX(iens,ibas)%QFUTURE_IRF(:) = 0._dp
        !          	  irfsize(iens,ibas)=nUH_DATA_MAX
        !          	end if
        !
        !            end do
        !          end do
        !
        !          print*, "Start routing"
        !          ! define flags
        !          LAKEFLAG=0  ! no lakes in the river network
        !
        !          ! define time
        !          T0 = 0._dp
        !          T1 = dt
