MODULE write_simoutput

! Moudle wide external modules
USE nrtype
USE var_lookup,only: ixRFLX, nVarsRFLX
USE public_var,only: iulog
USE public_var,only: integerMissing
USE globalData,only: meta_rflx
USE globalData,only: simout_nc
USE globalData,only: idxSUM, idxIRF, idxKWT, idxKW, idxMC, idxDW
USE io_netcdf, only: ncd_int
USE io_netcdf, only: ncd_float, ncd_double
USE io_netcdf, only: ncd_unlimited
USE io_netcdf, only: def_nc                 ! define netcdf
USE io_netcdf, only: def_var                ! define netcdf variable
USE io_netcdf, only: def_dim                ! define netcdf dimension
USE io_netcdf, only: put_global_attr        ! write global attributes
USE io_netcdf, only: end_def                ! end defining netcdf
USE io_netcdf, only: open_nc                ! open netcdf
USE io_netcdf, only: close_nc               ! close netcdf
USE io_netcdf, only: write_nc               ! write a variable to the NetCDF file

implicit none

! The following variables used only in this module
integer(i4b),         save        :: jTime            ! time step in output netCDF

private

public::prep_output
public::output

CONTAINS

 ! *********************************************************************
 ! public subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE output(ierr, message)    ! out:   error control
  !Dependent modules
  USE globalData, ONLY: nHRU, nRch          ! number of ensembles, HRUs and river reaches
  USE globalData, ONLY: RCHFLX              ! Reach fluxes (ensembles, space [reaches])
  USE globalData, ONLY: runoff_data         ! runoff data for one time step for LSM HRUs and River network HRUs

  implicit none

  ! input variables: none
  ! output variables
  integer(i4b), intent(out)       :: ierr             ! error code
  character(*), intent(out)       :: message          ! error message
  ! local variables
  real(dp),    allocatable        :: array_temp(:)
  integer(i4b)                    :: ix               ! loop index
  integer(i4b)                    :: iens             ! temporal
  character(len=strLen)           :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='output/'

  iens = 1

  allocate(array_temp(nRch), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif

  ! write time -- note time is just carried across from the input
  call write_nc(simout_nc%ncid, 'time', (/runoff_data%time/), (/jTime/), (/1/), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (meta_rflx(ixRFLX%basRunoff)%varFile) then
   ! write the basin runoff at HRU (m/s)
   call write_nc(simout_nc%ncid, 'basRunoff', runoff_data%basinRunoff, (/1,jTime/), (/nHRU,1/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%instRunoff)%varFile) then
   ! write instataneous local runoff in each stream segment (m3/s)
   call write_nc(simout_nc%ncid, 'instRunoff', RCHFLX(iens,:)%BASIN_QI, (/1,jTime/), (/nRch,1/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
   ! write routed local runoff in each stream segment (m3/s)
   call write_nc(simout_nc%ncid, 'dlayRunoff', RCHFLX(iens,:)%BASIN_QR(1), (/1,jTime/), (/nRch,1/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile) then
    do ix=1,nRCH
      array_temp(ix) = RCHFLX(iens, ix)%ROUTE(idxSUM)%REACH_Q
    end do
    call write_nc(simout_nc%ncid, 'sumUpstreamRunoff', array_temp, (/1,jTime/), (/nRch,1/), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%KWTroutedRunoff)%varFile) then
    do ix=1,nRCH
      array_temp(ix) = RCHFLX(iens, ix)%ROUTE(idxKWT)%REACH_Q
    end do
    call write_nc(simout_nc%ncid, 'KWTroutedRunoff', array_temp, (/1,jTime/), (/nRch,1/), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

  if (meta_rflx(ixRFLX%IRFroutedRunoff)%varFile) then
    do ix=1,nRCH
      array_temp(ix) = RCHFLX(iens, ix)%ROUTE(idxIRF)%REACH_Q
    end do
    call write_nc(simout_nc%ncid, 'IRFroutedRunoff', array_temp, (/1,jTime/), (/nRch,1/), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

  if (meta_rflx(ixRFLX%KWroutedRunoff)%varFile) then
    do ix=1,nRCH
      array_temp(ix) = RCHFLX(iens, ix)%ROUTE(idxKW)%REACH_Q
    end do
    call write_nc(simout_nc%ncid, 'KWroutedRunoff', array_temp, (/1,jTime/), (/nRch,1/), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

  if (meta_rflx(ixRFLX%MCroutedRunoff)%varFile) then
    do ix=1,nRCH
      array_temp(ix) = RCHFLX(iens, ix)%ROUTE(idxMC)%REACH_Q
    end do
    call write_nc(simout_nc%ncid, 'MCroutedRunoff', array_temp, (/1,jTime/), (/nRch,1/), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

  if (meta_rflx(ixRFLX%DWroutedRunoff)%varFile) then
    do ix=1,nRCH
      array_temp(ix) = RCHFLX(iens, ix)%ROUTE(idxDW)%REACH_Q
    end do
    call write_nc(simout_nc%ncid, 'DWroutedRunoff', array_temp(1:nRch), (/1,jTime/), (/nRch,1/), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

 END SUBROUTINE output


 ! *********************************************************************
 ! public subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE prep_output(ierr, message)    ! out:   error control

 USE ascii_util_module,   ONLY: lower
 ! saved public variables (usually parameters, or values not modified)
 USE public_var,          only : output_dir        ! output directory
 USE public_var,          only : case_name         ! simulation name ==> output filename head
 USE public_var,          only : calendar          ! calendar name
 USE public_var,          only : newFileFrequency  ! frequency for new output files (day, month, annual)
 USE public_var,          only : time_units        ! time units (seconds, hours, or days)
 ! saved global data
 USE globalData,          only : basinID,reachID   ! HRU and reach ID in network
 USE globalData,          only : modTime           ! previous and current model time
 USE globalData,          only : nEns, nHRU, nRch  ! number of ensembles, HRUs and river reaches

 implicit none

 ! input variables: none
 ! output variables
 integer(i4b), intent(out)       :: ierr             ! error code
 character(*), intent(out)       :: message          ! error message
 ! local variables
 character(len=strLen)           :: cmessage         ! error message of downwind routine
 integer(i4b)                    :: sec_in_day       ! second within day
 logical(lgt)                    :: defnewoutputfile ! flag to define new output file
 character(len=50),parameter     :: fmtYMDS='(a,I0.4,a,I0.2,a,I0.2,a,I0.5,a)'

 ierr=0; message='prep_output/'

 ! print progress
 write(iulog,'(a,I4,4(x,I4))') new_line('a'), modTime(1)%year(), modTime(1)%month(), modTime(1)%day(), modTime(1)%hour(), modTime(1)%minute()

 ! check need for the new file
 select case(lower(trim(newFileFrequency)))
   case('single'); defNewOutputFile=(modTime(0)%year() ==integerMissing)
   case('yearly'); defNewOutputFile=(modTime(1)%year() /=modTime(0)%year())
   case('monthly');  defNewOutputFile=(modTime(1)%month()/=modTime(0)%month())
   case('daily');    defNewOutputFile=(modTime(1)%day()  /=modTime(0)%day())
   case default; ierr=20; message=trim(message)//'Accepted <newFileFrequency> options (case-insensitive): single yearly, monthly, or daily '; return
 end select

 ! define new file
 if(defNewOutputFile)then

   if (simout_nc%status == 2) then
     call close_nc(simout_nc%ncid, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     simout_nc%ncname = 'empty'
     simout_nc%ncid   = integerMissing
     simout_nc%status = integerMissing
   endif

   ! initialize time
   jTime=1

   ! update filename
   sec_in_day = modTime(1)%hour()*60*60+modTime(1)%minute()*60+nint(modTime(1)%sec())
   write(simout_nc%ncname, fmtYMDS) trim(output_dir)//trim(case_name)//'.h.', &
                                     modTime(1)%year(), '-', modTime(1)%month(), '-', modTime(1)%day(), '-',sec_in_day,'.nc'

   call defineFile(simout_nc%ncname,                      &  ! input: file name
                   nEns,                                  &  ! input: number of ensembles
                   nHRU,                                  &  ! input: number of HRUs
                   nRch,                                  &  ! input: number of stream segments
                   time_units,                            &  ! input: time units
                   calendar,                              &  ! input: calendar
                   ierr,cmessage)                            ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call open_nc(simout_nc%ncname, 'w', simout_nc%ncid, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   simout_nc%status = 2

   if (meta_rflx(ixRFLX%basRunoff)%varFile) then
     call write_nc(simout_nc%ncid, 'basinID', int(basinID,kind(i4b)), (/1/), (/nHRU/), ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   call write_nc(simout_nc%ncid, 'reachID', reachID, (/1/), (/nRch/), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! no new file requested: increment time
 else

   jTime = jTime+1

 endif

 END SUBROUTINE prep_output


 ! *********************************************************************
 ! private subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE defineFile(fname,           &  ! input: filename
                       nEns_in,         &  ! input: number of ensembles
                       nHRU_in,         &  ! input: number of HRUs
                       nRch_in,         &  ! input: number of stream segments
                       units_time,      &  ! input: time units
                       calendar,        &  ! input: calendar
                       ierr, message)      ! output: error control
 !Dependent modules
 USE public_var, ONLY: mizuRouteVersion
 USE public_var, ONLY: netcdf_format
 USE globalData, ONLY: meta_qDims
 USE var_lookup, ONLY: ixQdims, nQdims

 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 integer(i4b), intent(in)        :: nEns_in      ! number of ensembles
 integer(i4b), intent(in)        :: nHRU_in      ! number of HRUs
 integer(i4b), intent(in)        :: nRch_in      ! number of stream segments
 character(*), intent(in)        :: units_time   ! time units
 character(*), intent(in)        :: calendar     ! calendar
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 character(len=strLen),allocatable :: dim_array(:)
 integer(i4b)                      :: nDims
 integer(i4b)                      :: ixDim
 integer(i4b)                      :: ncid         ! NetCDF file ID
 integer(i4b)                      :: jDim, iVar   ! dimension, and variable index
 character(len=strLen)             :: cmessage     ! error message of downwind routine

 ! initialize error control
 ierr=0; message='defineFile/'

! populate q dimension meta (not sure if this should be done here...)
 meta_qDims(ixQdims%seg)%dimLength = nRch_in
 meta_qDims(ixQdims%hru)%dimLength = nHRU_in
 meta_qDims(ixQdims%ens)%dimLength = nEns_in

 ! --------------------
 ! define file
 ! --------------------
 call def_nc(trim(fname), ncid, ierr, cmessage, nctype=netcdf_format)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 do jDim =1,nQdims
   if (jDim ==ixQdims%time) then ! time dimension (unlimited)
     call def_dim(ncid, trim(meta_qDims(jDim)%dimName), ncd_unlimited, meta_qDims(jDim)%dimId, ierr, cmessage)
   else if (jDim==ixQdims%hru) then
     if (meta_rflx(ixRFLX%basRunoff)%varFile) then
       call def_dim(ncid, trim(meta_qDims(jDim)%dimName), meta_qDims(jDim)%dimLength ,meta_qDims(jDim)%dimId, ierr, cmessage)
     end if
   else
     call def_dim(ncid, trim(meta_qDims(jDim)%dimName), meta_qDims(jDim)%dimLength ,meta_qDims(jDim)%dimId, ierr, cmessage)
   endif
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do

 ! Define coordinate variable for time
 call def_var(ncid, 'time', (/meta_qDims(ixQdims%time)%dimName/), ncd_double, ierr, cmessage, vdesc='time', vunit=trim(units_time), vcal=calendar)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! Define ID variable for time
 call def_var(ncid, 'reachID', (/meta_qDims(ixQdims%seg)%dimName/), ncd_int, ierr, cmessage, vdesc='reach ID', vunit='-')
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 if (meta_rflx(ixRFLX%basRunoff)%varFile) then
   call def_var(ncid, 'basinID', (/meta_qDims(ixQdims%hru)%dimName/), ncd_int, ierr, cmessage, vdesc='basin ID', vunit='-')
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! define variables
 do iVar=1, nVarsRFLX

  if (.not.meta_rflx(iVar)%varFile) cycle

  ! define dimension ID array
  nDims = size(meta_rflx(iVar)%varDim)
  if (allocated(dim_array)) then
    deallocate(dim_array)
  endif
  allocate(dim_array(nDims))
  do ixDim = 1, nDims
    dim_array(ixDim) = meta_qDims(meta_rflx(iVar)%varDim(ixDim))%dimName
  end do

  call def_var(ncid, meta_rflx(iVar)%varName, dim_array, meta_rflx(iVar)%varType, ierr, cmessage, vdesc=meta_rflx(iVar)%varDesc, vunit=meta_rflx(iVar)%varUnit )
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end do

 call put_global_attr(ncid, 'version', trim(mizuRouteVersion) ,ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call end_def(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call close_nc(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE defineFile


END MODULE write_simoutput
