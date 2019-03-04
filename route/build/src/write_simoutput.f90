MODULE write_simoutput

! Moudle wide external modules
USE nrtype
USE netcdf
USE public_var
implicit none

private

public::prep_output

CONTAINS

 ! *********************************************************************
 ! new subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE prep_output(nEns,           & ! in:    number of ensembles
                        nHRU,           & ! in:    number of HRUs
                        nRch,           & ! in:    number of reaches
                        iTime,          & ! in:    time index
                        modJulday,      & ! out:   current simulation julian day
                        jTime,          & ! out:   time step in output netCDF
                        ierr, message)    ! out:   error control
 !Dependent modules
 USE dataTypes,           only : time              ! time data type
 USE public_var,          only : refJulday         ! julian day: reference
 USE public_var,          only : startJulday       ! julian day: start
 USE public_var,          only : endJulday         ! julian day: end
 USE public_var,          only : calendar          ! calendar name
 USE public_var,          only : newFileFrequency  ! frequency for new output files (day, month, annual)
 USE public_var,          only : time_units        ! time units (seconds, hours, or days)
 USE public_var,          only : annual,month,day  ! time frequency named variable for output files
 USE globalData,          only : basinID,reachID   ! HRU and reach ID in network
 USE globalData,          only : timeVar           ! time variable
 USE globalData,          only : modTime           ! previous and current model time
 USE time_utils_module,   only : compCalday        ! compute calendar day
 USE time_utils_module,   only : compCalday_noleap ! compute calendar day
 USE write_netcdf,        only : write_nc          ! write a variable to the NetCDF file

 implicit none

 ! input variables
 integer(i4b), intent(in)        :: nEns             ! number of ensembles
 integer(i4b), intent(in)        :: nHRU             ! number of HRUs
 integer(i4b), intent(in)        :: nRch             ! number of stream segments
 integer(i4b), intent(in)        :: iTime            ! simulation time index
 ! output variables
 real(dp),     intent(out)       :: modJulday        ! julian day: model simulation
 integer(i4b), intent(out)       :: jTime            ! time step in output netCDF
 integer(i4b), intent(out)       :: ierr             ! error code
 character(*), intent(out)       :: message          ! error message
 ! local variables
 character(len=strLen)           :: fileout          ! name of the output file
 real(dp)                        :: convTime2Days    ! conversion factor to convert time to units of days
 logical(lgt)                    :: defnewoutputfile ! flag to define new output file
 character(len=strLen)           :: cmessage         ! error message of downwind routine

 ! initialize error control
 ierr=0; message='prep_output/'

 ! get the time multiplier needed to convert time to units of days
 select case( trim( time_units(1:index(time_units,' ')) ) )
  case('seconds'); convTime2Days=86400._dp
  case('hours');   convTime2Days=24._dp
  case('days');    convTime2Days=1._dp
  case default;    ierr=20; message=trim(message)//'unable to identify time units'; return
 end select

 ! get the julian day of the model simulation
 modJulday = refJulday + timeVar(iTime)/convTime2Days

 if(modJulday < startJulday .or. modJulday > endJulday) then
   return
 endif

 ! get the time
 select case(trim(calendar))
  case('noleap')
   call compCalday_noleap(modJulday,modTime(1)%iy,modTime(1)%im,modTime(1)%id,modTime(1)%ih,modTime(1)%imin,modTime(1)%dsec,ierr,cmessage)
  case ('standard','gregorian','proleptic_gregorian')
   call compCalday(modJulday,modTime(1)%iy,modTime(1)%im,modTime(1)%id,modTime(1)%ih,modTime(1)%imin,modTime(1)%dsec,ierr,cmessage)
  case default;    ierr=20; message=trim(message)//'calendar name: '//trim(calendar)//' invalid'; return
 end select
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! print progress
  !print*, modTime(1)%iy,modTime(1)%im,modTime(1)%id,modTime(1)%ih,modTime(1)%imin

  ! *****
  ! *** Define model output file...
  ! *******************************

  ! check need for the new file
  select case(newFileFrequency)
   case(annual); defNewOutputFile=(modTime(1)%iy/=modTime(0)%iy)
   case(month);  defNewOutputFile=(modTime(1)%im/=modTime(0)%im)
   case(day);    defNewOutputFile=(modTime(1)%id/=modTime(0)%id)
   case default; ierr=20; message=trim(message)//'unable to identify the option to define new output files'; return
  end select

  ! define new file
  if(defNewOutputFile)then

    ! initialize time
    jTime=1

    ! update filename
    write(fileout,'(a,3(i0,a))') trim(output_dir)//trim(fname_output)//'_', modTime(1)%iy, '-', modTime(1)%im, '-', modTime(1)%id, '.nc'

    ! define output file
    call defineFile(trim(fileout),                         &  ! input: file name
                    nEns,                                  &  ! input: number of ensembles
                    nHRU,                                  &  ! input: number of HRUs
                    nRch,                                  &  ! input: number of stream segments
                    time_units,                            &  ! input: time units
                    calendar,                              &  ! input: calendar
                    ierr,cmessage)                            ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! define basin ID
    call write_nc(trim(fileout), 'basinID', basinID, (/1/), (/nHRU/), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! define reach ID
    call write_nc(trim(fileout), 'reachID', reachID, (/1/), (/nRch/), ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! no new file requested: increment time
  else
    jTime = jTime+1
  endif

  modTime(0) = modTime(1)

 END SUBROUTINE prep_output


 ! *********************************************************************
 ! new subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE defineFile(fname,           &  ! input: filename
                       nEns,            &  ! input: number of ensembles
                       nHRU,            &  ! input: number of HRUs
                       nSeg,            &  ! input: number of stream segments
                       units_time,      &  ! input: time units
                       calendar,        &  ! input: calendar
                       ierr, message)      ! output: error control
 !Dependent modules
 USE globalData, ONLY: meta_qDims
 USE var_lookup, ONLY: ixQdims, nQdims

 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 integer(i4b), intent(in)        :: nEns         ! number of ensembles
 integer(i4b), intent(in)        :: nHRU         ! number of HRUs
 integer(i4b), intent(in)        :: nSeg         ! number of stream segments
 character(*), intent(in)        :: units_time   ! time units
 character(*), intent(in)        :: calendar     ! calendar
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: jDim, iVar   ! dimension, and variable index
 integer(i4b),parameter          :: nVars=8      ! number of variables
 character(len=strLen)           :: cmessage     ! error message of downwind routine

 ! initialize error control
 ierr=0; message='defineFile/'

 associate (dim_seg  => meta_qDims(ixQdims%seg)%dimName,    &
            dim_hru  => meta_qDims(ixQdims%hru)%dimName,    &
            dim_ens  => meta_qDims(ixQdims%ens)%dimName,    &
            dim_time => meta_qDims(ixQdims%time)%dimName)

! populate q dimension meta (not sure if this should be done here...)
 meta_qDims(ixQdims%seg)%dimLength = nSeg
 meta_qDims(ixQdims%hru)%dimLength = nHRU
 meta_qDims(ixQdims%ens)%dimLength = nEns

 ! --------------------
 ! define file
 ! --------------------
 ierr = nf90_create(trim(fname),nf90_classic_model,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 do jDim =1,nQdims
   if (jDim ==ixQdims%time) then ! time dimension (unlimited)
    ierr = nf90_def_dim(ncid, trim(meta_qDims(jDim)%dimName), nf90_unlimited, meta_qDims(jDim)%dimId)
   else
    ierr = nf90_def_dim(ncid, trim(meta_qDims(jDim)%dimName), meta_qDims(jDim)%dimLength ,meta_qDims(jDim)%dimId)
   endif
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 end do

 ! define coordinate variable for time
 call defvar(ncid, trim(dim_time), (/dim_time/), nf90_double, ierr, cmessage, vdesc=trim(dim_time), vunit=trim(units_time), vcal=calendar)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! define variables
 do iVar=1,nVars
  ! define variable
  select case(iVar)
   ! define network topology (integers)
   case( 1); call defvar(ncid, 'basinID',           (/dim_hru/),          nf90_int,   ierr,cmessage, vdesc='basin ID',                            vunit='-'   )
   case( 2); call defvar(ncid, 'reachID',           (/dim_seg/),          nf90_int,   ierr,cmessage, vdesc='reach ID',                            vunit='-'   )
   ! define runoff variables (double precision)
   case( 3); call defvar(ncid, 'basRunoff',         (/dim_hru,dim_time/), nf90_float, ierr,cmessage, vdesc='basin runoff',                        vunit='m/s' )
   case( 4); call defvar(ncid, 'instRunoff',        (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='instantaneous runoff in each reach',  vunit='m3/s')
   case( 5); call defvar(ncid, 'dlayRunoff',        (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='delayed runoff in each reach',        vunit='m3/s')
   case( 6); call defvar(ncid, 'sumUpstreamRunoff', (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='sum of upstream runoff in each reach',vunit='m3/s')
   case( 7); call defvar(ncid, 'KWTroutedRunoff',   (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='KWT routed runoff in each reach',     vunit='m3/s')
   case( 8); call defvar(ncid, 'IRFroutedRunoff',   (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='IRF routed runoff in each reach',     vunit='m3/s')
   case default; ierr=20; message=trim(message)//'unable to identify variable index'; return
  end select
  ! check errors
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 end do

 end associate

 ! end definitions
 ierr = nf90_enddef(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 END SUBROUTINE defineFile

 ! *********************************************************************
 ! private subroutine: define variable attributes NetCDF file
 ! *********************************************************************
 SUBROUTINE defvar(ncid, vname, dimNames, ivtype, ierr, message, vdesc, vunit, vcal)
  ! input
  integer(i4b), intent(in)             :: ncid                   ! Input: netcdf fine ID
  character(*), intent(in)             :: vname                  ! Input: variable name
  character(*), intent(in)             :: dimNames(:)            ! Input: variable dimension names
  integer(i4b), intent(in)             :: ivtype                 ! Input: variable type
  character(*), intent(in), optional   :: vdesc                  ! Input: variable description
  character(*), intent(in), optional   :: vunit                  ! Input: variable units
  character(*), intent(in), optional   :: vcal                   ! Input: calendar (if time variable)
  ! output
  integer(i4b), intent(out)            :: ierr                   ! error code
  character(*), intent(out)            :: message                ! error message
  ! local
  character(len=strLen)                :: calendar_str           ! calendar string
  character(len=strLen)                :: unit_str               ! unit string
  character(len=strLen)                :: desc_str               ! long_name string
  integer(i4b)                         :: id                     ! loop through dimensions
  integer(i4b)                         :: dimIDs(size(dimNames)) ! vector of dimension IDs
  integer(i4b)                         :: iVarId                 ! variable ID

  ! define dimension IDs
  do id=1,size(dimNames)
   ierr=nf90_inq_dimid(ncid,trim(dimNames(id)),dimIDs(id))
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end do

  ! define variable
  ierr = nf90_def_var(ncid,trim(vname),ivtype,dimIds,iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  if (present(vdesc)) then ! add long_name
    desc_str = trim(vdesc)
    ierr = nf90_put_att(ncid,iVarId,'long_name',trim(desc_str))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end if

  if (present(vunit)) then ! add variable unit
    unit_str = trim(vunit)
    ierr = nf90_put_att(ncid,iVarId,'units',trim(unit_str))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end if

  if (present(vcal)) then ! add time calendar
    calendar_str = trim(vcal)
    ierr = nf90_put_att(ncid,iVarId,'calendar',trim(calendar_str))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end if

 END SUBROUTINE defvar


END MODULE write_simoutput
