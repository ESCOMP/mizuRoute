MODULE write_restart

! Moudle wide external modules
USE nrtype, ONLY: i4b, dp, lgt, strLen
USE public_var
USE dataTypes, ONLY: time
USE io_netcdf, ONLY: ncd_int
USE io_netcdf, ONLY: ncd_float, ncd_double
USE io_netcdf, ONLY: ncd_unlimited
USE io_netcdf, only: def_nc                 ! define netcdf
USE io_netcdf, ONLY: def_var                ! define netcdf variable
USE io_netcdf, ONLY: def_dim                ! define netcdf dimension
USE io_netcdf, ONLY: end_def                ! end defining netcdf
USE io_netcdf, only: open_nc                ! open netcdf
USE io_netcdf, ONLY: close_nc               ! close netcdf
USE io_netcdf, ONLY: write_nc

implicit none

integer(i4b),  parameter :: currTimeStep = 1
integer(i4b),  parameter :: nextTimeStep = 2

private

public::main_restart

CONTAINS

 ! *********************************************************************
 ! public subroutine: restart write main routine
 ! *********************************************************************
 SUBROUTINE main_restart(ierr, message)

  USE globalData, ONLY: restartAlarm   ! logical to make alarm for restart writing

  implicit none
  ! output variables
  integer(i4b),   intent(out)          :: ierr             ! error code
  character(*),   intent(out)          :: message          ! error message
  ! local variables
  character(len=strLen)                :: cmessage         ! error message of downwind routine

  ierr=0; message='main_restart/'

  call restart_alarm(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (restartAlarm) then
    call restart_output(ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

 END SUBROUTINE main_restart


 ! *********************************************************************
 ! private subroutine: restart alarming
 ! *********************************************************************
 SUBROUTINE restart_alarm(ierr, message)

   USE public_var,        ONLY: calendar
   USE public_var,        ONLY: restart_write  ! restart write options
   USE public_var,        ONLY: restart_day
   USE globalData,        ONLY: restartAlarm   ! logical to make alarm for restart writing
   USE globalData,        ONLY: restCal        ! restart Calendar time
   USE globalData,        ONLY: dropCal        ! restart drop off Calendar time
   USE globalData,        ONLY: modTime        ! previous and current model time
   ! external routine
   USE time_utils_module, ONLY: ndays_month    ! compute number of days in a month

   implicit none

   ! output
   integer(i4b),   intent(out)          :: ierr             ! error code
   character(*),   intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)                :: cmessage         ! error message of downwind routine
   integer(i4b)                         :: nDays            ! number of days in a month

   ierr=0; message='restart_alarm/'

   ! adjust restart dropoff day if the dropoff day is outside number of days in particular month
   dropCal%id=restart_day
   call ndays_month(modTime(1)%iy, modTime(1)%im, calendar, nDays, ierr, cmessage)
   if (dropCal%id > nDays) then
     dropCal%id=nDays
   end if

   ! adjust dropoff day further if restart day is actually outside number of days in a particular month
   if (restCal%id > nDays) then
     dropCal%id=dropCal%id-1
   end if

   select case(trim(restart_write))
     case('Specified','specified','Last','last')
       restartAlarm = (dropCal%iy==modTime(1)%iy .and. dropCal%im==modTime(1)%im .and. dropCal%id==modTime(1)%id .and. &
                       dropCal%ih==modTime(1)%ih .and. dropCal%imin==modTime(1)%imin .and. nint(dropCal%dsec)==nint(modTime(1)%dsec))
     case('Annual','annual')
       restartAlarm = (dropCal%im==modTime(1)%im .and. dropCal%id==modTime(1)%id .and. &
                       dropCal%ih==modTime(1)%ih .and. dropCal%imin==modTime(1)%imin .and. nint(dropCal%dsec)==nint(modTime(1)%dsec))
     case('Monthly','monthly')
       restartAlarm = (dropCal%id==modTime(1)%id .and. &
                       dropCal%ih==modTime(1)%ih .and. dropCal%imin==modTime(1)%imin .and. nint(dropCal%dsec)==nint(modTime(1)%dsec))
     case('Daily','daily')
       restartAlarm = (dropCal%ih==modTime(1)%ih .and. dropCal%imin==modTime(1)%imin .and. nint(dropCal%dsec)==nint(modTime(1)%dsec))
     case('Never','never')
       restartAlarm = .false.
     case default
       ierr=20; message=trim(message)//'Current accepted <restart_write> options: L[l]ast, N[n]ever, S[s]pecified, Annual, Monthly, or Daily '; return
   end select

 END SUBROUTINE restart_alarm


 ! *********************************************************************
 ! private subroutine: write restart netCDF
 ! *********************************************************************
 SUBROUTINE restart_output(ierr, message)

  USE public_var, ONLY: routOpt
  USE public_var, ONLY: time_units
  USE public_var, ONLY: dt
  USE globalData, ONLY: runoff_data    ! runoff data for one time step for LSM HRUs and River network HRUs
  USE globalData, ONLY: TSEC
  USE globalData, ONLY: reachID
  USE globalData, ONLY: modTime           ! previous and current model time

  implicit none

  ! output variables
  integer(i4b),   intent(out)          :: ierr             ! error code
  character(*),   intent(out)          :: message          ! error message
  ! local variables
  real(dp)                             :: TSEC1, TSEC2
  character(len=strLen)                :: cmessage         ! error message of downwind routine
  character(len=strLen)                :: fnameRestart     ! name of the restart file name
  character(len=50),parameter          :: fmtYMDHMS='(2a,I0.4,a,I0.2,a,I0.2,x,I0.2,a,I0.2,a,I0.2)'

  ierr=0; message='restart_output/'

  write(iulog,fmtYMDHMS) new_line('a'),'Write restart file at ', &
                         modTime(1)%iy,'-',modTime(1)%im, '-', modTime(1)%id, modTime(1)%ih,':',modTime(1)%imin,':',nint(modTime(1)%dsec)

  call restart_fname(fnameRestart, nextTimeStep, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call define_state_nc(fnameRestart, time_units, routOpt, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! update model time step bound
  TSEC1 = TSEC(0) + dt
  TSEC2 = TSEC1   + dt

  call write_state_nc(fnameRestart,                            &  ! Input: state netcdf name
                      routOpt,                                 &  ! input: which routing options
                      runoff_data%time, 1, TSEC1, TSEC2,       &  ! Input: time, time step, start and end time [sec]
                      reachID,                                 &  ! Input: segment id vector
                      ierr, message)                              ! Output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE restart_output


 ! *********************************************************************
 ! private subroutine: define restart NetCDF file name
 ! *********************************************************************
 SUBROUTINE restart_fname(fnameRestart, timeStamp, ierr, message)

   USE public_var,          ONLY: restart_dir
   USE public_var,          ONLY: case_name        ! simulation name ==> output filename head
   USE public_var,          ONLY: calendar
   USE public_var,          ONLY: secprday
   USE public_var,          ONLY: dt
   USE globalData,          ONLY: modJulday        ! current model Julian day
   USE process_time_module, ONLY: conv_julian2cal  ! compute data and time from julian day

   implicit none

   ! input
   integer(i4b),   intent(in)           :: timeStamp        ! optional:
   ! output
   character(*),   intent(out)          :: fnameRestart     ! name of the restart file name
   integer(i4b),   intent(out)          :: ierr             ! error code
   character(*),   intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)                :: cmessage         ! error message of downwind routine
   real(dp)                             :: timeStampJulday  ! Julidan days corresponding to file name time stamp
   integer(i4b)                         :: sec_in_day       ! second within day
   type(time)                           :: timeStampCal     ! calendar date at next time step (for restart file name)
   character(len=50),parameter          :: fmtYMDS='(a,I0.4,a,I0.2,a,I0.2,a,I0.5,a)'

   ierr=0; message='restart_fname/'

   select case(timeStamp)
     case(currTimeStep); timeStampJulday = modJulday
     case(nextTimeStep); timeStampJulday = modJulday + dt/secprday
     case default;       ierr=20; message=trim(message)//'time stamp option in restart filename: invalid -> 1: current time Step or 2: next time step'; return
   end select

   call conv_julian2cal(timeStampJulday, calendar, timeStampCal, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   sec_in_day = timeStampCal%ih*60*60+timeStampCal%imin*60+nint(timeStampCal%dsec)

   write(fnameRestart, fmtYMDS) trim(restart_dir)//trim(case_name)//'.r.', &
                                timeStampCal%iy, '-', timeStampCal%im, '-', timeStampCal%id, '-',sec_in_day,'.nc'

 END SUBROUTINE restart_fname


 ! *********************************************************************
 ! subroutine: define restart NetCDF file
 ! *********************************************************************
 SUBROUTINE define_state_nc(fname,           &  ! input: filename
                            units_time,      &  ! input: time units
                            opt,             &  ! input: which routing options (state variables depends on routing options)
                            ierr, message)      ! output: error control
 ! External modules
 USE globalData, ONLY: meta_stateDims
 USE var_lookup, ONLY: ixStateDims, nStateDims
 implicit none
 ! input variables
 character(*),   intent(in)           :: fname            ! filename
 integer(i4b),   intent(in)           :: opt              ! routing option 0=all, 1=kwt, 2=irf
 character(*),   intent(in)           :: units_time       ! time units
 ! output variables
 integer(i4b),   intent(out)          :: ierr             ! error code
 character(*),   intent(out)          :: message          ! error message
 ! local variables
 integer(i4b)                         :: jDim             ! loop index for dimension
 integer(i4b)                         :: ncid             ! NetCDF file ID
 integer(i4b)                         :: ixDim_common(4)  ! custom dimension ID array
 character(len=strLen)                :: cmessage         ! error message of downwind routine

 ! initialize error control
 ierr=0; message='define_state_nc/'

 associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimName,     &
           dim_ens     => meta_stateDims(ixStateDims%ens)%dimName,     &
           dim_time    => meta_stateDims(ixStateDims%time)%dimName,    &
           dim_tbound  => meta_stateDims(ixStateDims%tbound)%dimName)

 ! Create file
 call def_nc(trim(fname), ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! For common dimension/variables - seg id, time, time-bound -----------
 ixDim_common = (/ixStateDims%seg, ixStateDims%ens, ixStateDims%time, ixStateDims%tbound/)

 ! Define dimensions
 do jDim = 1,size(ixDim_common)
  associate(ixDim_tmp => ixDim_common(jDim))
  if (meta_stateDims(ixDim_tmp)%dimLength == integerMissing) then
   call set_dim_len(ixDim_tmp, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//' for '//trim(meta_stateDims(ixDim_tmp)%dimName); return; endif
  endif
  call def_dim(ncid, meta_stateDims(ixDim_tmp)%dimName, meta_stateDims(ixDim_tmp)%dimLength, meta_stateDims(ixDim_tmp)%dimId, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end associate
 end do

 ! Define variable
 call def_var(ncid, 'reachID', (/dim_seg/), ncd_int, ierr, cmessage, vdesc='reach ID', vunit='-')
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call def_var(ncid, 'time ', (/dim_time/), ncd_double, ierr, cmessage, vdesc='time', vunit=units_time)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call def_var(ncid,'time_bound', (/dim_tbound, dim_time/), ncd_double, ierr, cmessage, vdesc='time bound at last time step', vunit='sec')
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end associate

 ! Routing specific variables --------------

 ! basin IRF
 if (doesBasinRoute == 1) then
  call define_IRFbas_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! KWT routing
 if (opt==allRoutingMethods .or. opt==kinematicWave) then
  call define_KWT_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! IRF routing
 if (opt==allRoutingMethods .or. opt==impulseResponseFunc) then
  call define_IRF_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! Finishing up definition -------
 ! end definitions
 call end_def(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! close NetCDF file
 call close_nc(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 CONTAINS

  SUBROUTINE set_dim_len(ixDim, ierr, message1)
   ! State/flux data structures
   USE globalData,   ONLY: meta_stateDims  ! states dimension meta
   USE globalData,   ONLY: nRch
   USE globalData,   ONLY: FRAC_FUTURE     ! To get size of q future for basin IRF
   implicit none
   ! input
   integer(i4b), intent(in)   :: ixDim    ! ixDim
   ! output
   integer(i4b), intent(out)  :: ierr     ! error code
   character(*), intent(out)  :: message1  ! error message

   ! initialize error control
   ierr=0; message1='set_dim_len/'

   select case(ixDim)
    case(ixStateDims%time);    meta_stateDims(ixStateDims%time)%dimLength    = ncd_unlimited
    case(ixStateDims%seg);     meta_stateDims(ixStateDims%seg)%dimLength     = nRch
    case(ixStateDims%ens);     meta_stateDims(ixStateDims%ens)%dimLength     = 1
    case(ixStateDims%tbound);  meta_stateDims(ixStateDims%tbound)%dimLength  = 2
    case(ixStateDims%tdh);     meta_stateDims(ixStateDims%tdh)%dimLength     = size(FRAC_FUTURE)
    case(ixStateDims%tdh_irf); meta_stateDims(ixStateDims%tdh_irf)%dimLength = 20   !just temporarily
    case(ixStateDims%wave);    meta_stateDims(ixStateDims%wave)%dimLength    = MAXQPAR
    case default; ierr=20; message1=trim(message1)//'unable to identify dimension variable index'; return
   end select

  END SUBROUTINE

  SUBROUTINE define_IRFbas_state(ierr, message1)
   ! External modules
   USE globalData, ONLY: meta_irf_bas
   USE var_lookup, ONLY: ixIRFbas, nVarsIRFbas
   implicit none
   ! output
   integer(i4b), intent(out)         :: ierr          ! error code
   character(*), intent(out)         :: message1      ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim    ! index loop for variables
   integer(i4b)                      :: nDims         ! number of dimensions
   character(len=strLen),allocatable :: dim_IRFbas(:) ! dimensions combination case 4

   ! initialize error control
   ierr=0; message1='define_IRFbas_state/'

   associate(dim_seg    => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens    => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_time   => meta_stateDims(ixStateDims%time)%dimName,    &
             dim_tdh    => meta_stateDims(ixStateDims%tdh)%dimName)

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%tdh)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh)%dimName); return; endif
   end if

   call def_dim(ncid, meta_stateDims(ixStateDims%tdh)%dimName, meta_stateDims(ixStateDims%tdh)%dimLength, meta_stateDims(ixStateDims%tdh)%dimId, ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsIRFbas

     nDims = size(meta_irf_bas(iVar)%varDim)
     if (allocated(dim_IRFbas)) then
       deallocate(dim_IRFbas)
     end if
     allocate(dim_IRFbas(nDims))
     do ixDim = 1, nDims
       dim_IRFbas(ixDim) = meta_stateDims(meta_irf_bas(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_irf_bas(iVar)%varName, dim_IRFbas, meta_irf_bas(iVar)%varType, ierr, cmessage, vdesc=meta_irf_bas(iVar)%varDesc, vunit=meta_irf_bas(iVar)%varUnit )
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   end do
  end associate

  END SUBROUTINE define_IRFbas_state

  SUBROUTINE define_KWT_state(ierr, message1)
   ! External modules
   USE globalData, ONLY: meta_kwt
   USE var_lookup, ONLY: ixKWT, nVarsKWT
   implicit none
   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   character(len=strLen),allocatable :: dim_kwt(:)  ! dimensions combination case 4

   ! initialize error control
   ierr=0; message1='define_KWT_state/'

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_time    => meta_stateDims(ixStateDims%time)%dimName,    &
             dim_wave    => meta_stateDims(ixStateDims%wave)%dimName)

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%wave)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%wave, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%wave)%dimName); return; endif
   end if

   ! Define dimension needed for this routing specific state variables
   call def_dim(ncid, meta_stateDims(ixStateDims%wave)%dimName, meta_stateDims(ixStateDims%wave)%dimLength, meta_stateDims(ixStateDims%wave)%dimId, ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   call def_var(ncid, 'numWaves', (/dim_seg,dim_ens,dim_time/), ncd_int, ierr, cmessage, vdesc='number of waves in a reach', vunit='-')
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsKWT

     nDims = size(meta_kwt(iVar)%varDim)
     if (allocated(dim_kwt)) then
       deallocate(dim_kwt)
     endif
     allocate(dim_kwt(nDims))
     do ixDim = 1, nDims
       dim_kwt(ixDim) = meta_stateDims(meta_kwt(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_kwt(iVar)%varName, dim_kwt, meta_kwt(iVar)%varType, ierr, cmessage, vdesc=meta_kwt(iVar)%varDesc, vunit=meta_kwt(iVar)%varUnit)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   end do

   end associate

  END SUBROUTINE define_KWT_state

  SUBROUTINE define_IRF_state(ierr, message1)
   ! External modules
   USE globalData, ONLY: meta_irf
   USE var_lookup, ONLY: ixIRF, nVarsIRF
   implicit none
   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   character(len=strLen),allocatable :: dim_irf(:)  ! dimensions combination case 4
   ! initialize error control
   ierr=0; message1='define_IRF_state/'

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_time    => meta_stateDims(ixStateDims%time)%dimName,    &
             dim_tdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimName)

   ! define dimension ID array

   if (meta_stateDims(ixStateDims%tdh_irf)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh_irf, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh_irf)%dimName); return; endif
   endif

   ! Define dimension needed for this routing specific state variables
   call def_dim(ncid, meta_stateDims(ixStateDims%tdh_irf)%dimName, meta_stateDims(ixStateDims%tdh_irf)%dimLength, meta_stateDims(ixStateDims%tdh_irf)%dimId, ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   call def_var(ncid, 'numQF', (/dim_seg,dim_ens,dim_time/), ncd_int, ierr, cmessage, vdesc='number of future q time steps in a reach', vunit='-')
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsIRF

     nDims = size(meta_irf(iVar)%varDim)
     if (allocated(dim_irf)) then
       deallocate(dim_irf)
     endif
     allocate(dim_irf(nDims))
     do ixDim = 1, nDims
       dim_irf(ixDim) = meta_stateDims(meta_irf(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_irf(iVar)%varName, dim_irf, meta_irf(iVar)%varType, ierr, cmessage, vdesc=meta_irf(iVar)%varDesc, vunit=meta_irf(iVar)%varUnit)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   end do

   end associate

  END SUBROUTINE define_IRF_state

 END SUBROUTINE define_state_nc


 ! *********************************************************************
 ! public subroutine: writing routing state NetCDF file
 ! *********************************************************************
 SUBROUTINE write_state_nc(fname,                &   ! Input: state netcdf name
                           opt,                  &   ! input: which routing options
                           time, iTime, T0, T1,  &   ! Input: time, time step, start and end time [sec]
                           seg_id,               &   ! Input: segment id vector
                           ierr, message)            ! Output: error control
 ! External module
 USE dataTypes,    ONLY: states
 ! meta data
 USE globalData,   ONLY: meta_stateDims  ! dimension for state variables
 ! Named variables
 USE var_lookup,   ONLY: ixStateDims, nStateDims
 implicit none
 ! input variables
 character(*), intent(in)        :: fname           ! filename
 integer(i4b), intent(in)        :: opt             ! routing option 0=all, 1=kwt, 2=irf
 real(dp),     intent(in)        :: time            ! calendar time
 integer(i4b), intent(in)        :: iTime           ! ith Time step
 real(dp),     intent(in)        :: T0              ! beginning time [sec] of ith time step - lapse time from the beginning of the simulation
 real(dp),     intent(in)        :: T1              ! ending time [sec] ith time step - lapse time from the beginning of the simulation
 integer(i4b), intent(in)        :: seg_id(:)       ! segment id vector
 ! output variables
 integer(i4b), intent(out)       :: ierr            ! error code
 character(*), intent(out)       :: message         ! error message
 ! local variables
 integer(i4b)                    :: ncid            ! netCDF ID
 type(states)                    :: state(0:2)      ! temporal state data structures -currently 2 river routing scheme + basin IRF routing
 character(len=strLen)           :: cmessage        ! error message of downwind routine

 ! initialize error control
 ierr=0; message='write_state_nc/'

 ! -- open netCDF
 call open_nc(fname, 'w', ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! -- Write out to netCDF
 ! Miscellaneous variables - seg id, time etc
 call write_nc(ncid,'reachID', seg_id, (/1/), (/size(seg_id)/), ierr, cmessage);
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_nc(ncid,'time', (/time/), (/iTime/), (/1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_nc(ncid,'time_bound', (/T0,T1/), (/1,iTime/), (/2,1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 if (doesBasinRoute == 1) then
  call write_IRFbas_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (opt==allRoutingMethods .or. opt==impulseResponseFunc)then
  call write_IRF_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (opt==allRoutingMethods .or. opt==kinematicWave)then
  call write_KWT_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! -- close netCDF
 call close_nc(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 CONTAINS

  ! Basin IRF writing procedures
  SUBROUTINE write_IRFbas_state(ierr, message1)
  ! external module
  USE globalData,   ONLY: meta_irf_bas
  USE var_lookup,   ONLY: ixIRFbas, nVarsIRFbas
  USE globalData,   ONLY: RCHFLX          ! To get q future for basin IRF and IRF (these should not be in this data strucuture)
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

  ! initialize error control
  ierr=0; message1='write_IRFbas_state/'

  associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,  &
            nens     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            ntdh     => meta_stateDims(ixStateDims%tdh)%dimLength)      ! maximum future q time steps among basins

  allocate(state(0)%var(nVarsIRFbas), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRFbas

   select case(iVar)
    case(ixIRFbas%qfuture); allocate(state(0)%var(iVar)%array_3d_dp(nSeg, ntdh, nens), stat=ierr)
    case(ixIRFbas%q);       allocate(state(0)%var(iVar)%array_2d_dp(nSeg, nens),       stat=ierr)
    case default; ierr=20; message1=trim(message1)//'unable to identify basin routing variable index'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//'problem allocating space for basin IRF routing state '//trim(meta_irf_bas(iVar)%varName); return; endif

  end do

 ! --Convert data structures to arrays
  do iens=1,nens
   do iSeg=1,nSeg
     do iVar=1,nVarsIRFbas

      select case(iVar)
       case(ixIRFbas%qfuture); state(0)%var(iVar)%array_3d_dp(iSeg,:,iens) = RCHFLX(iens,iSeg)%QFUTURE
       case(ixIRFbas%q);       state(0)%var(iVar)%array_2d_dp(iSeg,iens)   = RCHFLX(iens,iSeg)%BASIN_QR(1)
       case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF state variable index'; return
      end select

    enddo
   enddo
  enddo

  do iVar=1,nVarsIRFbas

   select case(iVar)
    case(ixIRFbas%q);       call write_nc(ncid, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_2d_dp, (/1,1,iTime/), (/nSeg,nens,1/), ierr, cmessage)
    case(ixIRFbas%qfuture); call write_nc(ncid, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_3d_dp, (/1,1,1,iTime/), (/nSeg,ntdh,nens,1/), ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF variable index for nc writing'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  enddo

  end associate

  END SUBROUTINE write_IRFbas_state

  ! KWT writing procedures
  SUBROUTINE write_KWT_state(ierr, message1)
  ! External module
  USE globalData,   ONLY: meta_kwt
  USE var_lookup,   ONLY: ixKWT, nVarsKWT
  USE globalData,   ONLY: KROUTE
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively
  integer(i4b), allocatable  :: RFvec(:)        ! temporal vector
  integer(i4b), allocatable  :: numWaves(:,:)   ! number of waves for each ensemble and segment
  ! initialize error control
  ierr=0; message1='write_KWT_state/'

  associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,  &
            nens     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            nwave    => meta_stateDims(ixStateDims%wave)%dimLength)     ! maximum waves allowed in a reach

  allocate(state(kinematicWave)%var(nVarsKWT), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! array to store number of wave per segment and ensemble
  allocate(numWaves(nens,nSeg), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT
    select case(iVar)
     case(ixKWT%routed); allocate(state(kinematicWave)%var(iVar)%array_3d_int(nSeg, nwave, nens), stat=ierr)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      allocate(state(kinematicWave)%var(iVar)%array_3d_dp(nSeg, nwave, nens), stat=ierr)
     case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWT routing state '//trim(meta_kwt(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nens
   do iSeg=1,nSeg

    numWaves(iens,iseg) = size(KROUTE(iens,iseg)%KWAVE)

    do iVar=1,nVarsKWT

     select case(iVar)
      case(ixKWT%tentry)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = KROUTE(iens,iSeg)%KWAVE(:)%TI
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%texit)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = KROUTE(iens,iSeg)%KWAVE(:)%TR
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%qwave)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = KROUTE(iens,iSeg)%KWAVE(:)%QF
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%qwave_mod)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = KROUTE(iens,iSeg)%KWAVE(:)%QM
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%routed) ! this is suppposed to be logical variable, but put it as 0 or 1 in double now
       if (allocated(RFvec)) deallocate(RFvec, stat=ierr)
       allocate(RFvec(numWaves(iens,iSeg)),stat=ierr); RFvec=0_i4b
       where (KROUTE(iens,iSeg)%KWAVE(:)%RF) RFvec=1_i4b
       state(kinematicWave)%var(iVar)%array_3d_int(iSeg,1:numWaves(iens,iSeg),iens) = RFvec
       state(kinematicWave)%var(iVar)%array_3d_int(iSeg,numWaves(iens,iSeg)+1:,iens) = integerMissing
      case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
     end select

    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! Writing netCDF
  call write_nc(ncid, 'numWaves', numWaves, (/1,1,iTime/), (/nSeg,nens,1/), ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT

    select case(iVar)
     case(ixKWT%routed)
       call write_nc(ncid, trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_int, (/1,1,1,iTime/), (/nSeg,nwave,nens,1/), ierr, cmessage)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      call write_nc(ncid, trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_dp, (/1,1,1,iTime/), (/nSeg,nwave,nens,1/), ierr, cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc writing'; return
    end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  end do

  end associate

  END SUBROUTINE write_KWT_state


  ! IRF writing procedures
  SUBROUTINE write_IRF_state(ierr, message1)
  ! external module
  USE globalData,   ONLY: meta_irf        ! IRF routing
  USE var_lookup,   ONLY: ixIRF, nVarsIRF
  USE globalData,   ONLY: RCHFLX          ! To get q future for basin IRF and IRF (these should not be in this data strucuture)
  USE globalData,   ONLY: NETOPO          ! To get UH (this should not be in this data strucuture)
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively
  integer(i4b), allocatable  :: numQF(:,:)      ! number of future Q time steps for each ensemble and segment
  ! initialize error control
  ierr=0; message1='write_IRF_state/'

  associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,  &
            nens     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            ntdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimLength)    ! maximum future q time steps among reaches

  allocate(state(impulseResponseFunc)%var(nVarsIRF), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! array to store number of wave per segment and ensemble
  allocate(numQF(nens,nSeg), stat=ierr)
  if(ierr/=0)then; message1=trim(message1)//'problem allocating space for numQF'; return; endif

  do iVar=1,nVarsIRF
   select case(iVar)
    case(ixIRF%qfuture); allocate(state(impulseResponseFunc)%var(iVar)%array_3d_dp(nSeg, ntdh_irf, nens), stat=ierr)
    case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//'problem allocating space for IRF routing state '//trim(meta_irf(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nens
   do iSeg=1,nSeg

    numQF(iens,iseg) = size(NETOPO(iSeg)%UH)

    do iVar=1,nVarsIRF

     select case(iVar)
      case(ixIRF%qfuture)
       state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,1:numQF(iens,iSeg),iens) = RCHFLX(iens,iSeg)%QFUTURE_IRF
       state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,numQF(iens,iSeg)+1:ntdh_irf,iens) = realMissing
      case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
     end select

    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! writing netcdf
  call write_nc(ncid, 'numQF', numQF, (/1,1,iTime/), (/nSeg,nens,1/), ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRF

   select case(iVar)
    case(ixIRF%qfuture)
     call write_nc(ncid, trim(meta_irf(iVar)%varName), state(impulseResponseFunc)%var(iVar)%array_3d_dp, (/1,1,1,iTime/), (/nSeg,ntdh_irf,nens,1/), ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc writing'; return
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end select

  end do

  end associate

  END SUBROUTINE write_IRF_state

 END SUBROUTINE write_state_nc

END MODULE write_restart
