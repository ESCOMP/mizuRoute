MODULE write_restart

USE nrtype
USE date_time, ONLY: datetime
USE io_netcdf, ONLY: ncd_int
USE io_netcdf, ONLY: ncd_float, ncd_double
USE io_netcdf, ONLY: ncd_unlimited
USE io_netcdf, ONLY: def_nc                 ! define netcdf
USE io_netcdf, ONLY: def_var                ! define netcdf variable
USE io_netcdf, ONLY: put_global_attr        ! write global attribute
USE io_netcdf, ONLY: def_dim                ! define netcdf dimension
USE io_netcdf, ONLY: end_def                ! end defining netcdf
USE io_netcdf, ONLY: open_nc                ! open netcdf
USE io_netcdf, ONLY: close_nc               ! close netcdf
USE io_netcdf, ONLY: write_nc
USE globalData, ONLY: onRoute               ! logical to indicate which routing method(s) is on
USE public_var, ONLY: iulog                 ! i/o logical unit number
USE public_var, ONLY: integerMissing
USE public_var, ONLY: realMissing
USE public_var, ONLY: dt
USE public_var, ONLY: doesBasinRoute
USE public_var, ONLY: impulseResponseFunc
USE public_var, ONLY: kinematicWaveTracking
USE public_var, ONLY: kinematicWave
USE public_var, ONLY: muskingumCunge
USE public_var, ONLY: diffusiveWave

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

   USE ascii_util_module, ONLY: lower
   USE public_var,        ONLY: calendar
   USE public_var,        ONLY: restart_write  ! restart write options
   USE public_var,        ONLY: restart_day
   USE globalData,        ONLY: restartAlarm   ! logical to make alarm for restart writing
   USE globalData,        ONLY: restCal        ! restart Calendar time
   USE globalData,        ONLY: dropCal        ! restart drop off Calendar time
   USE globalData,        ONLY: modTime        ! previous and current model time

   implicit none

   ! output
   integer(i4b),   intent(out)          :: ierr             ! error code
   character(*),   intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)                :: cmessage         ! error message of downwind routine
   integer(i4b)                         :: nDays            ! number of days in a month

   ierr=0; message='restart_alarm/'

   ! adjust restart dropoff day if the dropoff day is outside number of days in particular month
   call dropCal%set_datetime(dropCal%year(), dropCal%month(), restart_day, dropCal%hour(), dropCal%minute(), dropCal%sec())
   nDays = modTime(1)%ndays_month(calendar, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   if (dropCal%day() > nDays) then
     call dropCal%set_datetime(dropCal%year(), dropCal%month(), nDays, dropCal%hour(), dropCal%minute(), dropCal%sec())
   end if

   ! adjust dropoff day further if restart day is actually outside number of days in a particular month
   if (restCal%day() > nDays) then
     dropCal = dropCal%add_day(-1, calendar, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   select case(lower(trim(restart_write)))
     case('specified','last')
       restartAlarm = (dropCal==modTime(1))
     case('yearly')
       restartAlarm = (dropCal%is_equal_mon(modTime(1)) .and. dropCal%is_equal_day(modTime(1)) .and. dropCal%is_equal_time(modTime(1)))
     case('monthly')
       restartAlarm = (dropCal%is_equal_day(modTime(1)) .and. dropCal%is_equal_time(modTime(1)))
     case('daily')
       restartAlarm = dropCal%is_equal_time(modTime(1))
     case('never')
       restartAlarm = .false.
     case default
       ierr=20; message=trim(message)//'Accepted <restart_write> options (case insensitive): last, never, specified, yearly, monthly, or daily '; return
   end select

 END SUBROUTINE restart_alarm


 ! *********************************************************************
 ! private subroutine: write restart netCDF
 ! *********************************************************************
 SUBROUTINE restart_output(ierr, message)

  USE globalData, ONLY: TSEC
  USE globalData, ONLY: reachID

  implicit none

  ! output variables
  integer(i4b),   intent(out)          :: ierr             ! error code
  character(*),   intent(out)          :: message          ! error message
  ! local variables
  real(dp)                             :: TSEC1, TSEC2
  character(len=strLen)                :: cmessage         ! error message of downwind routine
  character(len=strLen)                :: fnameRestart     ! name of the restart file name

  ierr=0; message='restart_output/'

  call restart_fname(fnameRestart, nextTimeStep, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call define_state_nc(fnameRestart, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! update model time step bound
  TSEC1 = TSEC(0) + dt
  TSEC2 = TSEC1   + dt

  call write_state_nc(fnameRestart,                            &  ! Input: state netcdf name
                      TSEC1, TSEC2,                            &  ! Input: time, time step, start and end time [sec]
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
   USE globalData,          ONLY: modTime          ! current model datetime

   implicit none

   ! input
   integer(i4b),   intent(in)           :: timeStamp        ! optional:
   ! output
   character(*),   intent(out)          :: fnameRestart     ! name of the restart file name
   integer(i4b),   intent(out)          :: ierr             ! error code
   character(*),   intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)                :: cmessage         ! error message of downwind routine
   type(datetime)                       :: timeStampCal     ! datetime corresponding to file name time stamp
   integer(i4b)                         :: sec_in_day       ! second within day
   character(len=50),parameter          :: fmtYMDS='(a,I0.4,a,I0.2,a,I0.2,a,I0.5,a)'
   character(len=50),parameter          :: fmtYMDHMS='(2a,I0.4,a,I0.2,a,I0.2,x,I0.2,a,I0.2,a,I0.2)'

   ierr=0; message='restart_fname/'

   select case(timeStamp)
     case(currTimeStep); timeStampCal = modTime(1)
     case(nextTimeStep); timeStampCal = modTime(1)%add_sec(dt, calendar, ierr, cmessage)
     case default;       ierr=20; message=trim(message)//'time stamp option in restart filename: invalid -> 1: current time Step or 2: next time step'; return
   end select

   write(iulog,fmtYMDHMS) new_line('a'),'Write restart file for ', &
                          timeStampCal%year(),'-',timeStampCal%month(),'-',timeStampCal%day(),timeStampCal%hour(),':',timeStampCal%minute(),':',nint(timeStampCal%sec())

   sec_in_day = timeStampCal%hour()*60*60+timeStampCal%minute()*60+nint(timeStampCal%sec())

   write(fnameRestart, fmtYMDS) trim(restart_dir)//trim(case_name)//'.r.', &
                                timeStampCal%year(), '-', timeStampCal%month(), '-', timeStampCal%day(), '-',sec_in_day,'.nc'

 END SUBROUTINE restart_fname


 ! *********************************************************************
 ! subroutine: define restart NetCDF file
 ! *********************************************************************
 SUBROUTINE define_state_nc(fname,           &  ! input: filename
                            ierr, message)      ! output: error control

 USE globalData, ONLY: meta_stateDims
 USE globalData, ONLY: modTime                 ! current model datetime
 USE public_var, ONLY: calendar
 USE var_lookup, ONLY: ixStateDims, nStateDims

 implicit none

 ! input variables
 character(*),   intent(in)           :: fname            ! filename
 ! output variables
 integer(i4b),   intent(out)          :: ierr             ! error code
 character(*),   intent(out)          :: message          ! error message
 ! local variables
 type(datetime)                       :: timeStampCal     ! datetime corresponding to file name time stamp
 character(len=50),parameter          :: fmtYMDHMS='(I0.4,a,I0.2,a,I0.2,x,I0.2,a,I0.2,a,I0.2)'
 character(len=strLen)                :: globalDesc       ! global attributes: description
 integer(i4b)                         :: jDim             ! loop index for dimension
 integer(i4b)                         :: ncid             ! NetCDF file ID
 integer(i4b)                         :: ixDim_common(3)  ! custom dimension ID array
 character(len=strLen)                :: cmessage         ! error message of downwind routine

 ierr=0; message='define_state_nc/'

 associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimName,     &
           dim_ens     => meta_stateDims(ixStateDims%ens)%dimName,     &
           dim_tbound  => meta_stateDims(ixStateDims%tbound)%dimName)

 ! Create file
 call def_nc(trim(fname), ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! For common dimension/variables - seg id, time, time-bound -----------
 ixDim_common = (/ixStateDims%seg, ixStateDims%ens, ixStateDims%tbound/)

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

 call def_var(ncid,'time_bound', (/dim_tbound/), ncd_double, ierr, cmessage, vdesc='time bound at last time step', vunit='sec')
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end associate

 ! Write global attribute
 call put_global_attr(ncid, 'Title', 'mizuRoute restart file', ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get which time is for restarting
 timeStampCal = modTime(1)%add_sec(dt, calendar, ierr, cmessage)
 write(globalDesc, fmtYMDHMS) timeStampCal%year(),'-',timeStampCal%month(),'-',timeStampCal%day(),timeStampCal%hour(),':',timeStampCal%minute(),':',nint(timeStampCal%sec())

 call put_global_attr(ncid, 'Restart time', trim(globalDesc), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! previous-time step hru inflow into reach
 call define_basinQ_state(ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! Routing specific variables --------------
 if (doesBasinRoute == 1) then
  call define_IRFbas_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(kinematicWaveTracking)) then
   call define_KWT_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(impulseResponseFunc))then
   call define_IRF_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(kinematicWave)) then
   call define_KW_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(muskingumCunge)) then
   call define_MC_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(diffusiveWave)) then
   call define_DW_state(ierr, cmessage)
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

   USE globalData,   ONLY: meta_stateDims  ! states dimension meta
   USE globalData,   ONLY: nRch
   USE globalData,   ONLY: nMolecule
   USE public_var,   ONLY: MAXQPAR
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
     case(ixStateDims%tdh_irf); meta_stateDims(ixStateDims%tdh_irf)%dimLength = 50   !just temporarily
     case(ixStateDims%wave);    meta_stateDims(ixStateDims%wave)%dimLength    = MAXQPAR
     case(ixStateDims%mol_kw);  meta_stateDims(ixStateDims%mol_kw)%dimLength  = nMolecule%KW_ROUTE
     case(ixStateDims%mol_mc);  meta_stateDims(ixStateDims%mol_mc)%dimLength  = nMolecule%MC_ROUTE
     case(ixStateDims%mol_dw);  meta_stateDims(ixStateDims%mol_dw)%dimLength  = nMolecule%DW_ROUTE
     case default; ierr=20; message1=trim(message1)//'unable to identify dimension variable index'; return
   end select

  END SUBROUTINE


  SUBROUTINE define_basinQ_state(ierr, message1)

   USE globalData, ONLY: meta_basinQ
   USE var_lookup, ONLY: ixBasinQ, nVarsBasinQ

   implicit none

   ! output
   integer(i4b), intent(out)         :: ierr          ! error code
   character(*), intent(out)         :: message1      ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim    ! index loop for variables
   integer(i4b)                      :: nDims         ! number of dimensions
   character(len=strLen),allocatable :: dim_set(:)    ! dimensions combination

   ! initialize error control
   ierr=0; message1='define_basinQ_state/'

   associate(dim_seg    => meta_stateDims(ixStateDims%seg)%dimName, &
             dim_ens    => meta_stateDims(ixStateDims%ens)%dimName)

   do iVar=1,nVarsBasinQ
     nDims = size(meta_basinQ(iVar)%varDim)
     if (allocated(dim_set)) deallocate(dim_set)
     allocate(dim_set(nDims))

     do ixDim = 1, nDims
       dim_set(ixDim) = meta_stateDims(meta_basinQ(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_basinQ(iVar)%varName, dim_set, meta_basinQ(iVar)%varType, ierr, cmessage, vdesc=meta_basinQ(iVar)%varDesc, vunit=meta_basinQ(iVar)%varUnit )
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do
  end associate

  END SUBROUTINE define_basinQ_state


  SUBROUTINE define_IRFbas_state(ierr, message1)

   USE globalData, ONLY: meta_irf_bas
   USE var_lookup, ONLY: ixIRFbas, nVarsIRFbas

   implicit none

   ! output
   integer(i4b), intent(out)         :: ierr          ! error code
   character(*), intent(out)         :: message1      ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim    ! index loop for variables
   integer(i4b)                      :: nDims         ! number of dimensions
   character(len=strLen),allocatable :: dim_set(:)    ! dimensions combination case 4

   ierr=0; message1='define_IRFbas_state/'

   associate(dim_seg    => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens    => meta_stateDims(ixStateDims%ens)%dimName,     &
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
     if (allocated(dim_set)) deallocate(dim_set)
     allocate(dim_set(nDims))

     do ixDim = 1, nDims
       dim_set(ixDim) = meta_stateDims(meta_irf_bas(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_irf_bas(iVar)%varName, dim_set, meta_irf_bas(iVar)%varType, ierr, cmessage, vdesc=meta_irf_bas(iVar)%varDesc, vunit=meta_irf_bas(iVar)%varUnit )
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do
  end associate

  END SUBROUTINE define_IRFbas_state


  SUBROUTINE define_KWT_state(ierr, message1)

   USE globalData, ONLY: meta_kwt
   USE var_lookup, ONLY: ixKWT, nVarsKWT

   implicit none

   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   character(len=strLen),allocatable :: dim_set(:)  ! dimensions combination case 4

   ierr=0; message1='define_KWT_state/'

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_wave    => meta_stateDims(ixStateDims%wave)%dimName)

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%wave)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%wave, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%wave)%dimName); return; endif
   end if

   ! Define dimension needed for this routing specific state variables
   call def_dim(ncid, meta_stateDims(ixStateDims%wave)%dimName, meta_stateDims(ixStateDims%wave)%dimLength, meta_stateDims(ixStateDims%wave)%dimId, ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   call def_var(ncid, 'numWaves', (/dim_seg,dim_ens/), ncd_int, ierr, cmessage, vdesc='number of waves in a reach', vunit='-')
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsKWT
     nDims = size(meta_kwt(iVar)%varDim)
     if (allocated(dim_set)) deallocate(dim_set)
     allocate(dim_set(nDims))

     do ixDim = 1, nDims
       dim_set(ixDim) = meta_stateDims(meta_kwt(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_kwt(iVar)%varName, dim_set, meta_kwt(iVar)%varType, ierr, cmessage, vdesc=meta_kwt(iVar)%varDesc, vunit=meta_kwt(iVar)%varUnit)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do

   end associate

  END SUBROUTINE define_KWT_state


  SUBROUTINE define_IRF_state(ierr, message1)

   USE globalData, ONLY: meta_irf
   USE var_lookup, ONLY: ixIRF, nVarsIRF

   implicit none

   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   character(len=strLen),allocatable :: dim_set(:)  ! dimensions combination case 4

   ierr=0; message1='define_IRF_state/'

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_tdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimName)

   ! define dimension ID array
   if (meta_stateDims(ixStateDims%tdh_irf)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh_irf, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh_irf)%dimName); return; endif
   endif

   ! Define dimension needed for this routing specific state variables
   call def_dim(ncid, meta_stateDims(ixStateDims%tdh_irf)%dimName, meta_stateDims(ixStateDims%tdh_irf)%dimLength, meta_stateDims(ixStateDims%tdh_irf)%dimId, ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   call def_var(ncid, 'numQF', (/dim_seg,dim_ens/), ncd_int, ierr, cmessage, vdesc='number of future q time steps in a reach', vunit='-')
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsIRF
     nDims = size(meta_irf(iVar)%varDim)
     if (allocated(dim_set)) deallocate(dim_set)
     allocate(dim_set(nDims))

     do ixDim = 1, nDims
       dim_set(ixDim) = meta_stateDims(meta_irf(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_irf(iVar)%varName, dim_set, meta_irf(iVar)%varType, ierr, cmessage, vdesc=meta_irf(iVar)%varDesc, vunit=meta_irf(iVar)%varUnit)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do

   end associate

  END SUBROUTINE define_IRF_state


  SUBROUTINE define_KW_state(ierr, message1)

   USE globalData, ONLY: meta_kw
   USE var_lookup, ONLY: ixKW, nVarsKW

   implicit none

   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   character(len=strLen),allocatable :: dim_set(:)  ! dimensions combination

   ierr=0; message1='define_KW_state/'

   associate(dim_seg  => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens  => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_mesh => meta_stateDims(ixStateDims%mol_kw)%dimName)

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%mol_kw)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%mol_kw, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%mol_kw)%dimName); return; endif
   end if

   ! Define dimension needed for this routing specific state variables
   call def_dim(ncid, meta_stateDims(ixStateDims%mol_kw)%dimName,   &
                      meta_stateDims(ixStateDims%mol_kw)%dimLength, &
                      meta_stateDims(ixStateDims%mol_kw)%dimId, ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsKW
     nDims = size(meta_kw(iVar)%varDim)
     if (allocated(dim_set)) deallocate(dim_set)
     allocate(dim_set(nDims))

     do ixDim = 1, nDims
       dim_set(ixDim) = meta_stateDims(meta_kw(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_kw(iVar)%varName, dim_set, meta_kw(iVar)%varType, ierr, cmessage, vdesc=meta_kw(iVar)%varDesc, vunit=meta_kw(iVar)%varUnit)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do

   end associate

  END SUBROUTINE define_KW_state


  SUBROUTINE define_MC_state(ierr, message1)

   USE globalData, ONLY: meta_mc
   USE var_lookup, ONLY: ixMC, nVarsMC

   implicit none

   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   character(len=strLen),allocatable :: dim_set(:)   ! dimensions combination case 4

   ierr=0; message1='define_MC_state/'

   associate(dim_seg  => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens  => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_mesh => meta_stateDims(ixStateDims%mol_mc)%dimName)

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%mol_mc)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%mol_mc, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%mol_mc)%dimName); return; endif
   end if

   ! Define dimension needed for this routing specific state variables
   call def_dim(ncid, meta_stateDims(ixStateDims%mol_mc)%dimName,   &
                      meta_stateDims(ixStateDims%mol_mc)%dimLength, &
                      meta_stateDims(ixStateDims%mol_mc)%dimId, ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsMC
     nDims = size(meta_mc(iVar)%varDim)
     if (allocated(dim_set)) deallocate(dim_set)
     allocate(dim_set(nDims))

     do ixDim = 1, nDims
       dim_set(ixDim) = meta_stateDims(meta_mc(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_mc(iVar)%varName, dim_set, meta_mc(iVar)%varType, ierr, cmessage, vdesc=meta_mc(iVar)%varDesc, vunit=meta_mc(iVar)%varUnit)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do

   end associate

  END SUBROUTINE define_MC_state


  SUBROUTINE define_DW_state(ierr, message1)

   USE globalData, ONLY: meta_dw
   USE var_lookup, ONLY: ixDW, nVarsDW

   implicit none

   ! output
   integer(i4b), intent(out)         :: ierr        ! error code
   character(*), intent(out)         :: message1    ! error message
   ! local
   integer(i4b)                      :: iVar,ixDim  ! index loop for variables
   integer(i4b)                      :: nDims       ! number of dimensions
   character(len=strLen),allocatable :: dim_set(:)   ! dimensions combination case 4

   ierr=0; message1='define_DW_state/'

   associate(dim_seg  => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens  => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_mesh => meta_stateDims(ixStateDims%mol_dw)%dimName)

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%mol_dw)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%mol_dw, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%mol_dw)%dimName); return; endif
   end if

   ! Define dimension needed for this routing specific state variables
   call def_dim(ncid, meta_stateDims(ixStateDims%mol_dw)%dimName,   &
                      meta_stateDims(ixStateDims%mol_dw)%dimLength, &
                      meta_stateDims(ixStateDims%mol_dw)%dimId, ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsDW
     nDims = size(meta_dw(iVar)%varDim)
     if (allocated(dim_set)) deallocate(dim_set)
     allocate(dim_set(nDims))

     do ixDim = 1, nDims
       dim_set(ixDim) = meta_stateDims(meta_dw(iVar)%varDim(ixDim))%dimName
     end do

     call def_var(ncid, meta_dw(iVar)%varName, dim_set, meta_dw(iVar)%varType, ierr, cmessage, vdesc=meta_dw(iVar)%varDesc, vunit=meta_dw(iVar)%varUnit)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do

   end associate

  END SUBROUTINE define_DW_state

 END SUBROUTINE define_state_nc


 ! *********************************************************************
 ! public subroutine: writing routing state NetCDF file
 ! *********************************************************************
 SUBROUTINE write_state_nc(fname,                &   ! Input: state netcdf name
                           T0, T1,               &   ! Input: time, time step, start and end time [sec]
                           seg_id,               &   ! Input: segment id vector
                           ierr, message)            ! Output: error control

 USE dataTypes,    ONLY: states
 USE globalData,   ONLY: RCHFLX
 USE globalData,   ONLY: RCHSTA
 USE globalData,   ONLY: meta_stateDims  ! dimension meta for state variables
 USE var_lookup,   ONLY: ixStateDims, nStateDims

 implicit none

 ! input variables
 character(*), intent(in)        :: fname           ! filename
 real(dp),     intent(in)        :: T0              ! beginning time [sec] of ith time step - lapse time from the beginning of the simulation
 real(dp),     intent(in)        :: T1              ! ending time [sec] ith time step - lapse time from the beginning of the simulation
 integer(i4b), intent(in)        :: seg_id(:)       ! segment id vector
 ! output variables
 integer(i4b), intent(out)       :: ierr            ! error code
 character(*), intent(out)       :: message         ! error message
 ! local variables
 integer(i4b)                    :: ncid            ! netCDF ID
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

 call write_nc(ncid,'time_bound', (/T0,T1/), (/1/), (/2/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_basinQ_state(ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 if (doesBasinRoute == 1) then
   call write_IRFbas_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(impulseResponseFunc)) then
   call write_IRF_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(kinematicWaveTracking)) then
   call write_KWT_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(kinematicWave)) then
   call write_KW_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(muskingumCunge)) then
   call write_MC_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (onRoute(diffusiveWave)) then
   call write_DW_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! -- close netCDF
 call close_nc(ncid, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 CONTAINS

  ! reach inflow writing procedure
  SUBROUTINE write_basinQ_state(ierr, message1)

  USE globalData,   ONLY: meta_basinQ
  USE var_lookup,   ONLY: ixBasinQ, nVarsBasinQ

  implicit none

  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  type(states)               :: state           ! temporal state data structures
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

  ! initialize error control
  ierr=0; message1='write_basinQ_state/'

  associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,  &
            nens     => meta_stateDims(ixStateDims%ens)%dimLength)

  allocate(state%var(nVarsBasinQ), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsBasinQ
    select case(iVar)
      case(ixBasinQ%q); allocate(state%var(iVar)%array_2d_dp(nSeg, nens), stat=ierr)
      case default; ierr=20; message1=trim(message1)//'unable to identify basin routing variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for basin IRF routing state '//trim(meta_basinQ(iVar)%varName); return; endif
  end do

 ! --Convert data structures to arrays
  do iens=1,nens
    do iSeg=1,nSeg
      do iVar=1,nVarsBasinQ
        select case(iVar)
          case(ixBasinQ%q); state%var(iVar)%array_2d_dp(iSeg,iens) = RCHFLX(iens,iSeg)%BASIN_QR(1)
          case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF state variable index'; return
        end select
      end do
    end do
  end do

  do iVar=1,nVarsBasinQ
    select case(iVar)
      case(ixBasinQ%q); call write_nc(ncid, meta_basinQ(iVar)%varName, state%var(iVar)%array_2d_dp, (/1,1/), (/nSeg,nens/), ierr, cmessage)
      case default; ierr=20; message1=trim(message1)//'unable to identify reach inflow variable index for nc writing'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
  end do

  end associate

  END SUBROUTINE write_basinQ_state

  ! Basin IRF writing procedures
  SUBROUTINE write_IRFbas_state(ierr, message1)

  USE globalData,   ONLY: meta_irf_bas
  USE var_lookup,   ONLY: ixIRFbas, nVarsIRFbas

  implicit none

  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  type(states)               :: state           ! temporal state data structures
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

  ! initialize error control
  ierr=0; message1='write_IRFbas_state/'

  associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,  &
            nens     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            ntdh     => meta_stateDims(ixStateDims%tdh)%dimLength)      ! maximum future q time steps among basins

  allocate(state%var(nVarsIRFbas), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRFbas
    select case(iVar)
      case(ixIRFbas%qfuture); allocate(state%var(iVar)%array_3d_dp(nSeg, ntdh, nens), stat=ierr)
      case default; ierr=20; message1=trim(message1)//'unable to identify basin routing variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for basin IRF routing state '//trim(meta_irf_bas(iVar)%varName); return; endif
  end do

 ! --Convert data structures to arrays
  do iens=1,nens
    do iSeg=1,nSeg
      do iVar=1,nVarsIRFbas
        select case(iVar)
          case(ixIRFbas%qfuture); state%var(iVar)%array_3d_dp(iSeg,:,iens) = RCHFLX(iens,iSeg)%QFUTURE
          case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF state variable index'; return
        end select
      end do
    end do
  end do

  do iVar=1,nVarsIRFbas
   select case(iVar)
    case(ixIRFbas%qfuture); call write_nc(ncid, meta_irf_bas(iVar)%varName, state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,ntdh,nens/), ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF variable index for nc writing'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
  enddo

  end associate

  END SUBROUTINE write_IRFbas_state

  ! KWT writing procedures
  SUBROUTINE write_KWT_state(ierr, message1)

  USE globalData,   ONLY: meta_kwt
  USE var_lookup,   ONLY: ixKWT, nVarsKWT

  implicit none

  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  type(states)               :: state           ! temporal state data structures
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively
  integer(i4b), allocatable  :: RFvec(:)        ! temporal vector
  integer(i4b), allocatable  :: numWaves(:,:)   ! number of waves for each ensemble and segment
  ! initialize error control
  ierr=0; message1='write_KWT_state/'

  associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,  &
            nens     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            nwave    => meta_stateDims(ixStateDims%wave)%dimLength)     ! maximum waves allowed in a reach

  allocate(state%var(nVarsKWT), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! array to store number of wave per segment and ensemble
  allocate(numWaves(nens,nSeg), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT
    select case(iVar)
     case(ixKWT%routed); allocate(state%var(iVar)%array_3d_int(nSeg, nwave, nens), stat=ierr)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      allocate(state%var(iVar)%array_3d_dp(nSeg, nwave, nens), stat=ierr)
     case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWT routing state '//trim(meta_kwt(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nens
    do iSeg=1,nSeg
      numWaves(iens,iseg) = size(RCHSTA(iens, iseg)%LKW_ROUTE%KWAVE)
      do iVar=1,nVarsKWT
        select case(iVar)
          case(ixKWT%tentry)
            state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA(iens, iSeg)%LKW_ROUTE%KWAVE(:)%TI
            state%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
          case(ixKWT%texit)
            state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA(iens, iSeg)%LKW_ROUTE%KWAVE(:)%TR
            state%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
          case(ixKWT%qwave)
            state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA(iens, iSeg)%LKW_ROUTE%KWAVE(:)%QF
            state%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
          case(ixKWT%qwave_mod)
            state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA(iens, iSeg)%LKW_ROUTE%KWAVE(:)%QM
            state%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
          case(ixKWT%routed) ! this is suppposed to be logical variable, but put it as 0 or 1 in double now
            if (allocated(RFvec)) deallocate(RFvec, stat=ierr)
            allocate(RFvec(numWaves(iens,iSeg)),stat=ierr); RFvec=0_i4b
            where (RCHSTA(iens, iSeg)%LKW_ROUTE%KWAVE(:)%RF) RFvec=1_i4b
            state%var(iVar)%array_3d_int(iSeg,1:numWaves(iens,iSeg),iens) = RFvec
            state%var(iVar)%array_3d_int(iSeg,numWaves(iens,iSeg)+1:,iens) = integerMissing
          case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
        end select
      enddo ! variable loop
    enddo ! seg loop
  enddo ! ensemble loop

  ! Writing netCDF
  call write_nc(ncid, 'numWaves', numWaves, (/1,1/), (/nSeg,nens/), ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT
    select case(iVar)
     case(ixKWT%routed)
       call write_nc(ncid, trim(meta_kwt(iVar)%varName), state%var(iVar)%array_3d_int, (/1,1,1/), (/nSeg,nwave,nens/), ierr, cmessage)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      call write_nc(ncid, trim(meta_kwt(iVar)%varName), state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nwave,nens/), ierr, cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc writing'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
  end do

  end associate

  END SUBROUTINE write_KWT_state


  ! IRF writing procedures
  SUBROUTINE write_IRF_state(ierr, message1)

  USE globalData,   ONLY: meta_irf
  USE globalData,   ONLY: idxIRF
  USE var_lookup,   ONLY: ixIRF, nVarsIRF
  USE globalData,   ONLY: NETOPO          ! To get UH (this should not be in this data strucuture)

  implicit none

  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  type(states)               :: state           ! temporal state data structures
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively
  integer(i4b), allocatable  :: numQF(:,:)      ! number of future Q time steps for each ensemble and segment

  ierr=0; message1='write_IRF_state/'

  associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,  &
            nens     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            ntdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimLength)    ! maximum future q time steps among reaches

  allocate(state%var(nVarsIRF), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! array to store number of wave per segment and ensemble
  allocate(numQF(nens,nSeg), stat=ierr)
  if(ierr/=0)then; message1=trim(message1)//'problem allocating space for numQF'; return; endif

  do iVar=1,nVarsIRF
    select case(iVar)
      case(ixIRF%qfuture); allocate(state%var(iVar)%array_3d_dp(nSeg, ntdh_irf, nens), stat=ierr)
      case(ixIRF%irfVol);  allocate(state%var(iVar)%array_2d_dp(nSeg, nens), stat=ierr)
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
            state%var(iVar)%array_3d_dp(iSeg,1:numQF(iens,iSeg),iens) = RCHFLX(iens,iSeg)%QFUTURE_IRF
            state%var(iVar)%array_3d_dp(iSeg,numQF(iens,iSeg)+1:ntdh_irf,iens) = realMissing
          case(ixIRF%irfVol)
            state%var(iVar)%array_2d_dp(iSeg,iens) = RCHFLX(iens,iSeg)%ROUTE(idxIRF)%REACH_VOL(1)
          case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
        end select
      enddo ! variable loop
    enddo ! seg loop
  enddo ! ensemble loop

  ! writing netcdf
  call write_nc(ncid, 'numQF', numQF, (/1,1/), (/nSeg,nens/), ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRF
    select case(iVar)
      case(ixIRF%qfuture)
        call write_nc(ncid, trim(meta_irf(iVar)%varName), state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,ntdh_irf,nens/), ierr, cmessage)
      case(ixIRF%irfVol)
        call write_nc(ncid, trim(meta_irf(iVar)%varName), state%var(iVar)%array_2d_dp, (/1,1/), (/nSeg,nens/), ierr, cmessage)
      case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc writing'; return
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
    end select
  end do

  end associate

  END SUBROUTINE write_IRF_state

  ! KW writing procedures
  SUBROUTINE write_KW_state(ierr, message1)

    USE globalData,   ONLY: meta_kw
    USE var_lookup,   ONLY: ixKW, nVarsKW

    implicit none

    ! output
    integer(i4b), intent(out)  :: ierr            ! error code
    character(*), intent(out)  :: message1        ! error message
    ! local variables
    type(states)               :: state           ! temporal state data structures -currently 2 river routing scheme + basin IRF routing
    integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

    ierr=0; message1='write_KW_state/'

    associate(nSeg  => size(RCHFLX),                         &
              nEns  => meta_stateDims(ixStateDims%ens)%dimLength,  &
              nMesh => meta_stateDims(ixStateDims%mol_kw)%dimLength) ! number of computing molecule used for finite difference

    allocate(state%var(nVarsKW), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

    do iVar=1,nVarsKW
      select case(iVar)
       case(ixKW%qsub)
        allocate(state%var(iVar)%array_3d_dp(nSeg, nMesh, nEns), stat=ierr)
       case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KW routing state '//trim(meta_kw(iVar)%varName); return; endif
    end do

    ! --Convert data structures to arrays
    do iens=1,nEns
      do iSeg=1,nSeg
        do iVar=1,nVarsKW
          select case(iVar)
            case(ixKW%qsub)
              state%var(iVar)%array_3d_dp(iSeg,1:nMesh,iens) = RCHSTA(iens, iSeg)%KW_ROUTE%molecule%Q(1:nMesh)
            case default; ierr=20; message1=trim(message1)//'unable to identify KW routing state variable index'; return
          end select
        enddo ! variable loop
      enddo ! seg loop
    enddo ! ensemble loop

    ! Writing netCDF
    do iVar=1,nVarsKW
      select case(iVar)
       case(ixKW%qsub)
         call write_nc(ncid, trim(meta_kw(iVar)%varName), state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nMesh,nens/), ierr, cmessage)
       case default; ierr=20; message1=trim(message1)//'unable to identify KW variable index for nc writing'; return
      end select
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

    end do

    end associate

  END SUBROUTINE write_KW_state

  ! MC writing procedures
  SUBROUTINE write_MC_state(ierr, message1)

    USE globalData,   ONLY: meta_mc
    USE var_lookup,   ONLY: ixMC, nVarsMC

    implicit none

    ! output
    integer(i4b), intent(out)  :: ierr            ! error code
    character(*), intent(out)  :: message1        ! error message
    ! local variables
    type(states)               :: state           ! temporal state data structures
    integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

    ierr=0; message1='write_MC_state/'

    associate(nSeg  => size(RCHFLX),                         &
              nEns  => meta_stateDims(ixStateDims%ens)%dimLength,  &
              nMesh => meta_stateDims(ixStateDims%mol_mc)%dimLength) ! number of computing molecule used for finite difference

    allocate(state%var(nVarsMC), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

    do iVar=1,nVarsMC
      select case(iVar)
       case(ixMC%qsub)
        allocate(state%var(iVar)%array_3d_dp(nSeg, nMesh, nEns), stat=ierr)
       case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for MC routing state '//trim(meta_mc(iVar)%varName); return; endif
    end do

    ! --Convert data structures to arrays
    do iens=1,nEns
      do iSeg=1,nSeg
        do iVar=1,nVarsMC
          select case(iVar)
            case(ixMC%qsub)
              state%var(iVar)%array_3d_dp(iSeg,1:nMesh,iens) = RCHSTA(iens, iSeg)%MC_ROUTE%molecule%Q(1:nMesh)
            case default; ierr=20; message1=trim(message1)//'unable to identify MC routing state variable index'; return
          end select
        enddo ! variable loop
      enddo ! seg loop
    enddo ! ensemble loop

    ! Writing netCDF
    do iVar=1,nVarsMC
      select case(iVar)
       case(ixMC%qsub)
         call write_nc(ncid, trim(meta_mc(iVar)%varName), state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nMesh,nens/), ierr, cmessage)
       case default; ierr=20; message1=trim(message1)//'unable to identify MC variable index for nc writing'; return
      end select
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

    end do

    end associate

  END SUBROUTINE write_MC_state


  ! DW writing procedures
  SUBROUTINE write_DW_state(ierr, message1)

    USE globalData,   ONLY: meta_dw
    USE var_lookup,   ONLY: ixDW, nVarsDW

    implicit none

    ! output
    integer(i4b), intent(out)  :: ierr            ! error code
    character(*), intent(out)  :: message1        ! error message
    ! local variables
    type(states)               :: state           ! temporal state data structures
    integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

    ierr=0; message1='write_DW_state/'

    associate(nSeg  => size(RCHFLX),                         &
              nEns  => meta_stateDims(ixStateDims%ens)%dimLength,  &
              nMesh => meta_stateDims(ixStateDims%mol_dw)%dimLength) ! number of computing molecule used for finite difference

    allocate(state%var(nVarsDW), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

    do iVar=1,nVarsDW
      select case(iVar)
       case(ixDW%qsub)
        allocate(state%var(iVar)%array_3d_dp(nSeg, nMesh, nEns), stat=ierr)
       case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for DW routing state '//trim(meta_dw(iVar)%varName); return; endif
    end do

    ! --Convert data structures to arrays
    do iens=1,nEns
      do iSeg=1,nSeg
        do iVar=1,nVarsDW
          select case(iVar)
            case(ixDW%qsub)
              state%var(iVar)%array_3d_dp(iSeg,1:nMesh,iens) = RCHSTA(iens, iSeg)%DW_ROUTE%molecule%Q(1:nMesh)
            case default; ierr=20; message1=trim(message1)//'unable to identify DW routing state variable index'; return
          end select
        enddo ! variable loop
      enddo ! seg loop
    enddo ! ensemble loop

    ! Writing netCDF
    do iVar=1,nVarsDW
      select case(iVar)
       case(ixDW%qsub)
         call write_nc(ncid, trim(meta_dw(iVar)%varName), state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nMesh,nens/), ierr, cmessage)
       case default; ierr=20; message1=trim(message1)//'unable to identify DW variable index for nc writing'; return
      end select
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

    end do

    end associate

  END SUBROUTINE write_DW_state

 END SUBROUTINE write_state_nc

END MODULE write_restart
