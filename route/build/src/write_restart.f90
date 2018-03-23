MODULE write_restart

! Moudle wide external modules
USE nrtype, only: i4b, dp, &
                  strLen,  &
                  integerMissing
USE netcdf
USE public_var

implicit none

private

public::define_state_nc
public::write_state_nc

CONTAINS

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
 ierr = nf90_create(trim(fname),nf90_classic_model,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! For common dimension/variables - seg id, time, time-bound -----------
 ixDim_common = (/ixStateDims%seg, ixStateDims%ens, ixStateDims%time, ixStateDims%tbound/)

 ! Define dimensions
 ! time dimension
 ierr = nf90_def_dim(ncid, trim(dim_time), nf90_unlimited, meta_stateDims(ixStateDims%time)%dimId)

 do jDim = 1,size(ixDim_common)
  associate(ixDim_tmp => ixDim_common(jDim))
  if (meta_stateDims(ixDim_tmp)%dimLength == integerMissing) then
   call set_dim_len(ixDim_tmp, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//' for '//trim(meta_stateDims(ixDim_tmp)%dimName); return; endif
  endif
  ierr = nf90_def_dim(ncid, trim(meta_stateDims(ixDim_tmp)%dimName), meta_stateDims(ixDim_tmp)%dimLength ,meta_stateDims(jDim)%dimId)
  end associate
 end do

 ! Define variable
 call defvar(ncid, 'reachID', 'reach ID', '-', (/dim_seg/), nf90_int, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 call defvar(ncid, 'time ', 'time', trim(units_time), (/dim_time/), nf90_double, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 call defvar(ncid,'time_bound','time bound at last time step','sec', (/dim_tbound, dim_time/),nf90_double, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end associate

 ! Routing specific variables --------------

 ! basin IRF
 if (.true.) then
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
 ierr = nf90_enddef(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 CONTAINS

  SUBROUTINE set_dim_len(ixDim, ierr, message1)
   ! State/flux data structures
   USE globalData,   ONLY: meta_stateDims  ! states dimension meta
   USE globalData,   ONLY: NETOPO          ! To get segment size
   USE globalData,   ONLY: RCHFLX          ! To get size of q future for IRF
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
    case(ixStateDims%time);    meta_stateDims(ixStateDims%time)%dimLength    = 1
    case(ixStateDims%seg);     meta_stateDims(ixStateDims%seg)%dimLength     = size(NETOPO)
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
   integer(i4b), intent(out)  :: ierr          ! error code
   character(*), intent(out)  :: message1       ! error message
   ! local
   integer(i4b)               :: iVar          ! index loop for variables
   character(len=strLen)      :: dim_IRFbas(3) ! dimensions combination case 4
   character(len=strLen)      :: dim_q(2)      ! dimensions combination case 4

   ! initialize error control
   ierr=0; message1='define_IRFbas_state/'

   associate(dim_seg    => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens    => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_tdh    => meta_stateDims(ixStateDims%tdh)%dimName)

   ! Array dimension sets
   dim_IRFbas(:) = (/dim_seg, dim_tdh, dim_ens/)
   dim_q(:)      = (/dim_seg, dim_ens/)

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%tdh)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh)%dimName); return; endif
   end if

   ierr = nf90_def_dim(ncid, trim(meta_stateDims(ixStateDims%tdh)%dimName), meta_stateDims(ixStateDims%tdh)%dimLength ,meta_stateDims(ixStateDims%tdh)%dimId)
   if(ierr/=0)then; message1=trim(message1)//trim(nf90_strerror(ierr)); return; endif

   do iVar=1,nVarsIRFbas
    select case(iVar)
     case(ixIRFbas%qfuture); call defvar(ncid,trim(meta_irf_bas(iVar)%varName),trim(meta_irf_bas(iVar)%varDesc),trim(meta_irf_bas(iVar)%varUnit), dim_IRFbas, nf90_double, ierr, cmessage)
     case(ixIRFbas%q);       call defvar(ncid,trim(meta_irf_bas(iVar)%varName),trim(meta_irf_bas(iVar)%varDesc),trim(meta_irf_bas(iVar)%varUnit), dim_q, nf90_double, ierr, cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify hill-slope routing state variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//trim(nf90_strerror(ierr)); return; endif
   end do
  end associate

  END SUBROUTINE define_IRFbas_state

  SUBROUTINE define_KWT_state(ierr, message1)
   ! External modules
   USE globalData, ONLY: meta_kwt
   USE var_lookup, ONLY: ixKWT, nVarsKWT
   implicit none
   ! output
   integer(i4b), intent(out)  :: ierr        ! error code
   character(*), intent(out)  :: message1     ! error message
   ! local
   integer(i4b)               :: iVar        ! index loop for variables
   character(len=strLen)      :: dim_kwt(3)  ! dimensions combination case 4

   ! initialize error control
   ierr=0; message1='define_KWT_state/'

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_wave    => meta_stateDims(ixStateDims%wave)%dimName)

   dim_kwt(:) = (/dim_seg, dim_wave, dim_ens/)

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%wave)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%wave, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%wave)%dimName); return; endif
   end if

   ! Define dimension needed for this routing specific state variables
   ierr = nf90_def_dim(ncid, trim(meta_stateDims(ixStateDims%wave)%dimName), meta_stateDims(ixStateDims%wave)%dimLength ,meta_stateDims(ixStateDims%wave)%dimId)
   if(ierr/=0)then; message1=trim(message1)//trim(nf90_strerror(ierr)); return; endif

   call defvar(ncid,'numWaves', 'number of waves in a reach', '-', (/dim_seg, dim_ens/), nf90_int,ierr,cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(nf90_strerror(ierr)); return; endif

   do iVar=1,nVarsKWT

    if (iVar==ixKWT%q) cycle  !skip writing final discharge

    select case(iVar)
     case(ixKWT%tentry);    call defvar(ncid,trim(meta_kwt(iVar)%varName),trim(meta_kwt(iVar)%varDesc),trim(meta_kwt(iVar)%varUnit),dim_kwt,nf90_double,ierr,cmessage)
     case(ixKWT%texit);     call defvar(ncid,trim(meta_kwt(iVar)%varName),trim(meta_kwt(iVar)%varDesc),trim(meta_kwt(iVar)%varUnit),dim_kwt,nf90_double,ierr,cmessage)
     case(ixKWT%qwave);     call defvar(ncid,trim(meta_kwt(iVar)%varName),trim(meta_kwt(iVar)%varDesc),trim(meta_kwt(iVar)%varUnit),dim_kwt,nf90_double,ierr,cmessage)
     case(ixKWT%qwave_mod); call defvar(ncid,trim(meta_kwt(iVar)%varName),trim(meta_kwt(iVar)%varDesc),trim(meta_kwt(iVar)%varUnit),dim_kwt,nf90_double,ierr,cmessage)
     case(ixKWT%routed);    call defvar(ncid,trim(meta_kwt(iVar)%varName),trim(meta_kwt(iVar)%varDesc),trim(meta_kwt(iVar)%varUnit),dim_kwt,nf90_int,   ierr,cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//trim(nf90_strerror(ierr)); return; endif

   end do

   end associate

  END SUBROUTINE define_KWT_state

  SUBROUTINE define_IRF_state(ierr, message1)
   ! External modules
   USE globalData, ONLY: meta_irf
   USE var_lookup, ONLY: ixIRF, nVarsIRF
   implicit none
   ! output
   integer(i4b), intent(out)  :: ierr        ! error code
   character(*), intent(out)  :: message1     ! error message
   ! local
   integer(i4b)               :: iVar        ! index loop for variables
   character(len=strLen)      :: dim_irf(3)  ! dimensions combination case 4
   ! initialize error control
   ierr=0; message1='define_IRF_state/'

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimName,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimName,     &
             dim_tdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimName)

   dim_irf(:) = (/dim_seg, dim_tdh_irf, dim_ens/)

   if (meta_stateDims(ixStateDims%tdh_irf)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh_irf, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh_irf)%dimName); return; endif
   endif

   ! Define dimension needed for this routing specific state variables
   ierr = nf90_def_dim(ncid, trim(meta_stateDims(ixStateDims%tdh_irf)%dimName), meta_stateDims(ixStateDims%tdh_irf)%dimLength ,meta_stateDims(ixStateDims%tdh_irf)%dimId)
   if(ierr/=0)then; message1=trim(message1)//trim(nf90_strerror(ierr)); return; endif

   call defvar(ncid,'numQF','number of future q time steps in a reach','-',(/dim_seg, dim_ens/), nf90_int, ierr,cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(nf90_strerror(ierr)); return; endif

   do iVar=1,nVarsIRF

    if (iVar==ixIRF%q) cycle  ! skip writing final discharge

    select case(iVar)
     case(ixIRF%qfuture); call defvar(ncid,trim(meta_irf(iVar)%varName),trim(meta_irf(iVar)%varDesc),trim(meta_irf(iVar)%varUnit), dim_irf, nf90_double,ierr,cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify IRF routing state variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//trim(nf90_strerror(ierr)); return; endif

   end do

   end associate

  END SUBROUTINE define_IRF_state

 END SUBROUTINE define_state_nc


 ! *********************************************************************
 ! public subroutine: writing routing state NetCDF file
 ! *********************************************************************
 SUBROUTINE write_state_nc(&
                           fname,                &   ! Input: state netcdf name
                           opt,                  &   ! input: which routing options
                           time, iTime, T0, T1,  &   ! Input: time, time step, start and end time [sec]
                           seg_id,               &   ! Input: segment id vector
                           ierr, message)            ! Output: error control
 ! External module
 USE write_netcdf, ONLY: write_nc
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
 type(states)                    :: state(0:2)      ! temporal state data structures -currently 2 river routing scheme + basin IRF routing
 character(len=strLen)           :: cmessage        ! error message of downwind routine

 ! initialize error control
 ierr=0; message='write_state_nc/'

 ! -- Write out to netCDF
 ! Miscellaneous variables - seg id, time etc
 call write_nc(fname,'reachID', seg_id, (/1/), (/size(seg_id)/), ierr, cmessage);
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_nc(fname,'time', (/time/), (/iTime/), (/1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_nc(fname,'time_bound', (/T0,T1/), (/1,iTime/), (/2,1/), ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 if (.true.)then
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
    case(ixIRFbas%q);       call write_nc(fname, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_2d_dp, (/1,1/), (/nSeg,nens/), ierr, cmessage)
    case(ixIRFbas%qfuture); call write_nc(fname, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,ntdh,nens/), ierr, cmessage)
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
    if (iVar==ixKWT%q) cycle  ! not writing out river flow in state file
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

     if (iVar==ixKWT%q) cycle ! not writing out KWT routed flow

     select case(iVar)
      case(ixKWT%tentry);    state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = KROUTE(iens,iSeg)%KWAVE(:)%TI
      case(ixKWT%texit);     state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = KROUTE(iens,iSeg)%KWAVE(:)%TR
      case(ixKWT%qwave);     state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = KROUTE(iens,iSeg)%KWAVE(:)%QF
      case(ixKWT%qwave_mod); state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = KROUTE(iens,iSeg)%KWAVE(:)%QM
      case(ixKWT%routed) ! this is suppposed to be logical variable, but put it as 0 or 1 in double now
       if (allocated(RFvec)) deallocate(RFvec, stat=ierr)
       allocate(RFvec(numWaves(iens,iSeg)),stat=ierr); RFvec=0_i4b
       where (KROUTE(iens,iSeg)%KWAVE(:)%RF) RFvec=1_i4b
       state(kinematicWave)%var(iVar)%array_3d_int(iSeg,1:numWaves(iens,iSeg),iens) = RFvec
      case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
     end select

    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! Writing netCDF
  call write_nc(fname, 'numWaves', numWaves, (/1,1/), (/nSeg,nens/), ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT

    if (iVar==ixKWT%q) cycle  ! not writing out river flow in state file

    select case(iVar)
     case(ixKWT%routed)
       call write_nc(fname, trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_int, (/1,1,1/), (/nSeg,nwave,nens/), ierr, cmessage)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      call write_nc(fname, trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nwave,nens/), ierr, cmessage)
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
   if (iVar==ixIRF%q) cycle
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

     if (iVar==ixIRF%q) cycle ! not writing out IRF routed flow

     select case(iVar)
      case(ixIRF%qfuture); state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,1:numQF(iens,iSeg),iens) = RCHFLX(iens,iSeg)%QFUTURE_IRF
      case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
     end select

    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! writing netcdf
  call write_nc(fname, 'numQF', numQF, (/1,1/), (/nSeg,nens/), ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRF

   call write_nc(fname, trim(meta_irf(iVar)%varName), state(impulseResponseFunc)%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,ntdh_irf,nens/), ierr, cmessage)
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  end do

  end associate

  END SUBROUTINE write_IRF_state

 END SUBROUTINE write_state_nc


 ! *********************************************************************
 ! private subroutine: define variable attributes NetCDF file
 ! *********************************************************************
 SUBROUTINE defvar(ncid, vname, vdesc, vunit, dimNames, ivtype, ierr, message)
  implicit none
  ! input
  integer(i4b), intent(in)   :: ncid        ! Input: netcdf fine ID
  character(*), intent(in)   :: vname       ! Input: variable name
  character(*), intent(in)   :: vdesc       ! Input: variable description
  character(*), intent(in)   :: vunit       ! Input: variable units
  character(*), intent(in)   :: dimNames(:) ! Input: variable dimension names
  integer(i4b), intent(in)   :: ivtype      ! Input: variable type
  ! output
  integer(i4b), intent(out)  :: ierr        ! error code
  character(*), intent(out)  :: message     ! error message
  ! local
  integer(i4b)               :: id          ! loop through dimensions
  integer(i4b)               :: dimIDs(size(dimNames))  ! vector of dimension IDs
  integer(i4b)               :: iVarId      ! variable ID

  ! define dimension IDs
  do id=1,size(dimNames)
   ierr=nf90_inq_dimid(ncid,trim(dimNames(id)),dimIDs(id))
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end do

  ! define variable
  ierr = nf90_def_var(ncid,trim(vname),ivtype,dimIds,iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! add variable description
  ierr = nf90_put_att(ncid,iVarId,'long_name',trim(vdesc))
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! add variable units
  ierr = nf90_put_att(ncid,iVarId,'units',trim(vunit))
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 END SUBROUTINE defvar


END MODULE write_restart
