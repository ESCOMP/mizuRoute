MODULE read_restart

! Moudle wide external modules
USE nrtype, ONLY: i4b, dp, &
                  strLen
USE public_var

implicit none

private

public::read_state_nc

CONTAINS

 ! *********************************************************************
 ! public subroutine: read routing state NetCDF file
 ! *********************************************************************
 SUBROUTINE read_state_nc(&
                          fname,           &   ! Input:  state netcdf name
                          opt,             &   ! input:  which routing options
                          T0, T1,          &   ! output: start and end time [sec]
                          ierr, message)       ! Output: error control
 ! External module
 USE io_netcdf, ONLY: get_nc, &
                      get_nc_dim_len
 USE dataTypes, ONLY: states
 ! meta data
 USE var_lookup, ONLY: ixStateDims, nStateDims
 USE globalData, ONLY: meta_stateDims            ! dimension for state variables
 USE globalData, ONLY: RCHFLX                    ! reach flux data structure for the entire domain
 USE globalData, ONLY: RCHSTA                    ! reach state data structure for the entire domain
 USE globalData, ONLY: ixRch_order

 implicit none
 ! input variables
 character(*), intent(in)      :: fname                ! filename
 integer(i4b), intent(in)      :: opt                  ! routing option 0=all, 1=kwt, 2=irf
 ! output variables
 real(dp),     intent(out)     :: T0                   ! beginning time [sec] of ith time step - lapse time from the beginning of the simulation
 real(dp),     intent(out)     :: T1                   ! ending time [sec] ith time step - lapse time from the beginning of the simulation
 integer(i4b), intent(out)     :: ierr                 ! error code
 character(*), intent(out)     :: message              ! error message
 ! local variables
 real(dp)                      :: TB(2)                ! 2 element-time bound vector
 type(states)                  :: state(0:2)           ! temporal state data structures -currently 2 river routing scheme + basin IRF routing
 integer(i4b)                  :: nSeg,nens            ! dimenion sizes
 integer(i4b)                  :: nTime,ntbound        ! dimenion sizes
 integer(i4b)                  :: ixDim_common(4)      ! custom dimension ID array
 integer(i4b)                  :: jDim                 ! index loops for dimension
 character(len=strLen)         :: cmessage             ! error message of downwind routine

 ! initialize error control
 ierr=0; message='read_state_nc/'

 ! get Dimension sizes
 ! For common dimension/variables - seg id, time, time-bound -----------
 ixDim_common = (/ixStateDims%seg, ixStateDims%ens, ixStateDims%time, ixStateDims%tbound/)

 do jDim=1,size(ixDim_common)
   associate (ixDim_tmp => ixDim_common(jDim))
   select case(ixDim_tmp)
    case(ixStateDims%seg);     call get_nc_dim_len(fname, trim(meta_stateDims(ixDim_tmp)%dimName), nSeg,    ierr, cmessage)
    case(ixStateDims%ens);     call get_nc_dim_len(fname, trim(meta_stateDims(ixDim_tmp)%dimName), nens,    ierr, cmessage)
    case(ixStateDims%time);    call get_nc_dim_len(fname, trim(meta_stateDims(ixDim_tmp)%dimName), nTime,   ierr, cmessage)
    case(ixStateDims%tbound);  call get_nc_dim_len(fname, trim(meta_stateDims(ixDim_tmp)%dimName), ntbound, ierr, cmessage)
    case default; ierr=20; message=trim(message)//'unable to identify dimension name index'; return
   end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end associate
 enddo

 allocate(RCHFLX(nens,nSeg), RCHSTA(nens,nSeg), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating [RCHFLX, RCHSTA]'; return; endif

 ! Read variables
 ! time bound
 call get_nc(fname,'time_bound',TB(:), 1, 2, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 T0=TB(1); T1=TB(2)

 ! routing specific variables
 call read_IRFbas_state(ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return;endif

 if (opt==allRoutingMethods .or. opt==kinematicWave) then
  call read_KWT_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif
 endif

 if (opt==allRoutingMethods .or. opt==impulseResponseFunc) then
  call read_IRF_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif
 end if

 CONTAINS

  SUBROUTINE read_IRFbas_state(ierr, message1)
  ! meta data
  USE globalData, ONLY: meta_irf_bas              ! basin IRF routing
  ! Named variables
  USE var_lookup, ONLY: ixIRFbas, nVarsIRFbas
  implicit none
  ! output
  integer(i4b), intent(out)     :: ierr           ! error code
  character(*), intent(out)     :: message1       ! error message
  ! local variables
  integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
  integer(i4b)                  :: jSeg           ! index loops for reaches respectively
  integer(i4b)                  :: ntdh           ! dimension size

  ! initialize error control
  ierr=0; message1='read_IRFbas_state/'

  call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%tdh)%dimName), ntdh, ierr, cmessage)
  if(ierr/=0)then;  message1=trim(message1)//trim(cmessage); return; endif

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

  do iVar=1,nVarsIRFbas

   select case(iVar)
    case(ixIRFbas%q);       call get_nc(fname, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_2d_dp, (/1,1/), (/nSeg,nens/), ierr, cmessage)
    case(ixIRFbas%qfuture); call get_nc(fname, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,ntdh,nens/), ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF variable index for nc writing'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  enddo

  do iens=1,nens
   do iSeg=1,nSeg

    jSeg = ixRch_order(iSeg)

    allocate(RCHFLX(iens,jSeg)%QFUTURE(ntdh), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

    do iVar=1,nVarsIRFbas

     select case(iVar)
      case(ixIRFbas%q);       RCHFLX(iens,jSeg)%BASIN_QR(1) = state(0)%var(iVar)%array_2d_dp(iSeg,iens)
      case(ixIRFbas%qfuture); RCHFLX(iens,jSeg)%QFUTURE(:)  = state(0)%var(iVar)%array_3d_dp(iSeg,:,iens)
      case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF state variable index'; return
     end select

    enddo
   enddo
  enddo

  END SUBROUTINE read_IRFbas_state


  SUBROUTINE read_IRF_state(ierr, message1)
  ! meta data
  USE globalData,  ONLY: meta_irf               ! IRF routing
  ! Named variables
  USE var_lookup,  ONLY: ixIRF, nVarsIRF
  implicit none
  integer(i4b), intent(out)     :: ierr           ! error code
  character(*), intent(out)     :: message1       ! error message
  ! local variables
  integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
  integer(i4b)                  :: jSeg           ! index loops for reaches respectively
  integer(i4b), allocatable     :: numQF(:,:)     ! number of future Q time steps for each ensemble and segment
  integer(i4b)                  :: ntdh_irf       ! dimenion sizes
  ! initialize error control
  ierr=0; message1='read_IRF_state/'

  call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%tdh_irf)%dimName), ntdh_irf, ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  allocate(state(impulseResponseFunc)%var(nVarsIRF), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  allocate(numQF(nens,nSeg), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRF

   select case(iVar)
    case(ixIRF%qfuture); allocate(state(impulseResponseFunc)%var(iVar)%array_3d_dp(nSeg, ntdh_irf, nens), stat=ierr)
    case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//'problem allocating space for IRF routing state '//trim(meta_irf(iVar)%varName); return; endif

  end do

  call get_nc(fname,'numQF',numQF,(/1,1/),(/nSeg,nens/),ierr,cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRF

   select case(iVar)
    case(ixIRF%qfuture); call get_nc(fname, meta_irf(iVar)%varName, state(impulseResponseFunc)%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,ntdh_irf,nens/), ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc reading'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  end do

  do iens=1,nens
   do iSeg=1,nSeg

    jSeg = ixRch_order(iSeg)

    allocate(RCHFLX(iens,jSeg)%QFUTURE_IRF(numQF(iens,iSeg)), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

    do iVar=1,nVarsIRF

     select case(iVar)
      case(ixIRF%qfuture); RCHFLX(iens,jSeg)%QFUTURE_IRF = state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,1:numQF(iens,iSeg),iens)
      case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
     end select

    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  END SUBROUTINE read_IRF_state


  SUBROUTINE read_KWT_state(ierr, message1)
  ! meta data
  USE globalData, ONLY: meta_kwt                  ! kwt routing
  ! Named variables
  USE var_lookup, ONLY: ixKWT, nVarsKWT
  implicit none
  integer(i4b), intent(out)     :: ierr           ! error code
  character(*), intent(out)     :: message1       ! error message
  ! output
  integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
  integer(i4b)                  :: jSeg           ! index loops for reaches respectively
  integer(i4b)                  :: nwave          ! dimenion sizes
  integer(i4b), allocatable     :: RFvec(:)       ! temporal vector
  integer(i4b), allocatable     :: numWaves(:,:)  ! number of waves for each ensemble and segment
  ! initialize error control
  ierr=0; message1='read_KWT_state/'

  allocate(state(kinematicWave)%var(nVarsKWT), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  allocate(numWaves(nens,nSeg), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! get Dimension sizes
  call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%wave)%dimName), nwave, ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT

    select case(iVar)
     case(ixKWT%routed); allocate(state(kinematicWave)%var(iVar)%array_3d_dp(nSeg, nwave, nens), stat=ierr)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      allocate(state(kinematicWave)%var(iVar)%array_3d_dp(nSeg, nwave, nens), stat=ierr)
     case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWT routing state '//trim(meta_kwt(iVar)%varName); return; endif
  end do

  call get_nc(fname,'numWaves',numWaves, (/1,1/), (/nSeg,nens/), ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT

    select case(iVar)
     case(ixKWT%routed)
      call get_nc(fname,trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nwave,nens/), ierr, cmessage)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      call get_nc(fname,trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nwave,nens/), ierr, cmessage)
     case default; ierr=20; message1=trim(message)//'unable to identify KWT variable index for nc reading'; return
    end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
  end do

  do iens=1,nens
   do iSeg=1,nSeg

    jSeg = ixRch_order(iSeg)

    allocate(RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1), stat=ierr)

    do iVar=1,nVarsKWT

     select case(iVar)
      case(ixKWT%tentry);    RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%TI = state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
      case(ixKWT%texit);     RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%TR = state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
      case(ixKWT%qwave);     RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%QF = state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
      case(ixKWT%qwave_mod); RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%QM = state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
      case(ixKWT%routed) ! this is suppposed to be logical variable, but put it as 0 or 1 in double now
       if (allocated(RFvec)) deallocate(RFvec, stat=ierr)
       allocate(RFvec(0:numWaves(iens,iSeg)-1),stat=ierr)
       RFvec = nint(state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens))
       RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%RF=.False.
       where (RFvec==1_i4b) RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%RF=.True.
      case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
     end select

    enddo
   enddo
  enddo

  END SUBROUTINE read_KWT_state

 END SUBROUTINE read_state_nc


END MODULE read_restart
