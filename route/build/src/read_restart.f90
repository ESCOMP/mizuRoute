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
 SUBROUTINE read_state_nc(fname,           &   ! input:  state netcdf name
                          T0, T1,          &   ! output: start and end time [sec]
                          ierr, message)       ! output: error control
 ! External module
 USE ncio_utils, ONLY: get_nc, &
                       get_nc_dim_len
 USE ncio_utils, ONLY: open_nc, close_nc
 USE ncio_utils, ONLY: check_variable
 USE dataTypes,  ONLY: states
 ! meta data
 USE var_lookup, ONLY: ixStateDims, nStateDims
 USE globalData, ONLY: meta_stateDims            ! dimension for state variables
 USE globalData, ONLY: RCHFLX                    ! reach flux data structure for the entire domain
 USE globalData, ONLY: RCHSTA                    ! reach state data structure for the entire domain
 USE globalData, ONLY: ixRch_order
 USE globalData, ONLY: nNodes                    ! number of MPI tasks
 USE globalData, ONLY: onRoute                   ! logical to indicate which routing method(s) is on
 USE globalData, ONLY: nRoutes                   ! number of active routing methods
 USE public_var, ONLY: impulseResponseFunc
 USE public_var, ONLY: kinematicWaveTracking
 USE public_var, ONLY: kinematicWave
 USE public_var, ONLY: muskingumCunge
 USE public_var, ONLY: diffusiveWave

 implicit none
 ! argument variables
 character(*), intent(in)      :: fname                ! filename
 real(dp),     intent(out)     :: T0                   ! beginning time [sec] of ith time step - lapse time from the beginning of the simulation
 real(dp),     intent(out)     :: T1                   ! ending time [sec] ith time step - lapse time from the beginning of the simulation
 integer(i4b), intent(out)     :: ierr                 ! error code
 character(*), intent(out)     :: message              ! error message
 ! local variables
 real(dp)                      :: TB(2)                ! 2 element-time bound vector
 integer(i4b)                  :: nSeg,nens            ! dimenion sizes
 integer(i4b)                  :: ntbound              ! dimenion sizes
 integer(i4b)                  :: ixDim_common(3)      ! custom dimension ID array
 integer(i4b)                  :: jDim                 ! index loops for dimension
 integer(i4b)                  :: iSeg                 ! index loops for reach
 integer(i4b)                  :: nNodes_in            ! number of MPI tasks for restart file
 character(len=strLen)         :: cmessage             ! error message of downwind routine

 ierr=0; message='read_state_nc/'

 ! get Dimension sizes
 ! For common dimension/variables - seg id, time-bound -----------
 ixDim_common = (/ixStateDims%seg, ixStateDims%ens, ixStateDims%tbound/)

 do jDim=1,size(ixDim_common)
   associate (ixDim_tmp => ixDim_common(jDim))
   select case(ixDim_tmp)
    case(ixStateDims%seg);     call get_nc_dim_len(fname, trim(meta_stateDims(ixDim_tmp)%dimName), nSeg,    ierr, cmessage)
    case(ixStateDims%ens);     call get_nc_dim_len(fname, trim(meta_stateDims(ixDim_tmp)%dimName), nens,    ierr, cmessage)
    case(ixStateDims%tbound);  call get_nc_dim_len(fname, trim(meta_stateDims(ixDim_tmp)%dimName), ntbound, ierr, cmessage)
    case default; ierr=20; message=trim(message)//'unable to identify dimension name index'; return
   end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end associate
 enddo

 allocate(RCHFLX(nens,nSeg), RCHSTA(nens,nSeg), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating [RCHFLX, RCHSTA]'; return; endif

 do iSeg=1,nSeg
   allocate(RCHFLX(nens,iSeg)%ROUTE(nRoutes))
 end do

 ! Read variables
 call get_nc(fname,'nNodes',nNodes_in, 1, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 if ( nNodes /= nNodes_in )then
    write(cmessage,'(i4)') nNodes_in
    ierr = 100
    message=trim(message)//'Number of MPI tasks on restart file is different from number using now' // &
            ', they must be the same, nNodes_in='//trim(cmessage)
    return
 end if

 ! time bound
 call get_nc(fname,'time_bound',TB(:), 1, 2, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 T0=TB(1); T1=TB(2)

 call read_basinQ_state(ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return;endif

 ! routing specific variables
 if (doesBasinRoute == 1) then
   call read_IRFbas_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return;endif
 endif

 if (onRoute(impulseResponseFunc)) then
   call read_IRF_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif
 end if

 if (onRoute(kinematicWaveTracking)) then
   call read_KWT_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif
 end if

 if (onRoute(kinematicWave)) then
   call read_KW_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif
 end if

 if (onRoute(muskingumCunge)) then
   call read_MC_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif
 end if

 if (onRoute(diffusiveWave)) then
   call read_DW_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif
 end if

 CONTAINS

  SUBROUTINE read_basinQ_state(ierr, message1)

    USE globalData, ONLY: meta_basinQ              ! reach inflow from basin at previous time step
    USE globalData, ONLY: RCHFLX                    ! To get q future for basin IRF and IRF (these should not be in this data strucuture)
    USE var_lookup, ONLY: ixBasinQ, nVarsBasinQ
    implicit none
    ! output
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! loop indices for variables, ensembles, reaches respectively
    integer(i4b)                  :: jSeg           ! sorted index for reaches

    ! initialize error control
    ierr=0; message1='read_basinQ_state/'

    allocate(state%var(nVarsBasinQ), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsBasinQ
      select case(iVar)
        case(ixBasinQ%q); allocate(state%var(iVar)%array_2d_dp(nSeg, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify basin routing variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for reach inflow:'//trim(meta_basinQ(iVar)%varName); return; endif
    end do

    do iVar=1,nVarsBasinQ
      select case(iVar)
        case(ixBasinQ%q); call get_nc(fname, meta_basinQ(iVar)%varName, state%var(iVar)%array_2d_dp, [1,1], [nSeg,nens], ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify previous time step reach inflow variable index for nc writing'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_basinQ(iVar)%varName); return; endif
    enddo

    do iens=1,nens
      do iSeg=1,nSeg
        jSeg = ixRch_order(iSeg)
        do iVar=1,nVarsBasinQ
          select case(iVar)
            case(ixBasinQ%q); RCHFLX(iens,jSeg)%BASIN_QR(1) = state%var(iVar)%array_2d_dp(iSeg,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify previous time step reach inflow variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_basinQ_state

  SUBROUTINE read_IRFbas_state(ierr, message1)

    USE globalData, ONLY: meta_irf_bas              ! basin IRF routing
    USE var_lookup, ONLY: ixIRFbas, nVarsIRFbas
    implicit none
    ! output
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
    integer(i4b)                  :: jSeg           ! index loops for reaches respectively
    integer(i4b)                  :: ntdh           ! dimension size

    ierr=0; message1='read_IRFbas_state/'

    call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%tdh)%dimName), ntdh, ierr, cmessage1)
    if(ierr/=0)then;  message1=trim(message1)//trim(cmessage1); return; endif

    allocate(state%var(nVarsIRFbas), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsIRFbas
      select case(iVar)
        case(ixIRFbas%qfuture); allocate(state%var(iVar)%array_3d_dp(nSeg, ntdh, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify basin routing variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for basin IRF routing state:'//trim(meta_irf_bas(iVar)%varName); return; endif
    end do

    do iVar=1,nVarsIRFbas
      select case(iVar)
        case(ixIRFbas%qfuture); call get_nc(fname, meta_irf_bas(iVar)%varName, state%var(iVar)%array_3d_dp, [1,1,1], [nSeg,ntdh,nens], ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF variable index for nc writing'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_irf_bas(iVar)%varName); return; endif
    enddo

    do iens=1,nens
      do iSeg=1,nSeg
        jSeg = ixRch_order(iSeg)
        allocate(RCHFLX(iens,jSeg)%QFUTURE(ntdh), stat=ierr, errmsg=cmessage1)
        if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif
        do iVar=1,nVarsIRFbas
          select case(iVar)
            case(ixIRFbas%qfuture); RCHFLX(iens,jSeg)%QFUTURE(:)  = state%var(iVar)%array_3d_dp(iSeg,:,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF state variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_IRFbas_state


  SUBROUTINE read_IRF_state(ierr, message1)

    USE globalData,  ONLY: meta_irf               ! IRF routing
    USE globalData,  ONLY: idxIRF
    USE var_lookup,  ONLY: ixIRF, nVarsIRF
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: ncidRestart    ! restart netcdf id
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
    integer(i4b)                  :: jSeg           ! index loops for reaches respectively
    integer(i4b), allocatable     :: numQF(:,:)     ! number of future Q time steps for each ensemble and segment
    integer(i4b)                  :: ntdh_irf       ! dimenion sizes
    integer(i4b)                  :: nTbound=2      ! dimenion sizes

    ierr=0; message1='read_IRF_state/'

    call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%tdh_irf)%dimName), ntdh_irf, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    allocate(state%var(nVarsIRF), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    allocate(numQF(nens,nSeg), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsIRF
      select case(iVar)
        case(ixIRF%qfuture); allocate(state%var(iVar)%array_3d_dp(nSeg, ntdh_irf, nens), stat=ierr)
        case(ixIRF%vol : ixIRF%qerror); allocate(state%var(iVar)%array_3d_dp(nSeg, nTbound, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for IRF routing state:'//trim(meta_irf(iVar)%varName); return; endif
    end do

    call get_nc(fname,'numQF',numQF,(/1,1/),(/nSeg,nens/),ierr,cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':numQF'; return; endif

    do iVar=1,nVarsIRF
      select case(iVar)
        case(ixIRF%qfuture); call get_nc(fname, meta_irf(iVar)%varName, state%var(iVar)%array_3d_dp, [1,1,1], [nSeg,ntdh_irf,nens], ierr, cmessage1)
        case(ixIRF%vol);     call get_nc(fname, meta_irf(iVar)%varName, state%var(iVar)%array_3d_dp, [1,1,1], [nSeg,nTbound, nens], ierr, cmessage1)
        case(ixIRF%qerror)
          call open_nc(fname, 'r', ncidRestart, ierr, cmessage1)
          if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif
          if (check_variable(ncidRestart,meta_irf(iVar)%varName)) then
            call get_nc(fname, meta_irf(iVar)%varName, state%var(iVar)%array_2d_dp, (/1,1/), (/nSeg, nens/), ierr, cmessage1)
          else
            state%var(iVar)%array_2d_dp = 0._dp
          end if
          call close_nc(ncidRestart, ierr, cmessage1)
          if(ierr/=0)then; message=trim(message1)//trim(cmessage1); return; endif
        case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc reading'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_irf(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        jSeg = ixRch_order(iSeg)
        allocate(RCHFLX(iens,jSeg)%QFUTURE_IRF(numQF(iens,iSeg)), stat=ierr, errmsg=cmessage1)
        if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif
        do iVar=1,nVarsIRF
          select case(iVar)
            case(ixIRF%qfuture); RCHFLX(iens,jSeg)%QFUTURE_IRF    = state%var(iVar)%array_3d_dp(iSeg,1:numQF(iens,iSeg),iens)
            case(ixIRF%vol);     RCHFLX(iens,jSeg)%ROUTE(idxIRF)%REACH_VOL(0:1) = state%var(iVar)%array_3d_dp(iSeg,1:2,iens)
            case(ixIRF%qerror);  RCHFLX(iens,iSeg)%ROUTE(idxIRF)%Qerror = state%var(iVar)%array_2d_dp(iSeg,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
          end select
        enddo ! variable loop
      enddo ! seg loop
    enddo ! ensemble loop

  END SUBROUTINE read_IRF_state

  SUBROUTINE read_KWT_state(ierr, message1)

    USE globalData, ONLY: meta_kwt                  ! kwt routing
    USE globalData, ONLY: idxKWT
    USE var_lookup, ONLY: ixKWT, nVarsKWT
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
    integer(i4b)                  :: jSeg           ! index loops for reaches respectively
    integer(i4b)                  :: nwave          ! dimenion sizes
    integer(i4b), allocatable     :: RFvec(:)       ! temporal vector
    integer(i4b), allocatable     :: numWaves(:,:)  ! number of waves for each ensemble and segment

    ierr=0; message1='read_KWT_state/'

    allocate(state%var(nVarsKWT), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    allocate(numWaves(nens,nSeg), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    ! get Dimension sizes
    call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%wave)%dimName), nwave, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsKWT
      select case(iVar)
        case(ixKWT%routed); allocate(state%var(iVar)%array_3d_dp(nSeg, nwave, nens), stat=ierr)
        case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
          allocate(state%var(iVar)%array_3d_dp(nSeg, nwave, nens), stat=ierr)
        case(ixKWT%vol);  allocate(state%var(iVar)%array_2d_dp(nSeg, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWT routing state:'//trim(meta_kwt(iVar)%varName); return; endif
    end do

    call get_nc(fname,'numWaves',numWaves, (/1,1/), (/nSeg,nens/), ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//'numWaves'; return; endif

    do iVar=1,nVarsKWT
      select case(iVar)
        case(ixKWT%routed)
          call get_nc(fname,trim(meta_kwt(iVar)%varName), state%var(iVar)%array_3d_dp, [1,1,1], [nSeg,nwave,nens], ierr, cmessage1)
        case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
          call get_nc(fname,trim(meta_kwt(iVar)%varName), state%var(iVar)%array_3d_dp, [1,1,1], [nSeg,nwave,nens], ierr, cmessage1)
        case(ixKWT%vol)
          call get_nc(fname, meta_kwt(iVar)%varName, state%var(iVar)%array_2d_dp, [1,1], [nSeg, nens], ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify KWT variable index for nc reading'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_kwt(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        jSeg = ixRch_order(iSeg)
        allocate(RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1), stat=ierr)
        do iVar=1,nVarsKWT
          select case(iVar)
            case(ixKWT%tentry);    RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%TI = state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
            case(ixKWT%texit);     RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%TR = state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
            case(ixKWT%qwave);     RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%QF = state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
            case(ixKWT%qwave_mod); RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%QM = state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
            case(ixKWT%vol);       RCHFLX(iens,jSeg)%ROUTE(idxKWT)%REACH_VOL(1) = state%var(iVar)%array_2d_dp(iSeg,iens)
            case(ixKWT%routed) ! this is suppposed to be logical variable, but put it as 0 or 1 in double now
              if (allocated(RFvec)) deallocate(RFvec, stat=ierr)
              allocate(RFvec(0:numWaves(iens,iSeg)-1),stat=ierr)
              RFvec = nint(state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens))
              RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%RF=.False.
              where (RFvec==1_i4b) RCHSTA(iens,jSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%RF=.True.
            case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_KWT_state

  SUBROUTINE read_KW_state(ierr, message1)

    USE globalData, ONLY: meta_kw
    USE globalData, ONLY: idxKW
    USE globalData, ONLY: nMolecule
    USE var_lookup, ONLY: ixKW, nVarsKW
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    integer(i4b)                  :: ncidRestart    ! restart netcdf id
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
    integer(i4b)                  :: jSeg           ! index loops for reaches respectively

    ierr=0; message1='read_KW_state/'

    allocate(state%var(nVarsKW), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    ! get Dimension sizes
    call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%mol_kw)%dimName), nMolecule%KW_ROUTE, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsKW
      select case(iVar)
        case(ixKW%qsub); allocate(state%var(iVar)%array_3d_dp(nSeg, nMolecule%KW_ROUTE, nens), stat=ierr)
        case(ixKW%vol : ixKW%qerror);  allocate(state%var(iVar)%array_2d_dp(nSeg, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KW routing state:'//trim(meta_kw(iVar)%varName); return; endif
    end do

    do iVar=1,nVarsKW
      select case(iVar)
        case(ixKW%qsub); call get_nc(fname, meta_kw(iVar)%varName, state%var(iVar)%array_3d_dp, [1,1,1], [nSeg,nMolecule%KW_ROUTE,nens], ierr, cmessage1)
        case(ixKW%vol);  call get_nc(fname, meta_kw(iVar)%varName, state%var(iVar)%array_2d_dp, [1,1], [nSeg, nens], ierr, cmessage1)
        case(ixKW%qerror)
          call open_nc(fname, 'r', ncidRestart, ierr, cmessage1)
          if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif
          if (check_variable(ncidRestart,meta_kw(iVar)%varName)) then
            call get_nc(fname, meta_kw(iVar)%varName, state%var(iVar)%array_2d_dp, (/1,1/), (/nSeg, nens/), ierr, cmessage1)
          else
            state%var(iVar)%array_2d_dp = 0._dp
          end if
          call close_nc(ncidRestart, ierr, cmessage1)
          if(ierr/=0)then; message=trim(message1)//trim(cmessage1); return; endif
        case default; ierr=20; message1=trim(message1)//'unable to identify KW variable index for nc reading'; return
      end select
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_kw(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        jSeg = ixRch_order(iSeg)
        allocate(RCHSTA(iens,jSeg)%KW_ROUTE%molecule%Q(nMolecule%KW_ROUTE), stat=ierr, errmsg=cmessage1)
        do iVar=1,nVarsKW
          select case(iVar)
            case(ixKW%qsub); RCHSTA(iens,jSeg)%KW_ROUTE%molecule%Q(1:nMolecule%KW_ROUTE) = state%var(iVar)%array_3d_dp(iSeg,1:nMolecule%KW_ROUTE,iens)
            case(ixKW%vol);  RCHFLX(iens,jSeg)%ROUTE(idxKW)%REACH_VOL(1) = state%var(iVar)%array_2d_dp(iSeg,iens)
            case(ixKW%qerror); RCHFLX(iens,iSeg)%ROUTE(idxKW)%Qerror = state%var(iVar)%array_2d_dp(iSeg,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify KW routing state variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_KW_state

  SUBROUTINE read_MC_state(ierr, message1)

    USE globalData, ONLY: meta_mc
    USE globalData, ONLY: idxMC
    USE globalData, ONLY: nMolecule
    USE var_lookup, ONLY: ixMC, nVarsMC
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    integer(i4b)                  :: ncidRestart    ! restart netcdf id
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
    integer(i4b)                  :: jSeg           ! index loops for reaches respectively

    ierr=0; message1='read_MC_state/'

    allocate(state%var(nVarsMC), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    ! get Dimension sizes
    call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%mol_mc)%dimName), nMolecule%MC_ROUTE, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsMC
      select case(iVar)
        case(ixMC%qsub); allocate(state%var(iVar)%array_3d_dp(nSeg, nMolecule%MC_ROUTE, nens), stat=ierr)
        case(ixMC%vol : ixMC%qerror); allocate(state%var(iVar)%array_2d_dp(nSeg, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for MC routing state:'//trim(meta_mc(iVar)%varName); return; endif
    end do

    do iVar=1,nVarsMC
      select case(iVar)
        case(ixMC%qsub); call get_nc(fname, meta_mc(iVar)%varName, state%var(iVar)%array_3d_dp, [1,1,1], [nSeg,nMolecule%MC_ROUTE,nens], ierr, cmessage1)
        case(ixMC%vol);  call get_nc(fname, meta_mc(iVar)%varName, state%var(iVar)%array_2d_dp, [1,1], [nSeg, nens], ierr, cmessage1)
        case(ixMC%qerror)
          call open_nc(fname, 'r', ncidRestart, ierr, cmessage1)
          if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif
          if (check_variable(ncidRestart,meta_mc(iVar)%varName)) then
            call get_nc(fname, meta_mc(iVar)%varName, state%var(iVar)%array_2d_dp, (/1,1/), (/nSeg, nens/), ierr, cmessage1)
          else
            state%var(iVar)%array_2d_dp = 0._dp
          end if
          call close_nc(ncidRestart, ierr, cmessage1)
          if(ierr/=0)then; message=trim(message1)//trim(cmessage1); return; endif
        case default; ierr=20; message1=trim(message1)//'unable to identify MC variable index for nc reading'; return
      end select
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_mc(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        jSeg = ixRch_order(iSeg)
        allocate(RCHSTA(iens,jSeg)%MC_ROUTE%molecule%Q(nMolecule%MC_ROUTE), stat=ierr, errmsg=cmessage1)
        do iVar=1,nVarsMC
          select case(iVar)
            case(ixMC%qsub); RCHSTA(iens,jSeg)%MC_ROUTE%molecule%Q(1:nMolecule%MC_ROUTE) = state%var(iVar)%array_3d_dp(iSeg,1:nMolecule%MC_ROUTE,iens)
            case(ixMC%vol);  RCHFLX(iens,jSeg)%ROUTE(idxMC)%REACH_VOL(1) = state%var(iVar)%array_2d_dp(iSeg,iens)
            case(ixMC%qerror); RCHFLX(iens,iSeg)%ROUTE(idxMC)%Qerror = state%var(iVar)%array_2d_dp(iSeg,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify MC routing state variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_MC_state

  SUBROUTINE read_DW_state(ierr, message1)

    USE globalData, ONLY: meta_dw
    USE globalData, ONLY: idxDW
    USE globalData, ONLY: nMolecule
    USE var_lookup, ONLY: ixDW, nVarsDW
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    integer(i4b)                  :: ncidRestart    ! restart netcdf id
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
    integer(i4b)                  :: jSeg           ! index loops for reaches respectively

    ierr=0; message1='read_DW_state/'

    allocate(state%var(nVarsDW), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    ! get Dimension sizes
    call get_nc_dim_len(fname, trim(meta_stateDims(ixStateDims%mol_dw)%dimName), nMolecule%DW_ROUTE, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsDW
      select case(iVar)
        case(ixDW%qsub); allocate(state%var(iVar)%array_3d_dp(nSeg, nMolecule%DW_ROUTE, nens), stat=ierr)
        case(ixDW%vol : ixDW%qerror); allocate(state%var(iVar)%array_2d_dp(nSeg, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for DW routing state:'//trim(meta_dw(iVar)%varName); return; endif
    end do

    do iVar=1,nVarsDW
      select case(iVar)
        case(ixDW%qsub); call get_nc(fname, meta_dw(iVar)%varName, state%var(iVar)%array_3d_dp, [1,1,1], [nSeg,nMolecule%DW_ROUTE,nens], ierr, cmessage1)
        case(ixDW%vol);  call get_nc(fname, meta_dw(iVar)%varName, state%var(iVar)%array_2d_dp, [1,1], [nSeg, nens], ierr, cmessage1)
        case(ixDW%qerror)
          call open_nc(fname, 'r', ncidRestart, ierr, cmessage1)
          if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif
          if (check_variable(ncidRestart,meta_dw(iVar)%varName)) then
            call get_nc(fname, meta_dw(iVar)%varName, state%var(iVar)%array_2d_dp, (/1,1/), (/nSeg, nens/), ierr, cmessage1)
          else
            state%var(iVar)%array_2d_dp = 0._dp
          end if
          call close_nc(ncidRestart, ierr, cmessage1)
          if(ierr/=0)then; message=trim(message1)//trim(cmessage1); return; endif
        case default; ierr=20; message1=trim(message1)//'unable to identify DW variable index for nc reading'; return
      end select
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_dw(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        jSeg = ixRch_order(iSeg)
        allocate(RCHSTA(iens,jSeg)%DW_ROUTE%molecule%Q(nMolecule%DW_ROUTE), stat=ierr, errmsg=cmessage1)
        do iVar=1,nVarsDW
          select case(iVar)
            case(ixDW%qsub); RCHSTA(iens,jSeg)%DW_ROUTE%molecule%Q(1:nMolecule%DW_ROUTE) = state%var(iVar)%array_3d_dp(iSeg,1:nMolecule%DW_ROUTE,iens)
            case(ixDW%vol);  RCHFLX(iens,jSeg)%ROUTE(idxDW)%REACH_VOL(1) = state%var(iVar)%array_2d_dp(iSeg,iens)
            case(ixDW%qerror); RCHFLX(iens,iSeg)%ROUTE(idxDW)%Qerror = state%var(iVar)%array_2d_dp(iSeg,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify DW routing state variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_DW_state

 END SUBROUTINE read_state_nc


END MODULE read_restart
