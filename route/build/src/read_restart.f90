MODULE read_restart

! Moudle wide external modules
USE nrtype, only: i4b, dp, strLen
USE public_var
USE io_netcdf, ONLY: open_nc
USE io_netcdf, ONLY: close_nc
USE io_netcdf, ONLY: get_nc
USE io_netcdf, ONLY: get_nc_dim_len

implicit none

private

public::read_state_nc

CONTAINS

 ! *********************************************************************
 ! public subroutine: read routing state NetCDF file
 ! *********************************************************************
 SUBROUTINE read_state_nc(&
                          fname,           &   ! Input:  state netcdf name
                          T0, T1,          &   ! output: start and end time [sec]
                          ierr, message)       ! Output: error control

 USE globalData, ONLY: RCHFLX                    ! To get q future for basin IRF and IRF (these should not be in this data strucuture)
 USE globalData, ONLY: RCHSTA                    ! restart state data structure
 USE dataTypes,  ONLY: states
 USE globalData, ONLY: meta_stateDims            ! dimension for state variables
 USE globalData, ONLY: onRoute                   ! logical to indicate which routing method(s) is on
 USE public_var, ONLY: impulseResponseFunc
 USE public_var, ONLY: kinematicWaveTracking
 USE public_var, ONLY: kinematicWave
 USE public_var, ONLY: muskingumCunge
 USE public_var, ONLY: diffusiveWave
 USE var_lookup, ONLY: ixStateDims, nStateDims

 implicit none
 ! input variables
 character(*), intent(in)      :: fname                ! filename
 ! output variables
 real(dp),     intent(out)     :: T0                   ! beginning time [sec] of ith time step - lapse time from the beginning of the simulation
 real(dp),     intent(out)     :: T1                   ! ending time [sec] ith time step - lapse time from the beginning of the simulation
 integer(i4b), intent(out)     :: ierr                 ! error code
 character(*), intent(out)     :: message              ! error message
 ! local variables
 integer(i4b)                  :: ncidRestart          ! restart netcdf id
 real(dp)                      :: TB(2)                ! 2 element-time bound vector
 integer(i4b)                  :: nSeg,nens            ! dimenion sizes
 integer(i4b)                  :: ntbound              ! dimenion sizes
 integer(i4b)                  :: ixDim_common(3)      ! custom dimension ID array
 integer(i4b)                  :: jDim                 ! index loops for dimension
 character(len=strLen)         :: cmessage             ! error message of downwind routine

 ierr=0; message='read_state_nc/'

 ! get Dimension sizes
 ! For common dimension/variables - seg id, time-bound -----------
 ixDim_common = (/ixStateDims%seg, ixStateDims%ens, ixStateDims%tbound/)

 call open_nc(fname, 'r', ncidRestart, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 do jDim=1,size(ixDim_common)
   associate (ixDim_tmp => ixDim_common(jDim))
   select case(ixDim_tmp)
    case(ixStateDims%seg);     call get_nc_dim_len(ncidRestart, meta_stateDims(ixDim_tmp)%dimName, nSeg,    ierr, cmessage)
    case(ixStateDims%ens);     call get_nc_dim_len(ncidRestart, meta_stateDims(ixDim_tmp)%dimName, nens,    ierr, cmessage)
    case(ixStateDims%tbound);  call get_nc_dim_len(ncidRestart, meta_stateDims(ixDim_tmp)%dimName, ntbound, ierr, cmessage)
    case default; ierr=20; message=trim(message)//'unable to identify dimension name index'; return
   end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end associate
 enddo

 ! Read variables
 ! time bound
 call get_nc(ncidRestart,'time_bound',TB(:), 1, 2, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 T0=TB(1); T1=TB(2)

 ! routing specific variables
 if (doesBasinRoute == 1) then
   call read_IRFbas_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return;endif
 endif

 call read_basinQ_state(ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return;endif

 if (onRoute(kinematicWaveTracking)) then
   call read_KWT_state(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage);return; endif
 end if

 if (onRoute(impulseResponseFunc)) then
   call read_IRF_state(ierr, cmessage)
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

 call close_nc(ncidRestart, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 CONTAINS

  SUBROUTINE read_basinQ_state(ierr, message1)
    USE globalData, ONLY: meta_basinQ               ! reach inflow from basin at previous time step
    USE var_lookup, ONLY: ixBasinQ, nVarsBasinQ
    implicit none
    ! output
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively

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
       case(ixBasinQ%q); call get_nc(ncidRestart, meta_basinQ(iVar)%varName, state%var(iVar)%array_2d_dp, (/1,1/), (/nSeg,nens/), ierr, cmessage1)
       case default; ierr=20; message1=trim(message1)//'unable to identify previous time step reach inflow variable index for nc writing'; return
     end select
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_basinQ(iVar)%varName); return; endif

    enddo

    do iens=1,nens
      do iSeg=1,nSeg
        do iVar=1,nVarsBasinQ
          select case(iVar)
            case(ixBasinQ%q); RCHFLX(iens,iSeg)%BASIN_QR(1) = state%var(iVar)%array_2d_dp(iSeg,iens)
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
    integer(i4b)                  :: ntdh           ! dimension size

    ! initialize error control
    ierr=0; message1='read_IRFbas_state/'

    call get_nc_dim_len(ncidRestart, meta_stateDims(ixStateDims%tdh)%dimName, ntdh, ierr, cmessage1)
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
        case(ixIRFbas%qfuture); call get_nc(ncidRestart, meta_irf_bas(iVar)%varName, state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,ntdh,nens/), ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF variable index for nc writing'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_irf_bas(iVar)%varName); return; endif
    enddo

    do iens=1,nens
      do iSeg=1,nSeg
        allocate(RCHFLX(iens,iSeg)%QFUTURE(ntdh), stat=ierr, errmsg=cmessage1)
        if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif
        do iVar=1,nVarsIRFbas
          select case(iVar)
            case(ixIRFbas%qfuture); RCHFLX(iens,iSeg)%QFUTURE(:)  = state%var(iVar)%array_3d_dp(iSeg,:,iens)
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
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
    integer(i4b), allocatable     :: numQF(:,:)     ! number of future Q time steps for each ensemble and segment
    integer(i4b)                  :: ntdh_irf       ! dimenion sizes

    ierr=0; message1='read_IRF_state/'

    call get_nc_dim_len(ncidRestart, meta_stateDims(ixStateDims%tdh_irf)%dimName, ntdh_irf, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    allocate(state%var(nVarsIRF), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    allocate(numQF(nens,nSeg), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsIRF
      select case(iVar)
        case(ixIRF%qfuture); allocate(state%var(iVar)%array_3d_dp(nSeg, ntdh_irf, nens), stat=ierr)
        case(ixIRF%irfVol);  allocate(state%var(iVar)%array_2d_dp(nSeg, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for IRF routing state:'//trim(meta_irf(iVar)%varName); return; endif
    end do

    call get_nc(ncidRestart,'numQF',numQF,(/1,1/),(/nSeg,nens/),ierr,cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_irf(iVar)%varName); return; endif

    do iVar=1,nVarsIRF
      select case(iVar)
        case(ixIRF%qfuture); call get_nc(ncidRestart, meta_irf(iVar)%varName, state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,ntdh_irf,nens/), ierr, cmessage1)
        case(ixIRF%irfVol);  call get_nc(ncidRestart, meta_irf(iVar)%varName, state%var(iVar)%array_2d_dp, (/1,1/), (/nSeg, nens/), ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc reading'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_irf(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        allocate(RCHFLX(iens,iSeg)%QFUTURE_IRF(numQF(iens,iSeg)), stat=ierr, errmsg=cmessage1)
        if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif
        do iVar=1,nVarsIRF
          select case(iVar)
            case(ixIRF%qfuture); RCHFLX(iens,iSeg)%QFUTURE_IRF  = state%var(iVar)%array_3d_dp(iSeg,1:numQF(iens,iSeg),iens)
            case(ixIRF%irfVol);  RCHFLX(iens,iSeg)%ROUTE(idxIRF)%REACH_VOL(1) = state%var(iVar)%array_2d_dp(iSeg,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
          end select
        enddo ! variable loop
      enddo ! seg loop
    enddo ! ensemble loop

  END SUBROUTINE read_IRF_state


  SUBROUTINE read_KWT_state(ierr, message1)
    USE globalData, ONLY: meta_kwt
    USE var_lookup, ONLY: ixKWT, nVarsKWT
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively
    integer(i4b)                  :: nwave          ! dimenion sizes
    integer(i4b), allocatable     :: RFvec(:)       ! temporal vector
    integer(i4b), allocatable     :: numWaves(:,:)  ! number of waves for each ensemble and segment
    ! initialize error control
    ierr=0; message1='read_KWT_state/'

    allocate(state%var(nVarsKWT), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    allocate(numWaves(nens,nSeg), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    ! get Dimension sizes
    call get_nc_dim_len(ncidRestart, meta_stateDims(ixStateDims%wave)%dimName, nwave, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsKWT

      select case(iVar)
       case(ixKWT%routed); allocate(state%var(iVar)%array_3d_dp(nSeg, nwave, nens), stat=ierr)
       case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
        allocate(state%var(iVar)%array_3d_dp(nSeg, nwave, nens), stat=ierr)
       case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWT routing state:'//trim(meta_kwt(iVar)%varName); return; endif
    end do

    call get_nc(ncidRestart,'numWaves',numWaves, (/1,1/), (/nSeg,nens/), ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsKWT
      select case(iVar)
        case(ixKWT%routed)
          call get_nc(ncidRestart, meta_kwt(iVar)%varName, state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nwave,nens/), ierr, cmessage1)
        case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
          call get_nc(ncidRestart, meta_kwt(iVar)%varName, state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nwave,nens/), ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify KWT variable index for nc reading'; return
      end select
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_kwt(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        allocate(RCHSTA(iens,iSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1), stat=ierr)
        do iVar=1,nVarsKWT
          select case(iVar)
            case(ixKWT%tentry);    RCHSTA(iens,iSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%TI = state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
            case(ixKWT%texit);     RCHSTA(iens,iSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%TR = state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
            case(ixKWT%qwave);     RCHSTA(iens,iSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%QF = state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
            case(ixKWT%qwave_mod); RCHSTA(iens,iSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%QM = state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens)
            case(ixKWT%routed) ! this is suppposed to be logical variable, but put it as 0 or 1 in double now
              if (allocated(RFvec)) deallocate(RFvec, stat=ierr)
              allocate(RFvec(0:numWaves(iens,iSeg)-1),stat=ierr)
              RFvec = nint(state%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens))
              RCHSTA(iens,iSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%RF=.False.
              where (RFvec==1_i4b) RCHSTA(iens,iSeg)%LKW_ROUTE%KWAVE(0:numWaves(iens,iSeg)-1)%RF=.True.
            case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
          end select
        end do
      end do
    end do

  END SUBROUTINE read_KWT_state


  SUBROUTINE read_KW_state(ierr, message1)
    USE globalData, ONLY: meta_kw
    USE globalData, ONLY: nMolecule
    USE var_lookup, ONLY: ixKW, nVarsKW
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively

    ierr=0; message1='read_KW_state/'

    allocate(state%var(nVarsKW), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    ! get Dimension sizes
    call get_nc_dim_len(ncidRestart, trim(meta_stateDims(ixStateDims%mol_kw)%dimName), nMolecule%KW_ROUTE, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsKW
      select case(iVar)
        case(ixKW%qsub); allocate(state%var(iVar)%array_3d_dp(nSeg, nMolecule%KW_ROUTE, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KW routing state:'//trim(meta_kw(iVar)%varName); return; endif
    end do

    do iVar=1,nVarsKW
      select case(iVar)
        case(ixKW%qsub)
          call get_nc(ncidRestart, trim(meta_kw(iVar)%varName), state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nMolecule%KW_ROUTE,nens/), ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify KW variable index for nc reading'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_kw(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        allocate(RCHSTA(iens,iSeg)%KW_ROUTE%molecule%Q(nMolecule%KW_ROUTE), stat=ierr)
        do iVar=1,nVarsKW
          select case(iVar)
            case(ixKW%qsub); RCHSTA(iens,iSeg)%KW_ROUTE%molecule%Q(1:nMolecule%KW_ROUTE) = state%var(iVar)%array_3d_dp(iSeg,1:nMolecule%KW_ROUTE,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify KW routing state variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_KW_state


  SUBROUTINE read_MC_state(ierr, message1)
    USE globalData, ONLY: meta_mc
    USE globalData, ONLY: nMolecule
    USE var_lookup, ONLY: ixMC, nVarsMC
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively

    ierr=0; message1='read_MC_state/'

    allocate(state%var(nVarsMC), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    ! get Dimension sizes
    call get_nc_dim_len(ncidRestart, trim(meta_stateDims(ixStateDims%mol_mc)%dimName), nMolecule%MC_ROUTE, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsMC
      select case(iVar)
        case(ixMC%qsub); allocate(state%var(iVar)%array_3d_dp(nSeg, nMolecule%MC_ROUTE, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for MC routing state:'//trim(meta_mc(iVar)%varName); return; endif
    end do

    do iVar=1,nVarsMC
      select case(iVar)
        case(ixMC%qsub)
          call get_nc(ncidRestart, trim(meta_mc(iVar)%varName), state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nMolecule%MC_ROUTE,nens/), ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify MC variable index for nc reading'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_mc(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        allocate(RCHSTA(iens,iSeg)%MC_ROUTE%molecule%Q(nMolecule%MC_ROUTE), stat=ierr)
        do iVar=1,nVarsMC
          select case(iVar)
            case(ixMC%qsub); RCHSTA(iens,iSeg)%MC_ROUTE%molecule%Q(1:nMolecule%MC_ROUTE) = state%var(iVar)%array_3d_dp(iSeg,1:nMolecule%MC_ROUTE,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify MC routing state variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_MC_state


  SUBROUTINE read_DW_state(ierr, message1)
    USE globalData, ONLY: meta_dw
    USE globalData, ONLY: nMolecule
    USE var_lookup, ONLY: ixDW, nVarsDW
    implicit none
    integer(i4b), intent(out)     :: ierr           ! error code
    character(*), intent(out)     :: message1       ! error message
    ! local variables
    character(len=strLen)         :: cmessage1      ! error message of downwind routine
    type(states)                  :: state          ! temporal state data structures
    integer(i4b)                  :: iVar,iens,iSeg ! index loops for variables, ensembles, reaches respectively

    ierr=0; message1='read_DW_state/'

    allocate(state%var(nVarsDW), stat=ierr, errmsg=cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    ! get Dimension sizes
    call get_nc_dim_len(ncidRestart, trim(meta_stateDims(ixStateDims%mol_dw)%dimName), nMolecule%DW_ROUTE, ierr, cmessage1)
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage1); return; endif

    do iVar=1,nVarsDW
      select case(iVar)
        case(ixDW%qsub); allocate(state%var(iVar)%array_3d_dp(nSeg, nMolecule%DW_ROUTE, nens), stat=ierr)
        case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//'problem allocating space for DW routing state:'//trim(meta_dw(iVar)%varName); return; endif
    end do

    do iVar=1,nVarsDW
      select case(iVar)
        case(ixDW%qsub)
          call get_nc(ncidRestart, trim(meta_dw(iVar)%varName), state%var(iVar)%array_3d_dp, (/1,1,1/), (/nSeg,nMolecule%DW_ROUTE,nens/), ierr, cmessage1)
        case default; ierr=20; message1=trim(message1)//'unable to identify DW variable index for nc reading'; return
      end select
      if(ierr/=0)then; message1=trim(message1)//trim(cmessage1)//':'//trim(meta_dw(iVar)%varName); return; endif
    end do

    do iens=1,nens
      do iSeg=1,nSeg
        allocate(RCHSTA(iens,iSeg)%DW_ROUTE%molecule%Q(nMolecule%DW_ROUTE), stat=ierr)
        do iVar=1,nVarsDW
          select case(iVar)
            case(ixDW%qsub); RCHSTA(iens,iSeg)%DW_ROUTE%molecule%Q(1:nMolecule%DW_ROUTE) = state%var(iVar)%array_3d_dp(iSeg,1:nMolecule%DW_ROUTE,iens)
            case default; ierr=20; message1=trim(message1)//'unable to identify DW routing state variable index'; return
          end select
        enddo
      enddo
    enddo

  END SUBROUTINE read_DW_state

 END SUBROUTINE read_state_nc

END MODULE read_restart
