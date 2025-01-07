MODULE pio_decomp_data

  USE nrtype

  USE var_lookup, ONLY: ixStateDims
  USE var_lookup, ONLY: ixHFLX
  USE globalData, ONLY: masterproc
  USE globalData, ONLY: pid
  USE globalData, ONLY: ioDesc_hru_float
  USE globalData, ONLY: ioDesc_rch_float
  USE globalData, ONLY: ioDesc_gauge_float
  USE globalData, ONLY: ioDesc_rch_int
  USE globalData, ONLY: ioDesc_rch_double
  USE globalData, ONLY: ioDesc_hru_double
  USE globalData, ONLY: ioDesc_wave_int
  USE globalData, ONLY: ioDesc_wave_double
  USE globalData, ONLY: ioDesc_mesh_kw_double
  USE globalData, ONLY: ioDesc_mesh_mc_double
  USE globalData, ONLY: ioDesc_mesh_dw_double
  USE globalData, ONLY: ioDesc_irf_double
  USE globalData, ONLY: ioDesc_vol_double
  USE globalData, ONLY: ioDesc_irf_bas_double
  USE globalData, ONLY: index_write_gage

  private

  public::set_pio_decomp
  public::get_compdof_all_network ! used in write_restart_pio

CONTAINS

  SUBROUTINE set_pio_decomp(ierr, message)

    USE public_var, ONLY: doesBasinRoute
    USE public_var, ONLY: impulseResponseFunc
    USE public_var, ONLY: kinematicWaveTracking
    USE public_var, ONLY: kinematicWave
    USE public_var, ONLY: muskingumCunge
    USE public_var, ONLY: diffusiveWave
    USE globalData, ONLY: meta_stateDims         !
    USE globalData, ONLY: meta_hflx         !
    USE globalData, ONLY: nRch                   !
    USE globalData, ONLY: nHRU                   !
    USE globalData, ONLY: onRoute                !
    USE globalData, ONLY: gage_data
    USE globalData, ONLY: pioSystem              !
    USE pio_utils,  ONLY: pio_decomp
    USE pio_utils,  ONLY: ncd_int
    USE pio_utils,  ONLY: ncd_float
    USE pio_utils,  ONLY: ncd_double

    ! populate state netCDF dimension size
    USE public_var, ONLY: MAXQPAR, outputAtGage
    USE globalData, ONLY: nMolecule
    USE globalData, ONLY: maxtdh           ! maximum unit-hydrogrph future time
    USE globalData, ONLY: FRAC_FUTURE      ! To get size of q future for basin IRF

    implicit none
    ! argument variables
    integer(i4b),   intent(out)     :: ierr                ! error code
    character(*),   intent(out)     :: message             ! error message
    ! local variables
    integer(i4b), allocatable       :: compdof_rch(:)      !
    integer(i4b), allocatable       :: compdof_rch_gage(:) !
    integer(i4b), allocatable       :: compdof_hru(:)      !
    character(len=strLen)           :: cmessage            ! error message of downwind routine

    ierr=0; message='set_pio_decomp/'

    associate(ndim_seg      => meta_stateDims(ixStateDims%seg)%dimLength,      & !
              ndim_hru      => meta_stateDims(ixStateDims%hru)%dimLength,      & !
              ndim_hist_fil => meta_stateDims(ixStateDims%hist_fil)%dimLength, & !
              ndim_tdh      => meta_stateDims(ixStateDims%tdh)%dimLength,      & ! maximum future q time steps among basins
              ndim_tdh_irf  => meta_stateDims(ixStateDims%tdh_irf)%dimLength,  & ! maximum future q time steps among reaches
              ndim_tbound   => meta_stateDims(ixStateDims%tbound)%dimLength,   & ! time bound
              ndim_Mesh_kw  => meta_stateDims(ixStateDims%mol_kw)%dimLength,   & ! kw_finite difference mesh points
              ndim_Mesh_mc  => meta_stateDims(ixStateDims%mol_mc)%dimLength,   & ! mc_finite difference mesh points
              ndim_Mesh_dw  => meta_stateDims(ixStateDims%mol_dw)%dimLength,   & ! dw_finite difference mesh points
              ndim_Wave     => meta_stateDims(ixStateDims%wave)%dimLength)      ! maximum waves allowed in a reach

    ! set state file dimension length
    ndim_seg = nRch
    ndim_hru = nHru
    ndim_tbound = 2
    if (doesBasinRoute==1)              ndim_tdh = size(FRAC_FUTURE)
    if (onRoute(impulseResponseFunc))   ndim_tdh_irf = maxtdh
    if (onRoute(kinematicWave))         ndim_Mesh_kw = nMolecule%KW_ROUTE
    if (onRoute(muskingumCunge))        ndim_Mesh_mc = nMolecule%MC_ROUTE
    if (onRoute(diffusiveWave))         ndim_Mesh_dw = nMolecule%DW_ROUTE
    if (onRoute(kinematicWaveTracking)) ndim_wave = MAXQPAR
    if (outputAtGage) then
      ndim_hist_fil= 2
    else
      ndim_hist_fil= 1
    end if

    ! compute PIO domain decomposition for entire network
    call get_compdof_all_network(compdof_rch, compdof_hru)

    ! For gauge history files
    if (outputAtGage) then
      call get_compdof_gage(compdof_rch_gage, ierr, cmessage)

      call pio_decomp(pioSystem,         & ! inout: pio system descriptor
                      ncd_float,         & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [gage_data%nGage], & ! input: dimension length == global array size
                      compdof_rch_gage,  & ! input: local->global mapping
                      ioDesc_gauge_float)
    end if

    ! For history files for entire network ---------
    call pio_decomp(pioSystem,         & ! inout: pio system descriptor
                    ncd_float,         & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                    [ndim_seg],        & ! input: dimension length == global array size
                    compdof_rch,       & ! input: local->global mapping
                    ioDesc_rch_float)

    call pio_decomp(pioSystem,         & ! inout: pio system descriptor
                    ncd_float,         & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                    [ndim_hru],        & ! input: dimension length == global array size
                    compdof_hru,       & ! input: local->global mapping
                    ioDesc_hru_float)

    ! For restart file ----------
    ! type: float  dim: [dim_seg]  -- channel runoff coming from hru
    call pio_decomp(pioSystem,              & ! inout: pio system descriptor
                    ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                    [ndim_seg],             & ! input: dimension length == global array size
                    compdof_rch,            & ! input:
                    ioDesc_rch_double)

    ! type: int  dim: [dim_seg] -- number of wave or uh future time steps
    call pio_decomp(pioSystem,              & ! inout: pio system descriptor
                    ncd_int,                & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                    [ndim_seg],             & ! input: dimension length == global array size
                    compdof_rch,            & ! input:
                    ioDesc_rch_int)

    if (meta_hflx(ixHFLX%basRunoff)%varFile) then
      ! type: single precision float  dim: [dim_hru,dim_ens]  --
      call pio_decomp(pioSystem,              & ! inout: pio system descriptor
                      ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_hru],             & ! input: dimension length == global array size
                      compdof_hru,            & ! input:
                      ioDesc_hru_double)
    end if

    if (doesBasinRoute==1) then
      ! type: float dim: [dim_seg, dim_tdh_irf]
      call pio_decomp(pioSystem,                     & ! inout: pio system descriptor
                      ncd_double,                    & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_seg,ndim_tdh],           & ! input: dimension length == global array size
                      compdof_rch,                   & ! input:
                      ioDesc_irf_bas_double)
    end if

    if (onRoute(impulseResponseFunc))then
      ! type: float dim: [dim_seg, dim_tdh_irf]
      call pio_decomp(pioSystem,                         & ! inout: pio system descriptor
                      ncd_double,                        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_seg,ndim_tdh_irf],           & ! input: dimension length == global array size
                      compdof_rch,                       & ! input:
                      ioDesc_irf_double)

      ! type: float dim: [dim_seg]
      call pio_decomp(pioSystem,                         & ! inout: pio system descriptor
                      ncd_double,                        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_seg,ndim_tbound],            & ! input: dimension length == global array size
                      compdof_rch,                       & ! input:
                      ioDesc_vol_double)
    end if

    if (onRoute(kinematicWaveTracking)) then
      ! type: int, dim: [dim_seg, dim_wave]
      call pio_decomp(pioSystem,                         & ! inout: pio system descriptor
                      ncd_int,                           & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_seg,ndim_wave],              & ! input: dimension length == global array size
                      compdof_rch,                       & ! input:
                      ioDesc_wave_int)

      ! type: float, dim: [dim_seg, dim_wave]
      call pio_decomp(pioSystem,                         & ! inout: pio system descriptor
                      ncd_double,                        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_seg,ndim_wave],              & ! input: dimension length == global array size
                      compdof_rch,                       & ! input:
                      ioDesc_wave_double)
    end if

    if (onRoute(kinematicWave)) then
      ! type: float, dim: [dim_seg, dim_mesh]
      call pio_decomp(pioSystem,                         & ! inout: pio system descriptor
                      ncd_double,                        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_seg,ndim_Mesh_kw],           & ! input: dimension length == global array size
                      compdof_rch,                       & ! input:
                      ioDesc_mesh_kw_double)
    end if

    if (onRoute(muskingumCunge)) then
      ! type: float, dim: [dim_seg, dim_mesh]
      call pio_decomp(pioSystem,                         & ! inout: pio system descriptor
                      ncd_double,                        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_seg,ndim_Mesh_mc],           & ! input: dimension length == global array size
                      compdof_rch,                       & ! input:
                      ioDesc_mesh_mc_double)
    end if

    if (onRoute(diffusiveWave)) then
      ! type: float, dim: [dim_seg, dim_mesh]
      call pio_decomp(pioSystem,                         & ! inout: pio system descriptor
                      ncd_double,                        & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                      [ndim_seg,ndim_Mesh_dw],           & ! input: dimension length == global array size
                      compdof_rch,                       & ! input:
                      ioDesc_mesh_dw_double)
    end if

    end associate

  END SUBROUTINE set_pio_decomp

  ! *********************************************************************
  ! private subroutine: get pio local-global mapping data
  ! *********************************************************************
  SUBROUTINE get_compdof_all_network(compdof_rch, compdof_hru)

    USE globalData,        ONLY: hru_per_proc      ! number of hrus assigned to each proc (size = num of procs+1)
    USE globalData,        ONLY: rch_per_proc      ! number of reaches assigned to each proc (size = num of procs+1)
    USE nr_utils,          ONLY: arth

    implicit none
    ! Argument variables
    integer(i4b), allocatable, intent(out) :: compdof_rch(:) !
    integer(i4b), allocatable, intent(out) :: compdof_hru(:) !
    ! Local variables
    integer(i4b)                :: ix1, ix2               ! frst and last indices of global array for local array chunk

    if (masterproc) then
      allocate(compdof_rch(sum(rch_per_proc(-1:pid))))
      allocate(compdof_hru(sum(hru_per_proc(-1:pid))))
    else
      allocate(compdof_rch(rch_per_proc(pid)))
      allocate(compdof_hru(hru_per_proc(pid)))
    end if

    ! For reach flux/volume
    if (masterproc) then
      ix1 = 1
    else
      ix1 = sum(rch_per_proc(-1:pid-1))+1
    endif
    ix2 = sum(rch_per_proc(-1:pid))
    compdof_rch = arth(ix1, 1, ix2-ix1+1)

    ! For HRU flux/volume
    if (masterproc) then
      ix1 = 1
    else
      ix1 = sum(hru_per_proc(-1:pid-1))+1
    endif
    ix2 = sum(hru_per_proc(-1:pid))
    compdof_hru = arth(ix1, 1, ix2-ix1+1)

  END SUBROUTINE get_compdof_all_network

  ! *********************************************************************
  ! private subroutine: get pio local-global mapping data
  ! *********************************************************************
  SUBROUTINE get_compdof_gage(compdof_rch, ierr, message)

    USE globalData, ONLY: rch_per_proc        ! number of reaches assigned to each proc (size = num of procs+1)
    USE globalData, ONLY: nRch_mainstem       ! number of mainstem reaches
    USE globalData, ONLY: nTribOutlet         ! number of tributary outlets flowing to the mainstem
    USE globalData, ONLY: NETOPO_main         ! mainstem Reach neteork
    USE globalData, ONLY: NETOPO_trib         ! tributary Reach network
    USE globalData, ONLY: gage_data
    USE process_gage_meta, ONLY: reach_subset

    implicit none
    ! Argument variables
    integer(i4b), allocatable, intent(out) :: compdof_rch(:)   !
    integer(i4b), intent(out)              :: ierr             ! error code
    character(*), intent(out)              :: message          ! error message
    ! Local variables
    integer(i4b)                           :: ix
    integer(i4b)                           :: nRch_local
    integer(i4b),allocatable               :: reachID_local(:) !
    character(len=strLen)                  :: cmessage         ! error message of downwind routine

    ierr=0; message='get_compdof_gage/'

    ! Only For reach flux/volume
    if (masterproc) then
      nRch_local = nRch_mainstem+rch_per_proc(0)
      allocate(reachID_local(nRch_local), stat=ierr, errmsg=cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [reachID_local]'; return; endif

      if (nRch_mainstem>0)   reachID_local(1:nRch_mainstem) = NETOPO_main(1:nRch_mainstem)%REACHID
      if (rch_per_proc(0)>0) reachID_local(nRch_mainstem+1:nRch_local) = NETOPO_trib(:)%REACHID
    else
      nRch_local = rch_per_proc(pid)
      allocate(reachID_local(nRch_local), stat=ierr, errmsg=cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [reachID_local]'; return; endif
      reachID_local = NETOPO_trib(:)%REACHID
    endif

    call reach_subset(reachID_local, gage_data, ierr, cmessage, compdof=compdof_rch, index2=index_write_gage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Need to adjust tributary indices in root processor
    ! This is because RCHFLX has three components in the order: mainstem, halo, tributary
    ! mainstem (1,2,..,nRch_mainstem),
    ! halo(nRch_mainstem+1,..,nRch_mainstem+nTribOutlet), and
    ! tributary(nRch_mainstem+nTribOutlet+1,...,nRch_total)
    ! index_write_gage is computed based on reachID consisting of mainstem and tributary
    if (masterproc) then
      do ix=1,size(index_write_gage)
        if (index_write_gage(ix)<=nRch_mainstem) cycle
        index_write_gage(ix) = index_write_gage(ix) + nTribOutlet
      end do
    end if

  END SUBROUTINE get_compdof_gage

END MODULE pio_decomp_data
