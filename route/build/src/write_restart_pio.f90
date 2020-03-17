MODULE write_restart_pio

! Moudle wide shared data
USE nrtype
USE dataTypes,         ONLY: STRFLX            ! fluxes in each reach
USE dataTypes,         ONLY: STRSTA            ! state in each reach
USE dataTypes,         ONLY: RCHTOPO           ! Network topology
USE dataTypes,         ONLY: states
USE public_var,        ONLY: iulog             ! i/o logical unit number
USE public_var,        ONLY: integerMissing
USE public_var,        ONLY: realMissing
USE public_var,        ONLY: rpntfil           ! ascii containing last restart file (used in coupled mode)
USE globalData,        ONLY: pid, nNodes
USE globalData,        ONLY: masterproc
USE globalData,        ONLY: mpicom_route
USE globalData,        ONLY: pio_netcdf_format
USE globalData,        ONLY: pio_typename
USE globalData,        ONLY: pio_numiotasks
USE globalData,        ONLY: pio_rearranger
USE globalData,        ONLY: pio_root
USE globalData,        ONLY: pio_stride
! Moudle wide external modules
USE nr_utility_module, ONLY: arth
USE pio_utils

implicit none

! The following variables used only in this module
type(iosystem_desc_t),save :: pioSystemState
type(file_desc_t),    save :: pioFileDescState     ! contains data identifying the file
type(io_desc_t),      save :: iodesc_state_int
type(io_desc_t),      save :: iodesc_state_double
type(io_desc_t),      save :: iodesc_wave_int
type(io_desc_t),      save :: iodesc_wave_double
type(io_desc_t),      save :: iodesc_mesh_double
type(io_desc_t),      save :: iodesc_irf_float
type(io_desc_t),      save :: iodesc_irf_bas_double

integer(i4b),parameter     :: recordDim=-999       ! record dimension Indicator

integer(i4b),         save :: kTime                ! Time index in restart netCDF

private

public::define_state_nc
public::output_state

CONTAINS

 ! *********************************************************************
 ! Public subroutine: output restart netcdf
 ! *********************************************************************
 subroutine output_state(ierr, message)

  USE public_var, ONLY: output_dir
  USE public_var, ONLY: case_name         ! simulation name ==> output filename head
  USE globalData, ONLY: modTime           ! previous and current model time

  implicit none
  ! output variables
  integer(i4b),   intent(out)          :: ierr             ! error code
  character(*),   intent(out)          :: message          ! error message
  ! local variables
  character(len=strLen)                :: fileout_state    ! name of the output file
  character(*),parameter               :: fmtYMDS='(a,I0.4,a,I0.2,a,I0.2,a,I0.5,a)'
  integer(i4b)                         :: sec_in_day      ! second within day
  character(len=strLen)                :: cmessage         ! error message of downwind routine

  ! Define filename
  sec_in_day = 0
  write(fileout_state, fmtYMDS) trim(output_dir)//trim(case_name)//'.mizuRoute.r.', &
                          modTime(0)%iy, '-', modTime(0)%im, '-', modTime(0)%id, '-',sec_in_day,'.nc'

  ! Define output state netCDF
  call define_state_nc(trim(fileout_state), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call write_state_nc(fileout_state, ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  open(1, file = trim(output_dir)//trim(rpntfil), status='replace', action='write')
  write(1,*) fileout_state

 end subroutine output_state

 ! *********************************************************************
 ! subroutine: define restart NetCDF file
 ! *********************************************************************
 SUBROUTINE define_state_nc(fname,           &  ! input: filename
                            ierr, message)      ! output: error control
 ! External modules
 USE var_lookup, ONLY: ixStateDims, nStateDims  ! name index and number for state dimensions
 USE public_var, ONLY: time_units               ! time units (seconds, hours, or days)
 USE public_var, ONLY: calendar                 ! calendar name
 USE public_var, ONLY: routOpt
 USE public_var, ONLY: allRoutingMethods
 USE public_var, ONLY: doesBasinRoute
 USE public_var, ONLY: kinematicWave
 USE public_var, ONLY: kinematicWaveEuler
 USE public_var, ONLY: impulseResponseFunc
 USE globalData, ONLY: meta_stateDims           ! dimension meta for state
 USE globalData, ONLY: rch_per_proc             ! number of reaches assigned to each proc (size = num of procs+1)
 USE globalData, ONLY: nRch                     ! number of ensembles and river reaches

 implicit none
 ! input variables
 character(*),   intent(in)      :: fname            ! filename
 ! output variables
 integer(i4b),   intent(out)     :: ierr             ! error code
 character(*),   intent(out)     :: message          ! error message
 ! local variables
 integer(i4b)                    :: ix1, ix2         ! frst and last indices of global array for local array chunk
 integer(i4b)                    :: ixRch(nRch)      ! global reach index used for output
 integer(i4b)                    :: jDim             ! loop index for dimension
 integer(i4b)                    :: ixDim_common(4)  ! custom dimension ID array
 character(len=strLen)           :: cmessage         ! error message of downwind routine

 ! Initialize error control
 ierr=0; message='define_state_nc/'

 ! Initialize netCDF time index
 kTime = 0

 associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimId,     &
           dim_ens     => meta_stateDims(ixStateDims%ens)%dimId,     &
           dim_time    => meta_stateDims(ixStateDims%time)%dimId,    &
           dim_tbound  => meta_stateDims(ixStateDims%tbound)%dimId)

 ! ----------------------------------
 ! pio initialization for restart netCDF
 ! ----------------------------------
 pio_numiotasks = nNodes/pio_stride
 call pio_sys_init(pid, mpicom_route,          & ! input: MPI related parameters
                   pio_stride, pio_numiotasks, & ! input: PIO related parameters
                   pio_rearranger, pio_root,   & ! input: PIO related parameters
                   pioSystemState)               ! output: PIO system descriptors

 ! ----------------------------------
 ! Create file
 ! ----------------------------------
 call createFile(pioSystemState, trim(fname), pio_typename, pio_netcdf_format, pioFileDescState, ierr, cmessage)
 if(ierr/=0)then; message=trim(cmessage)//'cannot create state netCDF'; return; endif

 ! For common dimension/variables - seg id, time, time-bound -----------
 ixDim_common = [ixStateDims%seg, ixStateDims%ens, ixStateDims%time, ixStateDims%tbound]

 ! ----------------------------------
 ! Define dimensions
 ! ----------------------------------
 do jDim = 1,size(ixDim_common)
  associate(ixDim => ixDim_common(jDim))
  if (meta_stateDims(ixDim)%dimLength == integerMissing) then
   call set_dim_len(ixDim, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//' for '//trim(meta_stateDims(ixDim)%dimName); return; endif
  endif
  call defdim(pioFileDescState, trim(meta_stateDims(ixDim)%dimName), meta_stateDims(ixDim)%dimLength, meta_stateDims(ixDim)%dimId)
  end associate
 end do

 ! ----------------------------------
 ! Define variable
 ! ----------------------------------
 call defvar(pioFileDescState, 'reachID', [dim_seg], ncd_int, ierr, cmessage, vdesc='reach ID',  vunit='-' )
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call defVar(pioFileDescState, 'time', [dim_time], ncd_float, ierr, cmessage, vdesc='time', vunit=trim(time_units), vcal=calendar)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call defVar(pioFileDescState, 'time_bound', [dim_tbound, dim_time], ncd_float, ierr, cmessage, vdesc='time bound at last time step', vunit='sec')
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end associate

 ! Routing specific variables --------------
 ! basin IRF
 if (doesBasinRoute == 1) then
  call define_IRFbas_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! KWT routing
 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
  call define_KWT_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! KWE routing
 if (routOpt==kinematicWaveEuler) then
  call define_KWE_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! IRF routing
 if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
  call define_IRF_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! ----------------------------------
 ! pio initialization of decomposition
 ! ----------------------------------
 associate(nSeg     => meta_stateDims(ixStateDims%seg)%dimLength,     &
           nEns     => meta_stateDims(ixStateDims%ens)%dimLength,     &
           ntdh     => meta_stateDims(ixStateDims%tdh)%dimLength,     & ! maximum future q time steps among basins
           ntdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimLength, & ! maximum future q time steps among reaches
           nFdmesh  => meta_stateDims(ixStateDims%fdmesh)%dimLength,  & ! finite difference mesh points
           nWave    => meta_stateDims(ixStateDims%wave)%dimLength)      ! maximum waves allowed in a reach

 if (masterproc) then
   ix1 = 1
 else
   ix1 = sum(rch_per_proc(-1:pid-1))+1
 endif
 ix2 = sum(rch_per_proc(-1:pid))
 ixRch = arth(1,1,nSeg)

 ! type: float  dim: [dim_seg, dim_ens, dim_time]  -- channel runoff coming from hru
 call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                 ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nSeg,nEns],            & ! input: dimension length == global array size
                 ixRch(ix1:ix2),         & ! input:
                 iodesc_state_double)

 ! type: int  dim: [dim_seg, dim_ens, dim_time]  -- number of wave or uh future time steps
 call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                 ncd_int,                & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nSeg,nEns],            & ! input: dimension length == global array size
                 ixRch(ix1:ix2),         & ! input:
                 iodesc_state_int)

 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
   ! type: int, dim: [dim_seg, dim_wave, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_int,                & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,nWave,nEns],      & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_wave_int)

   ! type: float, dim: [dim_seg, dim_wave, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,nWave,nEns],      & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_wave_double)
 end if

 if (routOpt==kinematicWaveEuler) then
   ! type: float, dim: [dim_seg, dim_fdmesh, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,nFdmesh,nEns],    & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_mesh_double)
 end if

 if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
   ! type: float dim: [dim_seg, dim_tdh_irf, dim_ens, dim_time]
   call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                   ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                   [nSeg,ntdh_irf,nEns],   & ! input: dimension length == global array size
                   ixRch(ix1:ix2),         & ! input:
                   iodesc_irf_float)

 end if
 ! type: float dim: [dim_seg, dim_tdh_irf, dim_ens, dim_time]
 call pio_decomp(pioSystemState,         & ! input: pio system descriptor
                 ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nSeg,ntdh,nEns],       & ! input: dimension length == global array size
                 ixRch(ix1:ix2),         & ! input:
                 iodesc_irf_bas_double)

 end associate

 ! Finishing up definition -------
 call endDef(pioFileDescState, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 CONTAINS

  SUBROUTINE set_dim_len(ixDim, ierr, message1)
   ! populate state netCDF dimension size
   USE public_var, ONLY: MAXQPAR
   USE globalData, ONLY: meta_stateDims  ! states dimension meta
   USE globalData, ONLY: FRAC_FUTURE     ! To get size of q future for basin IRF
   USE globalData, ONLY: nEns, nRch      ! number of ensembles and river reaches
   implicit none
   ! input
   integer(i4b), intent(in)   :: ixDim    ! ixDim
   ! output
   integer(i4b), intent(out)  :: ierr     ! error code
   character(*), intent(out)  :: message1  ! error message
   ! initialize error control
   ierr=0; message1='set_dim_len/'

   select case(ixDim)
    case(ixStateDims%time);    meta_stateDims(ixStateDims%time)%dimLength    = recordDim
    case(ixStateDims%seg);     meta_stateDims(ixStateDims%seg)%dimLength     = nRch
    case(ixStateDims%ens);     meta_stateDims(ixStateDims%ens)%dimLength     = nEns
    case(ixStateDims%tbound);  meta_stateDims(ixStateDims%tbound)%dimLength  = 2
    case(ixStateDims%tdh);     meta_stateDims(ixStateDims%tdh)%dimLength     = size(FRAC_FUTURE)
    case(ixStateDims%tdh_irf); meta_stateDims(ixStateDims%tdh_irf)%dimLength = 20   !just temporarily
    case(ixStateDims%fdmesh);  meta_stateDims(ixStateDims%fdmesh)%dimLength  = 4
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
   integer(i4b), intent(out) :: ierr          ! error code
   character(*), intent(out) :: message1       ! error message
   ! local
   integer(i4b)              :: iVar          ! index loop for variables
   integer(i4b)              :: dim_IRFbas(4) ! dimensions combination case 4
   integer(i4b)              :: dim_q(3)      ! dimensions combination case 4

   ! initialize error control
   ierr=0; message1='define_IRFbas_state/'

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%tdh)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh)%dimName); return; endif
   end if

   call defdim(pioFileDescState, trim(meta_stateDims(ixStateDims%tdh)%dimName), meta_stateDims(ixStateDims%tdh)%dimLength, meta_stateDims(ixStateDims%tdh)%dimId)
   if(ierr/=0)then; ierr=20; message1=trim(message1)//'cannot define dimension'; return; endif

   associate(dim_seg    => meta_stateDims(ixStateDims%seg)%dimId,     &
             dim_ens    => meta_stateDims(ixStateDims%ens)%dimId,     &
             dim_time   => meta_stateDims(ixStateDims%time)%dimId,    &
             dim_tdh    => meta_stateDims(ixStateDims%tdh)%dimId)

   ! Array dimension sets
   dim_IRFbas(:) = [dim_seg, dim_tdh, dim_ens, dim_time]
   dim_q(:)      = [dim_seg, dim_ens, dim_time]

   do iVar=1,nVarsIRFbas
    select case(iVar)
     case(ixIRFbas%qfuture); call defVar(pioFileDescState, trim(meta_irf_bas(iVar)%varName), dim_IRFbas, ncd_double, ierr, cmessage, vdesc=trim(meta_irf_bas(iVar)%varDesc), vunit=trim(meta_irf_bas(iVar)%varUnit))
     case(ixIRFbas%q);       call defVar(pioFileDescState, trim(meta_irf_bas(iVar)%varName), dim_q,      ncd_double, ierr, cmessage, vdesc=trim(meta_irf_bas(iVar)%varDesc), vunit=trim(meta_irf_bas(iVar)%varUnit))
     case default; ierr=20; message1=trim(message1)//'unable to identify hill-slope routing state variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do

  end associate

  END SUBROUTINE define_IRFbas_state


  SUBROUTINE define_KWE_state(ierr, message1)
   ! External modules
   USE globalData, ONLY: meta_KWE
   USE var_lookup, ONLY: ixKWE, nVarsKWE
   implicit none
   ! output
   integer(i4b), intent(out) :: ierr          ! error code
   character(*), intent(out) :: message1       ! error message
   ! local
   integer(i4b)              :: iVar          ! index loop for variables
   integer(i4b)              :: dim_a(4)      ! dimensions combination case 4
   integer(i4b)              :: dim_q(4)      ! dimensions combination case 4

   ! initialize error control
   ierr=0; message1='define_KWE_state/'

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%fdmesh)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%fdmesh, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%fdmesh)%dimName); return; endif
   end if

   call defdim(pioFileDescState, trim(meta_stateDims(ixStateDims%fdmesh)%dimName), meta_stateDims(ixStateDims%fdmesh)%dimLength, meta_stateDims(ixStateDims%fdmesh)%dimId)
   if(ierr/=0)then; ierr=20; message1=trim(message1)//'cannot define dimension'; return; endif

   associate(dim_seg    => meta_stateDims(ixStateDims%seg)%dimId,     &
             dim_ens    => meta_stateDims(ixStateDims%ens)%dimId,     &
             dim_time   => meta_stateDims(ixStateDims%time)%dimId,    &
             dim_fdmesh => meta_stateDims(ixStateDims%fdmesh)%dimId)

   ! Array dimension sets
   dim_a(:) = [dim_seg, dim_fdmesh, dim_ens, dim_time]
   dim_q(:) = [dim_seg, dim_fdmesh, dim_ens, dim_time]

   do iVar=1,nVarsKWE
    select case(iVar)
     case(ixKWE%a); call defVar(pioFileDescState, trim(meta_KWE(iVar)%varName), dim_a, ncd_double, ierr, cmessage, vdesc=trim(meta_KWE(iVar)%varDesc), vunit=trim(meta_KWE(iVar)%varUnit))
     case(ixKWE%q); call defVar(pioFileDescState, trim(meta_KWE(iVar)%varName), dim_q, ncd_double, ierr, cmessage, vdesc=trim(meta_KWE(iVar)%varDesc), vunit=trim(meta_KWE(iVar)%varUnit))
     case default; ierr=20; message1=trim(message1)//'unable to identify Euler kinematic wave routing state variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do

  end associate

  END SUBROUTINE define_KWE_state


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
   integer(i4b)               :: dim_kwt(4)  ! dimensions combination case 4

   ! initialize error control
   ierr=0; message1='define_KWT_state/'

   ! Check dimension length is populated
   if (meta_stateDims(ixStateDims%wave)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%wave, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%wave)%dimName); return; endif
   end if

   ! Define dimension needed for this routing specific state variables
   call defdim(pioFileDescState, trim(meta_stateDims(ixStateDims%wave)%dimName), meta_stateDims(ixStateDims%wave)%dimLength, meta_stateDims(ixStateDims%wave)%dimId)
   if(ierr/=0)then; ierr=20; message1=trim(message1)//'cannot define dimension'; return; endif

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimId,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimId,     &
             dim_time    => meta_stateDims(ixStateDims%time)%dimId,    &
             dim_wave    => meta_stateDims(ixStateDims%wave)%dimId)

   dim_kwt(:) = [dim_seg, dim_wave, dim_ens, dim_time]

   call defVar(pioFileDescState, 'numWaves', [dim_seg,dim_ens,dim_time], ncd_int, ierr, cmessage, vdesc='number of waves in a reach', vunit='-')
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsKWT
    if (iVar==ixKWT%q) cycle  !skip writing final discharge
    select case(iVar)
     case(ixKWT%tentry);    call defVar(pioFileDescState, trim(meta_kwt(iVar)%varName), dim_kwt, ncd_double, ierr, cmessage, vdesc=trim(meta_kwt(iVar)%varDesc), vunit=trim(meta_kwt(iVar)%varUnit))
     case(ixKWT%texit);     call defVar(pioFileDescState, trim(meta_kwt(iVar)%varName), dim_kwt, ncd_double, ierr, cmessage, vdesc=trim(meta_kwt(iVar)%varDesc), vunit=trim(meta_kwt(iVar)%varUnit))
     case(ixKWT%qwave);     call defVar(pioFileDescState, trim(meta_kwt(iVar)%varName), dim_kwt, ncd_double, ierr, cmessage, vdesc=trim(meta_kwt(iVar)%varDesc), vunit=trim(meta_kwt(iVar)%varUnit))
     case(ixKWT%qwave_mod); call defVar(pioFileDescState, trim(meta_kwt(iVar)%varName), dim_kwt, ncd_double, ierr, cmessage, vdesc=trim(meta_kwt(iVar)%varDesc), vunit=trim(meta_kwt(iVar)%varUnit))
     case(ixKWT%routed);    call defVar(pioFileDescState, trim(meta_kwt(iVar)%varName), dim_kwt, ncd_int,   ierr, cmessage, vdesc=trim(meta_kwt(iVar)%varDesc), vunit=trim(meta_kwt(iVar)%varUnit))
     case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
    end select
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
   integer(i4b), intent(out)  :: ierr        ! error code
   character(*), intent(out)  :: message1    ! error message
   ! local
   integer(i4b)               :: iVar        ! index loop for variables
   integer(i4b)               :: dim_irf(4)  ! dimensions combination case 4
   ! initialize error control
   ierr=0; message1='define_IRF_state/'

   if (meta_stateDims(ixStateDims%tdh_irf)%dimLength == integerMissing) then
     call set_dim_len(ixStateDims%tdh_irf, ierr, cmessage)
     if(ierr/=0)then; message1=trim(message1)//trim(cmessage)//' for '//trim(meta_stateDims(ixStateDims%tdh_irf)%dimName); return; endif
   endif

   ! Define dimension needed for this routing specific state variables
   call defdim(pioFileDescState, trim(meta_stateDims(ixStateDims%tdh_irf)%dimName), meta_stateDims(ixStateDims%tdh_irf)%dimLength, meta_stateDims(ixStateDims%tdh_irf)%dimId)
   if(ierr/=0)then; ierr=20; message1=trim(message1)//'cannot define dimension'; return; endif

   associate(dim_seg     => meta_stateDims(ixStateDims%seg)%dimId,     &
             dim_ens     => meta_stateDims(ixStateDims%ens)%dimId,     &
             dim_time    => meta_stateDims(ixStateDims%time)%dimId,    &
             dim_tdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimId)

   dim_irf(:) = [dim_seg, dim_tdh_irf, dim_ens, dim_time]

   call defVar(pioFileDescState, 'numQF', [dim_seg,dim_ens,dim_time], ncd_int, ierr, cmessage, vdesc='number of future q time steps in a reach', vunit='-')
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

   do iVar=1,nVarsIRF
    if (iVar==ixIRF%q) cycle  ! skip writing final discharge
    select case(iVar)
     case(ixIRF%qfuture); call defVar(pioFileDescState, trim(meta_irf(iVar)%varName), dim_irf, ncd_double, ierr, cmessage, vdesc=trim(meta_irf(iVar)%varDesc), vunit=trim(meta_irf(iVar)%varUnit))
     case default; ierr=20; message1=trim(message1)//'unable to identify IRF routing state variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end do

   end associate

  END SUBROUTINE define_IRF_state

 END SUBROUTINE define_state_nc


 ! *********************************************************************
 ! public subroutine: writing routing state NetCDF file
 ! *********************************************************************
 SUBROUTINE write_state_nc(fname,                &   ! Input: state netcdf name
                           ierr, message)            ! Output: error control
 ! External module
 ! meta data
 USE globalData, ONLY: meta_stateDims  ! dimension for state variables
 ! Named variables
 USE var_lookup, ONLY: ixStateDims, nStateDims
 ! global/public data
 USE public_var, ONLY: routOpt
 USE public_var, ONLY: doesBasinRoute
 USE public_var, ONLY: allRoutingMethods
 USE public_var, ONLY: kinematicWave
 USE public_var, ONLY: kinematicWaveEuler
 USE public_var, ONLY: impulseResponseFunc
 USE globalData, ONLY: RCHFLX_main         ! mainstem reach fluxes (ensembles, reaches)
 USE globalData, ONLY: RCHFLX_trib         ! tributary reach fluxes (ensembles, reaches)
 USE globalData, ONLY: NETOPO_main         ! mainstem reach topology
 USE globalData, ONLY: NETOPO_trib         ! tributary reach topology
 USE globalData, ONLY: RCHSTA_main         ! mainstem reach state (ensembles, reaches)
 USE globalData, ONLY: RCHSTA_trib         ! tributary reach state (ensembles, reaches)
 USE globalData, ONLY: rch_per_proc        ! number of reaches assigned to each proc (size = num of procs+1)
 USE globalData, ONLY: nRch_mainstem       ! number of mainstem reaches
 USE globalData, ONLY: reachID         ! reach ID in network
 USE globalData, ONLY: nRch            ! number of reaches in network
 USE globalData, ONLY: TSEC            ! beginning/ending of simulation time step [sec]
 USE globalData, ONLY: timeVar         ! time variable
 USE globalData, ONLY: iTime           ! time index
 implicit none
 ! input variables
 character(*), intent(in)        :: fname           ! filename
 ! output variables
 integer(i4b), intent(out)       :: ierr            ! error code
 character(*), intent(out)       :: message         ! error message
 ! local variables
 integer(i4b)                    :: iens            ! temporal
 type(states)                    :: state(0:3)      ! temporal state data structures -currently 2 river routing scheme + basin IRF routing
 type(STRFLX), allocatable       :: RCHFLX_local(:) ! reordered reach flux data structure
 type(RCHTOPO),allocatable       :: NETOPO_local(:) ! reordered topology data structure
 type(STRSTA), allocatable       :: RCHSTA_local(:) ! reordered statedata structure
 character(len=strLen)           :: cmessage        ! error message of downwind routine

 ! initialize error control
 ierr=0; message='write_state_nc/'

 iens = 1
 kTime = kTime + 1

 if (masterproc) then
  associate(nRch_trib => rch_per_proc(0))
  allocate(RCHFLX_local(nRch_mainstem+nRch_trib), &
           NETOPO_local(nRch_mainstem+nRch_trib), &
           RCHSTA_local(nRch_mainstem+nRch_trib), stat=ierr)
  if (nRch_mainstem>0) then
    RCHFLX_local(1:nRch_mainstem) = RCHFLX_main(iens,1:nRch_mainstem)
    NETOPO_local(1:nRch_mainstem) = NETOPO_main(1:nRch_mainstem)
    RCHSTA_local(1:nRch_mainstem) = RCHSTA_main(iens,1:nRch_mainstem)
  end if
   if (nRch_trib>0) then
     RCHFLX_local(nRch_mainstem+1:nRch_mainstem+nRch_trib) = RCHFLX_trib(iens,:)
     NETOPO_local(nRch_mainstem+1:nRch_mainstem+nRch_trib) = NETOPO_trib(:)
     RCHSTA_local(nRch_mainstem+1:nRch_mainstem+nRch_trib) = RCHSTA_trib(iens,:)
   endif
   end associate
 else
  allocate(RCHFLX_local(rch_per_proc(pid)), &
           NETOPO_local(rch_per_proc(pid)), &
           RCHSTA_local(rch_per_proc(pid)),stat=ierr)
  RCHFLX_local = RCHFLX_trib(iens,:)
  NETOPO_local = NETOPO_trib(:)
  RCHSTA_local = RCHSTA_trib(iens,:)
 endif

 ! -- Write out to netCDF

 call openFile(pioSystemState, pioFileDescState, trim(fname),pio_typename, ncd_write, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! Miscellaneous variables - seg id, time etc
 call write_netcdf(pioFileDescState, 'reachID', reachID, [1], [nRch], ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_netcdf(pioFileDescState, 'time', [timeVar(iTime)], [kTime], [1], ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call write_netcdf(pioFileDescState, 'time_bound', [TSEC(0),TSEC(1)], [1,kTime], [2,1], ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 if (doesBasinRoute == 1) then
  call write_IRFbas_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
  call write_IRF_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (routOpt==kinematicWaveEuler) then
  call write_KWE_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
  call write_KWT_state(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 call closeFile(pioFileDescState)

 CONTAINS

  ! Basin IRF writing procedures
  SUBROUTINE write_IRFbas_state(ierr, message1)
  ! external module
  USE globalData, ONLY: meta_irf_bas
  USE var_lookup, ONLY: ixIRFbas, nVarsIRFbas
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively

  ! initialize error control
  ierr=0; message1='write_IRFbas_state/'

  associate(nSeg     => size(RCHFLX_local),                         &
            nEns     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            ntdh     => meta_stateDims(ixStateDims%tdh)%dimLength)      ! maximum future q time steps among basins

  allocate(state(0)%var(nVarsIRFbas), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRFbas

   select case(iVar)
    case(ixIRFbas%qfuture); allocate(state(0)%var(iVar)%array_3d_dp(nSeg, ntdh, nEns), stat=ierr)
    case(ixIRFbas%q);       allocate(state(0)%var(iVar)%array_2d_dp(nSeg, nEns),       stat=ierr)
    case default; ierr=20; message1=trim(message1)//'unable to identify basin routing variable index'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//'problem allocating space for basin IRF routing state '//trim(meta_irf_bas(iVar)%varName); return; endif

  end do

 ! --Convert data structures to arrays
  do iens=1,nEns
   do iSeg=1,nSeg
     do iVar=1,nVarsIRFbas
      select case(iVar)
       case(ixIRFbas%qfuture); state(0)%var(iVar)%array_3d_dp(iSeg,:,iens) = RCHFLX_local(iSeg)%QFUTURE
       case(ixIRFbas%q);       state(0)%var(iVar)%array_2d_dp(iSeg,iens)   = RCHFLX_local(iSeg)%BASIN_QR(1)
       case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF state variable index'; return
      end select
    enddo
   enddo
  enddo

  do iVar=1,nVarsIRFbas

   select case(iVar)
    case(ixIRFbas%qfuture); call write_pnetcdf_recdim(pioFileDescState, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_3d_dp, iodesc_irf_bas_double, kTime, ierr, cmessage)
    case(ixIRFbas%q);       call write_pnetcdf_recdim(pioFileDescState, meta_irf_bas(iVar)%varName, state(0)%var(iVar)%array_2d_dp, iodesc_state_double, kTime, ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify basin IRF variable index for nc writing'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  enddo

  end associate

  END SUBROUTINE write_IRFbas_state

  ! KWE writing procedures
  SUBROUTINE write_KWE_state(ierr, message1)
  ! External module
  USE globalData, ONLY: meta_kwe
  USE var_lookup, ONLY: ixKWE, nVarsKWE
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively
  ! initialize error control
  ierr=0; message1='write_KWE_state/'

  associate(nSeg     => size(RCHFLX_local),                         &
            nEns     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            nMesh    => meta_stateDims(ixStateDims%fdmesh)%dimLength)     ! maximum waves allowed in a reach

  allocate(state(kinematicWaveEuler)%var(nVarsKWE), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWE
    select case(iVar)
     case(ixKWE%a, ixKWE%q)
      allocate(state(kinematicWaveEuler)%var(iVar)%array_3d_dp(nSeg, nMesh, nEns), stat=ierr)
     case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWE routing state '//trim(meta_kwe(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nEns
   do iSeg=1,nSeg

    do iVar=1,nVarsKWE
     select case(iVar)
      case(ixKWE%a)
       state(kinematicWaveEuler)%var(iVar)%array_3d_dp(iSeg,1:4,iens) = RCHSTA_local(iSeg)%EKW_ROUTE%A(:)
      case(ixKWE%q)
       state(kinematicWaveEuler)%var(iVar)%array_3d_dp(iSeg,1:4,iens) = RCHSTA_local(iSeg)%EKW_ROUTE%Q(:)
      case default; ierr=20; message1=trim(message1)//'unable to identify KWE routing state variable index'; return
     end select
    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! Writing netCDF
  do iVar=1,nVarsKWE
    select case(iVar)
     case(ixKWE%a, ixKWE%q)
       call write_pnetcdf_recdim(pioFileDescState, trim(meta_kwe(iVar)%varName), state(kinematicWaveEuler)%var(iVar)%array_3d_dp, iodesc_mesh_double, kTime, ierr, cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify KWE variable index for nc writing'; return
    end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  end do

  end associate

  END SUBROUTINE write_KWE_state


  ! KWT writing procedures
  SUBROUTINE write_KWT_state(ierr, message1)
  ! External module
  USE globalData, ONLY: meta_kwt
  USE var_lookup, ONLY: ixKWT, nVarsKWT
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

  associate(nSeg     => size(RCHFLX_local),                         &
            nEns     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            nWave    => meta_stateDims(ixStateDims%wave)%dimLength)     ! maximum waves allowed in a reach

  allocate(state(kinematicWave)%var(nVarsKWT), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! array to store number of wave per segment and ensemble
  allocate(numWaves(nEns,nSeg), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT
    if (iVar==ixKWT%q) cycle  ! not writing out river flow in state file
    select case(iVar)
     case(ixKWT%routed); allocate(state(kinematicWave)%var(iVar)%array_3d_int(nSeg, nWave, nEns), stat=ierr)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
      allocate(state(kinematicWave)%var(iVar)%array_3d_dp(nSeg, nWave, nEns), stat=ierr)
     case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
    end select
    if(ierr/=0)then; message1=trim(message1)//'problem allocating space for KWT routing state '//trim(meta_kwt(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nEns
   do iSeg=1,nSeg

    numWaves(iens,iseg) = size(RCHSTA_local(iseg)%LKW_ROUTE%KWAVE)

    do iVar=1,nVarsKWT

     if (iVar==ixKWT%q) cycle ! not writing out KWT routed flow

     select case(iVar)
      case(ixKWT%tentry)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%TI
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%texit)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%TR
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%qwave)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%QF
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%qwave_mod)
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,1:numWaves(iens,iSeg),iens) = RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%QM
       state(kinematicWave)%var(iVar)%array_3d_dp(iSeg,numWaves(iens,iSeg)+1:,iens) = realMissing
      case(ixKWT%routed) ! this is suppposed to be logical variable, but put it as 0 or 1 in double now
       if (allocated(RFvec)) deallocate(RFvec, stat=ierr)
       allocate(RFvec(numWaves(iens,iSeg)),stat=ierr); RFvec=0_i4b
       where (RCHSTA_local(iSeg)%LKW_ROUTE%KWAVE(:)%RF) RFvec=1_i4b
       state(kinematicWave)%var(iVar)%array_3d_int(iSeg,1:numWaves(iens,iSeg),iens) = RFvec
       state(kinematicWave)%var(iVar)%array_3d_int(iSeg,numWaves(iens,iSeg)+1:,iens) = integerMissing
      case default; ierr=20; message1=trim(message1)//'unable to identify KWT routing state variable index'; return
     end select
    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! Writing netCDF
  call write_pnetcdf_recdim(pioFileDescState,'numWaves', numWaves, iodesc_state_int, kTime, ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsKWT

    if (iVar==ixKWT%q) cycle  ! not writing out river flow in state file

    select case(iVar)
     case(ixKWT%routed)
       call write_pnetcdf_recdim(pioFileDescState, trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_int, iodesc_wave_int, kTime, ierr, cmessage)
     case(ixKWT%tentry, ixKWT%texit, ixKWT%qwave, ixKWT%qwave_mod)
       call write_pnetcdf_recdim(pioFileDescState, trim(meta_kwt(iVar)%varName), state(kinematicWave)%var(iVar)%array_3d_dp, iodesc_wave_double, kTime, ierr, cmessage)
     case default; ierr=20; message1=trim(message1)//'unable to identify KWT variable index for nc writing'; return
    end select
   if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  end do

  end associate

  END SUBROUTINE write_KWT_state


  ! IRF writing procedures
  SUBROUTINE write_IRF_state(ierr, message1)
  ! external module
  USE globalData, ONLY: meta_irf        ! IRF routing
  USE var_lookup, ONLY: ixIRF, nVarsIRF
  implicit none
  ! output
  integer(i4b), intent(out)  :: ierr            ! error code
  character(*), intent(out)  :: message1        ! error message
  ! local variables
  integer(i4b)               :: iVar,iens,iSeg  ! index loops for variables, ensembles and segments respectively
  integer(i4b), allocatable  :: numQF(:,:)      ! number of future Q time steps for each ensemble and segment
  ! initialize error control
  ierr=0; message1='write_IRF_state/'

  associate(nSeg     => size(RCHFLX_local),                         &
            nEns     => meta_stateDims(ixStateDims%ens)%dimLength,  &
            ntdh_irf => meta_stateDims(ixStateDims%tdh_irf)%dimLength)    ! maximum future q time steps among reaches

  allocate(state(impulseResponseFunc)%var(nVarsIRF), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  ! array to store number of wave per segment and ensemble
  allocate(numQF(nEns,nSeg), stat=ierr)
  if(ierr/=0)then; message1=trim(message1)//'problem allocating space for numQF'; return; endif

  do iVar=1,nVarsIRF
   if (iVar==ixIRF%q) cycle
   select case(iVar)
    case(ixIRF%qfuture); allocate(state(impulseResponseFunc)%var(iVar)%array_3d_dp(nSeg, ntdh_irf, nEns), stat=ierr)
    case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
   end select
   if(ierr/=0)then; message1=trim(message1)//'problem allocating space for IRF routing state '//trim(meta_irf(iVar)%varName); return; endif
  end do

  ! --Convert data structures to arrays
  do iens=1,nEns
   do iSeg=1,nSeg

    numQF(iens,iseg) = size(NETOPO_local(iSeg)%UH)

    do iVar=1,nVarsIRF

     if (iVar==ixIRF%q) cycle ! not writing out IRF routed flow

     select case(iVar)
      case(ixIRF%qfuture)
       state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,1:numQF(iens,iSeg),iens) = RCHFLX_local(iSeg)%QFUTURE_IRF
       state(impulseResponseFunc)%var(iVar)%array_3d_dp(iSeg,numQF(iens,iSeg)+1:ntdh_irf,iens) = realMissing
      case default; ierr=20; message1=trim(message1)//'unable to identify variable index'; return
     end select

    enddo ! variable loop
   enddo ! seg loop
  enddo ! ensemble loop

  ! writing netcdf
  call write_pnetcdf_recdim(pioFileDescState, 'numQF', numQF, iodesc_state_int, kTime, ierr, cmessage)
  if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif

  do iVar=1,nVarsIRF

   if (iVar==ixIRF%q) cycle  ! not writing out river flow in state file

   select case(iVar)
    case(ixIRF%qfuture)
     call write_pnetcdf_recdim(pioFileDescState, trim(meta_irf(iVar)%varName), state(impulseResponseFunc)%var(iVar)%array_3d_dp, iodesc_irf_float, kTime, ierr, cmessage)
    case default; ierr=20; message1=trim(message1)//'unable to identify IRF variable index for nc writing'; return
    if(ierr/=0)then; message1=trim(message1)//trim(cmessage); return; endif
   end select

  end do

  end associate

  END SUBROUTINE write_IRF_state

 END SUBROUTINE write_state_nc

END MODULE write_restart_pio
