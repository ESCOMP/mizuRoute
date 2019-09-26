MODULE write_simoutput_pio

! Moudle wide external modules
USE nrtype
USE dataTypes,         ONLY: STRFLX            ! fluxes in each reach
USE public_var,        ONLY: root
USE public_var,        ONLY: integerMissing
USE globalData,        ONLY: pid, nNodes
USE nr_utility_module, ONLY: arth
USE pio_utils

implicit none

! The following variables used only in this module
character(len=strLen),save :: fileout              ! name of the output file
integer(i4b),         save :: jTime                ! time step in output netCDF
type(iosystem_desc_t),save :: pioSystem            ! PIO I/O system data
type(file_desc_t),    save :: pioFileDesc          ! PIO data identifying the file
type(io_desc_t),      save :: iodesc_rch_flx       ! PIO domain decomposition data for reach flux [nRch]
type(io_desc_t),      save :: iodesc_hru_ro        ! PIO domain decomposition data for hru runoff [nHRU]

integer(i4b),parameter     :: recordDim=-999       ! record dimension Indicator

private

public::prep_output
public::output
public::close_output_nc

contains

 ! *********************************************************************
 ! public subroutine: define routing output NetCDF file
 ! *********************************************************************
 subroutine output(ierr, message)

  !Dependent modules
  USE public_var,          only : doesBasinRoute      ! basin routing options   0-> no, 1->IRF, otherwise error
  USE public_var,          only : doesAccumRunoff     ! option to delayed runoff accumulation over all the upstream reaches. 0->no, 1->yes
  USE public_var,          only : routOpt             ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error
  USE public_var,          only : kinematicWave       ! kinematic wave
  USE public_var,          only : impulseResponseFunc ! impulse response function
  USE public_var,          only : allRoutingMethods   ! all routing methods
  USE globalData,          only : nHRU                ! number of ensembles, HRUs and river reaches
  USE globalData,          only : RCHFLX              ! global Reach fluxes (ensembles, space [reaches])
  USE globalData,          only : RCHFLX_trib         ! tributary Reach fluxes (ensembles, space [reaches])
  USE globalData,          only : iTime               ! time index at simulation time step
  USE globalData,          only : timeVar             ! time variables (unit given by runoff data)
  USE globalData,          only : runoff_data         ! runoff data for one time step for LSM HRUs and River network HRUs
  USE globalData,          only : ixRch_order         ! global reach index in the order of proc assignment (size = total number of reaches in the entire network)
  USE globalData,          only : rch_per_proc        ! number of reaches assigned to each proc (size = num of procs+1)

  implicit none

  ! input variables: none
  ! output variables
  integer(i4b), intent(out)       :: ierr             ! error code
  character(*), intent(out)       :: message          ! error message
  ! local variables
  integer(i4b)                    :: iens             ! temporal
  character(len=strLen)           :: cmessage         ! error message of downwind routine
  type(STRFLX),allocatable        :: RCHFLX_local(:)
  real(sp),    allocatable        :: tmp_array(:)
  real(sp),    allocatable        :: basinRunoff(:)
  integer(i4b)                    :: ix               ! error code

  ! initialize error control
  ierr=0; message='output/'

  iens = 1

  ! Need to combine mainstem RCHFLX and tributary RCHFLX into RCHFLX_local for root node
  if (pid==root) then
   associate(nRch_main => rch_per_proc(-1), nRch_trib => rch_per_proc(0))
   allocate(RCHFLX_local(nRch_main+nRch_trib), tmp_array(nRch_main+nRch_trib), stat=ierr)
   if (nRch_main/=0) then
     do ix = 1,nRch_main
      RCHFLX_local(ix) = RCHFLX(iens,ixRch_order(ix))
     enddo
   end if
   RCHFLX_local(nRch_main+1:nRch_main+nRch_trib) = RCHFLX_trib(iens,:)
   end associate
  else
   allocate(RCHFLX_local(rch_per_proc(pid)),tmp_array(rch_per_proc(pid)), stat=ierr)
   RCHFLX_local = RCHFLX_trib(iens,:)
  endif

  if (pid==root) then
   allocate(basinRunoff(nHRU))
   basinRunoff = real(runoff_data%basinRunoff, kind=sp)
  else
   allocate(basinRunoff(1))
  endif

  ! write time -- note time is just carried across from the input
  call write_netcdf(pioFileDesc, 'time', [timeVar(iTime)], [jTime], [1], ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write the basin runoff to the netcdf file
  call write_pnetcdf_recdim(pioFileDesc, 'basRunoff',basinRunoff, iodesc_hru_ro, jTime, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (doesBasinRoute == 1) then
   ! write instataneous local runoff in each stream segment (m3/s)
   tmp_array = real(RCHFLX_local(:)%BASIN_QI,kind=sp)
   call write_pnetcdf_recdim(pioFileDesc, 'instRunoff', tmp_array, iodesc_rch_flx, jTime, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! write routed local runoff in each stream segment (m3/s)
  tmp_array = real(RCHFLX_local(:)%BASIN_QR(1),kind=sp)
  call write_pnetcdf_recdim(pioFileDesc, 'dlayRunoff', tmp_array, iodesc_rch_flx, jTime, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write accumulated runoff (m3/s)
  if (doesAccumRunoff == 1) then
   tmp_array = real(RCHFLX_local(:)%UPSTREAM_QI,kind=sp)
   call write_pnetcdf_recdim(pioFileDesc, 'sumUpstreamRunoff', tmp_array, iodesc_rch_flx, jTime, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
   ! write routed runoff (m3/s)
   tmp_array = real(RCHFLX_local(:)%REACH_Q,kind=sp)
   call write_pnetcdf_recdim(pioFileDesc, 'KWTroutedRunoff', tmp_array, iodesc_rch_flx, jTime, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
   ! write routed runoff (m3/s)
   tmp_array = real(RCHFLX_local(:)%REACH_Q_IRF,kind=sp)
   call write_pnetcdf_recdim(pioFileDesc, 'IRFroutedRunoff', tmp_array, iodesc_rch_flx, jTime, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

 end subroutine output


 ! *********************************************************************
 ! public subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE prep_output(ierr, message)

 ! saved public variables (usually parameters, or values not modified)
 USE public_var,          only : output_dir        ! output directory
 USE public_var,          only : fname_output      ! output file name head
 USE public_var,          only : calendar          ! calendar name
 USE public_var,          only : newFileFrequency  ! frequency for new output files (day, month, annual, single)
 USE public_var,          only : time_units        ! time units (seconds, hours, or days)
 ! saved global data
 USE globalData,          only : basinID,reachID   ! HRU and reach ID in network
 USE globalData,          only : modJulday         ! julian day: at model time step
 USE globalData,          only : modTime           ! previous and current model time
 USE globalData,          only : nHRU, nRch        ! number of ensembles, HRUs and river reaches
 USE globalData,          only : isFileOpen        ! file open/close status
 ! subroutines
 USE time_utils_module,   only : compCalday        ! compute calendar day
 USE time_utils_module,   only : compCalday_noleap ! compute calendar day

 implicit none

 ! input variables: none
 ! output variables
 integer(i4b), intent(out)       :: ierr             ! error code
 character(*), intent(out)       :: message          ! error message
 ! local variables
 logical(lgt)                    :: defNewOutputFile ! flag to define new output file
 character(len=strLen)           :: cmessage         ! error message of downwind routine

 ! initialize error control
 ierr=0; message='prep_output/'

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
  if (pid==root) then
    print*, modTime(1)%iy,modTime(1)%im,modTime(1)%id,modTime(1)%ih,modTime(1)%imin
  endif

  ! check need for the new file
  select case(newFileFrequency)
   case('single'); defNewOutputFile=(modTime(0)%iy==integerMissing)
   case('annual'); defNewOutputFile=(modTime(1)%iy/=modTime(0)%iy)
   case('month');  defNewOutputFile=(modTime(1)%im/=modTime(0)%im)
   case('day');    defNewOutputFile=(modTime(1)%id/=modTime(0)%id)
   case default; ierr=20; message=trim(message)//'unable to identify the option to define new output files'; return
  end select

  ! define new file
  if(defNewOutputFile)then

   ! close netcdf only if is is open
   call close_output_nc()

   ! initialize time
   jTime=1

   ! Define filename
   write(fileout,'(a,3(i0,a))') trim(output_dir)//trim(fname_output)//'_', modTime(1)%iy, '-', modTime(1)%im, '-', modTime(1)%id, '.nc'

   ! define output file
   call defineFile(trim(fileout),                         &  ! input: file name
                   time_units,                            &  ! input: time units
                   calendar,                              &  ! input: calendar
                   ierr,cmessage)                            ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call openFile(pioSystem, pioFileDesc, trim(fileout), ncd_write, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! define basin ID
   call write_netcdf(pioFileDesc, 'basinID', basinID, [1], [nHRU], ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! define reach ID
   call write_netcdf(pioFileDesc, 'reachID', reachID, [1], [nRch], ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   isFileOpen = .True.

  ! no new file requested: increment time
  else

   jTime = jTime+1

  endif

  modTime(0) = modTime(1)

 END SUBROUTINE prep_output

 SUBROUTINE close_output_nc()
  USE globalData, only : isFileOpen   ! file open/close status
  implicit none
  if (isFileOpen) then
   call closeFile(pioFileDesc)
   isFileOpen=.false.
  endif
 END SUBROUTINE close_output_nc

 ! *********************************************************************
 ! private subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE defineFile(fname,           &  ! input: filename
                       units_time,      &  ! input: time units
                       calendar,        &  ! input: calendar
                       ierr, message)      ! output: error control
 !Dependent modules
 USE var_lookup, ONLY: ixQdims, nQdims
 USE globalData, ONLY: meta_qDims
 USE globalData, ONLY: rch_per_proc             ! number of reaches assigned to each proc (size = num of procs+1)
 USE globalData, ONLY: nEns, nHRU, nRch         ! number of ensembles, HRUs and river reaches

 implicit none
 ! input variables
 character(*), intent(in)    :: fname             ! filename
 character(*), intent(in)    :: units_time        ! time units
 character(*), intent(in)    :: calendar          ! calendar
 ! output variables
 integer(i4b), intent(out)   :: ierr              ! error code
 character(*), intent(out)   :: message           ! error message
 ! local variables
 character(len=strLen)       :: cmessage          ! error message of downwind routine
 integer(i4b)                :: jDim,iVar         ! dimension, and variable index
 integer(i4b)                :: ix1, ix2          ! frst and last indices of global array for local array chunk
 integer(i4b)                :: ixRch(nRch)        !
 integer(i4b)                :: nHRU_in
 integer(i4b),allocatable    :: dof_hru(:)        ! dof for basin runoff
 integer(i4b),parameter      :: nVars=8           ! number of variables

 ! initialize error control
 ierr=0; message='defineFile/'

 associate (dim_seg  => meta_qDims(ixQdims%seg)%dimId,    &
            dim_hru  => meta_qDims(ixQdims%hru)%dimId,    &
            dim_ens  => meta_qDims(ixQdims%ens)%dimId,    &
            dim_time => meta_qDims(ixQdims%time)%dimId)

! populate q dimension meta (not sure if this should be done here...)
 meta_qDims(ixQdims%seg)%dimLength = nRch
 meta_qDims(ixQdims%hru)%dimLength = nHRU
 meta_qDims(ixQdims%ens)%dimLength = nEns

 ! pio initialization
 call pio_sys_init(pid, nNodes, pioSystem)

 if (pid==root) then
   ix1 = 1_i4b
 else
   ix1 = sum(rch_per_proc(-1:pid-1))+1_i4b
 endif
 ix2 = sum(rch_per_proc(-1:pid))
 ixRch = arth(1,1,nRch)
 call pio_decomp(pioSystem,              & ! input: pio system descriptor
                 ncd_float,              & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nRch],                 & ! input: dimension length == global array size
                 ixRch(ix1:ix2),         & ! input:
                 iodesc_rch_flx)

! For runoff
 if (pid/=root) then
  nHRU_in = 1_i4b
 else
  nHRU_in = nHRU
 endif

 allocate(dof_hru(nHRU_in))

 if (pid==root) then
  dof_hru = arth(1,1,nHRU)
 else
  dof_hru = 0_i4b
 endif

 call pio_decomp(pioSystem,     & ! input: pio system descriptor
                 ncd_float,     & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nHRU],        & ! input: dimension length == global array size
                 dof_hru,       & ! input:
                 iodesc_hru_ro)

 call createFile(pioSystem, trim(fname), pioFileDesc, ierr, cmessage)
 if(ierr/=0)then; message=trim(cmessage)//'cannot create netCDF'; return; endif

 do jDim =1,nQdims
   if (jDim ==ixQdims%time) then ! time dimension (unlimited)
    call defdim(pioFileDesc, trim(meta_qDims(jDim)%dimName), recordDim, meta_qDims(jDim)%dimId)
   else
    call defdim(pioFileDesc, trim(meta_qDims(jDim)%dimName), meta_qDims(jDim)%dimLength, meta_qDims(jDim)%dimId)
   endif
  if(ierr/=0)then; message=trim(message)//'cannot define dimension'; return; endif
 end do

 ! define coordinate variable for time
 call defVar(pioFileDesc, trim(meta_qDims(ixQdims%time)%dimName), [dim_time], ncd_float, ierr, cmessage, vdesc=trim(meta_qDims(ixQdims%time)%dimName), vunit=trim(units_time), vcal=calendar)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! define variables
 do iVar=1,nVars
  ! define variable
  select case(iVar)
   ! define network topology (integers)
   case( 1); call defvar(pioFileDesc, 'basinID',           [dim_hru],          ncd_int,   ierr, cmessage, vdesc='basin ID',                            vunit='-'   )
   case( 2); call defvar(pioFileDesc, 'reachID',           [dim_seg],          ncd_int,   ierr, cmessage, vdesc='reach ID',                            vunit='-'   )
   ! define runoff variables (single precision)
   case( 3); call defvar(pioFileDesc, 'basRunoff',         [dim_hru,dim_time], ncd_float, ierr, cmessage, vdesc='basin runoff',                         vunit='m/s' )
   case( 4); call defVar(pioFileDesc, 'instRunoff',        [dim_seg,dim_time], ncd_float, ierr, cmessage, vdesc='instantaneous runoff in each reach',   vunit='m3/s')
   case( 5); call defVar(pioFileDesc, 'dlayRunoff',        [dim_seg,dim_time], ncd_float, ierr, cmessage, vdesc='delayed runoff in each reach',         vunit='m3/s')
   case( 6); call defVar(pioFileDesc, 'sumUpstreamRunoff', [dim_seg,dim_time], ncd_float, ierr, cmessage, vdesc='sum of upstream runoff in each reach', vunit='m3/s')
   case( 7); call defVar(pioFileDesc, 'KWTroutedRunoff',   [dim_seg,dim_time], ncd_float, ierr, cmessage, vdesc='KWT routed runoff in each reach',      vunit='m3/s')
   case( 8); call defVar(pioFileDesc, 'IRFroutedRunoff',   [dim_seg,dim_time], ncd_float, ierr, cmessage, vdesc='IRF routed runoff in each reach',      vunit='m3/s')
   case default; ierr=20; message=trim(message)//'unable to identify variable index'; return
  end select
  ! check errors
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do

 end associate

 ! end definitions
 call endDef(pioFileDesc, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE defineFile


END MODULE write_simoutput_pio
