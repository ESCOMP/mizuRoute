MODULE write_simoutput_pio

! Moudle wide shared data
USE nrtype
USE var_lookup,        ONLY: ixRFLX, nVarsRFLX
USE dataTypes,         ONLY: STRFLX            ! fluxes in each reach
USE public_var,        ONLY: iulog             ! i/o logical unit number
USE public_var,        ONLY: integerMissing
USE public_var,        ONLY: doesBasinRoute      ! basin routing options   0-> no, 1->IRF, otherwise error
USE public_var,        ONLY: doesAccumRunoff     ! option to delayed runoff accumulation over all the upstream reaches. 0->no, 1->yes
USE public_var,        ONLY: routOpt             ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error
USE public_var,        ONLY: kinematicWave       ! Lagrangian kinematic wave
USE public_var,        ONLY: kinematicWaveEuler  ! Eulerian kinematic wave
USE public_var,        ONLY: impulseResponseFunc ! impulse response function
USE public_var,        ONLY: allRoutingMethods   ! all routing methods
USE globalData,        ONLY: meta_rflx
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
character(300),       save :: fileout              ! name of the output file
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

CONTAINS

 ! *********************************************************************
 ! public subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE output(ierr, message)

  ! global data required only this routine
  USE globalData, ONLY: RCHFLX_main         ! mainstem Reach fluxes (ensembles, space [reaches])
  USE globalData, ONLY: RCHFLX_trib         ! tributary Reach fluxes (ensembles, space [reaches])
  USE globalData, ONLY: basinRunoff_main    ! mainstem only HRU runoff
  USE globalData, ONLY: basinRunoff_trib    ! tributary only HRU runoff
  USE globalData, ONLY: nHRU_mainstem       ! number of mainstem HRUs
  USE globalData, ONLY: iTime               ! time index at simulation time step
  USE globalData, ONLY: timeVar             ! time variables (unit given by runoff data)
  USE globalData, ONLY: nRch_mainstem       ! number of mainstem reaches
  USE globalData, ONLY: rch_per_proc        ! number of reaches assigned to each proc (size = num of procs+1)
  USE globalData, ONLY: hru_per_proc        ! number of hrus assigned to each proc (size = num of procs+1)

  implicit none

  ! input variables: none
  ! output variables
  integer(i4b), intent(out)       :: ierr             ! error code
  character(*), intent(out)       :: message          ! error message
  ! local variables
  integer(i4b)                    :: iens             ! temporal
  character(len=strLen)           :: cmessage         ! error message of downwind routine
  type(STRFLX),allocatable        :: RCHFLX_local(:)
  real(dp),    allocatable        :: basinRunoff(:)

  ! initialize error control
  ierr=0; message='output/'

  iens = 1

  ! Need to combine mainstem RCHFLX and tributary RCHFLX into RCHFLX_local for root node
  if (masterproc) then
   associate(nRch_trib => rch_per_proc(0))
   allocate(RCHFLX_local(nRch_mainstem+nRch_trib), stat=ierr)
   if (nRch_mainstem>0) then
     RCHFLX_local(1:nRch_mainstem) = RCHFLX_main(iens,1:nRch_mainstem)
   end if
   if (nRch_trib>0) then
     RCHFLX_local(nRch_mainstem+1:nRch_mainstem+nRch_trib) = RCHFLX_trib(iens,:)
   endif
   end associate
  else
   allocate(RCHFLX_local(rch_per_proc(pid)), stat=ierr)
   RCHFLX_local = RCHFLX_trib(iens,:)
  endif

  if (masterproc) then
    associate(nHRU_trib => hru_per_proc(0))
    allocate(basinRunoff(nHRU_mainstem+nHRU_trib), stat=ierr)
    if (nHRU_mainstem>0) then
      basinRunoff(1:nHRU_mainstem) = basinRunoff_main(1:nHRU_mainstem)
    end if
    if (nHRU_trib>0) then
      basinRunoff(nHRU_mainstem+1:nHRU_mainstem+nHRU_trib) = basinRunoff_trib(1:nHRU_trib)
    endif
    end associate
  else
    allocate(basinRunoff(hru_per_proc(pid)), stat=ierr)
    basinRunoff = basinRunoff_trib
  endif

  ! write time -- note time is just carried across from the input
  call write_netcdf(pioFileDesc, 'time', [timeVar(iTime)], [jTime], [1], ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (meta_rflx(ixRFLX%basRunoff)%varFile) then
    ! write the basin runoff at HRU (unit: the same as runoff input)
    call write_pnetcdf_recdim(pioFileDesc, 'basRunoff',basinRunoff, iodesc_hru_ro, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%instRunoff)%varFile) then
    ! write instataneous local runoff in each stream segment (m3/s)
    call write_pnetcdf_recdim(pioFileDesc, 'instRunoff', RCHFLX_local(:)%BASIN_QI, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
    ! write routed local runoff in each stream segment (m3/s)
    call write_pnetcdf_recdim(pioFileDesc, 'dlayRunoff', RCHFLX_local(:)%BASIN_QR(1), iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile) then
    ! write accumulated runoff (m3/s)
    call write_pnetcdf_recdim(pioFileDesc, 'sumUpstreamRunoff', RCHFLX_local(:)%UPSTREAM_QI, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%KWTroutedRunoff)%varFile) then
    ! write routed runoff (m3/s)
    call write_pnetcdf_recdim(pioFileDesc, 'KWTroutedRunoff', RCHFLX_local(:)%REACH_Q, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%KWEroutedRunoff)%varFile) then
    ! write routed runoff (m3/s)
    call write_pnetcdf_recdim(pioFileDesc, 'KWEroutedRunoff', RCHFLX_local(:)%REACH_Q, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%IRFroutedRunoff)%varFile) then
     ! write routed runoff (m3/s)
     call write_pnetcdf_recdim(pioFileDesc, 'IRFroutedRunoff', RCHFLX_local(:)%REACH_Q_IRF, iodesc_rch_flx, jTime, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

 END SUBROUTINE output


 ! *********************************************************************
 ! public subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE prep_output(ierr, message)

 ! saved public variables (usually parameters, or values not modified)
 USE public_var, ONLY: output_dir        ! output directory
 USE public_var, ONLY: case_name         ! simulation name ==> output filename head
 USE public_var, ONLY: calendar          ! calendar name
 USE public_var, ONLY: newFileFrequency  ! frequency for new output files (day, month, annual, single)
 USE public_var, ONLY: time_units        ! time units (seconds, hours, or days)
 ! saved global data
 USE globalData, ONLY: basinID,reachID   ! HRU and reach ID in network
 USE globalData, ONLY: modJulday         ! julian day: at model time step
 USE globalData, ONLY: modTime           ! previous and current model time
 USE globalData, ONLY: nEns, nHRU, nRch  ! number of ensembles, HRUs and river reaches
 USE globalData, ONLY: isFileOpen        ! file open/close status
 ! subroutines
 USE time_utils_module, ONLy: compCalday        ! compute calendar day
 USE time_utils_module, ONLy: compCalday_noleap ! compute calendar day

 implicit none

 ! input variables: none
 ! output variables
 integer(i4b), intent(out)       :: ierr             ! error code
 character(*), intent(out)       :: message          ! error message
 ! local variables
 logical(lgt)                    :: defNewOutputFile ! flag to define new output file
 integer(i4b)                    :: sec_in_day       ! second within day
 character(len=strLen)           :: cmessage         ! error message of downwind routine
 character(*),parameter          :: fmtYMDS='(a,I0.4,a,I0.2,a,I0.2,a,I0.5,a)'

 ! initialize error control
 ierr=0; message='prep_output/'

  ! get the time
  select case(trim(calendar))
   case('noleap','365_day')
    call compCalday_noleap(modJulday,modTime(1)%iy,modTime(1)%im,modTime(1)%id,modTime(1)%ih,modTime(1)%imin,modTime(1)%dsec,ierr,cmessage)
   case ('standard','gregorian','proleptic_gregorian')
    call compCalday(modJulday,modTime(1)%iy,modTime(1)%im,modTime(1)%id,modTime(1)%ih,modTime(1)%imin,modTime(1)%dsec,ierr,cmessage)
   case default;    ierr=20; message=trim(message)//'calendar name: '//trim(calendar)//' invalid'; return
  end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! print progress
  if (masterproc) then
    write(iulog,'(a,I4,4(x,I4))') new_line('a'), modTime(1)%iy, modTime(1)%im, modTime(1)%id, modTime(1)%ih, modTime(1)%imin
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
   sec_in_day = modTime(1)%ih*60*60+modTime(1)%imin*60+nint(modTime(1)%dsec)
   write(fileout, fmtYMDS) trim(output_dir)//trim(case_name)//'mizuRoute.h.', &
                           modTime(1)%iy, '-', modTime(1)%im, '-', modTime(1)%id, '-',sec_in_day,'.nc'
   write(iulog,*), trim(fileout)

   ! define output file
   call defineFile(trim(fileout),                         &  ! input: file name
                   nEns,                                  &  ! input: number of ensembles
                   nHRU,                                  &  ! input: number of HRUs
                   nRch,                                  &  ! input: number of stream segments
                   time_units,                            &  ! input: time units
                   calendar,                              &  ! input: calendar
                   ierr,cmessage)                            ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call openFile(pioSystem, pioFileDesc, trim(fileout), pio_typename, ncd_write, ierr, cmessage)
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
  USE globalData, ONLY: isFileOpen   ! file open/close status
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
                       nEns_in,         &  ! input: number of ensembles
                       nHRU_in,         &  ! input: number of HRUs
                       nRch_in,         &  ! input: number of stream segments
                       units_time,      &  ! input: time units
                       calendar,        &  ! input: calendar
                       ierr, message)      ! output: error control
 !Dependent modules
 USE var_lookup, ONLY: ixQdims, nQdims
 USE globalData, ONLY: meta_qDims
 USE globalData, ONLY: rch_per_proc             ! number of reaches assigned to each proc (size = num of procs+1)
 USE globalData, ONLY: hru_per_proc             ! number of hrus assigned to each proc (size = num of procs+1)

 implicit none
 ! input variables
 character(*), intent(in)          :: fname             ! filename
 integer(i4b), intent(in)          :: nEns_in           ! number of ensembles
 integer(i4b), intent(in)          :: nHRU_in           ! number of HRUs
 integer(i4b), intent(in)          :: nRch_in           ! number of stream segments
 character(*), intent(in)          :: units_time        ! time units
 character(*), intent(in)          :: calendar          ! calendar
 ! output variables
 integer(i4b), intent(out)         :: ierr              ! error code
 character(*), intent(out)         :: message           ! error message
 ! local variables
 integer(i4b),allocatable          :: dim_array(:)      ! dimension
 integer(i4b)                      :: jDim,iVar         ! dimension, and variable index
 integer(i4b)                      :: nDims             ! number of dimension
 integer(i4b)                      :: ix1, ix2          ! frst and last indices of global array for local array chunk
 integer(i4b)                      :: ixRch(nRch_in)    !
 integer(i4b)                      :: ixHRU(nHRU_in)    !
 integer(i4b)                      :: ixDim             !
 character(len=strLen)             :: cmessage          ! error message of downwind routine

 ! initialize error control
 ierr=0; message='defineFile/'

 ! populate q dimension meta (not sure if this should be done here...)
 meta_qDims(ixQdims%seg)%dimLength = nRch_in
 meta_qDims(ixQdims%hru)%dimLength = nHRU_in
 meta_qDims(ixQdims%ens)%dimLength = nEns_in

 ! Modify write option
 ! This is temporary
 if (routOpt==kinematicWave) then
   meta_rflx(ixRFLX%IRFroutedRunoff)%varFile = .false.
   meta_rflx(ixRFLX%KWEroutedRunoff)%varFile = .false.
 elseif (routOpt==kinematicWaveEuler) then
   meta_rflx(ixRFLX%IRFroutedRunoff)%varFile = .false.
   meta_rflx(ixRFLX%KWTroutedRunoff)%varFile = .false.
 elseif (routOpt==impulseResponseFunc) then
   meta_rflx(ixRFLX%KWTroutedRunoff)%varFile = .false.
   meta_rflx(ixRFLX%KWEroutedRunoff)%varFile = .false.
 elseif (routOpt==allRoutingMethods) then
   meta_rflx(ixRFLX%KWEroutedRunoff)%varFile = .false.
 end if
 ! runoff accumulation option
 if (doesAccumRunoff==0) then
   meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile = .false.
 end if
 ! basin runoff routing option
 if (doesBasinRoute==0) then
   meta_rflx(ixRFLX%instRunoff)%varFile = .false.
 end if

 ! pio initialization
 pio_numiotasks = nNodes/pio_stride
 call pio_sys_init(pid, mpicom_route,          & ! input: MPI related parameters
                   pio_stride, pio_numiotasks, & ! input: PIO related parameters
                   pio_rearranger, pio_root,   & ! input: PIO related parameters
                   pioSystem)                    ! output: PIO system descriptors

 ! For reach flux/volume
 if (masterproc) then
   ix1 = 1
 else
   ix1 = sum(rch_per_proc(-1:pid-1))+1
 endif
 ix2 = sum(rch_per_proc(-1:pid))
 ixRch = arth(1,1,nRch_in)

 call pio_decomp(pioSystem,              & ! input: pio system descriptor
                 ncd_double,             & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nRch_in],              & ! input: dimension length == global array size
                 ixRch(ix1:ix2),         & ! input:
                 iodesc_rch_flx)

! For HRU flux/volume
 if (masterproc) then
   ix1 = 1
 else
   ix1 = sum(hru_per_proc(-1:pid-1))+1
 endif
 ix2 = sum(hru_per_proc(-1:pid))
 ixHRU = arth(1,1,nHRU_in)

 call pio_decomp(pioSystem,     & ! input: pio system descriptor
                 ncd_double,    & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                 [nHRU_in],     & ! input: dimension length == global array size
                 ixHRU(ix1:ix2),& ! input:
                 iodesc_hru_ro)

 ! --------------------
 ! define file
 ! --------------------
 call createFile(pioSystem, trim(fname), pio_typename, pio_netcdf_format, pioFileDesc, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 do jDim =1,nQdims
   if (jDim ==ixQdims%time) then ! time dimension (unlimited)
    call def_dim(pioFileDesc, trim(meta_qDims(jDim)%dimName), recordDim, meta_qDims(jDim)%dimId)
   else
    call def_dim(pioFileDesc, trim(meta_qDims(jDim)%dimName), meta_qDims(jDim)%dimLength, meta_qDims(jDim)%dimId)
   endif
 end do

 ! define coordinate variable for time
 call def_var(pioFileDesc,                                &                                        ! pio file descriptor
             trim(meta_qDims(ixQdims%time)%dimName),      &                                        ! variable name
             [meta_qDims(ixQdims%time)%dimId], ncd_float, &                                        ! dimension array and type
             ierr, cmessage,                              &                                        ! error handle
             vdesc=trim(meta_qDims(ixQdims%time)%dimName), vunit=trim(units_time), vcal=calendar)  ! optional attributes
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! define hru ID and reach ID variables
 call def_var(pioFileDesc, 'basinID', [meta_qDims(ixQdims%hru)%dimId], ncd_int, ierr, cmessage, vdesc='basin ID', vunit='-')
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 call def_var(pioFileDesc, 'reachID', [meta_qDims(ixQdims%seg)%dimId], ncd_int, ierr, cmessage, vdesc='reach ID', vunit='-')
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! define flux variables
 do iVar=1,nVarsRFLX

  if (.not.meta_rflx(iVar)%varFile) cycle

  ! define dimension ID array
  nDims = size(meta_rflx(iVar)%varDim)
  if (allocated(dim_array)) then
    deallocate(dim_array)
  endif
  allocate(dim_array(nDims))
  do ixDim = 1, nDims
    dim_array(ixDim) = meta_qDims(meta_rflx(iVar)%varDim(ixDim))%dimId
  end do

  ! define variable
  call def_var(pioFileDesc,            &                 ! pio file descriptor
              meta_rflx(iVar)%varName, &                 ! variable name
              dim_array, ncd_float,    &                 ! dimension array and type
              ierr, cmessage,          &                 ! error handling
              vdesc=meta_rflx(iVar)%varDesc, vunit=meta_rflx(iVar)%varUnit)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end do

 ! end definitions
 call end_def(pioFileDesc, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE defineFile


END MODULE write_simoutput_pio
