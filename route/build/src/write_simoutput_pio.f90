MODULE write_simoutput_pio

! Moudle wide shared data
USE nrtype
USE var_lookup,        ONLY: ixRFLX, nVarsRFLX
USE dataTypes,         ONLY: STRFLX
USE public_var,        ONLY: iulog
USE public_var,        ONLY: integerMissing
USE globalData,        ONLY: idxSUM,idxIRF,idxKWT, &
                             idxKW,idxMC,idxDW
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
USE globalData,        ONLY: pioSystem
USE globalData,        ONLY: runMode
USE globalData,        ONLY: hfileout, rfileout
! Moudle wide external modules
USE nr_utility_module, ONLY: arth
USE pio_utils

implicit none

! module-wide variables
integer(i4b),         save :: jTime                ! time step in output netCDF
type(file_desc_t),    save :: pioFileDesc          ! PIO data identifying the file
type(io_desc_t),      save :: iodesc_rch_flx       ! PIO domain decomposition data for reach flux [nRch]
type(io_desc_t),      save :: iodesc_hru_ro        ! PIO domain decomposition data for hru runoff [nHRU]

integer(i4b),parameter     :: recordDim=-999       ! record dimension Indicator

private
public::main_new_file
public::init_histFile
public::prep_output
public::output
public::close_output_nc

CONTAINS

 ! *********************************************************************
 ! public subroutine: main routine to define new output file
 ! *********************************************************************
 SUBROUTINE main_new_file(ierr, message)

   USE globalData, ONLY: simDatetime   ! previous and current model time

  implicit none
  ! output variables
  integer(i4b),   intent(out)          :: ierr             ! error code
  character(*),   intent(out)          :: message          ! error message
  ! local variables
  logical(lgt)                         :: newFileAlarm     ! logical to make alarm for restart writing
  character(len=strLen)                :: cmessage         ! error message of downwind routine

  ierr=0; message='main_new_file/'

  call new_file_alarm(newFileAlarm, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (newFileAlarm) then
    call prep_output(ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

  simDatetime(0) = simDatetime(1)

 END SUBROUTINE main_new_file


 ! *********************************************************************
 ! private subroutine: restart alarming
 ! *********************************************************************
 SUBROUTINE new_file_alarm(newFileAlarm, ierr, message)

   USE public_var,        ONLY: newFileFrequency  ! frequency for new output files (day, month, annual, single)
   USE globalData,        ONLY: simDatetime       ! previous and current model time

   implicit none
   ! output
   logical(lgt),   intent(out)          :: newFileAlarm     ! logical to make alarm for creating new output file
   integer(i4b),   intent(out)          :: ierr             ! error code
   character(*),   intent(out)          :: message          ! error message

   ierr=0; message='new_file_alarm/'

   ! print progress
   if (masterproc) then
     write(iulog,'(a,I4,4(x,I4))') new_line('a'), simDatetime(1)%year(), simDatetime(1)%month(), simDatetime(1)%day(), simDatetime(1)%hour(), simDatetime(1)%minute()
   endif

   ! check need for the new file
   select case(newFileFrequency)
     case('single'); newFileAlarm=(simDatetime(0)%year() ==integerMissing)
     case('annual'); newFileAlarm=(simDatetime(1)%year() /=simDatetime(0)%year())
     case('month');  newFileAlarm=(simDatetime(1)%month()/=simDatetime(0)%month())
     case('day');    newFileAlarm=(simDatetime(1)%day()  /=simDatetime(0)%day())
     case default; ierr=20; message=trim(message)//'unable to identify the option to define new output files'; return
   end select

 END SUBROUTINE new_file_alarm


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
  integer(i4b)                    :: ix
  integer(i4b)                    :: nRch_local
  real(dp),    allocatable        :: basinRunoff(:)
  real(dp),    allocatable        :: array_temp(:)

  ierr=0; message='output/'

  iens = 1

  jTime = jTime+1

  ! Need to combine mainstem RCHFLX and tributary RCHFLX into RCHFLX_local for root node
  if (masterproc) then
    nRch_local = nRch_mainstem+rch_per_proc(0)
    allocate(RCHFLX_local(nRch_local), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHFLX_local]'; return; endif
    allocate(array_temp(nRch_local), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif

    if (nRch_mainstem>0) RCHFLX_local(1:nRch_mainstem) = RCHFLX_main(iens, 1:nRch_mainstem)
    if (rch_per_proc(0)>0) RCHFLX_local(nRch_mainstem+1:nRch_local) = RCHFLX_trib(iens,:)
  else
    nRch_local = rch_per_proc(pid)
    allocate(RCHFLX_local(nRch_local), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [RCHFLX_local]'; return; endif
    allocate(array_temp(nRch_local), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif

    RCHFLX_local = RCHFLX_trib(iens,:)
  endif

  if (masterproc) then
    associate(nHRU_trib => hru_per_proc(0))
    allocate(basinRunoff(nHRU_mainstem+nHRU_trib), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [basinRunoff]'; return; endif
    if (nHRU_mainstem>0) basinRunoff(1:nHRU_mainstem) = basinRunoff_main(1:nHRU_mainstem)
    if (nHRU_trib>0) basinRunoff(nHRU_mainstem+1:nHRU_mainstem+nHRU_trib) = basinRunoff_trib(1:nHRU_trib)
    end associate
  else
    allocate(basinRunoff(hru_per_proc(pid)), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [basinRunoff]'; return; endif
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

  ! write instataneous local runoff in each stream segment (m3/s)
  if (meta_rflx(ixRFLX%instRunoff)%varFile) then
    array_temp(1:nRch_local) = RCHFLX_local(:)%BASIN_QI
    call write_pnetcdf_recdim(pioFileDesc, 'instRunoff', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! write routed local runoff in each stream segment (m3/s)
  if (meta_rflx(ixRFLX%dlayRunoff)%varFile) then
    array_temp(1:nRch_local) = RCHFLX_local(:)%BASIN_QR(1)
    call write_pnetcdf_recdim(pioFileDesc, 'dlayRunoff', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile) then
    do ix=1,nRch_local
      array_temp(ix) = RCHFLX_local(ix)%ROUTE(idxSUM)%REACH_Q
    end do
    call write_pnetcdf_recdim(pioFileDesc, 'sumUpstreamRunoff', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%KWTroutedRunoff)%varFile) then
    do ix=1,nRch_local
      array_temp(ix) = RCHFLX_local(ix)%ROUTE(idxKWT)%REACH_Q
    end do
    call write_pnetcdf_recdim(pioFileDesc, 'KWTroutedRunoff', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%IRFroutedRunoff)%varFile) then
    do ix=1,nRch_local
      array_temp(ix) = RCHFLX_local(ix)%ROUTE(idxIRF)%REACH_Q
    end do
    call write_pnetcdf_recdim(pioFileDesc, 'IRFroutedRunoff', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%KWroutedRunoff)%varFile) then
    do ix=1,nRch_local
      array_temp(ix) = RCHFLX_local(ix)%ROUTE(idxKW)%REACH_Q
    end do
    call write_pnetcdf_recdim(pioFileDesc, 'KWroutedRunoff', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%MCroutedRunoff)%varFile) then
    do ix=1,nRch_local
      array_temp(ix) = RCHFLX_local(ix)%ROUTE(idxMC)%REACH_Q
    end do
    call write_pnetcdf_recdim(pioFileDesc, 'MCroutedRunoff', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%DWroutedRunoff)%varFile) then
    do ix=1,nRch_local
      array_temp(ix) = RCHFLX_local(ix)%ROUTE(idxDW)%REACH_Q
    end do
    call write_pnetcdf_recdim(pioFileDesc, 'DWroutedRunoff', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  if (meta_rflx(ixRFLX%volume)%varFile) then
    do ix=1,nRch_local
      array_temp(ix) = RCHFLX_local(ix)%ROUTE(idxIRF)%REACH_VOL(1)
    end do
    call write_pnetcdf_recdim(pioFileDesc, 'volume', array_temp, iodesc_rch_flx, jTime, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  call sync_file(pioFileDesc, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE output


 ! *********************************************************************
 ! public subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE prep_output(ierr, message)

 ! saved public variables (usually parameters, or values not modified)
 USE public_var, ONLY: output_dir        ! output directory
 USE public_var, ONLY: case_name         ! simulation name ==> output filename head
 USE public_var, ONLY: calendar          ! calendar name
 USE public_var, ONLY: time_units        ! time units (seconds, hours, or days)
 USE public_var, ONLY: restart_dir       ! directory for restart output
 USE public_var, ONLY: rpntfil           ! ascii containing last restart and history files
 ! saved global data
 USE globalData, ONLY: basinID, reachID  ! HRU and reach ID in network
 USE globalData, ONLY: simDatetime       ! previous and current model time
 USE globalData, ONLY: nEns, nHRU, nRch  ! number of ensembles, HRUs and river reaches
 USE globalData, ONLY: isHistFileOpen    ! history file open/close status
 ! subroutines
 USE time_utils_module, ONLy: compCalday        ! compute calendar day
 USE time_utils_module, ONLy: compCalday_noleap ! compute calendar day

 implicit none
 ! argument variables
 integer(i4b), intent(out)       :: ierr             ! error code
 character(*), intent(out)       :: message          ! error message
 ! local variables
 integer(i4b)                    :: sec_in_day       ! second within day
 character(len=strLen)           :: cmessage         ! error message of downwind routine
 character(*),parameter          :: fmtYMDS='(a,I0.4,a,I0.2,a,I0.2,a,I0.5,a)'
 logical(lgt),parameter          :: createNewFile=.true.

 ierr=0; message='prep_output/'

   ! close netcdf only if is is open
   call close_output_nc(isHistFileOpen)

   jTime=0

   ! Define filename
   sec_in_day = simDatetime(1)%hour()*60*60+simDatetime(1)%minute()*60+nint(simDatetime(1)%sec())
   select case(trim(runMode))
     case('cesm-coupling')
       write(hfileout, fmtYMDS) trim(output_dir)//trim(case_name)//'.mizuroute.h.', &
                               simDatetime(1)%year(),'-',simDatetime(1)%month(),'-',simDatetime(1)%day(),'-',sec_in_day,'.nc'
     case('standalone')
       write(hfileout, fmtYMDS) trim(output_dir)//trim(case_name)//'.h.', &
                               simDatetime(1)%year(),'-',simDatetime(1)%month(),'-',simDatetime(1)%day(),'-',sec_in_day,'.nc'
     case default; ierr=20; message=trim(message)//'unable to identify the run option. Avaliable options are standalone and cesm-coupling'; return
   end select

   call defineFile(trim(hfileout),                        &  ! input: file name
                   nEns,                                  &  ! input: number of ensembles
                   nHRU,                                  &  ! input: number of HRUs
                   nRch,                                  &  ! input: number of stream segments
                   time_units,                            &  ! input: time units
                   calendar,                              &  ! input: calendar
                   createNewFile,                         &
                   ierr,cmessage)                            ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call openFile(pioSystem, pioFileDesc, trim(hfileout), pio_typename, ncd_write, isHistFileOpen, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (meta_rflx(ixRFLX%basRunoff)%varFile) then
     call write_netcdf(pioFileDesc, 'basinID', basinID, [1], [nHRU], ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   call write_netcdf(pioFileDesc, 'reachID', reachID, [1], [nRch], ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (masterproc) then
     open(1, file = trim(restart_dir)//trim(rpntfil), status='replace', action='write')
     write(1,'(a)') trim(hfileout)
     write(1,'(a)') trim(rfileout)
     close(1)
   end if

 END SUBROUTINE prep_output


 SUBROUTINE close_output_nc(fileStatus)
  implicit none
  logical(lgt), intent(inout) :: fileStatus
  if (fileStatus) then
   call closeFile(pioFileDesc, fileStatus)
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
                       createNewFile,   &
                       ierr, message)      ! output: error control
 !Dependent modules
 USE var_lookup, ONLY: ixQdims, nQdims
 USE globalData, ONLY: meta_qDims
 USE globalData, ONLY: rch_per_proc             ! number of reaches assigned to each proc (size = num of procs+1)
 USE globalData, ONLY: hru_per_proc             ! number of hrus assigned to each proc (size = num of procs+1)

 implicit none
 ! argument variables
 character(*), intent(in)          :: fname             ! filename
 integer(i4b), intent(in)          :: nEns_in           ! number of ensembles
 integer(i4b), intent(in)          :: nHRU_in           ! number of HRUs
 integer(i4b), intent(in)          :: nRch_in           ! number of stream segments
 character(*), intent(in)          :: units_time        ! time units
 character(*), intent(in)          :: calendar          ! calendar
 logical(lgt), intent(in)          :: createNewFile     !
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

 ierr=0; message='defineFile/'

 ! populate q dimension meta (not sure if this should be done here...)
 meta_qDims(ixQdims%seg)%dimLength = nRch_in
 meta_qDims(ixQdims%hru)%dimLength = nHRU_in
 meta_qDims(ixQdims%ens)%dimLength = nEns_in

 ! pio initialization
 if (trim(runMode)=='standalone') then
   pio_numiotasks = nNodes/pio_stride
   call pio_sys_init(pid, mpicom_route,          & ! input: MPI related parameters
                     pio_stride, pio_numiotasks, & ! input: PIO related parameters
                     pio_rearranger, pio_root,   & ! input: PIO related parameters
                     pioSystem)                    ! output: PIO system descriptors
 endif

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

 if (.not. createNewFile) return

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
             ncd_float,                                   &                                        ! variable type
             ierr, cmessage,                              &                                        ! error handle
             pioDimId=[meta_qDims(ixQdims%time)%dimId],   &                                        ! dimension array
             vdesc=trim(meta_qDims(ixQdims%time)%dimName), vunit=trim(units_time), vcal=calendar)  ! optional attributes
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! define hru ID and reach ID variables
 if (meta_rflx(ixRFLX%basRunoff)%varFile) then
   call def_var(pioFileDesc, 'basinID', ncd_int, ierr, cmessage, pioDimId=[meta_qDims(ixQdims%hru)%dimId], vdesc='basin ID', vunit='-')
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 call def_var(pioFileDesc, 'reachID', ncd_int, ierr, cmessage, pioDimId=[meta_qDims(ixQdims%seg)%dimId], vdesc='reach ID', vunit='-')
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
   call def_var(pioFileDesc,             &                 ! pio file descriptor
                meta_rflx(iVar)%varName, &                 ! variable name
                ncd_float,               &                 ! dimension array and type
                ierr, cmessage,          &                 ! error handling
                pioDimId=dim_array,      &                 ! dimension id
                vdesc=meta_rflx(iVar)%varDesc, vunit=meta_rflx(iVar)%varUnit)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do

 ! end definitions
 call end_def(pioFileDesc, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE defineFile

 SUBROUTINE init_histFile(ierr, message)

   USE public_var,        ONLY: calendar          ! calendar name
   USE public_var,        ONLY: time_units        ! time units (seconds, hours, or days)
   USE public_var,        ONLY: restart_dir       ! restart directory
   USE public_var,        ONLY: rpntfil           !
   USE globalData,        ONLY: isHistFileOpen    ! history file open/close status
   USE globalData,        ONLY: nEns, nHRU, nRch  ! number of ensembles, HRUs and river reaches

   implicit none
   ! argument variables
   integer(i4b),   intent(out)          :: ierr             ! error code
   character(*),   intent(out)          :: message          ! error message
   ! local variables
   character(len=strLen)                :: cmessage         ! error message of downwind routine
   logical(lgt),parameter               :: createNewFile=.false.

   ierr=0; message='init_histFile/'

   open(1, file = trim(restart_dir)//trim(rpntfil), status='old', action='read')
   read(1, '(A)') hfileout
   read(1, '(A)') rfileout
   close(1)

   call defineFile(trim(hfileout),                        &  ! input: file name
                   nEns,                                  &  ! input: number of ensembles
                   nHRU,                                  &  ! input: number of HRUs
                   nRch,                                  &  ! input: number of stream segments
                   time_units,                            &  ! input: time units
                   calendar,                              &  ! input: calendar
                   createNewFile,                         &  ! input:
                   ierr,cmessage)                            ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! if file is not exist, get not-descriptive error
   call openFile(pioSystem, pioFileDesc, trim(hfileout), pio_typename, ncd_write, isHistFileOpen, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call inq_dim_len(pioFileDesc, 'time', jTime)

   if (masterproc) then
     open(1, file = trim(restart_dir)//trim(rpntfil), status='replace', action='write')
     write(1,'(a)') trim(hfileout)
     write(1,'(a)') trim(rfileout)
     close(1)
   end if

 END SUBROUTINE init_histFile

END MODULE write_simoutput_pio
