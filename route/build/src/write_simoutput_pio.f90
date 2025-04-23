MODULE write_simoutput_pio

! Moudle wide shared data/routines
USE nrtype
USE var_lookup,     ONLY: ixRFLX, nVarsRFLX
USE dataTypes,      ONLY: STRFLX
USE datetime_data,  ONLY: datetime
USE public_var,     ONLY: iulog
USE public_var,     ONLY: integerMissing
USE globalData,     ONLY: runMode           ! 'standalone' or 'cesm-coupling'
USE globalData,     ONLY: pid, nNodes
USE globalData,     ONLY: masterproc
USE globalData,     ONLY: hfileout, hfileout_gage, rfileout
USE globalData,     ONLY: initHvars
USE historyFile,    ONLY: histFile
USE histVars_data,  ONLY: histVars
USE io_rpointfile,  ONLY: io_rpfile
USE ascii_utils,    ONLY: lower
USE pio_utils

implicit none
type(histVars), save            :: hVars                 !
type(histFile), save            :: hist_all_network      !
type(histFile), save            :: hist_gage             !

private
public::hVars
public::main_new_file
public::output
public::close_all
public::init_histFile

CONTAINS

 ! *********************************************************************
 ! public subroutine: main routine to define new output file
 ! *********************************************************************
 SUBROUTINE main_new_file(ierr, message)

    USE public_var, ONLY: outputAtGage
    USE public_var, ONLY: newFileFrequency  ! frequency for new output files (day, month, annual, single)
    USE globalData, ONLY: simDatetime       ! previous and current model time
    USE globalData, ONLY: reachID           !
    USE globalData, ONLY: basinID           !
    USE globalData, ONLY: nRch              !
    USE globalData, ONLY: nHRU              !
    USE globalData, ONLY: gage_meta_data    !

    implicit none
    ! Argument variables
    integer(i4b),   intent(out)     :: ierr             ! error code
    character(*),   intent(out)     :: message          ! error message
    ! local variables
    logical(lgt)                    :: createNewFile       ! logical to make alarm for restart writing
    integer(i4b)                    :: nGage               ! number of gauge points
    integer(i4b), allocatable       :: reach_id_local(:)   ! logical to make alarm for restart writing
    character(len=strLen)           :: cmessage            ! error message of downwind routine

    ierr=0; message='main_new_file/'

    createNewFile = newFileAlarm(simDatetime(0:1), newFileFrequency, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    if (createNewFile) then
      call close_all(ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call get_hfilename(simDatetime(1), ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! initialize and create history netcdfs
      hist_all_network = histFile(hfileout)

      call hist_all_network%createNC(ierr, cmessage, nRch_in=nRch, nHRU_in=nHRU)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call hist_all_network%openNC(ierr, message)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call hist_all_network%write_loc(reachID, basinID, ierr, message)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      if (outputAtGage) then

        hist_gage = histFile(hfileout_gage, gageOutput=.true.)
        nGage = gage_meta_data%gage_num()
        allocate(reach_id_local(nGage))
        reach_id_local = gage_meta_data%reach_id()

        call hist_gage%createNC(ierr, cmessage, nRch_in= nGage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        call hist_gage%openNC(ierr, message)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        call hist_gage%write_loc(reach_id_local, ierr, message)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

    end if

 END SUBROUTINE main_new_file

 ! *********************************************************************
 ! private subroutine: restart alarming
 ! *********************************************************************
 logical(lgt) FUNCTION newFileAlarm(inDatetime, alarmFrequency, ierr, message)

   implicit none
   ! Argument variables
   type(datetime), intent(in)    :: inDatetime(0:1)    ! datetime at previous and current timestep
   character(*),   intent(in)    :: alarmFrequency   ! new file frequency
   integer(i4b),   intent(out)   :: ierr             ! error code
   character(*),   intent(out)   :: message          ! error message

   ierr=0; message='new_file_alarm/'

   if (masterproc) then
     write(iulog,'(a,I4,4(1x,I4))') new_line('a'), inDatetime(1)%year(), inDatetime(1)%month(), &
           inDatetime(1)%day(), inDatetime(1)%hour(), inDatetime(1)%minute()
   endif

   ! check need for the new file
   select case(lower(trim(alarmFrequency)))
     case('single'); newFileAlarm=(inDatetime(0)%year() ==integerMissing)
     case('yearly'); newFileAlarm=(inDatetime(1)%year() /=inDatetime(0)%year())
     case('monthly');newFileAlarm=(inDatetime(1)%month()/=inDatetime(0)%month())
     case('daily');  newFileAlarm=(inDatetime(1)%day()  /=inDatetime(0)%day())
     case default; ierr=20; message=trim(message)//'unable to identify the option to define new output files'; return
   end select

 END FUNCTION newFileAlarm

 ! *********************************************************************
 ! public subroutine: output history file
 ! *********************************************************************
 SUBROUTINE output(ierr, message)

   USE public_var, ONLY: outputAtGage                 ! ascii containing last restart and history files
   USE public_var, ONLY: nOutFreq                     ! integer output frequency, i.e, written at every "nOutFreq" of simulation step
   USE public_var, ONLY: outputFrequency              ! writing frequency
   USE public_var, ONLY: histTimeStamp_offset         ! history time stamp offset from the start of time step [sec]
   USE globalData, ONLY: simDatetime                  ! previous,current and next model datetime
   USE globalData, ONLY: timeVar                      ! current simulation time variable
   USE globalData, ONLY: sec2tunit                    ! seconds per time unit
   USE globalData, ONLY: RCHFLX_trib                  ! reach flux data structure containing current flux variables
   USE globalData, ONLY: rch_per_proc                 ! number of reaches assigned to each proc (size = num of procs+1)
   USE globalData, ONLY: nRch_mainstem                ! number of mainstem reach
   USE globalData, ONLY: nTribOutlet                  ! number of
   USE globalData, ONLY: index_write_gage             ! reach index (w.r.t. global domain) corresponding gauge location
   USE globalData, ONLY: ioDesc_hru_float             ! pio decomposition descriptor for hru
   USE globalData, ONLY: ioDesc_rch_float             ! pio decomposition descriptor for reaches
   USE globalData, ONLY: ioDesc_gauge_float           ! pio decomposition descriptor for gauges
   USE nr_utils,   ONLY: arth

   implicit none
   ! Argument variables:
   integer(i4b), intent(out)   :: ierr               ! error code
   character(*), intent(out)   :: message            ! error message
   ! local variables:
   logical(lgt)                :: writeHistory          !
   integer(i4b)                :: nRch_local         ! number of reaches per processors
   real(dp)                    :: timeVar_local(1:2)
   real(dp),     allocatable   :: basinRunoff(:)
   integer(i4b), allocatable   :: index_write_all(:) ! indices in RCHFLX_trib to be written in netcdf
   character(strLen)           :: cmessage           ! error message of downwind routine

   ierr=0; message='output/'

   ! initialize if  histVars data - hVar is not initalized
   if (.not. initHvars) then
     call init_histVar_data(ierr, message)
     initHvars = .true.
   end if

   ! get unified basinRunoff array
   call get_proc_flux(ierr, cmessage, basinRunoff=basinRunoff)

   ! accumulate history variables
   timeVar_local = timeVar/sec2tunit
   call hVars%aggregate(timeVar_local, basinRunoff, RCHFLX_trib, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! history output alarm
   if (nOutFreq/=integerMissing) then
     writeHistory = writeAlarmNumeric(nOutFreq)
   else
     writeHistory = writeAlarmLiteral(simDatetime(1:2), outputFrequency, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

   if (writeHistory) then

     ! compute final output values (compute mean)
     call hVars%finalize()

     ! compute index array for each processors (herer
     if (masterproc) then
       nRch_local = sum(rch_per_proc(-1:pid))
       allocate(index_write_all(nRch_local))
       if (nRch_mainstem>0) then
         index_write_all(1:nRch_mainstem) = arth(1,1,nRch_mainstem)
       end if
       index_write_all(nRch_mainstem+1:nRch_local) = arth(nRch_mainstem+nTribOutlet+1, 1, rch_per_proc(0))
     else
       nRch_local = rch_per_proc(pid)
       allocate(index_write_all(nRch_local))
       index_write_all = arth(1,1,nRch_local)
     end if

     ! For cesm coupling mode, the timestamp of the average history output is the midpoint of time bounds
     if (trim(runMode)=='cesm-coupling') then
       histTimeStamp_offset = (hVars%timeVar(2) - hVars%timeVar(1))*sec2tunit/2
     end if

     ! write time variables (time and time bounds)
     call hist_all_network%write_time(hVars, histTimeStamp_offset, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! write out output variables in history files
     ! fluxes in HRUs (catchments)
     call hist_all_network%write_flux_hru(hVars, ioDesc_hru_float, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! fluxes in RCHs (river reaches)
     call hist_all_network%write_flux_rch(hVars, ioDesc_rch_float, index_write_all, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     call hist_all_network%sync(ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     if (outputAtGage) then
       ! write time variables (time and time bounds)
       call hist_gage%write_time(hVars, histTimeStamp_offset, ierr, cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

       call hist_gage%write_flux_rch(hVars, ioDesc_gauge_float, index_write_gage, ierr, cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

       call hist_gage%sync(ierr, cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end if

     ! refresh history output buffers
     call hVars%refresh()

   end if

 END SUBROUTINE output

 ! *********************************************************************
 ! private subroutine: restart alarming
 ! *********************************************************************
 logical(lgt) FUNCTION writeAlarmNumeric(alarmFrequency)

   implicit none
   ! Argument variables
   integer(i4b),   intent(in)    :: alarmFrequency   ! numeric output frequency

   writeAlarmNumeric = .false.
   if (hVars%nt == alarmFrequency) writeAlarmNumeric = .true.

 END FUNCTION writeAlarmNumeric

 ! *********************************************************************
 ! private subroutine: restart alarming
 ! *********************************************************************
 logical(lgt) FUNCTION writeAlarmLiteral(inDatetime, alarmFrequency, ierr, message)

   implicit none
   ! Argument variables
   type(datetime), intent(in)    :: inDatetime(1:2)  ! datetime at previous and current timestep
   character(*),   intent(in)    :: alarmFrequency   ! literal output frequency
   integer(i4b),   intent(out)   :: ierr             ! error code
   character(*),   intent(out)   :: message          ! error message

   writeAlarmLiteral = .false.
   select case(lower(trim(alarmFrequency)))
     case('yearly'); writeAlarmLiteral=(inDatetime(1)%year() /=inDatetime(2)%year())
     case('monthly');writeAlarmLiteral=(inDatetime(1)%month()/=inDatetime(2)%month())
     case('daily');  writeAlarmLiteral=(inDatetime(1)%day()  /=inDatetime(2)%day())
     case default; ierr=20; message=trim(message)//'unable to identify the option to define new output files'; return
   end select

 END FUNCTION writeAlarmLiteral

 ! *********************************************************************
 ! private subroutine - organize flux data array or structure into per processor
 ! *********************************************************************
 SUBROUTINE get_proc_flux(ierr, message, basinRunoff)

   USE globalData, ONLY: nHRU_mainstem       ! number of mainstem HRUs
   USE globalData, ONLY: basinRunoff_main    ! mainstem only HRU runoff
   USE globalData, ONLY: basinRunoff_trib    ! tributary only HRU runoff
   USE globalData, ONLY: hru_per_proc        ! number of hrus assigned to each proc (size = num of procs+1)

   implicit none
   ! Argument variables
   real(dp),    allocatable, optional, intent(out) :: basinRunoff(:)
   integer(i4b),                       intent(out) :: ierr             ! error code
   character(*),                       intent(out) :: message          ! error message
   ! local variables
   integer(i4b)                          :: nHRU_local
   character(strLen)                     :: cmessage         ! error message of downwind routine

   if (present(basinRunoff)) then
     if (masterproc) then
       nHRU_local = nHRU_mainstem + hru_per_proc(0)
       allocate(basinRunoff(nHRU_local), stat=ierr, errmsg=cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [basinRunoff]'; return; endif
       if (nHRU_mainstem>0) basinRunoff(1:nHRU_mainstem) = basinRunoff_main(1:nHRU_mainstem)
       if (hru_per_proc(0)>0) basinRunoff(nHRU_mainstem+1:nHRU_local) = basinRunoff_trib(:)
     else
       nHRU_local = hru_per_proc(pid)
       allocate(basinRunoff(nHRU_local), stat=ierr, errmsg=cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [basinRunoff]'; return; endif
       basinRunoff = basinRunoff_trib
     endif
   end if

 END SUBROUTINE get_proc_flux

 ! *********************************************************************
 ! private subroutine: define history NetCDF file name
 ! *********************************************************************
 SUBROUTINE get_hfilename(inDatetime, ierr, message)

   USE public_var, ONLY: output_dir          ! output directory
   USE public_var, ONLY: case_name           ! simulation name ==> output filename head
   USE public_var, ONLY: newFileFrequency    ! new file frequency
   USE public_var, ONLY: outputAtGage        ! ascii containing last restart and history files
   USE public_var, ONLY: secprhour,secprmin  ! time conversion to sec
   USE globalData, ONLY: hfile_dayStamp      ! day stamp for history file name - period-start or period-end

   implicit none
   ! argument variables
   type(datetime), intent(in)  :: inDatetime     ! datetime at current timestep
   integer(i4b),   intent(out) :: ierr           ! error code
   character(*),   intent(out) :: message        ! error message
   ! local variables
   character(strLen)           :: cmessage        ! error message from subroutine
   type(datetime)              :: datetime_stamp  ! datetime used in history file name
   integer(i4b)                :: sec_in_day      ! second within day
   character(len=strLen)       :: timeString
   character(len=27),parameter :: fmtYMDS='(I0.4,a,I0.2,a,I0.2,a,I0.5)'
   character(len=20),parameter :: fmtYMD ='(I0.4,a,I0.2,a,I0.2)'
   character(len=13),parameter :: fmtYM  ='(I0.4,a,I0.2)'
   character(len=6),parameter  :: fmtY   ='(I0.4)'

   ierr=0; message='get_hfilename/'

   ! construct time stamp string in history file: timeString
   select case(trim(newFileFrequency))
     case('single')
       sec_in_day = inDatetime%hour()*int(secprhour)+inDatetime%minute()*int(secprmin)+nint(inDatetime%sec())
       write(timeString, fmtYMDS) inDatetime%year(),'-',inDatetime%month(),'-',inDatetime%day(),'-',sec_in_day
     case('daily')
       if (trim(hfile_dayStamp)=='period-start') then
         sec_in_day = inDatetime%hour()*int(secprhour)+inDatetime%minute()*int(secprmin)+nint(inDatetime%sec())
         write(timeString, fmtYMDS) inDatetime%year(),'-',inDatetime%month(),'-',inDatetime%day(),'-',sec_in_day
       else if (trim(hfile_dayStamp)=='period-end') then
         datetime_stamp = inDatetime%add_day(1, ierr, cmessage)
         sec_in_day = datetime_stamp%hour()*int(secprhour)+datetime_stamp%minute()*int(secprmin)+nint(datetime_stamp%sec())
         write(timeString, fmtYMDS) datetime_stamp%year(),'-',datetime_stamp%month(),'-',datetime_stamp%day(),'-',sec_in_day
       end if
     case('monthly')
       write(timeString, fmtYM) inDatetime%year(),'-',inDatetime%month()
     case('yearly')
       write(timeString, fmtY) inDatetime%year()
     case default; ierr=20; message=trim(message)//'Invalid file frequency'; return
   end select

   ! construct history file name: hfileout
   select case(trim(runMode))
     case('cesm-coupling')
       write(hfileout,'(a)') trim(output_dir)//trim(case_name)//'.mizuroute.h.'//trim(timeString)//'.nc'
       if (outputAtGage) then
         write(hfileout_gage,'(a)') trim(output_dir)//trim(case_name)//'.mizuroute.h_gauge.'//trim(timeString)//'.nc'
       end if
     case('standalone')
       write(hfileout, '(a)') trim(output_dir)//trim(case_name)//'.h.'//trim(timeString)//'.nc'
       if (outputAtGage) then
         write(hfileout_gage, '(a)') trim(output_dir)//trim(case_name)//'.h_gauge.'//trim(timeString)//'.nc'
       end if
     case default; ierr=20; message=trim(message)//'unable to identify the run option. Avaliable options are standalone and cesm-coupling'; return
   end select

 END SUBROUTINE get_hfilename

 ! *********************************************************************
 ! private subroutine: close all history files
 ! *********************************************************************
 SUBROUTINE close_all(ierr, message)

   USE public_var, ONLY: outputAtGage        ! ascii containing last restart and history files

   implicit none
   ! argument variables
   integer(i4b),   intent(out)     :: ierr                ! error code
   character(*),   intent(out)     :: message             ! error message

   ierr=0; message='close_all/'

   if (hist_all_network%fileOpen()) then
     call hist_all_network%closeNC()
   end if

   if (outputAtGage) then
     if (hist_gage%fileOpen()) then
       call hist_gage%closeNC()
     end if
   end if
 END SUBROUTINE

 ! *********************************************************************
 ! public subroutine: initialize history file for contiue run
 ! *********************************************************************
 SUBROUTINE init_histFile(ierr, message)

   USE public_var,  ONLY: outputAtGage      ! ascii containing last restart and history files

   implicit none
   ! argument variables
   integer(i4b),   intent(out)     :: ierr                ! error code
   character(*),   intent(out)     :: message             ! error message
   ! local variables
   character(len=strLen)           :: cmessage            ! error message of downwind routine

   ierr=0; message='init_histFile/'

   ! get history file names to append and assign it to hfileout
   call io_rpfile('r', ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! initialize and create history netcdfs
   hist_all_network = histFile(hfileout)

   call hist_all_network%openNc(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (outputAtGage) then
     hist_gage = histFile(hfileout_gage, gageOutput=.true.)

     call hist_gage%openNC(ierr, message)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

 END SUBROUTINE init_histFile

 ! *********************************************************************
 ! public subroutine: initialize history variable buffers
 ! *********************************************************************
 SUBROUTINE init_histVar_data(ierr, message)          ! output: error control

   USE globalData, ONLY: nRch_mainstem, nRch_trib
   USE globalData, ONLY: nHRU_mainstem, nHRU_trib
   USE globalData, ONLY: nTribOutlet

   implicit none
   ! Argument variables:
   integer(i4b),         intent(out) :: ierr             ! error code
   character(*),         intent(out) :: message          ! error message
   ! local variable
   character(len=strLen)             :: cmessage         ! error message from subroutine
   integer(i4b)                      :: nHru_local       ! nHRU in each processor
   integer(i4b)                      :: nRch_local       ! nRch in each processor

   ierr=0; message='init_histVar_data/'

   if (masterproc) then
     nHru_local=nHru_mainstem+nHru_trib
     nRch_local=nRch_mainstem+nRch_trib+nTribOutlet
   else
     nHru_local=nHru_trib
     nRch_local=nRch_trib
   end if
   hVars = histVars(nHru_local, nRch_local, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 END SUBROUTINE init_histVar_data

END MODULE write_simoutput_pio
