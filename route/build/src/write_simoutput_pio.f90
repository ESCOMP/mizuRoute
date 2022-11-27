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
USE historyFile,    ONLY: histFile
USE io_rpointfile,  ONLY: io_rpfile
USE ascii_utils,    ONLY: lower
USE pio_utils

implicit none
type(histFile), save            :: hist_all_network      !
type(histFile), save            :: hist_gage             !
integer(i4b), allocatable, save :: index_write_gage(:)   !

private
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
    USE globalData, ONLY: gage_data         !
    USE globalData, ONLY: pioSystem

    implicit none
    ! Argument variables
    integer(i4b),   intent(out)     :: ierr             ! error code
    character(*),   intent(out)     :: message          ! error message
    ! local variables
    integer(i4b), allocatable       :: compdof_rch(:)      !
    integer(i4b), allocatable       :: compdof_rch_gage(:) !
    integer(i4b), allocatable       :: compdof_hru(:)      !
    logical(lgt)                    :: createNewFile       ! logical to make alarm for restart writing
    character(len=strLen)           :: cmessage            ! error message of downwind routine

    ierr=0; message='main_new_file/'

    createNewFile = newFileAlarm(simDatetime, newFileFrequency, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    if (createNewFile) then
      call close_all()

      call get_hfilename(simDatetime(1), ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! initialize and create history netcdfs
      select case(trim(runMode))
        case('standalone');    hist_all_network = histFile(hfileout)
        case('cesm-coupling'); hist_all_network = histFile(hfileout, pioSys=pioSystem)
        case default; ierr=20; message=trim(message)//'rumMode: '//trim(runMode)//': must be standalone or cesm-coupling'; return
      end select

      call get_compdof_all_network(compdof_rch, compdof_hru)

      call hist_all_network%set_compdof(compdof_rch, compdof_hru, nRch, nHRU)

      call hist_all_network%createNC(nRch, nHRU, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call hist_all_network%openNC(ierr, message)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call hist_all_network%write_loc(reachID, basinID, ierr, message)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      if (outputAtGage) then

        ! initialize and create history netcdfs
        select case(trim(runMode))
          case('standalone');    hist_gage = histFile(hfileout_gage, gageOutput=.true.)
          case('cesm-coupling'); hist_gage = histFile(hfileout_gage, pioSys=pioSystem, gageOutput=.true.)
          case default; ierr=20; message=trim(message)//'rumMode: '//trim(runMode)//': must be standalone or cesm-coupling'; return
        end select

        call get_compdof_gage(compdof_rch_gage, ierr, cmessage)

        call hist_gage%set_compdof(compdof_rch_gage, gage_data%nGage)

        call hist_gage%createNC(gage_data%nGage, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        call hist_gage%openNC(ierr, message)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        call hist_gage%write_loc(gage_data%reachID, ierr, message)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

    end if

    simDatetime(0) = simDatetime(1)

    ! update history files
    call io_rpfile('w', ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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
     write(iulog,'(a,I4,4(x,I4))') new_line('a'), inDatetime(1)%year(), inDatetime(1)%month(), inDatetime(1)%day(), inDatetime(1)%hour(), inDatetime(1)%minute()
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

   USE public_var, ONLY: outputAtGage      ! ascii containing last restart and history files
   USE globalData, ONLY: iTime
   USE globalData, ONLY: timeVar
   USE globalData, ONLY: rch_per_proc      ! number of reaches assigned to each proc (size = num of procs+1)
   USE globalData, ONLY: nRch_mainstem     ! number of mainstem reach
   USE globalData, ONLY: nTribOutlet       ! number of
   USE nr_utils,   ONLY: arth

   implicit none
   ! Argument variables:
   integer(i4b), intent(out)   :: ierr               ! error code
   character(*), intent(out)   :: message            ! error message
   ! local variables:
   integer(i4b)                :: nRch_local         ! number of reaches per processors
   integer(i4b)                :: nRch_root          ! number of reaches in root processors (including halo reaches)
   integer(i4b), allocatable   :: index_write_all(:) ! indices in RCHFLX_trib to be written in netcdf
   character(strLen)           :: cmessage           ! error message of downwind routine

   ierr=0; message='output/'

   ! compute index array for each processors (herer
   if (masterproc) then
     nRch_local = sum(rch_per_proc(-1:pid))
     nRch_root = nRch_mainstem+nTribOutlet+rch_per_proc(0)
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

   call hist_all_network%write_flux(timeVar(iTime), index_write_all, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (outputAtGage) then
     call hist_gage%write_flux(timeVar(iTime), index_write_gage, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if

 END SUBROUTINE output

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

   call reach_subset(reachID_local, gage_data, compdof=compdof_rch, index2=index_write_gage)

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

 ! *********************************************************************
 ! private subroutine: define history NetCDF file name
 ! *********************************************************************
 SUBROUTINE get_hfilename(inDatetime, ierr, message)

   USE public_var, ONLY: output_dir        ! output directory
   USE public_var, ONLY: case_name         ! simulation name ==> output filename head
   USE public_var, ONLY: outputAtGage      ! ascii containing last restart and history files

   implicit none
   ! argument variables
   type(datetime), intent(in)  :: inDatetime     ! datetime at previous and current timestep
   integer(i4b),   intent(out) :: ierr           ! error code
   character(*),   intent(out) :: message        ! error message
   ! local variables
   integer(i4b)              :: sec_in_day       ! second within day
   character(*),parameter    :: fmtYMDS='(a,I0.4,a,I0.2,a,I0.2,a,I0.5,a)'

   ierr=0; message='get_hfilename/'

   ! Define history filename
   sec_in_day = inDatetime%hour()*60*60+inDatetime%minute()*60+nint(inDatetime%sec())
   select case(trim(runMode))
     case('cesm-coupling')
       write(hfileout, fmtYMDS) trim(output_dir)//trim(case_name)//'.mizuroute.h.', &
                             inDatetime%year(),'-',inDatetime%month(),'-',inDatetime%day(),'-',sec_in_day,'.nc'

       if (outputAtGage) then
         write(hfileout_gage, fmtYMDS) trim(output_dir)//trim(case_name)//'_gauge.mizuroute.h.', &
                                    inDatetime%year(),'-',inDatetime%month(),'-',inDatetime%day(),'-',sec_in_day,'.nc'
       end if
     case('standalone')
       write(hfileout, fmtYMDS) trim(output_dir)//trim(case_name)//'.h.', &
                             inDatetime%year(),'-',inDatetime%month(),'-',inDatetime%day(),'-',sec_in_day,'.nc'
       if (outputAtGage) then
         write(hfileout_gage, fmtYMDS) trim(output_dir)//trim(case_name)//'_gauge.h.', &
                                    inDatetime%year(),'-',inDatetime%month(),'-',inDatetime%day(),'-',sec_in_day,'.nc'
       end if
     case default; ierr=20; message=trim(message)//'unable to identify the run option. Avaliable options are standalone and cesm-coupling'; return
   end select

 END SUBROUTINE get_hfilename

 ! *********************************************************************
 ! private subroutine: close all history files
 ! *********************************************************************
 SUBROUTINE close_all()
   implicit none
   if (hist_all_network%fileOpen()) then
     call hist_all_network%cleanup()
     call hist_all_network%closeNC()
   end if
   if (hist_gage%fileOpen()) then
     call hist_gage%cleanup()
     call hist_gage%closeNC()
   end if
 END SUBROUTINE

 ! *********************************************************************
 ! public subroutine: initialize history file for contiue run
 ! *********************************************************************
 SUBROUTINE init_histFile(ierr, message)

   USE globalData,  ONLY: nHRU, nRch        ! number of HRUs and river reaches
   USE globalData,  ONLY: gage_data
   USE public_var,  ONLY: outputAtGage      ! ascii containing last restart and history files
   USE globalData,  ONLY: pioSystem

   implicit none
   ! argument variables
   integer(i4b),   intent(out)     :: ierr                ! error code
   character(*),   intent(out)     :: message             ! error message
   ! local variables
   integer(i4b), allocatable       :: compdof_rch(:)      !
   integer(i4b), allocatable       :: compdof_rch_gage(:) !
   integer(i4b), allocatable       :: compdof_hru(:)      !
   character(len=strLen)           :: cmessage            ! error message of downwind routine

   ierr=0; message='init_histFile/'

   ! get history file names to append
   call io_rpfile('r', ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! initialize and create history netcdfs
   select case(trim(runMode))
     case('standalone');    hist_all_network = histFile(hfileout)
     case('cesm-coupling'); hist_all_network = histFile(hfileout, pioSys=pioSystem)
     case default; ierr=20; message=trim(message)//'rumMode: '//trim(runMode)//': must be standalone or cesm-coupling'; return
   end select

   call hist_all_network%openNc(ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call get_compdof_all_network(compdof_rch, compdof_hru)

   call hist_all_network%set_compdof(compdof_rch, compdof_hru, nRch, nHRU)

   if (outputAtGage) then
     select case(trim(runMode))
       case('standalone');    hist_gage = histFile(hfileout_gage, gageOutput=.true.)
       case('cesm-coupling'); hist_gage = histFile(hfileout_gage, pioSys=pioSystem, gageOutput=.true.)
       case default; ierr=20; message=trim(message)//'rumMode: '//trim(runMode)//': must be standalone or cesm-coupling'; return
     end select

     call hist_gage%openNC(ierr, message)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     call get_compdof_gage(compdof_rch_gage, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     call hist_gage%set_compdof(compdof_rch_gage, gage_data%nGage)
   end if

 END SUBROUTINE init_histFile

END MODULE write_simoutput_pio
