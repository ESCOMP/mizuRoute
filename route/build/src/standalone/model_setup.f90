MODULE model_setup

USE nrtype,    ONLY: i4b,dp,lgt          ! variable types, etc.
USE nrtype,    ONLY: strLen              ! length of characters

USE public_var, ONLY: iulog              ! i/o logical unit number
USE public_var, ONLY: debug
USE public_var, ONLY: integerMissing
USE public_var, ONLY: realMissing
USE public_var, ONLY : charMissing

USE init_model_data,  ONLY: init_ntopo_data
USE init_model_data,  ONLY: init_state_data
USE init_model_data,  ONLY: get_mpi_omp
USE init_model_data,  ONLY: update_time

USE nr_utility_module, ONLY : unique  ! get unique element array
USE nr_utility_module, ONLY : indexx  ! get rank of data value

implicit none

! privacy -- everything private unless declared explicitly
private
public :: init_mpi
public :: init_data

CONTAINS

 ! *********************************************************************
 ! public subroutine: initialize MPI for stand-alone program
 ! *********************************************************************
 SUBROUTINE init_mpi()

  ! Initialize MPI and get OMP thread

  ! shared data used
  USE globalData, ONLY: mpicom_route ! communicator id
  ! subroutines: populate metadata
  USE mpi_mod, ONLY: shr_mpi_init

  implicit none

  ! input:  None
  ! output: None
  ! local variables
  character(len=strLen)       :: message             ! error message

  ! initialize error control
  message='init_mpi/'

  call shr_mpi_init(mpicom_route, message)

  call get_mpi_omp(mpicom_route)

 END SUBROUTINE init_mpi


 ! *********************************************************************
 ! public subroutine: initialize all the data - For stand-alone
 ! *********************************************************************
 SUBROUTINE init_data(pid,           & ! input: proc id
                      nNodes,        & ! input: number of procs
                      comm,          & ! input: communicator
                      ierr, message)   ! output: error control

   USE globalData,  ONLY: isFileOpen             ! file open/close status
   ! external subroutines
   USE mpi_routine, ONLY: pass_global_data       ! mpi globaldata copy to slave proc

   IMPLICIT NONE
   integer(i4b),              intent(in)    :: pid              ! proc id
   integer(i4b),              intent(in)    :: nNodes           ! number of procs
   integer(i4b),              intent(in)    :: comm             ! communicator
   ! output: error control
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local variable
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   ! initialize error control
   ierr=0; message='init_data/'

   ! network topology data initialization
   call init_ntopo_data(pid, nNodes, comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! runoff and mapping data initialization
   call init_runoff_data(pid, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! broadcast public and some global variables
   call pass_global_data(comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! channel state initialization
   call init_state_data(pid, nNodes, comm, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   isFileOpen=.false.

 END SUBROUTINE init_data


 ! ********************************************************************************
 ! public subroutine: initialize runoff, and runoff-mapping data - For stand-alone
 ! ********************************************************************************
 SUBROUTINE init_runoff_data(pid,           & ! input: proc id
                             ierr, message)   ! output: error control

   USE public_var, ONLY: is_remap             ! logical whether or not runnoff needs to be mapped to river network HRU
   USE globalData, ONLY: remap_data           ! runoff mapping data structure
   USE globalData, ONLY: runoff_data          ! runoff data structure

   implicit none
   ! input:
   integer(i4b),              intent(in)    :: pid              ! proc id
   ! output: error control
   integer(i4b),              intent(out)   :: ierr             ! error code
   character(*),              intent(out)   :: message          ! error message
   ! local:
   character(len=strLen)                    :: cmessage         ! error message of downwind routine

   if (pid==0) then
     ! runoff and remap data initialization (TO DO: split runoff and remap initialization)
     call init_runoff(is_remap,        & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                      remap_data,      & ! output: data structure to remap data
                      runoff_data,     & ! output: data structure for runoff
                      ierr, cmessage)    ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! DateTime initialization
     call init_time(runoff_data%ntime, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   end if  ! if processor=0 (root)

 END SUBROUTINE init_runoff_data


 ! *********************************************************************
 ! private subroutine: initialize time data
 ! *********************************************************************
 SUBROUTINE init_time(nTime,     &    ! input: number of time steps
                      ierr, message)  ! output

  ! subroutines:
  USE process_time_module, ONLY: process_time  ! process time information
  USE io_netcdf,           ONLY: get_nc        ! netcdf input
  ! derived datatype
  USE dataTypes, ONLY: time           ! time data type
  ! Shared data
  USE public_var, ONLY: input_dir     ! directory containing input data
  USE public_var, ONLY: fname_qsim    ! simulated runoff netCDF name
  USE public_var, ONLY: vname_time    ! variable name for time
  USE public_var, ONLY: time_units    ! time units (seconds, hours, or days)
  USE public_var, ONLY: simStart      ! date string defining the start of the simulation
  USE public_var, ONLY: simEnd        ! date string defining the end of the simulation
  USE public_var, ONLY: calendar      ! calendar name
  USE globalData, ONLY: timeVar       ! time variables (unit given by runoff data)
  USE globalData, ONLY: iTime         ! time index at simulation time step
  USE globalData, ONLY: convTime2Days ! conversion multipliers for time unit of runoff input to day
  USE globalData, ONLY: refJulday     ! julian day: reference
  USE globalData, ONLY: startJulday   ! julian day: start of routing simulation
  USE globalData, ONLY: endJulday     ! julian day: end of routing simulation
  USE globalData, ONLY: modJulday     ! julian day: at model time step
  USE globalData, ONLY: modTime       ! model time data (yyyy:mm:dd:hh:mm:ss)

  implicit none

  ! input:
  integer(i4b),              intent(in)    :: nTime
  ! output: error control
  integer(i4b),              intent(out)   :: ierr             ! error code
  character(*),              intent(out)   :: message          ! error message
  ! local variable
  character(len=7)                         :: t_unit
  integer(i4b)                             :: ix
  character(len=strLen)                    :: cmessage         ! error message of downwind routine

  ! initialize error control
  ierr=0; message='init_time/'

  ! time initialization
  allocate(timeVar(nTime), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time data
  call get_nc(trim(input_dir)//trim(fname_qsim), vname_time, timeVar, 1, nTime, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the time multiplier needed to convert time to units of days
  t_unit = trim( time_units(1:index(time_units,' ')) )
  select case( trim(t_unit) )
   case('seconds','second','sec','s'); convTime2Days=86400._dp
   case('minutes','minute','min');     convTime2Days=1440._dp
   case('hours','hour','hr','h');      convTime2Days=24._dp
   case('days','day','d');             convTime2Days=1._dp
   case default
     ierr=20; message=trim(message)//'<time_units>= '//trim(t_unit)//': <time_units> must be seconds, minutes, hours or days.'; return
  end select

  ! extract time information from the control information
  call process_time(time_units,    calendar, refJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refJulday]'; return; endif
  call process_time(trim(simStart),calendar, startJulday, ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startJulday]'; return; endif
  call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endJulday]'; return; endif

  ! check that the dates are aligned
  if(endJulday<startJulday) then
    write(cmessage,'(7a)') 'simulation end is before simulation start:', new_line('a'), '<sim_start>= ', trim(simStart), new_line('a'), '<sim_end>= ', trim(simEnd)
    ierr=20; message=trim(message)//trim(cmessage); return
  endif

  ! fast forward time to time index at simStart and save iTime and modJulday
  ! need to convert time unit in timeVar to day
  do ix = 1, nTime
    modJulday = refJulday + timeVar(ix)/convTime2Days
    if( modJulday < startJulday ) cycle
    exit
  enddo
  iTime = ix

  ! initialize previous model time
  !modTime(0:1) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)
  modTime(0) = time(integerMissing, integerMissing, integerMissing, integerMissing, integerMissing, realMissing)

 END SUBROUTINE init_time


 ! *****
 ! private subroutine: get mapping data between runoff hru and river network hru
 ! *********************************************************************
 SUBROUTINE init_runoff(remap_flag,      & ! input:  logical whether or not runnoff needs to be mapped to river network HRU
                        remap_data_in,   & ! output: data structure to remap data
                        runoff_data_in,  & ! output: data structure for runoff
                        ierr, message)     ! output: error control

 USE public_var,  ONLY: ancil_dir              ! name of the ancillary directory
 USE public_var,  ONLY: input_dir              ! name of the runoff input directory
 USE public_var,  ONLY: fname_qsim             ! name of simulated runoff netCDF
 USE public_var,  ONLY: fname_remap            ! name of runoff mapping netCDF name
 USE public_var,  ONLY: calendar               ! name of calendar
 USE public_var,  ONLY: time_units             ! time units
 USE globalData,  ONLY: basinID                ! basin ID
 USE dataTypes,   ONLY: remap                  ! remapping data type
 USE dataTypes,   ONLY: runoff                 ! runoff data type
 USE read_runoff, ONLY: read_runoff_metadata   ! read meta data from runoff data
 USE read_remap,  ONLY: get_remap_data         ! read remap data

 implicit none
 ! data structures
 logical(lgt),      intent(in)      :: remap_flag       ! logical whether or not runnoff needs to be mapped to river network HRU
 type(remap)  , intent(out)         :: remap_data_in    ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
 type(runoff) , intent(out)         :: runoff_data_in   ! runoff for one time step for all HRUs
 ! error control
 integer(i4b), intent(out)          :: ierr             ! error code
 character(*), intent(out)          :: message          ! error message
 ! local variables
 integer(i4b), allocatable          :: unq_qhru_id(:)
 integer(i4b), allocatable          :: unq_idx(:)
 character(len=strLen)              :: cmessage         ! error message from subroutine

 ! initialize error control
 ierr=0; message='init_runoff/'

 ! get runoff metadata
 call read_runoff_metadata(trim(input_dir)//trim(fname_qsim), & ! input: filename
                          runoff_data_in,                     & ! output: runoff data structure
                          time_units, calendar,               & ! output: number of time steps, time units, calendar
                          ierr, cmessage)                       ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 !write(*,*) 'runoff_data_in%nSpace, nTime, trim(time_units) = ', runoff_data_in%nSpace(:), runoff_data_in%nTime, trim(time_units)

 ! need to remap runoff to HRUs
 if (remap_flag) then

   ! get runoff mapping file
   call get_remap_data(trim(ancil_dir)//trim(fname_remap),     & ! input: file name
                       runoff_data_in%nSpace,                  & ! input: number of spatial elements
                       remap_data_in,                          & ! output: data structure to remap data from a polygon
                       ierr, cmessage)                           ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! get indices of the HRU ids in the mapping file in the routing layer
   call get_qix(remap_data_in%hru_id, &  ! input: vector of ids in mapping file
                basinID,              &  ! input: vector of ids in the routing layer
                remap_data_in%hru_ix, &  ! output: indices of hru ids in routing layer
                ierr, cmessage)          ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (debug) then
     write(iulog,'(2a)') new_line('a'), 'DEBUG: Corresponding between River-Network(RN) hru in mapping data and RN hru in river network data'
     write(iulog,'(2x,a,I15)') '(1) number of RN hru in river-network = ', size(basinID)
     write(iulog,'(2x,a,I15)') '(2) number of RN hru in mapping       = ', size(remap_data_in%hru_id)
     write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two  = ', count(remap_data_in%hru_ix/=integerMissing)
     if(count(remap_data_in%hru_ix/=integerMissing)/=size(basinID))then
       message=trim(message)//'(1) not equal (2)'
       ierr=20; return
     endif
   end if

   if ( runoff_data_in%nSpace(2) == integerMissing ) then
     ! get indices of the "overlap HRUs" (the runoff input) in the runoff vector
     call get_qix(remap_data_in%qhru_id, &  ! input: vector of ids in mapping file
                  runoff_data_in%hru_id, &  ! input: vector of ids in runoff file
                  remap_data_in%qhru_ix, &  ! output: indices of mapping ids in runoff file
                  ierr, cmessage)           ! output: error control
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     if (debug) then
       call unique(remap_data_in%qhru_id, unq_qhru_id, unq_idx)
       write(iulog,'(2a)') new_line('a'),'DEBUG: corresponding between Hydro-Model (HM) hru in mapping data and HM hru in runoff data'
       write(iulog,'(2x,a,I15)') '(1) number of HM hru in hyrdo-model  = ', size(runoff_data_in%hru_id)
       write(iulog,'(2x,a,I15)') '(2) number of HM hru in mapping      = ', size(unq_qhru_id)
       write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two = ', count(remap_data_in%qhru_ix(unq_idx)/=integerMissing)
     end if
   end if

 else ! if runoff given in RN_HRU

   allocate(runoff_data_in%hru_ix(size(runoff_data_in%hru_id)), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating runoff_data_in%hru_ix'; return; endif

   ! get indices of the HRU ids in the runoff file in the routing layer
   call get_qix(runoff_data_in%hru_id,  &    ! input: vector of ids in mapping file
                basinID,                &    ! input: vector of ids in the routing layer
                runoff_data_in%hru_ix,  &    ! output: indices of hru ids in routing layer
                ierr, cmessage)              ! output: error control
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   if (debug) then
     write(iulog,'(2a)') new_line('a'), 'DEBUG: corresponding between River-Network (RN) hru in runoff data and RN hru in river network data'
     write(iulog,'(2x,a,I15)') '(1) number of RN hru in river-network = ', size(basinID)
     write(iulog,'(2x,a,I15)') '(2) number of RN hru in hyrdo-model   = ', size(runoff_data_in%hru_id)
     write(iulog,'(2x,a,I15)') '(3) number of mapped hru between two  = ', count(runoff_data_in%hru_ix/=integerMissing)
     if(count(runoff_data_in%hru_ix/=integerMissing)/=size(basinID))then
       message=trim(message)//'(1) not equal (2)'
       ierr=20; return
     endif
   end if

 endif

 END SUBROUTINE init_runoff


 ! *****
 ! private subroutine: get indices of mapping points within runoff file...
 ! ***********************************************************************
 SUBROUTINE get_qix(qid,qidMaster,qix,ierr,message)

 implicit none
 ! input
 integer(i4b), intent(in)  :: qid(:)                       ! ID of input vector
 integer(i4b), intent(in)  :: qidMaster(:)                 ! ID of master vector
 ! output
 integer(i4b), intent(out) :: qix(:)                       ! index within master vector
 integer(i4b), intent(out) :: ierr                         ! error code
 character(*), intent(out) :: message                      ! error message
 ! local
 integer(i4b)             :: rankID( size(qid) )           ! rank of input vector
 integer(i4b)             :: rankMaster( size(qidMaster) ) ! rank of master vector
 integer(i4b)             :: ix,jx,ixMaster                ! array indices
 integer(i4b)             :: nx                            ! counter

 ! initialize error control
 ierr=0; message='get_qix/'

 ! sort the data vector from smallest to largest
 call indexx(qid,       rankID)
 call indexx(qidMaster, rankMaster)

 qix(1:size(qid)) = integerMissing
 nx=0
 jx=1
 ! loop through id vector
 do ix=1,size(qid)

  ! find match
  do ixMaster=jx,size(qidMaster) ! normally a very short loop

   ! keep track of trials
   nx=nx+1
   !write(*,*) 'qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) ) = ', qid( rankId(ix) ), qidMaster( rankMaster(ixMaster) )

   ! find match
   if( qid( rankId(ix) ) == qidMaster( rankMaster(ixMaster) ) )then
    qix( rankId(ix) ) = rankMaster(ixMaster)
    jx = ixMaster
    exit
   endif

   ! unable to find match
   if( qidMaster( rankMaster(ixMaster) ) > qid( rankId(ix) ) )then
    qix( rankId(ix) ) = integerMissing
    jx = ixMaster
    exit
   endif

  end do  ! ixMaster

  ! print progress
  if(qix( rankId(ix) )/=integerMissing .and. mod(ix,1000000)==0)then
   write(iulog,*) trim(message)//'matching ids: ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) ) = ', &
                                                ix, qix( rankId(ix) ), qid( rankId(ix) ), qidMaster( qix( rankId(ix) ) )
  endif

 end do  ! looping through the vector

 ! check again
 do ix=1,size(qid)
  if(qix(ix) /= integerMissing)then
   if(qid(ix) /= qidMaster( qix(ix) ) )then
    write(iulog,'(a,2(x,I10,x,I15))') 'ERROR Mapping: ix, qid(ix), qix(ix), qidMaster(qix(ix))=', ix, qid(ix), qix(ix), qidMaster(qix(ix))
    message=trim(message)//'unable to find the match'
    ierr=20; return
   endif
  endif
 end do

 END SUBROUTINE get_qix

END MODULE model_setup
