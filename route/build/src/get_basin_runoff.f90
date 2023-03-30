MODULE get_runoff

USE nrtype

implicit none

private
public::get_hru_runoff

CONTAINS

 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 SUBROUTINE get_hru_runoff(ierr, message)

 ! populate runoff_data with runoff values at LSM domain and at simulation time step

  USE dataTypes,   ONLY: map_time                ! data type for time-step mapping between two time series
  USE public_var,  ONLY: input_dir               ! directory containing input data
  USE public_var,  ONLY: fname_qsim              ! simulated runoff netCDF name
  USE public_var,  ONLY: is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE public_var,  ONLY: qmodOption              ! options for streamflow modification (DA)
  USE public_var,  ONLY: takeWater               ! switch for water abstraction/injection
  USE public_var,  ONLY: integerMissing          !
  USE globalData,  ONLY: runoff_data             ! data structure to hru runoff data
  USE globalData,  ONLY: remap_data              ! data structure to remap data
  USE globalData,  ONLY: iTime                   !
  USE globalData,  ONLY: simDatetime             ! current model time data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,  ONLY: begDatetime             ! forcing data start datetime data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,  ONLY: roBegDatetime           ! forcing data start datetime data (yyyy:mm:dd:hh:mm:sec)
  USE globalData,  ONLY: gage_obs_data           !
  USE globalData,  ONLY: rch_qtake_data          !
  USE globalData,  ONLY: RCHFLX                  !
  USE read_runoff, ONLY: read_runoff_data        ! read runoff value into runoff_data data strucuture
  USE remapping,   ONLY: remap_runoff            ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE remapping,   ONLY: sort_runoff             ! mapping HM runoff to river network HRU runoff (HM_HRU == RN_HRU)

  implicit none
  ! argument variables
  integer(i4b), intent(out)     :: ierr               ! error code
  character(*), intent(out)     :: message            ! error message
  ! local variables
  type(map_time)                :: tmap_sim_ro        ! time-steps mapping data
  character(len=strLen)         :: cmessage           ! error message from subroutine
  integer(i4b)                  :: iens=1
  integer(i4b)                  :: ix, jx
  real(dp)                      :: qobs               ! elapsed time for the process
  integer(i4b), allocatable     :: reach_ix(:)
  integer(i4b), parameter       :: no_mod=0
  integer(i4b), parameter       :: direct_insert=1
  ! timing
!  integer*8                     :: startTime,endTime,cr ! star and end time stamp, rate
!  real(dp)                      :: elapsedTime          ! elapsed time for the process

  ierr=0; message='get_hru_runoff/'

!  call system_clock(count_rate=cr)

  ! time step mapping from runoff time step to simulation time step
  call timeMap_sim_forc(tmap_sim_ro,        &
                        begDatetime,        &
                        roBegDatetime,      &
                        runoff_data%nTime,  &
                        iTime,              &
                        ierr, cmessage)
  if(ierr/=0) then; message=trim(message)//trim(cmessage); return; endif

!  call system_clock(startTime)
  ! get the simulated runoff for the current time step - runoff_data%qsim(:) or %qsim2D(:,:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        tmap_sim_ro,                       & ! input: ro-sim time mapping at current simulation step
                        runoff_data,                       & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!  call system_clock(endTime)
!  elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
!  write(*,"(A,1PG15.7,A)") '  elapsed-time [runoff_input/read] = ', elapsedTime, ' s'

!  call system_clock(startTime)
  ! Get river network HRU runoff into runoff_data data structure
  if (is_remap) then ! remap LSM simulated runoff to the HRUs in the river network
    call remap_runoff(runoff_data, remap_data, runoff_data%basinRunoff, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  else ! runoff is already remapped to river network HRUs
    call sort_runoff(runoff_data, runoff_data%basinRunoff, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
!  call system_clock(endTime)
!  elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
!  write(*,"(A,1PG15.7,A)") '  elapsed-time [runoff_input/remap] = ', elapsedTime, ' s'

  select case(qmodOption)
    case(no_mod) ! do nothing
    case(direct_insert)
      ! read gage observation [m3/s] at current time
      jx = gage_obs_data%time_ix(simDatetime(1))

      if (jx/=integerMissing) then ! there is observation at this time step
        call gage_obs_data%read_obs(ierr, cmessage, index_time=jx)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! put qmod at right reach
        reach_ix = gage_obs_data%link_ix()
        do ix=1,size(reach_ix)
          if (reach_ix(ix)==integerMissing) cycle

          qobs = gage_obs_data%get_obs(tix=1, six=ix)

          if (isnan(qobs) .or. qobs<0) cycle
          RCHFLX(iens,reach_ix(ix))%QOBS = qobs
          RCHFLX(iens,reach_ix(ix))%Qelapsed = 0
        end do
      else
        RCHFLX(iens,:)%Qelapsed = RCHFLX(iens,:)%Qelapsed + 1 ! change only gauge point
      end if
    case default
      ierr=1; message=trim(message)//"Error: qmodOption invalid"; return
  end select

  ! initialize TAKE for water abstract/injection
  RCHFLX(:,:)%TAKE = 0.0_dp
  if (takeWater) then
    ! read reach water take [m3/s] at current time
    jx = rch_qtake_data%time_ix(simDatetime(1))

    if (jx/=integerMissing) then
      call rch_qtake_data%read_obs(ierr, cmessage, index_time=jx)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! put qmod at right reach
      reach_ix = rch_qtake_data%link_ix()
      do ix=1,size(reach_ix)
        if (reach_ix(ix)==integerMissing) cycle
        RCHFLX(iens,reach_ix(ix))%TAKE = rch_qtake_data%get_obs(tix=1, six=ix)
      end do
    end if
  end if

 END SUBROUTINE get_hru_runoff

 ! *********************************************************************
 ! private subroutine: initialize simulation time
 ! *********************************************************************
 SUBROUTINE timeMap_sim_forc(tmap_sim_forc,     & ! inout: time-steps mapping data
                             startSim, startRo, & ! input: starting datetime for simulation and input
                             nRo,               & ! input: numbers of runoff time steps
                             ixTime,            & ! input: simulation time index
                             ierr, message)

   USE nr_utility_module, ONLY: arth
   USE datetime_data,     ONLY: datetime       ! time data type
   USE dataTypes,         ONLY: map_time       ! data type for time-step mapping between two time series
   USE public_var,        ONLY: verySmall      ! smallest real values
   USE public_var,        ONLY: secprday       ! day to second conversion factor
   USE public_var,        ONLY: calendar       ! calender
   USE public_var,        ONLY: dt_sim         ! simulation time step
   USE public_var,        ONLY: dt_ro          ! input time step

   implicit none
   ! Argument variables
   type(map_time),              intent(out)   :: tmap_sim_forc    ! time-steps mapping data
   type(datetime),              intent(in)    :: startSim         ! simulation start datetime
   type(datetime),              intent(in)    :: startRo          ! runoff data start datetime
   integer(i4b),                intent(in)    :: nRo              ! number of runoff data time steps
   integer(i4b),                intent(in)    :: ixTime           ! number of simulation time steps
   integer(i4b),                intent(out)   :: ierr             ! error code
   character(*),                intent(out)   :: message          ! error message
   ! Local variables
   character(len=strLen)                      :: cmessage         ! error message from subroutine
   real(dp)                                   :: frcLapse(nRo+1)  !
   real(dp)                                   :: simLapse(2)      !
   real(dp)                                   :: startRoSec       ! starting runoff time relative to starting simulation time [sec]
   real(dp)                                   :: juldayRo         ! starting julian day in runoff time step [day]
   real(dp)                                   :: juldaySim        ! starting julian day in simulation time step [day]
   integer(i4b)                               :: ctr              ! counter
   integer(i4b)                               :: nRoSub           !
   integer(i4b)                               :: iRo              ! loop index of runoff time step
   integer(i4b)                               :: idxFront         ! index of r of which top is within ith model layer (i=1..nLyr)
   integer(i4b)                               :: idxEnd           ! index of the lowest soil layer of which bottom is within ith model layer (i=1..nLyr)

   ierr=0; message='timeMap_sim_forc/'

   call startRo%julianday(calendar, juldayRo, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   call startSim%julianday(calendar,juldaySim,  ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   startRoSec = (juldaySim - juldayRo)*secprday ! day->sec

   simLapse(1) = startRoSec+dt_sim*(ixTime-1)
   simLapse(2) = simLapse(1) + dt_sim
   frcLapse    = arth(0._dp,      dt_ro, nRo+1)

   !-- Find index of runoff time period of which end is within simulation period
   ! condition: from beginning to end of runoff time, first index of time period whose end exceeds the beginning of simulation time step
   do iRo = 1,nRo
     if (simLapse(1)<frcLapse(iRo+1)) then
       idxFront = iRo; exit
     end if
   end do
   !-- Find index of runoff time period of which beginning is within simulation period
   ! condition: from beginning to end of runoff time, first index of time period whose end exceeds the beginning of simulation time step
   do iRo = 1,nRo
     if ( simLapse(2)<frcLapse(iRo+1) .or. abs(simLapse(2)-frcLapse(iRo+1))<verySmall ) then
       idxEnd = iRo; exit
     end if
   end do

   ! Error check
   if (idxFront-idxEnd>0)then;ierr=30;message=trim(message)//'index of idxFront lower than idxEnd';return;endif

   !-- Compute weight of soil layer contributing to each model layer and populate lyrmap variable
   ctr = 1
   nRoSub = idxEnd-idxFront + 1

   allocate(tmap_sim_forc%iTime(nRoSub), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//'tmap_sim_forc%iTime'; return; endif

   if (idxFront == idxEnd)then ! if simulation period is completely within runoff period - one runoff time period per simulation period
     tmap_sim_forc%iTime(ctr) = idxFront
   else
     allocate(tmap_sim_forc%frac(nRoSub), stat=ierr, errmsg=cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage)//'tmap_sim_forc%frac'; return; endif

     ! loop frm the ealiest runoff time index to the latest, and compute fraction
     do iRo=idxFront, idxEnd
       tmap_sim_forc%iTime(ctr) = iRo
       if (iRo == idxFront)then        ! front side of simulation time step
         tmap_sim_forc%frac(ctr)   = (frcLapse(iRo+1) - simLapse(1))/dt_sim
       else if ( iRo == idxEnd ) then  ! end side of simulation time step
         tmap_sim_forc%frac(ctr)   = (simLapse(2)-frcLapse(iRo))/dt_sim
       else                            ! simulation time step is completely with forcing time step
         tmap_sim_forc%frac(ctr)   = (frcLapse(iRo+1)-frcLapse(iRo))/dt_sim
       endif
       ctr = ctr+1
     end do
   end if

 END SUBROUTINE timeMap_sim_forc

END MODULE get_runoff
