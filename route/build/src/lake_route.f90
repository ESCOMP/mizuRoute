module lake_route_module

!numeric type
USE nrtype
! data type
USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO        ! Network topology
USE dataTypes, ONLY: RCHPRP         ! Network parameter
! global parameters
USE public_var, ONLY: iulog          ! i/o logical unit number
USE public_var, ONLY: realMissing    ! missing value for real number
USE public_var, ONLY: integerMissing ! missing value for integer number

! privary
implicit none
private

public::lake_route

contains

 ! *********************************************************************
 ! subroutine: perform one segment lake routing
 ! *********************************************************************
 subroutine lake_route(&
                        ! input
                        iEns,       &    ! input: index of runoff ensemble to be processed
                        segIndex,   &    ! input: index of runoff ensemble to be processed
                        ixDesire,   &    ! input: reachID to be checked by on-screen pringing (here reachID can be lake)
                        NETOPO_in,  &    ! input: reach topology data structure
                        RPARAM_in,  &    ! input: reach parameter data strcuture
                        ! inout
                        RCHFLX_out, &    ! inout: reach flux data structure
                        ! output
                        ierr, message)   ! output: error control

 USE public_var, ONLY: dt, lakeWBTol

 implicit none
 ! Input
 INTEGER(I4B), intent(IN)                 :: iEns           ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)                 :: segIndex       ! segment where routing is performed
 INTEGER(I4B), intent(IN)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO), intent(IN),   allocatable :: NETOPO_in(:)   ! River Network topology
 type(RCHPRP),  intent(IN),   allocatable :: RPARAM_in(:)   ! River Network topology
 ! inout
 TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! Output
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables to
 real(dp)                                 :: q_upstream     ! total discharge at top of the reach being processed
 real(dp)                                 :: WB             ! water balance component in the lake
 type(STRFLX), allocatable                :: fluxstate(:)   ! upstream Reach fluxes
 INTEGER(I4B)                             :: nUps           ! number of upstream segment
 INTEGER(I4B)                             :: iUps           ! upstream reach index
 INTEGER(I4B)                             :: iRch_ups       ! index of upstream reach in NETOPO
 INTEGER(I4B)                             :: ntdh           ! number of time steps in IRF
 character(len=strLen)                    :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='lake_route/'

  ! identify number of upstream segments of the lake being processed
  nUps = size(NETOPO_in(segIndex)%UREACHI)

  allocate(fluxstate(nUps), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': fluxstate'; return; endif

  ! loop to get the streamlow of the usptream discharges to the lake
  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      fluxstate(iUps) = RCHFLX_out(iens,iRch_ups)
    end do
  endif

  ! Find out total q draining into a lake
  q_upstream = 0.0_dp
  if(nUps>0)then
   do iUps = 1,nUps
    q_upstream = q_upstream + fluxstate(iUps)%REACH_Q_IRF
   end do
  endif

  ! perform lake routing based on a fixed storage discharge relationship Q=kS
  ! no runoff input is added to the lake; the input are only precipitation and evaporation to the lake

  ! if(NETOPO_in(segIndex)%REACHIX == ixDesire)then   ! uncommnet when the ixDesire is fixed and not -9999
   print*, '------lake-simulation-------- '
   print*, 'node id that is lake .......= ', NETOPO_in(segIndex)%REACHID ! to check the reach id of lake
   print*, 'lake param RATECVA .........= ', RPARAM_in(segIndex)%RATECVA
   print*, 'lake param RATECVB .........= ', RPARAM_in(segIndex)%RATECVB
   print*, 'lake param RATECVC .........= ', RPARAM_in(segIndex)%RATECVC
   print*, 'lake param RATECVD .........= ', RPARAM_in(segIndex)%RATECVD
   print*, 'lake param RATECVE .........= ', RPARAM_in(segIndex)%RATECVE
   print*, 'lake param RATECVF .........= ', RPARAM_in(segIndex)%RATECVF
   print*, 'volume before simulation m3.= ', RCHFLX_out(iens,segIndex)%REACH_VOL(0)
   print*, 'upstream streamflow m3/s ...= ', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
   print*, 'upstream precipitation m3/s.= ', RCHFLX_out(iens,segIndex)%basinprecip
   print*, 'upstream evaporation m3/s ..= ', RCHFLX_out(iens,segIndex)%basinevapo
  ! endif

  RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(0) ! updating storage for current time
  RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) + q_upstream * dt  ! input upstream discharge from m3/s to m3
  RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) + RCHFLX_out(iens,segIndex)%basinprecip * dt ! input lake precipitation
  RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%basinevapo * dt ! output lake evaporaiton
  if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) .LT. 0) then; ! set the lake volume as 0 if it goes negative actual evaporation is not calculated here in case low storage mean low evaporaiton...
    RCHFLX_out(iens,segIndex)%REACH_VOL(1)=0
  endif
  ! A has the dimnesion of 1/s; the flux leaving per second
  RCHFLX_out(iens,segIndex)%REACH_Q_IRF = RPARAM_in(segIndex)%RATECVA * (RCHFLX_out(iens,segIndex)%REACH_VOL(1) ** RPARAM_in(segIndex)%RATECVB)! simplified level pool liner reservoir Q=AS^B
  RCHFLX_out(iens,segIndex)%REACH_Q_IRF = (min(RCHFLX_out(iens,segIndex)%REACH_Q_IRF * dt, RCHFLX_out(iens,segIndex)%REACH_VOL(1)) )/dt! in case is the output volume is more than lake volume
  RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_Q_IRF * dt ! updating the storage
  if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) .LT. 0) then; ! set the lake volume as 0 if it goes negative
    RCHFLX_out(iens,segIndex)%REACH_VOL(1)=0
  endif

  ! calculate water balance (in this water balance we dont have the actual evaporation)
  WB = q_upstream * dt + RCHFLX_out(iens,segIndex)%basinprecip * dt - RCHFLX_out(iens,segIndex)%REACH_Q_IRF * dt &
  - RCHFLX_out(iens,segIndex)%basinevapo * dt - (RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_VOL(0))

  !if(NETOPO_in(segIndex)%REACHIX == ixDesire)then    ! uncommnet when the ixDesire is fixed and not -9999
   print*, 'lake simulated output m3/s .= ', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
   print*, 'volume after simulation m3 .= ', RCHFLX_out(iens,segIndex)%REACH_VOL(1)
   print*, 'water balance error ........= ', WB
  !endif

  if(WB.ge.lakeWBTol)then;
   print*, 'Water balance for lake ID = ', NETOPO_in(segIndex)%REACHID, 'excees the Tolerance'
   cmessage = 'Water balance for lake ID = excees the Tolerance'
   ierr = 1; message=trim(message)//trim(cmessage);
  endif

  ! set the routed flag as .True.
  RCHFLX_out(iEns,segIndex)%isRoute=.True.

  ! pass the current storage for the past time step for the next time step simulation
  RCHFLX_out(iens,segIndex)%REACH_VOL(0) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) !shift on time step back

  ! assign the zero value as lake do not have a QFUTURE_IRF
  if (.not.allocated(RCHFLX_out(iens,segIndex)%QFUTURE_IRF))then
   ntdh = size(NETOPO_in(segIndex)%UH)
   allocate(RCHFLX_out(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX_out(iens,segIndex)%QFUTURE_IRF'; return; endif
   RCHFLX_out(iens,segIndex)%QFUTURE_IRF(1:ntdh) = 0._dp
  end if

 end subroutine lake_route

end module lake_route_module

