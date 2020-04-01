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

 USE public_var, ONLY: dt
 
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
 type(STRFLX), allocatable                :: fluxstate(:)   ! upstream Reach fluxes
 INTEGER(I4B)                             :: nUps           ! number of upstream segment
 INTEGER(I4B)                             :: iUps           ! upstream reach index
 INTEGER(I4B)                             :: iRch_ups       ! index of upstream reach in NETOPO
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
  ! there is no in basin routing for lake; the input to the lake (precipitation) should be added to the storage of the lake
  ! in the current implementation there is no input to the system from nodes that are lake
  ! if the entire network topology is lake the system response will be all zero
  print*, 'volume before simulation = ', RCHFLX_out(iens,segIndex)%REACH_VOL(0)
  print*, 'upstream streamflow = ', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
  print*, 'RATECVA = ', RPARAM_in(segIndex)%RATECVA
  RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(0) ! updating storage for current time
  RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) + q_upstream * dt  ! input upstream discharge  
  RCHFLX_out(iens,segIndex)%REACH_Q_IRF = min ( RPARAM_in(segIndex)%RATECVA * RCHFLX_out(iens,segIndex)%REACH_VOL(1) ** RPARAM_in(segIndex)%RATECVB , RCHFLX_out(iens,segIndex)%REACH_VOL(1)) / dt ! calculate the discharge
  RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_Q_IRF * dt ! updating the storage 

  ! set the routed flag as .True.
  RCHFLX_out(iEns,segIndex)%isRoute=.True.

  ! pass the current storage for the past time step for the next time step simulation
  print*, 'volume after simulation = ', RCHFLX_out(iens,segIndex)%REACH_VOL(1)
  RCHFLX_out(iens,segIndex)%REACH_VOL(0) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) !shift on time step back
  
 
 end subroutine lake_route

end module lake_route_module

