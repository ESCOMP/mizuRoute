module lake_route_module

!numeric type
USE nrtype
! data type
USE dataTypes, ONLY: LAKSTR         ! storage in each lake/reservoir *** the above subroutines should be adjusted to padd lake storage***
USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO        ! Network topology
USE dataTypes, ONLY: LAKPRP         ! lake parameter *** the above subroutine should be ajusted to pass lake paraemters ***
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
                        LAKEPR, &        ! input: lake param, should be pass to the subroutine from higher level models
                        ! inout
                        RCHFLX_out, &    ! inout: reach flux data structure
                        Lake_storage, &  ! inout: storage of each lake (river segment with have reservoir storage of zero)
                        ! output
                        ierr, message)   ! output: error control

 implicit none
 ! Input
 INTEGER(I4B), intent(IN)                 :: iEns           ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)                 :: segIndex       ! segment where routing is performed
 INTEGER(I4B), intent(IN)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO), intent(in),   allocatable :: NETOPO_in(:)   ! River Network topology
 type(LAKPRP),  intent(in),   allocatable :: LAKEPR         ! lake parameters for the given lake
 ! inout
 TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 TYPE(LAKSTR), intent(inout), allocatable :: Lake_storage(:,:) ! lake storages (ensembles, space [reaches/which include lakes])
 ! Output
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables to
 real(dp)                                 :: q_upstream     ! total discharge at top of the reach being processed
 type(STRFLX), allocatable                :: uprflux(:)     ! upstream Reach fluxes
 INTEGER(I4B)                             :: nUps           ! number of upstream segment
 INTEGER(I4B)                             :: iUps           ! upstream reach index
 INTEGER(I4B)                             :: iRch_ups       ! index of upstream reach in NETOPO
 INTEGER(I4B)                             :: ntdh           ! number of time steps in IRF
 character(len=strLen)                    :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='lake_route/'

 ! route streamflow through the river network
  if (.not.allocated(RCHFLX_out(iens,segIndex)%QFUTURE_IRF))then

   ntdh = size(NETOPO_in(segIndex)%UH)

   allocate(RCHFLX_out(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX_out(iens,segIndex)%QFUTURE_IRF'; return; endif

   RCHFLX_out(iens,segIndex)%QFUTURE_IRF(:) = 0._dp

  end if

  ! identify number of upstream segments of the lake being processed
  nUps = size(NETOPO_in(segIndex)%UREACHI)

  allocate(uprflux(nUps), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': uprflux'; return; endif

  ! loop to get the streamlow of the usptream discharges to the lake
  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      uprflux(iUps) = RCHFLX_out(iens,iRch_ups)
    end do
  endif

  ! Find out total q draining into a lake
  q_upstream = 0.0_dp
  if(nUps>0)then
   do iUps = 1,nUps
    q_upstream = q_upstream + uprflux(iUps)%REACH_Q_IRF
   end do
  endif

  ! perform lake routing based on a fixed storage discharge relationship (can be changes later)
  ! updating storage of the lake
  Lake_storage(iens,segIndex) = Lake_storage(iens,segIndex)+q_upstream * dt ! assume this dt is global and in sec?

  ! perform storage discharge relationship for a fixed sotrage discharge relashinship for all the lakes in this version
  ! at this moment we have explicit but later we can do a an iteretive solution of the storage and re
  RCHFLX_out(iens,segIndex) = Lake_storage(iens,segIndex)*0.01 ! here it is a simplified version of level pool I guess basically liner reservoir
  Lake_storage(iens,segIndex) = Lake_storage(iens,segIndex) - RCHFLX_out(iens,segIndex) * dt ! here it is a simplified version of level pool I guess basically liner reservoir
  ! later we should be able to pass the lake parameter for level pool I assume LAKPRP%RATECVA and LAKPRP%RATECVB

  ! se the routed flag as .True.
  RCHFLX_out(iEns,segIndex)%isRoute=.True.

  ! check I DONT KNOW HOW TO RIGHT THIS OUTPUT TO THE MODEL SIMULATION
  if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   write(iulog,*) 'RCHFLX_out(iens,segIndex)%BASIN_QR(1),RCHFLX_out(iens,segIndex)%REACH_Q_IRF = ', &
                   RCHFLX_out(iens,segIndex)%BASIN_QR(1),RCHFLX_out(iens,segIndex)%REACH_Q_IRF
  endif

 end subroutine lake_route

end module lake_route_module

