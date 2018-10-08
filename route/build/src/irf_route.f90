module irf_route_module

! global parameters
USE nrtype
USE dataTypes,  only : STRFLX        ! fluxes in each reach

! privary
implicit none
private

public::irf_route
public::make_uh

contains

 ! *********************************************************************
 ! subroutine: perform network UH routing
 ! *********************************************************************
 subroutine irf_route(&
                      ! input
                      iEns,       &    ! input: index of runoff ensemble to be processed
                      nRch,       &    ! input: index of reach to be processed
                      river_basin,&    ! input: river basin information (mainstem, tributary outlet etc.)
                      ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                      ! output
                      ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Convolute routed basisn flow volume at top of each of the upstream segment at one time step and at each segment
 !
 ! ----------------------------------------------------------------------------------------

 ! global routing data
 USE globalData, only : RCHFLX ! routing fluxes
 USE globalData, only : NETOPO ! Network topology
 USE dataTypes,  only : basin  ! river basin data type

 implicit none
 ! Input
 INTEGER(I4B), intent(IN)               :: iEns           ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)               :: nRch           ! number of reach segments in the network
 type(basin),  intent(in), allocatable  :: river_basin(:) ! river basin information (mainstem, tributary outlet etc.)
 INTEGER(I4B), intent(IN)               :: ixDesire       ! index of the reach for verbose output
 ! Output
 integer(i4b), intent(out)              :: ierr           ! error code
 character(*), intent(out)              :: message        ! error message
 ! Local variables to
 INTEGER(I4B)                           :: iRch           ! reach segment index
 INTEGER(I4B)                           :: jRch           ! reach segment to be routed
 character(len=strLen)                  :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='irf_route/'

 ! Initialize CHEC_IRF to False.
 RCHFLX(iEns,:)%CHECK_IRF=.False.

 ! route streamflow through the river network
 do iRch=1,nRch

  jRch = NETOPO(iRch)%RHORDER

  call segment_irf(iEns, jRch, ixDesire, ierr, message)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end do

 end subroutine irf_route


 ! *********************************************************************
 ! subroutine: perform one segment route UH routing
 ! *********************************************************************
 subroutine segment_irf(&
                        ! input
                        iEns,       &    ! input: index of runoff ensemble to be processed
                        segIndex,   &    ! input: index of runoff ensemble to be processed
                        ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                        ! output
                        ierr, message)   ! output: error control

 ! global routing data
 USE globalData, only : RCHFLX ! routing fluxes
 USE globalData, only : NETOPO ! Network topology

 implicit none
 ! Input
 INTEGER(I4B), intent(IN)               :: iEns           ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)               :: segIndex       ! segment where routing is performed
 INTEGER(I4B), intent(IN)               :: ixDesire       ! index of the reach for verbose output
 ! Output
 integer(i4b), intent(out)              :: ierr           ! error code
 character(*), intent(out)              :: message        ! error message
 ! Local variables to
 type(STRFLX), allocatable              :: uprflux(:)     ! upstream Reach fluxes
 INTEGER(I4B)                           :: nUps           ! number of upstream segment
 INTEGER(I4B)                           :: iUps           ! upstream reach index
 INTEGER(I4B)                           :: iRch_ups       ! index of upstream reach in NETOPO
 INTEGER(I4B)                           :: ntdh           ! number of time steps in IRF
 character(len=strLen)                  :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='segment_irf/'

 ! route streamflow through the river network
  if (.not.allocated(RCHFLX(iens,segIndex)%QFUTURE_IRF))then

   ntdh = size(NETOPO(segIndex)%UH)

   allocate(RCHFLX(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX(iens,segIndex)%QFUTURE_IRF'; return; endif

   RCHFLX(iens,segIndex)%QFUTURE_IRF(:) = 0._dp

  end if

  ! identify number of upstream segments of the reach being processed
  nUps = size(NETOPO(segIndex)%UREACHI)

  allocate(uprflux(nUps), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': uprflux'; return; endif

  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      uprflux(iUps) = RCHFLX(iens,iRch_ups)
    end do
  endif

  ! perform river network UH routing
  call conv_upsbas_qr(NETOPO(segIndex)%UH,   &    ! input: reach unit hydrograph
                      uprflux,               &    ! input: upstream reach fluxes
                      RCHFLX(iens,segIndex), &    ! inout: updated fluxes at reach
                      ierr, message)              ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! Check True since now this reach now routed
  RCHFLX(iEns,segIndex)%CHECK_IRF=.True.

  ! check
  if(NETOPO(segIndex)%REACHIX == ixDesire)then
   print*, 'RCHFLX(iens,segIndex)%BASIN_QR(1),RCHFLX(iens,segIndex)%REACH_Q_IRF = ', &
            RCHFLX(iens,segIndex)%BASIN_QR(1),RCHFLX(iens,segIndex)%REACH_Q_IRF
  endif

 end subroutine segment_irf



! *********************************************************************
! subroutine: compute normalized UH from Saint-Venant Eq. at sim. time step
!             for all the upstream segment
! *********************************************************************
 subroutine make_uh(&
                    ! Input
                    length,     &       ! input: river segment array [meter]
                    dt,         &       ! input: time step interval [sec]
                    velo,       &       ! input: IRF parameter 1 - celerity C for each stream segment [m/s]
                    diff,       &       ! input: IRF parameter 2 - diffusivity D for each stream segment [m^2/s]
                    ! Output
                    seg_uh,     &       ! output: unit hydrograph ordinates for a given segment length array
                    ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Calculate UH at given simulation time step, C, and D for each stream segment
 !   based on Saint-Venant equation - see equations (13)-(15) in Lohmann, et al. (1996) Tellus
 !   IRF and UH are used interchangebly.
 !
 ! ----------------------------------------------------------------------------------------
  ! global variables
  USE globalData, only : pi      ! pi
  USE dataTypes,  only : dlength
  implicit none
  ! input variables
  real(dp),                      intent(in)   :: length(:)     ! river segment length
  real(dp),                      intent(in)   :: dt            ! Time step interval [sec]
  real(dp),                      intent(in)   :: velo          ! Wave velocity C for each segment [m/sec]
  real(dp),                      intent(in)   :: diff          ! Diffusivity D for each segment [m2/sec]
  ! output variables
  type(dlength),allocatable,     intent(out)  :: seg_uh(:)     ! unit hydrograph ordinates for each segment
  integer(i4b),                  intent(out)  :: ierr          ! error code
  character(*),                  intent(out)  :: message       ! error message
  ! local variables
  real(dp), allocatable                       :: UHM(:)        ! Unit hydrograph derived from S-V equation
  real(dp), allocatable                       :: UHQ(:)        !
  real(dp), allocatable                       :: fr(:)         ! Unit runoff depth evenly distributed over the simulation duration at hourly step
  real(dp)                                    :: seg_length    ! length of rech segment
  real(dp)                                    :: POT           ! Inside of exponential function of IRF
  real(dp)                                    :: H             ! IRF, h(x,t) function h function in eq.15 Lohmann et al.1996
  real(dp)                                    :: UHQ0          !
  real(dp)                                    :: INTE          ! Accumulation (integration) of UH
  real(dp)                                    :: sec           ! hr at each time step  [hr]
  real(dp),parameter                          :: dTUH=3600.0   ! UH time step [sec]
  integer(i4b)                                :: nSeg          ! Numer of segment
  integer(i4b)                                :: iHrStrt       ! index of UH time step where rising limb of UH start
  integer(i4b)                                :: iHrLast       ! index of UH time step where recession limb of UH become zero
  integer(i4b)                                :: nTSub         ! number of time steps where 1/nTsub [m] of runoff is inserted
  integer(i4b)                                :: iSeg          ! Loop index
  integer(i4b)                                :: iHr,jHr       ! Loop index of hour
  integer(i4b)                                :: iTagg         ! index for aggregated (i.e. simulation) time step
  integer(i4b),parameter                      :: nTMAX=240     ! Maximum hour of UH [hr] - 10 days times 24hrs
  integer(i4b),parameter                      :: nHr=240       ! Maximum hour of UH [hr] - 10 days times 24hrs
  character(len=strLen)                       :: cmessage      ! error message from subroutine

 ! initialize error control
 ierr=0; message='make_uh/'

 ! Dynamically assigned parameters
 nTsub=ceiling(dt/dTUH)
 !nTsub=floor(dt/dTUH)
 nSeg = size(length)

 ! Memory allocation
 allocate(seg_uh(nSeg), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': seg_uh'; return; endif
 allocate(fr(nTMAX), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': fr'; return; endif
 allocate(UHQ(nTMAX),stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': UHQ'; return; endif
 allocate(UHM(nHr),stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': UHM'; return; endif

 ! Generate Unit Runoff depth [1/nTsub, 1/nTsub, ..., 1/nTsub, 0, 0, ...0]
 fr=0._dp
 fr(1:nTsub) = 1.0_dp / nTsub

 ! starting segment loop
 do iSeg=1,nSeg

  !Compute S-V based UH, UHM
  seg_length = length(iSeg)
  !seg_length = RPARAM(iSeg)%RLENGTH

  INTE = 0._dp
  sec = 0._dp
  UHM(:) = 0._dp
  do iHr=1,nHr
    sec = sec + dTUH
    if (velo .GT. 0.0) then
      POT = ((velo*sec-seg_length)**2.0)/(4.0_dp*diff*sec)
      if (POT .GT. 69.0_dp) then
        H = 0.0_dp
      else
        H = 1.0_dp/(2.0_dp*sqrt(pi*diff*sec))*seg_length*exp(-POT)
      endif
    else
      H = 0.0_dp
    endif
    UHM(iHr) = H
    INTE = INTE + H
  enddo

  !Normalize ordinates of S-V UH by its sum
  if (INTE > 0.0) then
    UHM= UHM/INTE
  end if

  !Get hour indices of start of rising part of UHM and end of recession part of UHM
  INTE = 0._dp
  do iHr=1,nTMAX
    INTE= INTE+UHM(iHr)
    iHrLast=iHr
    if (INTE > 0.99999) exit
  enddo
  INTE = 0._dp
  do iHr=nTMAX,1,-1
    INTE= INTE+UHM(iHr)
    iHrStrt=iHr
    if (INTE > 0.99999) exit
  enddo

  !Convolute with unit runoff depth, fr
  INTE = 0._dp
  UHQ(:) = 0._dp
  do jHr = 1,nTMAX
    UHQ0=0._dp
    do iHr = iHrStrt,iHrLast
      if ((jHr-iHr) > 0) then
        if (jHr-iHr <= nTsub ) then
          UHQ0 = UHQ0 + fr(jHr-iHr)*UHM(iHr)
        endif
      else
        exit
      endif
    enddo
    UHQ(jHr) = UHQ0
    INTE = INTE+UHQ0
  end do

  ! Normalize ordinates of UHQ by its sum
  if (INTE > 0.0_dp) then
    UHQ = UHQ/INTE
  end if

  !Get hour indices of end of recession part of UH
  INTE = 0._dp
  do iHr=1,nTMAX
    INTE= INTE+UHQ(iHr)
    iHrLast=iHr
    if (INTE > 0.9999) exit
  enddo

  !Aggregate hourly unit hydrograph to simulation time step
  allocate(seg_uh(iSeg)%dat((iHrLast+nTsub-1)/nTsub),stat=ierr,errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': seg_uh%dat'; return; endif

  seg_uh(iSeg)%dat(:)=0._dp
  do jHr = 1,iHrLast
    iTagg = (jHr+nTsub-1)/nTsub
    seg_uh(iSeg)%dat(iTagg) = seg_uh(iSeg)%dat(iTagg) + UHQ(jHr)
  end do

 end do ! sSeg loop

 end subroutine make_uh


 ! *********************************************************************
 ! subroutine: Compute delayed runoff from the upstream segments
 ! *********************************************************************
 subroutine conv_upsbas_qr(&
                           ! input
                           reach_uh,   &    ! input: reach unit hydrograph
                           rflux_ups,  &    ! input: upstream reach fluxes
                           rflux,      &    ! input: input flux at reach
                           ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Convolute routed basisn flow volume at top of each of the upstream segment at one time step and at each segment
 !
 ! ----------------------------------------------------------------------------------------

 implicit none
 ! Input
 real(dp),     intent(in)               :: reach_uh(:)  ! reach unit hydrograph
 type(STRFLX), intent(in)               :: rflux_ups(:) ! upstream Reach fluxes
 ! inout
 type(STRFLX), intent(inout)            :: rflux        ! current Reach fluxes
 ! Output
 integer(i4b), intent(out)              :: ierr         ! error code
 character(*), intent(out)              :: message      ! error message
 ! Local variables to
 real(dp)                               :: q_upstream   ! total discharge at top of the reach being processed
 INTEGER(I4B)                           :: nTDH         ! number of UH data
 INTEGER(I4B)                           :: iTDH         ! index of UH data (i.e.,future time step)
 INTEGER(I4B)                           :: nUps         ! number of all upstream segment
 INTEGER(I4B)                           :: iUps         ! loop indices for u/s reaches

 ! initialize error control
 ierr=0; message='conv_upsbas_qr/'

 ! identify number of upstream segments of the reach being processed
 nUps = size(rflux_ups)

 q_upstream = 0.0_dp
 if(nUps>0)then
   do iUps = 1,nUps
     ! Find out total q at top of a segment
     q_upstream = q_upstream + rflux_ups(iUps)%REACH_Q_IRF
   end do
 endif

 ! place a fraction of runoff in future time steps
 nTDH = size(reach_uh) ! identify the number of future time steps of UH for a given segment
 do iTDH=1,nTDH
   rflux%QFUTURE_IRF(iTDH) = rflux%QFUTURE_IRF(iTDH) &
                             + reach_uh(iTDH)*q_upstream
 enddo

 ! Add local routed flow
 rflux%REACH_Q_IRF = rflux%QFUTURE_IRF(1) + rflux%BASIN_QR(1)

 ! move array back   use eoshift
 rflux%QFUTURE_IRF=eoshift(rflux%QFUTURE_IRF,shift=1)

 end subroutine conv_upsbas_qr

end module irf_route_module

