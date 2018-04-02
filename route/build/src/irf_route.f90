module irf_route_module

USE nrtype

! global parameters
!USE globalData, only : RPARAM ! Reach parameters

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
 implicit none
 ! Input
 INTEGER(I4B), intent(IN)               :: iEns        ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)               :: nRch        ! number of reach segments in the network
 INTEGER(I4B), intent(IN)               :: ixDesire    ! index of the reach for verbose output
 ! Output
 integer(i4b), intent(out)              :: ierr        ! error code
 character(*), intent(out)              :: message     ! error message
 ! Local variables to
 INTEGER(I4B)                           :: ntdh        ! number of time steps in IRF
 INTEGER(I4B)                           :: iRch        ! reach segment index
 INTEGER(I4B)                           :: jRch        ! reach segment to be routed
 character(len=strLen)                  :: cmessage    ! error message from subroutine

 ! initialize error control
 ierr=0; message='irf_route/'

 ! Initialize CHEC_IRF to False.
 RCHFLX(iEns,:)%CHECK_IRF=.False.

 ! route streamflow through the river network
 do iRch=1,nRch

  ! identify reach to process
  jRch = NETOPO(iRch)%RHORDER

  if (.not.allocated(RCHFLX(iens,jRch)%QFUTURE_IRF))then

    ntdh = size(NETOPO(jRch)%UH)

    allocate(RCHFLX(iens,jRch)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX(iens,ibas)%QFUTURE'; return; endif

    RCHFLX(iens,jRch)%QFUTURE_IRF(:) = 0._dp

   endif

  ! perform river network UH routing
  call conv_upsbas_qr(iEns,          & ! input: index of runoff ensemble to be processed
                      jRch,          & ! input: index of reach to be processed
                      ierr, cmessage)   ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! check
  if(NETOPO(jRch)%REACHIX == ixDesire)then
   print*, 'RCHFLX(iens,jRch)%BASIN_QR(1),RCHFLX(iens,jRch)%REACH_Q_IRF = ', &
            RCHFLX(iens,jRch)%BASIN_QR(1),RCHFLX(iens,jRch)%REACH_Q_IRF
  endif

 end do  ! (looping through stream segments)

 end subroutine irf_route


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
 ! subroutine: Compute delayed runoff from all the upstream segments
 ! *********************************************************************
 subroutine conv_upsbas_qr(&
                           ! input
                           iEns,       &    ! input: index of runoff ensemble to be processed
                           iRch,       &    ! input: index of reach to be processed
                           ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Convolute routed basisn flow volume at top of each of the upstream segment at one time step and at each segment
 !
 ! ----------------------------------------------------------------------------------------
 !   Updates structure RCHFLX%BASIN_QR_IRF
 !                     RCHFLX%REACH_Q_IRF
 !                     RCHFLX%QFUTURE_IRF
 ! ----------------------------------------------------------------------------------------
 ! global routing data
 USE globalData, only : RCHFLX ! routing fluxes
 USE globalData, only : NETOPO ! Network topology
 implicit none
 ! Input
 INTEGER(I4B), intent(IN)               :: iEns        ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)               :: iRch        ! reach segment to be routed
 ! Output
 integer(i4b), intent(out)              :: ierr        ! error code
 character(*), intent(out)              :: message     ! error message
 ! Local variables to
 real(dp)                               :: q_upstream  ! total discharge at top of the reach being processed
 INTEGER(I4B)                           :: iRch_up     ! index of upstream of the reach being processed
 INTEGER(I4B)                           :: nTDH        ! number of UH data
 INTEGER(I4B)                           :: iTDH        ! index of UH data (i.e.,future time step)
 INTEGER(I4B)                           :: nUps        ! number of all upstream segment
 INTEGER(I4B)                           :: iUps        ! loop indices for u/s reaches
 character(len=strLen)                  :: intStr      ! integer in string

 ! initialize error control
 ierr=0; message='conv_upsbas_qr/'

 ! identify number of upstream segments of the reach being processed
 nUps = size(NETOPO(iRch)%UREACHI)

 q_upstream = 0.0_dp
 if(nUps>0)then
   do iUps = 1,nUps
     iRch_up = NETOPO(iRch)%UREACHI(iUps)      !  index of upstream of jRch reach
     if (RCHFLX(iEns,iRch_up)%CHECK_IRF) then
       ! Find out total q at top of a segment
       q_upstream = q_upstream + RCHFLX(iEns,iRch_up)%REACH_Q_IRF
     else
       write(intStr,'(i10)') NETOPO(iRch)%REACHID
       ierr=20; message=trim(message)//'Upstream segment of '//intStr//' has not not routed yet'; return
     endif
   end do
 endif

 ! place a fraction of runoff in future time steps
 nTDH = size(NETOPO(iRch)%UH) ! identify the number of future time steps of UH for a given segment
 do iTDH=1,nTDH
   RCHFLX(iens,iRch)%QFUTURE_IRF(iTDH) = RCHFLX(iens,iRch)%QFUTURE_IRF(iTDH) &
                                       + NETOPO(iRch)%UH(iTDH)*q_upstream
 enddo

 ! Add local routed flow
 RCHFLX(iEns,iRch)%REACH_Q_IRF = RCHFLX(iEns,iRch)%QFUTURE_IRF(1)+RCHFLX(iEns,iRch)%BASIN_QR(1)

 ! move array back   use eoshift
 RCHFLX(iEns,iRch)%QFUTURE_IRF=eoshift(RCHFLX(iEns,iRch)%QFUTURE_IRF,shift=1)

 ! Check True since now this reach now routed
 RCHFLX(iEns,iRch)%CHECK_IRF=.True.

 end subroutine conv_upsbas_qr

end module irf_route_module


!      ! ********************************************************************************
!      ! subroutine: Sum of basin routed runoff from all the immediate upstream basin
!      ! *********************************************************************************
!      subroutine get_upsbas_qr(nSeg,iEns, &     ! input
!                               ierr, message)   ! output=error control
!      ! ----------------------------------------------------------------------------------------
!      ! Purpose:
!      !   Find sum of routed runoff volume [m3/s] from immediate upstream basins of the segment (at top of
!      !   the segment). Segment of Headwater basin does not have flow at top of the segment
!      ! ----------------------------------------------------------------------------------------
!
!      implicit none
!      ! Input
!      INTEGER(I4B), intent(IN)               :: iEns        ! ensemble member
!      INTEGER(I4B), intent(IN)               :: nSeg        ! reach to process
!      ! Output
!      integer(i4b), intent(out)              :: ierr        ! error code
!      character(*), intent(out)              :: message     ! error message
!      ! Local variables
!      INTEGER(I4B)                           :: iUpsBas     ! loop through u/s reaches
!      INTEGER(I4B)                           :: jUpsBas     ! loop through u/s reaches
!      INTEGER(I4B)                           :: nUpsBas     ! number of immediate upstream basins
!      INTEGER(I4B)                           :: iSeg        ! index of reach with the earliest time
!      real(DP)                               :: area        ! Area of one upstream basin
!      real(DP)                               :: qrTotal     ! total routed runoff volume of all the upstream basins
!
!      ! initialize error control
!      ierr=0; message='get_upsbas_qr/'
!
!      seg_loop: do iSeg=1,nSeg
!        nUpsBas = SIZE(NETOPO(iSeg)%UREACHI)      ! number of upstream basins
!        qrTotal=0._dp
!    !    print *,'--------------------------'
!    !    print *,'Working on iSeg, ID, # of ups. basins= ', iSeg, NETOPO(iSeg)%REACHID, nUpsBas
!        if (nUpsBas > 0) then ! if segment is not located in headwater basin
!          do iUpsbas=1,nUpsBas
!            jUpsBas = NETOPO(iSeg)%UREACHI(iUpsBas)    ! index of the imediate upstream upstream reach
!            ! Get routed flow from all immediate upstream basins [m3/s] and sum them up
!            area=RPARAM(jUpsBas)%BASAREA
!            if ( area > 0) then
!              qrTotal = qrTotal + RCHFLX(iEns,jUpsBas)%BASIN_QR(1)!/area
!            endif
!    !        write(*,'(a,1x,i4,1x,f20.2,1x,f20.4)') 'IupBas, area, routed flow= ',jUpsBas, area, RCHFLX(iEns,jUpsBas)%BASIN_QR(1)
!          enddo !end of upstream basins
!        endif
!        RCHFLX(iEns,iSeg)%UPSBASIN_QR = qrTotal
!    !    write(*,'(a,1x,es14.7)') 'Sum of routed flow = ',RCHFLX(iEns,iSeg)%UPSBASIN_QR
!      enddo seg_loop
!
!      end subroutine

