module irf_route

USE nrtype

implicit none

!Followings accessible outside module
public::make_uh
public::get_upsbas_qr
public::conv_upsbas_qr
! anything else
private

contains

!    ! *********************************************************************
!    ! subroutine: Update total stream network length from a each esgment
!    ! *********************************************************************
!     subroutine upstrm_length(nSeg, &          ! input
!                              ierr, message)   ! error control
!     ! ----------------------------------------------------------------------------------------
!     ! Purpose:
!     !
!     !   Calculate total length of upstream reach network from each segment
!     !
!     ! ----------------------------------------------------------------------------------------
!     USE reachparam
!
!     implicit none
!     ! input variables
!     integer(I4B), intent(in)               :: nSeg          ! number of stream segments
!     ! output variables
!     integer(i4b), intent(out)              :: ierr          ! error code
!     character(*), intent(out)              :: message       ! error message
!     ! local variables
!     integer(I4B)                           :: iSeg          ! index for segment loop
!     integer(I4B)                           :: iUps          ! index for upstream segment loop
!     integer(I4B)                           :: jUps          ! index for upstream segment loop
!     integer(I4B)                           :: nUps          ! number of upstream reaches
!     real(DP)                               :: xLocal        ! length of one segement
!     real(DP)                               :: xTotal        ! total length of upstream segment
!
!     ! initialize error control
!     ierr=0; message='strmlength/'
!
!     seg_loop: do iSeg=1,nSeg !Loop through each segment
!
!       nUps = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment
!       allocate(NETOPO(iSeg)%UPSLENG(nUps),stat=ierr)
!       !print *,'--------------------------------------------'
!       !print *,'Seg ID, Num of upstream', iSeg, nUps
!       !print *,'list of upstream index',NETOPO(iSeg)%RCHLIST
!       upstrms_loop: do iUps=1,nUps !Loop through upstream segments of current segment
!         jUps=NETOPO(iSeg)%RCHLIST(iUps) !index of upstreamf segment
!         xTotal = 0.0 !Initialize total length of upstream segments
!         do
!           xLocal=RPARAM(jUps)%RLENGTH ! Get a segment length
!           xTotal=xTotal+xLocal
!           if (jUps.eq.NETOPO(iSeg)%REACHIX) exit
!           jUps = NETOPO(jUps)%DREACHI ! Get index of immediate downstream segment
!         enddo
!         NETOPO(iSeg)%UPSLENG(iUps)=xTotal
!       enddo upstrms_loop
!     enddo seg_loop
!
!     end subroutine

! *********************************************************************
! subroutine: compute normalized UH from Saint-Venant Eq. at sim. time step
!             for all the upstream segment
! *********************************************************************
  subroutine make_uh(nSeg,dt,velo,diff, &        ! input
                     ierr, message)              ! output=error control
  ! ----------------------------------------------------------------------------------------
  ! Purpose:
  !
  !   Calculate UH at given simulation time step from impulse response function
  !   for upstream reach network from each segment
  !   Use Saint-Venant equation -see equations (13)-(15) in Lohmann, et al. (1996) Tellus
  ! ----------------------------------------------------------------------------------------
  USE globalData, only : RPARAM  ! reach properties
  USE globalData, only : NETOPO  ! network topology
  USE globalData, only : PI      ! pi

  ! Define variables
  implicit none

  ! input variables
  integer(i4b), intent(in)              :: nSeg          ! Numer of segment
  real(dp), intent(in)                  :: dt            ! Interval of each time step [sec]
  real(dp), intent(in)                  :: velo          ! Wave velocity C for each segment [m/sec]
  real(dp), intent(in)                  :: diff          ! Diffusivity D for each segment [m2/sec]
  ! output variables
  integer(i4b), intent(out)             :: ierr          ! error code
  character(*), intent(out)             :: message       ! error message
  ! local variable
  real(dp),dimension(:),allocatable     :: UHM           !
  real(dp),dimension(:),allocatable     :: UHQ           !
  real(dp),dimension(:),allocatable     :: fr            ! Unit runoff depth evenly distributed over the simulation duration at hourly step
  real(dp)                              :: seg_length    ! length of rech segment
  real(dp)                              :: POT           ! Inside of exponential function of IRF
  real(dp)                              :: H             ! IRF, h(x,t) function h function in eq.15 Lohmann et al.1996
  real(dp)                              :: UHQ0          !
  real(dp)                              :: INTE          ! Accumulation of hrly UH for ordinate normalization
  real(dp)                              :: sec           ! hr at each time step  [hr]
  real(dp),parameter                    :: dTUH=3600.0   ! UH time step [sec]
  integer(i4b)                          :: iHrStrt       ! index of UH time step where rising limb of UH start
  integer(i4b)                          :: iHrLast       ! index of UH time step where recession limb of UH become zero
  integer(i4b)                          :: nTSub         ! # of time step where 1/nTsub [m] of runoff is inserted
  integer(i4b)                          :: iSeg          ! Loop index
  integer(i4b)                          :: iHr,jHr       ! Loop index of hour
  integer(i4b)                          :: iTagg         ! index for aggregated (i.e. simulation) time step
  integer(i4b),parameter                :: nTMAX=240     ! Maximum hour of UH [hr] - 10 days times 24hrs
  integer(i4b),parameter                :: nHr=240       ! Maximum hour of UH [hr] - 10 days times 24hrs

  ! initialize error control
  ierr=0; message='make_uh/'

  ! Dynamically assigned parameters
  nTsub=ceiling(dt/dTUH)
  !nTsub=floor(dt/dTUH)

  ! Memory allocation
  allocate(fr(nTMAX),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for fr'; return; endif
  allocate(UHQ(nTMAX),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for UHQ'; return; endif
  allocate(UHM(nHr),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for UHM'; return; endif

  ! Unit Runoff depth [1/nTsub, 1/nTsub, ..., 1/nTsub, 0, 0, ...0]
  fr=0._dp
  fr(1:nTsub) = 1.0 / nTsub

  do iSeg=1,nSeg

   !Compute IRF
   seg_length = RPARAM(iSeg)%RLENGTH
   INTE = 0._dp
   sec = 0._dp
   UHM(:) = 0._dp
   do iHr=1,nHr
     sec = sec + dTUH
     if (velo .GT. 0.0) then
       POT = ((velo*sec-seg_length)**2.0)/(4.0*diff*sec)
       if (POT .GT. 69.0) then
         H = 0.0
       else
         H = 1.0/(2.0*sqrt(PI*diff*sec))*seg_length*exp(-POT)
       endif
     else
       H = 0.0
     endif
     UHM(iHr) = H
     INTE = INTE + H
   enddo

   !Normalize ordinates by sum of IRF
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
     INTE=INTE+UHQ0
   end do
   ! Normalize ordinates by sum of total flow
   if (INTE > 0.0) then
     UHQ= UHQ/INTE
   end if

   !Get hour indices of start of rising part of UH and end of recession part of UH
   INTE = 0._dp
   do iHr=1,nTMAX
     INTE= INTE+UHQ(iHr)
     iHrLast=iHr
     if (INTE > 0.9999) exit
   enddo

   !Aggregate hourly unit hydrograph to simulation time step
   allocate(NETOPO(iSeg)%UH((iHrLast+nTsub-1)/nTsub),stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate space for UH%UH_DATA'; return; endif
   NETOPO(iSeg)%UH(:)=0._dp
   do jHr = 1,iHrLast
     iTagg = (jHr+nTsub-1)/nTsub
     NETOPO(iSeg)%UH(iTagg) = NETOPO(iSeg)%UH(iTagg)+UHQ(jHr)
   end do

  end do ! sSeg loop

  end subroutine make_uh

  ! ********************************************************************************
  ! subroutine: Sum of basin routed runoff from all the immediate upstream basin
  ! *********************************************************************************
  subroutine get_upsbas_qr(nSeg,iEns, &     ! input
                           ierr, message)   ! output=error control
  ! ----------------------------------------------------------------------------------------
  ! Purpose:
  !   Find sum of routed runoff volume [m3/s] from immediate upstream basins of the segment (at top of
  !   the segment). Segment of Headwater basin does not have flow at top of the segment
  ! ----------------------------------------------------------------------------------------
  USE globalData, only : NETOPO  ! network topology
  USE globalData, only : RPARAM  ! reach properties
  USE globalData, only : RCHFLX  ! reach fluxes

  implicit none
  ! Input
  INTEGER(I4B), intent(IN)               :: iEns        ! ensemble member
  INTEGER(I4B), intent(IN)               :: nSeg        ! reach to process
  ! Output
  integer(i4b), intent(out)              :: ierr        ! error code
  character(*), intent(out)              :: message     ! error message
  ! Local variables
  INTEGER(I4B)                           :: iUpsBas     ! loop through u/s reaches
  INTEGER(I4B)                           :: jUpsBas     ! loop through u/s reaches
  INTEGER(I4B)                           :: nUpsBas     ! number of immediate upstream basins
  INTEGER(I4B)                           :: iSeg        ! index of reach with the earliest time
  real(DP)                               :: area        ! Area of one upstream basin
  real(DP)                               :: qrTotal     ! total routed runoff volume of all the upstream basins

  ! initialize error control
  ierr=0; message='get_upsbas_qr/'

  seg_loop: do iSeg=1,nSeg
    nUpsBas = SIZE(NETOPO(iSeg)%UREACHI)      ! number of upstream basins
    qrTotal=0._dp
!    print *,'--------------------------'
!    print *,'Working on iSeg, ID, # of ups. basins= ', iSeg, NETOPO(iSeg)%REACHID, nUpsBas
    if (nUpsBas > 0) then ! if segment is not located in headwater basin
      do iUpsbas=1,nUpsBas
        jUpsBas = NETOPO(iSeg)%UREACHI(iUpsBas)    ! index of the imediate upstream upstream reach
        ! Get routed flow from all immediate upstream basins [m3/s] and sum them up
        area=RPARAM(jUpsBas)%BASAREA
        if ( area > 0) then
          qrTotal = qrTotal + RCHFLX(iEns,jUpsBas)%BASIN_QR(1)!/area
        endif
!        write(*,'(a,1x,i4,1x,f20.2,1x,f20.4)') 'IupBas, area, routed flow= ',jUpsBas, area, RCHFLX(iEns,jUpsBas)%BASIN_QR(1)
      enddo !end of upstream basins
    endif
    RCHFLX(iEns,iSeg)%UPSBASIN_QR = qrTotal
!    write(*,'(a,1x,es14.7)') 'Sum of routed flow = ',RCHFLX(iEns,iSeg)%UPSBASIN_QR
  enddo seg_loop

  end subroutine

  ! *********************************************************************
  ! subroutine: Compute delayed runoff from all the upstream segments
  ! *********************************************************************
  subroutine conv_upsbas_qr(iEns, &          ! input:
                            nSeg, &          ! input:
                            ierr, message)   ! output: error control
  ! ----------------------------------------------------------------------------------------
  ! Purpose:
  !
  !   Convolute routed basisn flow volume at top of each of the upstream segment at one time step and at each segment
  !
  ! ----------------------------------------------------------------------------------------
  !   Updates structure RCHFLX%BASIN_QR_IRF (module reach_flux)
  !                     RCHFLX%REACH_Q_IRF
  !                     RCHFLX%QFUTURE_IRF
  ! ----------------------------------------------------------------------------------------
  USE globalData, only : NETOPO  ! network topology
  USE globalData, only : RCHFLX  ! reach fluxes
  implicit none
  ! Input
  INTEGER(I4B), intent(IN)               :: iens        ! ensemble member
  INTEGER(I4B), intent(IN)               :: nSeg        ! number of total segments
  ! Output
  integer(i4b), intent(out)              :: ierr        ! error code
  character(*), intent(out)              :: message     ! error message
  ! Local variables to
  real(dp)                               :: q_upstream  ! total discharge at top of the reach being processed
  INTEGER(I4B)                           :: irch        ! index of reach to be processed
  INTEGER(I4B)                           :: irch_up     ! index of upstream of the reach being processed
  INTEGER(I4B)                           :: iSeg        ! index of segment
  INTEGER(I4B)                           :: nTDH        ! number of UH data
  INTEGER(I4B)                           :: iTDH        ! index of UH data (i.e.,future time step)
  INTEGER(I4B)                           :: nUps        ! number of all upstream segment
  INTEGER(I4B)                           :: iUps        ! loop indices for u/s reaches

  ! initialize error control
  ierr=0; message='conv_upsbas_qr/'

  do iSeg=1,nSeg

    ! identify index of reach to process
    irch = NETOPO(iSeg)%RHORDER

    ! identify number of upstream segments of the reach being processed
    nUps = size(NETOPO(irch)%UREACHI)

    q_upstream = 0.0_dp
    if(nUps>0)then
      do iUps = 1,nUps
        irch_up = NETOPO(irch)%UREACHI(iUps)      !  index of upstream of irch reach
        ! Find out total q at top of a segment
        q_upstream = q_upstream + RCHFLX(iEns,irch_up)%REACH_Q_IRF
      end do
    endif

    ! place a fraction of runoff in future time steps
    nTDH = size(NETOPO(irch)%UH) ! identify the number of future time steps of UH for a given segment
    do iTDH=1,nTDH
      RCHFLX(iens,irch)%QFUTURE_IRF(iTDH) = RCHFLX(iens,irch)%QFUTURE_IRF(iTDH) &
                                          + NETOPO(irch)%UH(iTDH)*q_upstream
    enddo

    ! Add local routed flow
    RCHFLX(iEns,irch)%REACH_Q_IRF = RCHFLX(iEns,irch)%QFUTURE_IRF(1)+RCHFLX(iEns,irch)%BASIN_QR(1)

    ! move array back   use eoshift
    RCHFLX(iEns,irch)%QFUTURE_IRF=eoshift(RCHFLX(iEns,irch)%QFUTURE_IRF,shift=1)

  end do

     !  Old routine
     !  ! route one time step routed runoff volue from each up stream basins to segment
     !  seg_loop: do iSeg=1,nSeg
     !    nUps = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment
     !    !print *,'--------------------------'
     !    !print *,'Working on iSeg, ID, # of ups. basins= ', iSeg, NETOPO(iSeg)%RCHLIST, nUps
     !    !print *,'upbasinflow= ',RCHFLX(iens,iSeg)%UPSBASIN_QR
     !
     !    upstrms_loop: do iUps=1,nUps
     !      ! index of the upstream upstream basin
     !      jUps = NETOPO(iSeg)%RCHLIST(iUps)
     !
     !      ! identify the number of future time steps for a given segment
     !      nTDH = size(NETOPO(iSeg)%UH(iUps)%UH_DATA)
     !
     !      ! place a fraction of runoff in future time steps
     !      UH_data: do iTDH=1,nTDH
     !        RCHFLX(iens,iSeg)%QFUTURE_IRF(iTDH) = RCHFLX(iens,iSeg)%QFUTURE_IRF(iTDH) &
     !                                            + NETOPO(iSeg)%UH(iUps)%UH_DATA(iTDH)*RCHFLX(iEns,jUps)%UPSBASIN_QR
     !      enddo UH_data
     !    enddo upstrms_loop
     !
     !    ! save the routed runoff
     !    RCHFLX(iEns,iSeg)%BASIN_QR_IRF(0) = RCHFLX(iEns,iSeg)%BASIN_QR_IRF(1)  ! (save the runoff from the previous time step)
     !    RCHFLX(iEns,iSeg)%BASIN_QR_IRF(1) = RCHFLX(iEns,iSeg)%QFUTURE_IRF(1)
     !    ! Add local routed flow
     !    RCHFLX(iEns,iSeg)%REACH_Q_IRF    = RCHFLX(iEns,iSeg)%BASIN_QR_IRF(1)+RCHFLX(iEns,iSeg)%BASIN_QR(1)
     !    ! move array back   use eoshift
     !    RCHFLX(iEns,iSeg)%QFUTURE_IRF=eoshift(RCHFLX(iEns,iSeg)%QFUTURE_IRF,shift=1)

     !  enddo seg_loop

  end subroutine

end module irf_route
