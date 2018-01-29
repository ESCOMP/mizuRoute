module irf_route

USE nrtype

implicit none

!Followings accessible outside module
public::upstrm_length
public::make_uh
public::get_upsbas_qr
public::conv_upsbas_qr
! anything else
private

contains

 ! *********************************************************************
 ! subroutine: Update total stream network length from a each esgment
 ! *********************************************************************
  subroutine upstrm_length(nSeg, &          ! input
                           ierr, message)   ! error control
  ! ----------------------------------------------------------------------------------------
  ! Creator(s):
  !
  !   Naoki Mizukami
  !
  ! ----------------------------------------------------------------------------------------
  ! Purpose:
  !
  !   Calculate total length of upstream reach network from each segment 
  !
  ! ----------------------------------------------------------------------------------------
  ! I/O:
  !
  !   INPUTS:
  !    nSeg: Number of stream segments in the river network (reaches)
  !
  ! ----------------------------------------------------------------------------------------
  ! Structures modified:
  !
  !   Updates structure NETOPO%UPSLENG (module reachparam)
  !
  ! ----------------------------------------------------------------------------------------
  USE reachparam
  
  implicit none 
  ! input variables
  integer(I4B), intent(in)               :: nSeg          ! number of stream segments
  ! output variables
  integer(i4b), intent(out)              :: ierr          ! error code
  character(*), intent(out)              :: message       ! error message
  ! local variables
  integer(I4B)                           :: iSeg          ! index for segment loop
  integer(I4B)                           :: iUps          ! index for upstream segment loop 
  integer(I4B)                           :: jUps          ! index for upstream segment loop 
  integer(I4B)                           :: nUps          ! number of upstream reaches
  real(DP)                               :: xLocal        ! length of one segement
  real(DP)                               :: xTotal        ! total length of upstream segment
  
  ! initialize error control
  ierr=0; message='strmlength/'
  
  seg_loop: do iSeg=1,nSeg !Loop through each segment
  
    nUps = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment 
    allocate(NETOPO(iSeg)%UPSLENG(nUps),stat=ierr)
    !print *,'--------------------------------------------'
    !print *,'Seg ID, Num of upstream', iSeg, nUps
    !print *,'list of upstream index',NETOPO(iSeg)%RCHLIST
    upstrms_loop: do iUps=1,nUps !Loop through upstream segments of current segment
      jUps=NETOPO(iSeg)%RCHLIST(iUps) !index of upstreamf segment
      xTotal = 0.0 !Initialize total length of upstream segments
      do 
        xLocal=RPARAM(jUps)%RLENGTH ! Get a segment length
        xTotal=xTotal+xLocal         
        if (jUps.eq.NETOPO(iSeg)%REACHIX) exit 
        jUps = NETOPO(jUps)%DREACHI ! Get index of immediate downstream segment 
      enddo
      NETOPO(iSeg)%UPSLENG(iUps)=xTotal
    enddo upstrms_loop 
  enddo seg_loop
  
  end subroutine

! *********************************************************************
! subroutine: compute normalized UH from Saint-Venant Eq. at sim. time step
!             for all the upstream segment
! *********************************************************************
  subroutine make_uh(nSeg,dt,velo,diff, &        ! input
                     ierr, message)              ! output=error control
  ! ----------------------------------------------------------------------------------------
  ! Creator(s):
  !
  !   Naoki Mizukami
  !
  ! ----------------------------------------------------------------------------------------
  ! Purpose:
  !
  !   Calculate UH at given simulation time step from impulse response function 
  !   for upstream reach network from each segment
  !   Use Saint-Venant equation -see equations (13)-(15) in Lohmann, et al. (1996) Tellus
  !
  ! ----------------------------------------------------------------------------------------
  ! I/O:
  !
  !   INPUTS:
  !    nSeg: Number of stream segments in the river network (reaches)
  !    dt:   Simulation time step [sec] 
  !    velo: velocity in Saint-Venant eq. [m/sec] 
  !    diff: diffusivity in Ssaint-Venant eq. [m2/sec] 
  !
  ! ----------------------------------------------------------------------------------------
  ! Structures modified:
  !
  !   Updates structure NETOPO%UH%UH_DATA (module reachparam)
  !
  ! ----------------------------------------------------------------------------------------
  USE reachparam
  USE reach_flux 

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
  real(dp)                              :: POT           ! Inside of exponential function of IRF 
  real(dp)                              :: H             ! IRF, h(x,t) function h function in eq.15 Lohmann et al.1996
  real(dp)                              :: UHQ0          ! 
  real(dp)                              :: INTE          ! Accumulation of hrly UH for ordinate normalization
  real(dp)                              :: sec           ! hr at each time step  [hr]
  real(dp),parameter                    :: dTUH=3600.0   ! UH time step [sec] 
  integer(I4B)                          :: nUps          ! # of upstream reaches
  integer(i4b)                          :: iHrStrt       ! index of UH time step where rising limb of UH start
  integer(i4b)                          :: iHrLast       ! index of UH time step where recession limb of UH become zero
  integer(i4b)                          :: nTSub         ! # of time step where 1/nTsub [m] of runoff is inserted  
  integer(i4b)                          :: iSeg          ! Loop index
  integer(i4b)                          :: iUps          ! Loop index
  integer(i4b)                          :: iHr           ! Loop index
  integer(i4b)                          :: iHr1          ! Loop index
  integer(i4b)                          :: iTagg         ! index for aggregated (i.e. simulation) time step
  integer(i4b),parameter                :: nTMAX=4560    ! Maximum hour of UH [hr] - 200 days times 24hrs
  integer(i4b),parameter                :: nHr=4560      ! Maximum hour of UH [hr] - 200 days times 24hrs

  ! Dynamically assigned parameters
  nTsub=ceiling(dt/dTUH)
  !nTsub=floor(dt/dTUH)
  
  ! initialize error control
  ierr=0; message='make_uh/'

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

  seg_loop: do iSeg=1,nSeg  
    nUps = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment 
    allocate(NETOPO(iSeg)%UH(nUps),stat=ierr)
    if(ierr/=0)then; message=trim(message)//'unable to allocate space for NETOPO%UH'; return; endif

    upstrms_loop: do iUps=1,nUps
      !Compute IRF
      INTE = 0._dp
      sec = 0._dp
      UHM(:) = 0._dp
      do iHr=1,nHr
        sec = sec + dTUH 
        if (velo .GT. 0.0) then
          POT = ((velo*sec-NETOPO(iSeg)%UPSLENG(iUps))**2.0)/(4.0*diff*sec)
          if (POT .GT. 69.0) then 
            H = 0.0
          else 
            H = 1.0/(2.0*sqrt(PI*diff*sec))*NETOPO(iSeg)%UPSLENG(iUps)*exp(-POT)
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
      do iHr1 = 1,nTMAX 
        UHQ0=0._dp 
        do iHr = iHrStrt,iHrLast
          if ((iHr1-iHr) > 0) then 
            if (iHr1-iHr <= nTsub ) then
              UHQ0 = UHQ0 + fr(iHr1-iHr)*UHM(iHr)
            endif
          else 
            exit 
          endif 
        enddo 
        UHQ(iHr1) = UHQ0
        INTE=INTE+UHQ0
      end do 
     !Normalize ordinates by sum of total flow
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
      !INTE = 0._dp
      !do iHr=nTMAX,1,-1
      !  INTE= INTE+UHQ(iHr)
      !  iHrStrt=iHr
      !  if (INTE > 0.9999) exit
      !enddo
      !Aggregate hourly unit hydrograph to simulation time step
      allocate(NETOPO(iSeg)%UH(iUps)%UH_DATA((iHrLast+nTsub-1)/nTsub),stat=ierr)
      if(ierr/=0)then; message=trim(message)//'unable to allocate space for UH%UH_DATA'; return; endif
      
      NETOPO(iSeg)%UH(iUps)%UH_DATA(:)=0._dp
      do iHr1 = 1,iHrLast
        iTagg = (iHr1+nTsub-1)/nTsub
        NETOPO(iSeg)%UH(iUps)%UH_DATA(iTagg) = NETOPO(iSeg)%UH(iUps)%UH_DATA(iTagg)+UHQ(iHr1)
      end do 

    enddo upstrms_loop
  enddo seg_loop
    
  end subroutine make_uh

  ! ********************************************************************************
  ! subroutine: Sum of basin routed runoff from all the immediate upstream basin 
  ! *********************************************************************************
  subroutine get_upsbas_qr(nSeg,iEns, &     ! input
                           ierr, message)   ! output=error control
  ! ----------------------------------------------------------------------------------------
  ! Creator(s):
  !
  !   Naoki Mizukami
  !
  ! ----------------------------------------------------------------------------------------
  ! Purpose:
  !
  !   Find sum of routed runoff volume [m3/s] from immediate upstream basins of the segment (at top of
  !   the segment). Segment of Headwater basin does not have flow at top of the segment 
  !
  ! ----------------------------------------------------------------------------------------
  ! I/O:
  !
  !   INPUTS:
  !    nSeg: Number of stream segments in the river network (reaches)
  !    iEns: index of ensemble flow 
  !
  ! ----------------------------------------------------------------------------------------
  ! Structures modified:
  !  
  !   Updates structure RCHFLX%UPSBASIN_QR (module reach_flux)
  !  
  ! ----------------------------------------------------------------------------------------
  USE reachparam
  USE reach_flux
 
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
  subroutine conv_upsbas_qr(nSeg,iEns, &    ! input
                           ierr, message)   ! output=error control
  ! ----------------------------------------------------------------------------------------
  ! Creator(s):
  !
  !   Naoki Mizukami
  !
  ! ----------------------------------------------------------------------------------------
  ! Purpose:
  !
  !   Convolute routed basisn flow volume at top of each of the upstream segment at one time step and at each segment 
  !    
  ! ----------------------------------------------------------------------------------------
  ! I/O:
  !
  !   INPUTS:
  !    nSeg: Number of stream segments in the river network (reaches)
  !    iEns: index of ensemble flow 
  !
  ! ----------------------------------------------------------------------------------------
  ! Structures modified:
  !  
  !   Updates structure RCHFLX%BASIN_QR_IRF (module reach_flux)
  !                     RCHFLX%REACH_Q_IRF
  ! ----------------------------------------------------------------------------------------
  USE reachparam
  USE reach_flux
 
  implicit none
  ! Input
  INTEGER(I4B), intent(IN)               :: iens      ! ensemble member
  INTEGER(I4B), intent(IN)               :: nSeg      ! reach to process
  ! Output
  integer(i4b), intent(out)              :: ierr      ! error code
  character(*), intent(out)              :: message   ! error message
  ! Local variables to 
  INTEGER(I4B)                           :: iUps      ! loop through u/s reaches
  INTEGER(I4B)                           :: jUps      ! 
  INTEGER(I4B)                           :: nUps      ! number of all upstream segment 
  INTEGER(I4B)                           :: iSeg      ! index of segment
  INTEGER(I4B)                           :: nTDH      ! number of UH data 
  INTEGER(I4B)                           :: iTDH      ! index of UH data (i.e.,future time step)

  ! initialize error control
  ierr=0; message='conv_upsbas_qr/'

  ! route one time step routed runoff volue from each up stream basins to segment 
  seg_loop: do iSeg=1,nSeg
    nUps = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment 
    !print *,'--------------------------'
    !print *,'Working on iSeg, ID, # of ups. basins= ', iSeg, NETOPO(iSeg)%RCHLIST, nUps
    !print *,'upbasinflow= ',RCHFLX(iens,iSeg)%UPSBASIN_QR

    upstrms_loop: do iUps=1,nUps
      ! index of the upstream upstream basin
      jUps = NETOPO(iSeg)%RCHLIST(iUps)    

      ! identify the number of future time steps for a given segment
      nTDH = size(NETOPO(iSeg)%UH(iUps)%UH_DATA)

      ! place a fraction of runoff in future time steps
      UH_data: do iTDH=1,nTDH
        RCHFLX(iens,iSeg)%QFUTURE_IRF(iTDH) = RCHFLX(iens,iSeg)%QFUTURE_IRF(iTDH) &
                                            + NETOPO(iSeg)%UH(iUps)%UH_DATA(iTDH)*RCHFLX(iEns,jUps)%UPSBASIN_QR
      enddo UH_data 
    enddo upstrms_loop 

    ! save the routed runoff
    RCHFLX(iEns,iSeg)%BASIN_QR_IRF(0) = RCHFLX(iEns,iSeg)%BASIN_QR_IRF(1)  ! (save the runoff from the previous time step)
    RCHFLX(iEns,iSeg)%BASIN_QR_IRF(1) = RCHFLX(iEns,iSeg)%QFUTURE_IRF(1)
    ! Add local routed flow
    RCHFLX(iEns,iSeg)%REACH_Q_IRF    = RCHFLX(iEns,iSeg)%BASIN_QR_IRF(1)+RCHFLX(iEns,iSeg)%BASIN_QR(1)
    ! move array back   use eoshift
    RCHFLX(iEns,iSeg)%QFUTURE_IRF=eoshift(RCHFLX(iEns,iSeg)%QFUTURE_IRF,shift=1)

  enddo seg_loop

  end subroutine

end module irf_route
