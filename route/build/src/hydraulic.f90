MODULE hydraulic

! hydraulic computations based on channel cross section assumption
!
! main channel assumes trapezoidal, rectangle, or triangle
!   y: flow depth [m], b: bottom width, zc: side slope (1:zc = height:horizontal)
!   if zc = 0, rectangle. if b = 0, triangle,
!
! floodplain
!   zf: side slope (1:zf = height:horizontal)
!
! all the functions use two optional inputs - zf and bankDepth. if they are not provided,
! channel is assumed to be unlimited bank depth.

USE nrtype
USE globalData, ONLY: bankDepth0 => high_depth

implicit none

private
public::Btop
public::Pwet
public::flow_area
public::uniformFlow
public::flow_depth
public::water_height
public::storage
public::celerity
public::diffusivity

! module-wide parameters
logical(lgt), parameter :: useFrictionSlope = .true.    ! .false. -> approximate friction slope with channel slope
real(dp),     parameter :: const13=1._dp/3._dp          ! constant
real(dp),     parameter :: const23=2._dp/3._dp          ! constant
real(dp),     parameter :: const43=4._dp/3._dp          ! constant
real(dp),     parameter :: const53=5._dp/3._dp          ! constant
real(dp),     parameter :: const103=10._dp/3._dp        ! constant
real(dp),     parameter :: err_thresh=0.005_dp          ! newton method convergence threshold
real(dp),     parameter :: zf0=1000._dp                 ! default floodplain slope: horizontal:vertcal=zf:1 [-]

CONTAINS

 ! *********************************************************************
 ! public function: top channel width [m]
 ! *********************************************************************
 elemental FUNCTION Btop(yin, b, zc, zf, bankDepth)
   implicit none
   ! Argument variables
   real(dp), intent(in)             :: yin          ! flow depth [m]
   real(dp), intent(in)             :: b            ! channel bottom width [m2] 0-> triangular channel
   real(dp), intent(in)             :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), optional, intent(in)   :: zf           ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp), optional, intent(in)   :: bankDepth    ! bankfull depth [m]
   real(dp)                         :: Btop         ! channel top width [m]
   ! Local variables
   real(dp)                         :: zf1          ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp)                         :: bankDepth1   ! bankfull depth [m]

   if (present(bankDepth)) then
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   if (yin <= bankDepth1) then
     Btop = b+2*yin*zc
   else
     Btop = b+2*bankDepth1*zc
     Btop = Btop + zf1*(yin-bankDepth1)*2
   end if

 END FUNCTION Btop

 ! *********************************************************************
 ! public function: wetted perimeter [m]
 ! *********************************************************************
 elemental FUNCTION Pwet(yin, b ,zc, zf, bankDepth)
   implicit none
   ! Argument variables
   real(dp), intent(in)             :: yin          ! flow depth [m]
   real(dp), intent(in)             :: b            ! channel bottom width [m2] 0-> triangular channel
   real(dp), intent(in)             :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), optional, intent(in)   :: zf           ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp), optional, intent(in)   :: bankDepth    ! bankfull depth [m]
   real(dp)                         :: Pwet         ! wetted perimeter [m]
   ! Local variables
   real(dp)                         :: zf1          ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp)                         :: bankDepth1   ! bankfull depth [m]

   if (present(bankDepth)) then
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   if (yin <= bankDepth1) then
     Pwet = b + 2*yin *sqrt(1 +zc*zc)
   else
     Pwet = b + 2*bankDepth1 *sqrt(1 +zc*zc)
     Pwet = Pwet + 2*(yin-bankDepth1) *sqrt(1 +zf1*zf1)
   end if

 END FUNCTION Pwet

 ! *********************************************************************
 ! public function: flow area [m2]
 ! *********************************************************************
 elemental FUNCTION flow_area(yin, b, zc, zf, bankDepth)
   implicit none
   ! Argument variables
   real(dp), intent(in)              :: yin          ! flow depth [m]
   real(dp), intent(in)              :: b            ! channel bottom width [m2] 0-> triangular channel
   real(dp), intent(in)              :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), optional, intent(in)    :: zf           ! floodplain slope: horizontal:vertcal=1:zf [-]
   real(dp), optional, intent(in)    :: bankDepth    ! bankfull depth [m]
   real(dp)                          :: flow_area    ! flow area [m2]
   ! Local variables
   real(dp)                          :: Bt           ! channel width at flow top [m]
   real(dp)                          :: Bt_bank      ! channel width at bankfull depth [m]
   real(dp)                          :: zf1          ! floodplain slope: horizontal:vertcal=1:zf [-]
   real(dp)                          :: bankDepth1   ! bankfull depth [m]

   if (present(bankDepth)) then
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   if (yin <= bankDepth1) then
     flow_area = yin*(b +zc*yin)
   else
     flow_area = bankDepth1*(b +zc*bankDepth1)
     Bt = Btop(yin, b, zc, zf=zf1, bankDepth=bankDepth1)
     Bt_bank = Btop(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
     flow_area = flow_area + (yin-bankDepth1)*(Bt+Bt_bank)/2._dp
   end if

 END FUNCTION flow_area

 ! *********************************************************************
 ! public function: water height if you know flow area [m]
 ! *********************************************************************
 elemental FUNCTION water_height(flowArea, b, zc, zf, bankDepth)
   implicit none
   ! Argument variables
   real(dp), intent(in)              :: flowArea     ! flow area [m2]
   real(dp), intent(in)              :: b            ! channel bottom width [m2] 0-> triangular channel
   real(dp), intent(in)              :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), optional, intent(in)    :: zf           ! floodplain slope: horizontal:vertcal=1:zf [-]
   real(dp), optional, intent(in)    :: bankDepth    ! bankfull depth [m]
   real(dp)                          :: water_height ! water height based on flow area [m]
   ! Local variables
   real(dp)                          :: Bt_bank      ! channel width at bankfull depth [m]
   real(dp)                          :: A_bank       ! channel bankful area [m2]
   real(dp)                          :: zf1          ! floodplain slope: horizontal:vertcal=1:zf [-]
   real(dp)                          :: bankDepth1   ! bankfull depth [m]

   if (present(bankDepth)) then
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   A_bank = flow_area(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
   if (flowArea>A_bank) then
     Bt_bank = Btop(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
     water_height = bankDepth1 + (sqrt(Bt_bank*Bt_bank/zf1/zf1 + 4._dp*flowArea/zf1 - 4._dp*bankDepth1*(b+bankDepth1*zc)/zf1)-Bt_bank/zf1)/2._dp
   else
     if (zc==0) then
       water_height = flowArea/b
     else
       water_height = (sqrt(b*b + 4._dp*flowArea*zc)-b)/2._dp/zc
     end if
   end if

 END FUNCTION water_height

 ! *********************************************************************
 ! public function: channel storage [m3]
 ! *********************************************************************
 elemental FUNCTION storage(yin, l, b, zc, zf, bankDepth)
   ! compute full capacity given a depth, NOT actual storage based on surface water profile
   implicit none
   ! Argument variables
   real(dp), intent(in)              :: yin          ! depth [m]
   real(dp), intent(in)              :: l            ! channel length [m]
   real(dp), intent(in)              :: b            ! channel bottom width [m2] 0-> triangular channel
   real(dp), intent(in)              :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), optional, intent(in)    :: zf           ! floodplain slope: horizontal:vertcal=1:zf [-]
   real(dp), optional, intent(in)    :: bankDepth    ! bankfull depth [m]
   real(dp)                          :: storage      ! storage [m3]
   ! Local variables
   real(dp)                          :: a            ! channel area [m2]
   real(dp)                          :: zf1          ! floodplain slope: horizontal:vertcal=1:zf [-]
   real(dp)                          :: bankDepth1   ! bankfull depth [m]

   if (present(bankDepth)) then
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   a = flow_area(yin, b, zc, zf=zf1, bankDepth=bankDepth1)
   storage = a*l

 END FUNCTION storage

 ! *********************************************************************
 ! public function: uniform flow [m3/s]
 ! *********************************************************************
 elemental FUNCTION uniformFlow(yin, b, zc, S, n, zf, bankDepth)
   implicit none
   ! Argument variables
   real(dp), intent(in)              :: yin          ! flow depth [m]
   real(dp), intent(in)              :: b            ! channel bottom width [m2] 0-> triangular channel
   real(dp), intent(in)              :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), intent(in)              :: S            ! channel slopw [-]
   real(dp), intent(in)              :: n            ! manning n [-]
   real(dp), optional, intent(in)    :: zf           ! floodplain slope: horizontal:vertcal=1:zf [-]
   real(dp), optional, intent(in)    :: bankDepth    ! bankfull depth [m]
   real(dp)                          :: uniformFlow  ! uniform flow [m3/s]
   ! Local variables
   real(dp)                          :: Abf          ! flow area at bankful depth [m2]
   real(dp)                          :: Pbf          ! wetted perimeter at bankful depth [m2]
   real(dp)                          :: Bbf          ! top width of at bankful depth [m2]
   real(dp)                          :: A            ! flow area [m2]
   real(dp)                          :: P            ! wetted perimeter [m]
   real(dp)                          :: A_channel    ! flow area within a channel [m2]
   real(dp)                          :: P_channel    ! wetted perimeter within a channel [m]
   real(dp)                          :: Q_channel    ! uniform flow within a channel [m3/s]
   real(dp)                          :: A_fp         ! flow area within floodplain [m2]
   real(dp)                          :: P_fp         ! wetted perimeter within floodplain [m]
   real(dp)                          :: Q_fp         ! uniform flow within floodplain [m3/s]
   real(dp)                          :: y_excess     ! flow depth above bankfull depth [m]
   real(dp)                          :: zf1          ! floodplain slope: horizontal:vertcal=1:zf [-]
   real(dp)                          :: bankDepth1   ! bankfull depth [m]

   if (present(bankDepth)) then
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   if (yin <= bankDepth1) then
     A = flow_area(yin, b, zc, zf=zf1, bankDepth=bankDepth1)
     P = Pwet(yin, b, zc, zf=zf1, bankdepth=bankDepth1)
     uniformFlow = A*(A/P)**const23 *sqrt(S)/n
   else
     Abf = flow_area(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
     Pbf = Pwet(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
     Bbf = Btop(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
     y_excess = yin - bankDepth1
     A_channel = Abf + Bbf*y_excess ! imaginary trapizoidal
     P_channel = Pbf                ! imaginary trapizoidal
     Q_channel = A_channel*(A_channel/P_channel)**const23*sqrt(S)/n

     A_fp      = y_excess*zf1*y_excess/2
     P_fp      = y_excess *sqrt(1+ zf1*zf1)
     Q_fp      = 2*(A_fp*(A_fp/P_fp)**const23*sqrt(S)/n)
     uniformFlow = Q_channel+Q_fp
   end if

 END FUNCTION uniformFlow

 ! *********************************************************************
 ! public function: flow depth [m]
 ! *********************************************************************
 elemental FUNCTION flow_depth(Qin, b, zc, S, n, zf, bankDepth)

   ! Solve y for 0 = f(y)/g(y) - Q
   ! h(y) = f(y)/g(y)-Q
   ! use Newton–Raphson

   implicit none
   ! Argument variables
   real(dp), intent(in)              :: Qin          ! uniform flow [m3/s]
   real(dp), intent(in)              :: b            ! channel bottom width [m2] 0-> triangular channel
   real(dp), intent(in)              :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), intent(in)              :: S            ! channel slopw [-]
   real(dp), intent(in)              :: n            ! manning n [-]
   real(dp), optional, intent(in)    :: zf           ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp), optional, intent(in)    :: bankDepth    ! bankfull depth [m]
   real(dp)                          :: flow_depth   ! flow depth [m]
   ! Local variables
   logical(lgt)                     :: floodplain   ! logical to account for floodplain
   real(dp)                         :: Coef1        ! constant
   real(dp)                         :: Coef2        ! constant
   real(dp)                         :: y0           ! estimated normal depth [m]
   real(dp)                         :: error        ! convergence error in depth [fraction]
   real(dp)                         :: Abf          ! flow area at bankful depth [m2]
   real(dp)                         :: Pbf          ! wetted perimeter at bankful depth [m2]
   real(dp)                         :: Bbf          ! top width of at bankful depth [m2]
   real(dp)                         :: Qbf          ! uniform flow at bankful depth [m3/s]
   real(dp)                         :: A            ! estimated flow area [m2]
   real(dp)                         :: P            ! estimtated wetted perimeter [m2]
   real(dp)                         :: Bt           ! estimated top width [m2]
   real(dp)                         :: y_excess     ! flow depth above bankfull depth [m]
   real(dp)                         :: h            !
   real(dp)                         :: dhdy         !
   real(dp)                         :: zf1          ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp)                         :: bankDepth1   ! bankfull depth [m]

   error = 100._dp ! initial dummy error

   floodplain = .false.
   if (present(bankDepth)) then
     floodplain = .true.
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   flow_depth=0._dp
   if (Qin>0._dp) then
     if (floodplain) then ! if floodplain exists

       ! compute flow area, wetted perimeter, top width, and uniform flow at bankful depth
       Abf = flow_area(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
       Pbf = Pwet(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
       Bbf = Btop(bankDepth1, b, zc, zf=zf1, bankDepth=bankDepth1)
       Qbf = Abf*(Abf/Pbf)**const23*sqrt(S)/n

       if (Qin < Qbf) then ! if Q < than Bankfull uniform flow (only channel)

         Coef1 = (sqrt(S)/n/Qin)**3
         Coef2 = 2*sqrt(zc*zc+1._dp) ! derivative of P w.r.t. y

         y0 = (1._dp/Coef1/b**3)**(1._dp/5._dp) ! first guess - based on rectanglur and b>>y

         do while (error > err_thresh) ! Newton–Raphson for flow depth below bankfull depth
           A  = flow_area(y0, b, zc, zf=zf1, bankDepth=bankDepth1)
           Bt = Btop(y0, b, zc, zf=zf1, bankDepth=bankDepth1)
           P  = Pwet(y0, b ,zc, zf=zf1, bankDepth=bankDepth1)

           h    = Coef1*A**5/P**2 -1._dp
           dhdy = Coef1*(5*A**4*Bt*P - 2*Coef2*A**5)/P**3

           flow_depth = y0 - h/dhdy
           error = abs((flow_depth-y0)/flow_depth)

           y0 = flow_depth
         end do

       else ! if Q > than Bankfull uniform flow

         y0 = bankDepth1 + 2._dp ! first guess - add 2 meter to bankfull depth

         ! some constants in manning equation and its derivative given channel parameter
         Coef1 = sqrt(S)/n/Pbf**const23
         Coef2 = 2*(zf1/2)**const53*sqrt(S)/n/(zf1*zf1+1._dp)**const13

         do while (error > err_thresh)
           y_excess = y0-bankDepth1 ! Newton–Raphson for flow depth above bankfull depth
           h    = Coef1*(Abf+Bbf*y_excess)**const53+ Coef2*y_excess**const103/y_excess**const23 -Qin
           dhdy = Coef1*const53*Bbf*(Abf+Bbf*y_excess)**const23 + Coef2*(const103-const23)*y_excess**const53

           flow_depth = y0 - h/dhdy
           error = abs((flow_depth-y0)/flow_depth)

           y0=flow_depth
         end do

       end if

     else ! no floodplain (only chaneel)

       Coef1 = (sqrt(S)/n/Qin)**3
       Coef2 = 2*sqrt(zc*zc+1._dp) ! derivative of P w.r.t. y

       y0 = (1._dp/Coef1/b**3)**(1._dp/5._dp) ! first guess - based on rectanglur and b>>y

       do while (error > err_thresh) ! Newton–Raphson for flow depth below bankfull depth
         A  = flow_area(y0, b, zc)
         Bt = Btop(y0, b, zc)
         P  = Pwet(y0, b ,zc)

         h    = Coef1*A**5/P**2 -1._dp
         dhdy = Coef1*(5*A**4*Bt*P - 2*Coef2*A**5)/P**3

         flow_depth = y0 - h/dhdy
         error = abs((flow_depth-y0)/flow_depth)

         y0 = flow_depth
       end do

     end if
   end if

 END FUNCTION flow_depth

 ! *********************************************************************
 ! public function: wave celerity [m/s]
 ! *********************************************************************
 elemental FUNCTION celerity(Qin, y, b, zc, S, n, zf, bankDepth)
   implicit none
   ! Argument variables
   real(dp), intent(in)               :: Qin          ! uniform flow [m3/s]
   real(dp), intent(in)               :: y            ! flow depth [m]
   real(dp), intent(in)               :: b            ! channel bottom width [m2]
   real(dp), intent(in)               :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), intent(in)               :: S            ! channel slopw [-]
   real(dp), intent(in)               :: n            ! manning n [-]
   real(dp), optional, intent(in)     :: zf           ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp), optional, intent(in)     :: bankDepth    ! bankfull depth [m]
   real(dp)                           :: celerity     ! celerity [m/s]
   ! Local variables
   real(dp)                           :: A            ! estimated flow area [m2]
   real(dp)                           :: P            ! estimtated wetted perimeter [m2]
   real(dp)                           :: Bt           ! estimated top width [m2]
   real(dp)                           :: Sf           ! friction slope [-]
   real(dp)                           :: zf1          ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp)                           :: bankDepth1   ! bankfull depth [m]

   if (present(bankDepth)) then
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   if (y>0._dp) then
     Bt = Btop(y, b, zc, zf=zf1, bankDepth=bankDepth1)
     if (useFrictionSlope) then
       A  = flow_area(y, b, zc, zf=zf1, bankDepth=bankDepth1)
       P  = Pwet(y, b ,zc, zf=zf1, bankDepth=bankDepth1)
       Sf = (Qin*n/A/(A/P)**const23)**2
     else
       Sf = S ! friction slope is approximated with channel slope
     end if
     celerity = const53*Sf**(0.3_dp)*Qin**(0.4_dp)/Bt**(0.4_dp)/n**(0.6_dp)
   else
     celerity = 0.0_dp
   end if

 END FUNCTION celerity

 ! *********************************************************************
 ! public function: wave diffusivity [m2/s]
 ! *********************************************************************
 elemental FUNCTION diffusivity(Qin, y, b, zc, S, n, zf, bankDepth)
   implicit none
   ! Argument variables
   real(dp), intent(in)               :: Qin          ! uniform flow [m3/s]
   real(dp), intent(in)               :: y            ! flow depth [m]
   real(dp), intent(in)               :: b            ! channel bottom width [m2]
   real(dp), intent(in)               :: zc           ! main channel side slope: horizontal:vertcal=zc:1 [-] 0-> rectangular channel
   real(dp), intent(in)               :: S            ! channel slopw [-]
   real(dp), intent(in)               :: n            ! manning n [-]
   real(dp), optional, intent(in)     :: zf           ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp), optional, intent(in)     :: bankDepth    ! bankfull depth [m]
   real(dp)                           :: diffusivity  ! diffusivity [m2/s]
   ! Local variables
   real(dp)                           :: A            ! estimated flow area [m2]
   real(dp)                           :: P            ! estimtated wetted perimeter [m2]
   real(dp)                           :: Bt           ! estimated top width [m2]
   real(dp)                           :: Sf           ! friction slope [-]
   real(dp)                           :: zf1          ! floodplain slope: horizontal:vertcal=zf:1 [-]
   real(dp)                           :: bankDepth1   ! bankfull depth [m]

   if (present(bankDepth)) then
     bankDepth1 = bankDepth
   else
     bankDepth1 = bankDepth0
   end if

   if (present(zf)) then
     zf1 = zf
   else
     zf1 = zf0
   end if

   if (y>0._dp) then
     Bt = Btop(y, b, zc, zf=zf1, bankDepth=bankDepth1)
     if (useFrictionSlope) then
       A  = flow_area(y, b, zc, zf=zf1, bankDepth=bankDepth1)
       P  = Pwet(y, b ,zc, zf=zf1, bankDepth=bankDepth1)
       Sf = (Qin*n/A/(A/P)**const23)**2
     else
       Sf = S ! friction slope is approximated with channel slope
     end if
     diffusivity = abs(Qin)/Sf/Bt/2._dp
   else
     diffusivity = 0.0_dp
   end if

 END FUNCTION diffusivity


END MODULE hydraulic
