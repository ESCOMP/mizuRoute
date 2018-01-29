module gamma_func_module
USE nrtype
USE nr_utility_module,only:arth
! contains functions that should really be part of the fortran standard, but are not
implicit none
private
public::gammp
contains

 ! ******************************************************************************************************************************
 ! public function gammp: compute cumulative probability using the Gamma distribution
 ! ******************************************************************************************************************************
 FUNCTION gammp(a,x)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,x
 REAL(DP) :: gammp
 if (x<a+1.0_dp) then
  gammp=gser(a,x)
 else
  gammp=1.0_dp-gcf(a,x)
 end if
 END FUNCTION gammp

 ! ******************************************************************************************************************************
 ! private function gser: series development of the incomplete Gamma function
 ! ******************************************************************************************************************************
 FUNCTION gser(a,x,gln)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,x
 REAL(DP), OPTIONAL, INTENT(OUT) :: gln
 REAL(DP) :: gser
 INTEGER(I4B), PARAMETER :: ITMAX=100
 REAL(DP), PARAMETER :: EPS=epsilon(x)
 INTEGER(I4B) :: n
 REAL(DP) :: ap,del,summ
 if (x == 0.0) then
  gser=0.0
  RETURN
 end if
 ap=a
 summ=1.0_dp/a
 del=summ
 do n=1,ITMAX
  ap=ap+1.0_dp
  del=del*x/ap
  summ=summ+del
  if (abs(del) < abs(summ)*EPS) exit
 end do
 if (n > ITMAX) stop 'a too large, ITMAX too small in gser'
 if (present(gln)) then
  gln=gammln(a)
  gser=summ*exp(-x+a*log(x)-gln)
 else
  gser=summ*exp(-x+a*log(x)-gammln(a))
 end if
 END FUNCTION gser

 ! ******************************************************************************************************************************
 ! private function gcf: continued fraction development of the incomplete Gamma function
 ! ******************************************************************************************************************************
 FUNCTION gcf(a,x,gln)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,x
 REAL(DP), OPTIONAL, INTENT(OUT) :: gln
 REAL(DP) :: gcf
 INTEGER(I4B), PARAMETER :: ITMAX=100
 REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
 INTEGER(I4B) :: i
 REAL(DP) :: an,b,c,d,del,h
 if (x == 0.0) then
  gcf=1.0
  RETURN
 end if
 b=x+1.0_dp-a
 c=1.0_dp/FPMIN
 d=1.0_dp/b
 h=d
 do i=1,ITMAX
  an=-i*(i-a)
  b=b+2.0_dp
  d=an*d+b
  if (abs(d) < FPMIN) d=FPMIN
  c=b+an/c
  if (abs(c) < FPMIN) c=FPMIN
  d=1.0_dp/d
  del=d*c
  h=h*del
  if (abs(del-1.0_dp) <= EPS) exit
 end do
 if (i > ITMAX) stop 'a too large, ITMAX too small in gcf'
 if (present(gln)) then
  gln=gammln(a)
  gcf=exp(-x+a*log(x)-gln)*h
 else
  gcf=exp(-x+a*log(x)-gammln(a))*h
 end if
 END FUNCTION gcf

 ! ******************************************************************************************************************************
 ! private function gammln: gamma function
 ! ******************************************************************************************************************************
 FUNCTION gammln(xx)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: xx
 REAL(DP) :: gammln
 REAL(DP) :: tmp,x
 REAL(DP) :: stp = 2.5066282746310005_dp
 REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
  -86.50532032941677_dp,24.01409824083091_dp,&
  -1.231739572450155_dp,0.1208650973866179e-2_dp,&
  -0.5395239384953e-5_dp/)
 if(xx <= 0._dp) stop 'xx > 0 in gammln'
 x=xx
 tmp=x+5.5_dp
 tmp=(x+0.5_dp)*log(tmp)-tmp
 gammln=tmp+log(stp*(1.000000000190015_dp+&
  sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
 END FUNCTION gammln

end module gamma_func_module
