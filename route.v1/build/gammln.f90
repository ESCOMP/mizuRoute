  FUNCTION gammln_s(xx)
  USE nrtype; USE nrutil, ONLY : arth,assert
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: xx
  REAL(SP) :: gammln_s
  REAL(SP) :: tmp,x
  REAL(SP) :: stp = 2.5066282746310005_sp
  REAL(SP), DIMENSION(6) :: coef = (/76.18009172947146_sp,&
    -86.50532032941677_sp,24.01409824083091_sp,&
    -1.231739572450155_sp,0.1208650973866179e-2_sp,&
    -0.5395239384953e-5_sp/)
  call assert(xx > 0.0, 'gammln_s arg')
  x=xx
  tmp=x+5.5_sp
  tmp=(x+0.5_sp)*log(tmp)-tmp
  gammln_s=tmp+log(stp*(1.000000000190015_sp+&
    sum(coef(:)/arth(x+1.0_sp,1.0_sp,size(coef))))/x)
  END FUNCTION gammln_s


  FUNCTION gammln_v(xx)
  USE nrtype; USE nrutil, ONLY: assert
  IMPLICIT NONE
  INTEGER(I4B) :: i
  REAL(SP), DIMENSION(:), INTENT(IN) :: xx
  REAL(SP), DIMENSION(size(xx)) :: gammln_v
  REAL(SP), DIMENSION(size(xx)) :: ser,tmp,x,y
  REAL(SP) :: stp = 2.5066282746310005_sp
  REAL(SP), DIMENSION(6) :: coef = (/76.18009172947146_sp,&
    -86.50532032941677_sp,24.01409824083091_sp,&
    -1.231739572450155_sp,0.1208650973866179e-2_sp,&
    -0.5395239384953e-5_sp/)
  if (size(xx) == 0) RETURN
  call assert(all(xx > 0.0), 'gammln_v arg')
  x=xx
  tmp=x+5.5_sp
  tmp=(x+0.5_sp)*log(tmp)-tmp
  ser=1.000000000190015_sp
  y=x
  do i=1,size(coef)
    y=y+1.0_sp
    ser=ser+coef(i)/y
  end do
  gammln_v=tmp+log(stp*ser/x)
  END FUNCTION gammln_v
