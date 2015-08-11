SUBROUTINE QTIMEDELAY(dt, scale_time, IERR, MESSAGE)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007 (FUSE)
! Martyn Clark, 2012 (modified to use allocatable arrays, and to remove option of no routing)
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the fraction of runoff in future time steps
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE reach_flux -- runoff fractions stored in FRAC_FUTURE(:)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE nr, ONLY : gammp                                  ! interface for the incomplete gamma function
USE reach_flux                                        ! reach fluxes
IMPLICIT NONE
! input
REAL(DP), INTENT(IN)                   :: DT          ! model time step
REAL(DP), INTENT(IN)                   :: SCALE_TIME  ! time parameter
! output
INTEGER(I4B), INTENT(OUT)              :: IERR        ! error code
CHARACTER(*), INTENT(OUT)              :: MESSAGE     ! error message
! locals
REAL(SP)                               :: ALPHA       ! shape parameter
REAL(SP)                               :: ALAMB       ! scale parameter
REAL(DP)                               :: ntdh_min    ! minimum number of time delay points
REAL(DP)                               :: ntdh_max    ! maximum number of time delay points
REAL(DP)                               :: ntdh_try    ! trial number of time delay points
INTEGER(I4B)                           :: itry        ! index of trial value
INTEGER(I4B), PARAMETER                :: MAXTRY=100  ! maximum number of trial values
INTEGER(I4B)                           :: NTDH        ! number of values on the time delay histogram
INTEGER(I4B)                           :: JTIM        ! (loop through future time steps)
REAL(SP)                               :: TFUTURE     ! future time (units of dt)
REAL(SP)                               :: X_VALUE     ! xvalue to evaluate using gammp
REAL(SP)                               :: CUMPROB     ! cumulative probability at JTIM
REAL(SP)                               :: PSAVE       ! cumulative probability at JTIM-1
! ---------------------------------------------------------------------------------------
! use a Gamma distribution with shape parameter = 2.5, and time parameter input
ALPHA = 2.5_SP                                             ! shape parameter
ALAMB = ALPHA/REAL(SCALE_TIME, kind(sp))                   ! scale parameter
! find the desired number of future time steps
ntdh_min = 1._dp
ntdh_max = 1000._dp
ntdh_try = 0.5_dp*(ntdh_min + ntdh_max)
do itry=1,maxtry
 x_value = alamb*real(dt*ntdh_try, kind(sp))
 cumprob = gammp(alpha, x_value)
 !print*, scale_time, ntdh_try, cumprob
 if(cumprob < 0.99_dp)  ntdh_min = ntdh_try
 if(cumprob > 0.999_dp) ntdh_max = ntdh_try
 if(cumprob > 0.99_dp .and. cumprob < 0.999_dp) exit
 ntdh_try = 0.5_dp*(ntdh_min + ntdh_max)
 if(itry==maxtry)then; ierr=20; message=trim(message)//'cannot identify the maximum number of bins for the tdh'; return; endif
end do
ntdh = ceiling(ntdh_try)
! allocate space for the time-delay histogram
allocate(FRAC_FUTURE(ntdh), stat=ierr)
if(ierr/=0)then; message=trim(message)//'unable to allocate space for the time delay histogram'; return; endif
! loop through time steps and compute the fraction of runoff in future time steps
PSAVE = 0.                                                 ! cumulative probability at JTIM-1
DO JTIM=1,NTDH
 TFUTURE            = real(REAL(JTIM, kind(dp))*DT, kind(sp))       ! future time
 CUMPROB            = GAMMP(ALPHA,ALAMB*TFUTURE)    ! cumulative probability at JTIM
 FRAC_FUTURE(JTIM)  = MAX(0._DP, CUMPROB-PSAVE)     ! probability between JTIM-1 and JTIM
 PSAVE              = CUMPROB                       ! cumulative probability at JTIM-1
 !WRITE(*,'(I5,1X,F20.5,1X,2(F11.5))') JTIM, TFUTURE, FRAC_FUTURE(JTIM), CUMPROB
END DO
! ensure that the fractions sum to 1.0 (account for rounding errors, and not enough bins)
FRAC_FUTURE(:) = FRAC_FUTURE(:) / SUM(FRAC_FUTURE(:))
! ---------------------------------------------------------------------------------------
END SUBROUTINE QTIMEDELAY
