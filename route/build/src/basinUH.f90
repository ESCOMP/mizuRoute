MODULE basinUH_module

USE nrtype
USE public_var

implicit none

private

public::IRF_route_basin
public::basinUH

CONTAINS

 ! ---------------------------------------------------------------------------------------
 ! Public subroutine main driver for basin routing
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE IRF_route_basin(iens, nSeg, ierr, message)
 ! I created wrapper for basin routing because of potential other routing methods
 ! External modules
 USE globalData, ONLY : FRAC_FUTURE
 USE globalData, ONLY : RCHFLX
 implicit none
 ! input
 integer(i4b), intent(in)           :: iens        ! ith ensemble
 integer(i4b), intent(in)           :: nSeg        ! number of segments = number of basins
 ! output
 integer(I4B), intent(out)          :: ierr        ! error code
 character(*), intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                       :: ntdh        ! number of future time step
 integer(i4b)                       :: ibas        ! ith ensemble
 character(len=strLen)              :: cmessage    ! error message from subroutines

 ierr=0; message='IRF_route_basin/'

 ntdh = size(FRAC_FUTURE)

  do ibas=1,nSeg

   if (.not.allocated(RCHFLX(iens,ibas)%QFUTURE))then

    allocate(RCHFLX(iens,ibas)%QFUTURE(ntdh), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'unable to allocate space for RCHFLX(iens,ibas)%QFUTURE'; return; endif

    RCHFLX(iens,ibas)%QFUTURE(:) = 0._dp

   end if

   call irf_conv(FRAC_FUTURE,                 &  ! input: pre-computed normalized UH
                 RCHFLX(iens,ibas)%BASIN_QI,  &  ! input: basin average instantaneous runoff
                 RCHFLX(iens,ibas)%QFUTURE,   &  ! inout: Update convoluted runoff
                 RCHFLX(iens,ibas)%BASIN_QR,  &  ! inout: delayed runoff to segment at current and previous time step
                 ierr, message)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  end do

 CONTAINS

  ! ---------------------------------------------------------------------------------------
  ! Private subroutine: Perform UH convolutions
  ! ---------------------------------------------------------------------------------------
  SUBROUTINE irf_conv(&! output
                      uh,        &    ! input: normalized unit hydrograph
                      basin_iq,  &    ! input: basin instantaneous runoff
                      ! in-output
                      qfuture,   &    ! inout: convoluted runoff including future time step
                      basin_qr,  &    ! inout: delayed runoff to segment at a current and previous time step
                      ierr, message)
  implicit none
  ! input
  real(dp),             intent(in)     :: uh(:)         ! normalized unit hydrograph
  real(dp),             intent(in)     :: basin_iq      ! basin instantaneous runoff
  real(dp),             intent(inout)  :: qfuture(:)    ! convoluted runoff including future time steps
  real(dp),             intent(inout)  :: basin_qr(0:1) ! delayed runoff to a segment at a current and previous time step
  ! output
  integer(I4B),         intent(out)    :: ierr          ! error code
  character(*),         intent(out)    :: message       ! error message
  ! local variables
  integer(i4b)                         :: itdh          ! index loop for basin, time, respectively
  integer(i4b)                         :: ntdh          ! number of time step for future flow

  ! initialize error control
  ierr=0; message='irf_conv/'

  ! route streamflow through the basin
  ntdh = size(qfuture)

  ! place a fraction of runoff in future time steps and add to current state of q in the future
  do itdh=1,ntdh
   qfuture(itdh) = qfuture(itdh) + uh(itdh)*basin_iq
  end do
  ! save the routed runoff
  basin_qr(0) = basin_qr(1)  ! (save the runoff from the previous time step)
  basin_qr(1) = qfuture(1)
  ! move array back
  do itdh=2,ntdh
   qfuture(itdh-1) = qfuture(itdh)
  end do
  qfuture(ntdh)    = 0._dp

  END SUBROUTINE irf_conv

 END SUBROUTINE IRF_route_basin


  SUBROUTINE basinUH(dt, fshape, tscale, IERR, MESSAGE)
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
  USE gamma_func_module, ONLY : gammp                   ! interface for the incomplete gamma function
  USE globalData,        ONLY : FRAC_FUTURE             ! fraction of runoff in future time steps
  IMPLICIT NONE
  ! input
  REAL(DP), INTENT(IN)                   :: dt          ! model time step
  REAL(DP), INTENT(IN)                   :: fshape      ! shapef parameter in gamma distribution
  REAL(DP), INTENT(IN)                   :: tscale      ! time scale parameter
  ! output
  INTEGER(I4B), INTENT(OUT)              :: IERR        ! error code
  CHARACTER(*), INTENT(OUT)              :: MESSAGE     ! error message
  ! locals
  REAL(DP)                               :: alamb       ! scale parameter
  REAL(DP)                               :: ntdh_min    ! minimum number of time delay points
  REAL(DP)                               :: ntdh_max    ! maximum number of time delay points
  REAL(DP)                               :: ntdh_try    ! trial number of time delay points
  INTEGER(I4B)                           :: itry        ! index of trial value
  INTEGER(I4B), PARAMETER                :: MAXTRY=100  ! maximum number of trial values
  INTEGER(I4B)                           :: NTDH        ! number of values on the time delay histogram
  INTEGER(I4B)                           :: JTIM        ! (loop through future time steps)
  REAL(DP)                               :: TFUTURE     ! future time (units of dt)
  REAL(DP)                               :: X_VALUE     ! xvalue to evaluate using gammp
  REAL(DP)                               :: CUMPROB     ! cumulative probability at JTIM
  REAL(DP)                               :: PSAVE       ! cumulative probability at JTIM-1
  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  ierr=0; message='basinUH/'
  ! use a Gamma distribution with shape parameter, fshape = 2.5, and time parameter, tscale, input
  alamb = fshape/tscale                  ! scale parameter
  ! find the desired number of future time steps
  ntdh_min = 1._dp
  ntdh_max = 1000._dp
  ntdh_try = 0.5_dp*(ntdh_min + ntdh_max)
  do itry=1,maxtry
   x_value = alamb*dt*ntdh_try
   cumprob = gammp(fshape, x_value)
   !print*, tscale, ntdh_try, cumprob
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
   TFUTURE            = REAL(JTIM, kind(dp))*DT       ! future time
   CUMPROB            = gammp(fshape,alamb*TFUTURE)   ! cumulative probability at JTIM
   FRAC_FUTURE(JTIM)  = MAX(0._DP, CUMPROB-PSAVE)     ! probability between JTIM-1 and JTIM
   PSAVE              = CUMPROB                       ! cumulative probability at JTIM-1
   !WRITE(*,'(I5,1X,F20.5,1X,2(F11.5))') JTIM, TFUTURE, FRAC_FUTURE(JTIM), CUMPROB
  END DO
  ! ensure that the fractions sum to 1.0 (account for rounding errors, and not enough bins)
  FRAC_FUTURE(:) = FRAC_FUTURE(:) / SUM(FRAC_FUTURE(:))
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE basinUH


END MODULE basinUH_module
