MODULE basinUH_module

USE nrtype
USE public_var
USE dataTypes,          only : STRFLX         ! fluxes in each reach

implicit none

private

public::IRF_route_basin

CONTAINS

 ! ---------------------------------------------------------------------------------------
 ! Public subroutine main driver for basin routing
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE IRF_route_basin(&
                            iens,       &    ! input: ensemble index
                            nSeg,       &    ! input: number of reaches
                            ! inout
                            RCHFLX_out, &    ! inout: reach flux data structure
                            ! output
                            ierr, message)   ! output: error control
 ! I created wrapper for basin routing because of potential other routing methods
 ! External modules
 USE globalData, ONLY : FRAC_FUTURE
 implicit none
 ! input
 integer(i4b), intent(in)                 :: iens            ! ith ensemble
 integer(i4b), intent(in)                 :: nSeg            ! number of reach segments in the network
 ! inout
 type(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output
 integer(I4B), intent(out)                :: ierr            ! error code
 character(*), intent(out)                :: message         ! error message
 ! local variables
 integer(i4b)                             :: ntdh            ! number of future time step
 integer(i4b)                             :: ibas            ! ith ensemble
 character(len=strLen)                    :: cmessage        ! error message from subroutines

 ierr=0; message='IRF_route_basin/'

 ntdh = size(FRAC_FUTURE)

  do ibas=1,nSeg

   if (.not.allocated(RCHFLX_out(iens,ibas)%QFUTURE))then

    allocate(RCHFLX_out(iens,ibas)%QFUTURE(ntdh), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'unable to allocate space for RCHFLX_out(iens,ibas)%QFUTURE'; return; endif

    RCHFLX_out(iens,ibas)%QFUTURE(:) = 0._dp

   end if

   call irf_conv(FRAC_FUTURE,                     &  ! input: pre-computed normalized UH
                 RCHFLX_out(iens,ibas)%BASIN_QI,  &  ! input: basin average instantaneous runoff
                 RCHFLX_out(iens,ibas)%QFUTURE,   &  ! inout: Update convoluted runoff
                 RCHFLX_out(iens,ibas)%BASIN_QR,  &  ! inout: delayed runoff to segment at current and previous time step
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


END MODULE basinUH_module
