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
 SUBROUTINE IRF_route_basin(iens,          & ! input: ensemble index
                            RCHFLX_out,    & ! inout: reach flux data structure
                            ierr, message, & ! output: error control
                            ixSubRch)        ! optional input: subset of reach indices to be processed

 ! External modules
 USE globalData,        ONLY : FRAC_FUTURE
 USE nr_utility_module, ONLY : arth

 implicit none

 ! input
 integer(i4b), intent(in)                 :: iens            ! ith ensemble
 ! inout
 type(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output
 integer(I4B), intent(out)                :: ierr            ! error code
 character(*), intent(out)                :: message         ! error message
 ! input (optional)
 integer(i4b), intent(in),   optional     :: ixSubRch(:)     ! subset of reach indices to be processed
 ! local variables
 integer(i4b)                             :: ntdh            ! number of future time step
 integer(i4b)                             :: nSeg            ! number of reaches to be processed
 integer(i4b)                             :: iSeg, jSeg      ! reach loop indix
 integer(i4b), allocatable                :: ixRch(:)        ! a list of reach indices to be processed
 character(len=strLen)                    :: cmessage        ! error message from subroutines

 ierr=0; message='IRF_route_basin/'

 ntdh = size(FRAC_FUTURE)        ! number of UH future times

 ! if a subset of reaches is processed
 if (present(ixSubRch))then
  nSeg=size(ixSubRch)
  allocate(ixRch(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for [ixRch]'; return; endif
  ixRch = ixSubRch
 ! if all the reaches are processed
 else
  nSeg = size(RCHFLX_out(iens,:))
  allocate(ixRch(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for [ixRch]'; return; endif
  ixRch = arth(1,1,nSeg)
 endif

 do iSeg=1,nSeg

  jSeg = ixRch(iSeg)

  if (.not.allocated(RCHFLX_out(iens,jSeg)%QFUTURE))then

   allocate(RCHFLX_out(iens,jSeg)%QFUTURE(ntdh), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate space for RCHFLX_out(iens,iSeg)%QFUTURE'; return; endif

   RCHFLX_out(iens,jSeg)%QFUTURE(:) = 0._dp

  end if

  call irf_conv(FRAC_FUTURE,                     &  ! input: pre-computed normalized UH
                RCHFLX_out(iens,jSeg)%BASIN_QI,  &  ! input: basin average instantaneous runoff
                RCHFLX_out(iens,jSeg)%QFUTURE,   &  ! inout: Update convoluted runoff
                RCHFLX_out(iens,jSeg)%BASIN_QR,  &  ! inout: delayed runoff to segment at current and previous time step
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
