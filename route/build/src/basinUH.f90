MODULE basinUH_module

USE nrtype
USE public_var
USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO        ! network tiver topology
USE model_utils, ONLY: handle_err

implicit none

private
public::IRF_route_basin

CONTAINS

 ! ---------------------------------------------------------------------------------------
 ! Public subroutine main driver for basin routing
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE IRF_route_basin(NETOPO_in,     & ! input: reach topology
                            RCHFLX_out,    & ! inout: reach flux data structure
                            ierr, message, & ! output: error control
                            ixSubRch)        ! optional input: subset of reach indices to be processed
 implicit none
 ! Argument variables
 type(RCHTOPO),intent(in),   allocatable  :: NETOPO_in(:)    ! River Network topology
 type(STRFLX), intent(inout)              :: RCHFLX_out(:)   ! Reach fluxes (space [reaches]) for decomposed domains
 integer(I4B), intent(out)                :: ierr            ! error code
 character(*), intent(out)                :: message         ! error message
 integer(i4b), intent(in),   optional     :: ixSubRch(:)     ! subset of reach indices to be processed
 ! local variables
 integer(i4b)                             :: nSeg            ! number of reaches to be processed
 integer(i4b)                             :: iSeg            ! reach loop indix
 logical(lgt),               allocatable  :: doRoute(:)      ! logical to indicate which reaches are processed
 character(len=strLen)                    :: cmessage        ! error message from subroutines

 ierr=0; message='IRF_route_basin/'

 nSeg = size(RCHFLX_out)

 allocate(doRoute(nSeg), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'unable to allocate space for [doRoute]'; return; endif

 ! if a subset of reaches is processed
 if (present(ixSubRch))then
  doRoute(:) = .false.
  doRoute(ixSubRch) = .true. ! only subset of reaches are on
 else
  doRoute(:) = .true.
 endif

!$OMP PARALLEL DO schedule(dynamic,1)   &
!$OMP          private(iSeg)            & ! loop index
!$OMP          private(ierr, cmessage)  & ! private for a given thread
!$OMP          shared(doRoute)          & ! data array shared
!$OMP          shared(RCHFLX_out)       & ! data structure shared
!$OMP          shared(NETOPO_in)        & ! data structure shared
!$OMP          firstprivate(nSeg)
 do iSeg=1,nSeg
   if (.not. doRoute(iSeg)) cycle
   call hru_irf(iSeg, NETOPO_in, RCHFLX_out, ierr, cmessage)
   if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
 end do
!$OMP END PARALLEL DO

 END SUBROUTINE IRF_route_basin

 ! *********************************************************************
 ! subroutine: perform one basin UH routing to a iSeg reach at one time
 ! *********************************************************************
 SUBROUTINE hru_irf(iSeg,         &    ! input: index of segment to be processed
                    NETOPO_in,    &    ! input: reach topology
                    RCHFLX_out,   &    ! inout: reach flux data structure
                    ierr, message)     ! output: error control

 USE globalData, ONLY: FRAC_FUTURE     !
 USE public_var, ONLY: is_lake_sim     ! logical whether or not lake should be simulated
 USE public_var, ONLY: tracer          ! logical whether or not tracer is on

 implicit none
 ! Argument variables
 integer(i4b), intent(in)                 :: iSeg                  ! segment where routing is performed
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)          ! River Network topology
 type(STRFLX), intent(inout)              :: RCHFLX_out(:)         ! Reach fluxes (space [reaches]) for decomposed domains
 integer(i4b), intent(out)                :: ierr                  ! error code
 character(*), intent(out)                :: message               ! error message
 ! Local variables
 real(dp),     allocatable                :: FRAC_FUTURE_local(:)  ! local FRAC_FUTURE so that it can be changed for lakes to impulse
 integer(i4b)                             :: ntdh                  ! number of time steps in IRF
 character(len=strLen)                    :: cmessage              ! error message from subroutine

 ierr=0; message='hru_irf/'

 ! initialize the first time step q future
  if (.not.allocated(RCHFLX_out(iSeg)%QFUTURE))then
   ntdh = size(FRAC_FUTURE)

   allocate(RCHFLX_out(iSeg)%QFUTURE(ntdh), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate space for RCHFLX_out(segIndex)%QFUTURE'; return; endif

   RCHFLX_out(iSeg)%QFUTURE(:) = 0._dp
  end if

  if (.not.allocated(RCHFLX_out(iSeg)%solute_future) .and. tracer)then
    ntdh = size(FRAC_FUTURE)
    allocate(RCHFLX_out(iSeg)%solute_future(ntdh), source=0._dp, stat=ierr)
    if(ierr/=0)then; message=trim(message)//'unable to allocate space for RCHFLX_out(segIndex)%solute_future'; return; endif
  end if

  allocate(FRAC_FUTURE_local, source=FRAC_FUTURE, stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for FRAC_FUTURE_local'; return; endif

  ! if the segment is flaged as lake and is_lake is on then no lagged flow for lakes
  if ((NETOPO_in(iSeg)%islake).and.(is_lake_sim)) then;
    FRAC_FUTURE_local(:) = 0._dp
    FRAC_FUTURE_local(1) = 1._dp
  endif

  ! perform river network UH routing
  RCHFLX_out(iSeg)%BASIN_QR(0) = RCHFLX_out(iSeg)%BASIN_QR(1)
  call irf_conv(FRAC_FUTURE_local,            &    ! input: unit hydrograph
                RCHFLX_out(iSeg)%BASIN_QI,    &    ! input: upstream fluxes
                RCHFLX_out(iSeg)%QFUTURE,     &    ! inout: updated q future time series
                RCHFLX_out(iSeg)%BASIN_QR(1), &    ! inout: updated fluxes at reach
               ierr, message)                            ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if(tracer) then
    ! perform river network UH routing
    call irf_conv(FRAC_FUTURE_local,                       &  ! input: unit hydrograph
                  RCHFLX_out(iSeg)%BASIN_solute_inst, &  ! input: upstream fluxes
                  RCHFLX_out(iSeg)%solute_future,     &  ! inout: updated solute future time series
                  RCHFLX_out(iSeg)%BASIN_solute,      &  ! inout: updated fluxes at reach
                 ierr, message)                               ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if

 END SUBROUTINE hru_irf


 ! ---------------------------------------------------------------------------------------
 ! Private subroutine: Perform UH convolutions
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE irf_conv(uh,        &    ! input: normalized unit hydrograph
                     inq,       &    ! input: instantaneous runoff
                     qfuture,   &    ! inout: convoluted runoff including future time step
                     delayq,    &    ! inout: delayed runoff to segment at a current and previous time step
                     ierr, message)
  implicit none
  ! Argument variables
  real(dp),             intent(in)     :: uh(:)         ! normalized unit hydrograph
  real(dp),             intent(in)     :: inq           ! basin instantaneous runoff
  real(dp),             intent(inout)  :: qfuture(:)    ! convoluted runoff including future time steps
  real(dp),             intent(inout)  :: delayq        ! delayed runoff to a segment at a current and previous time step
  integer(I4B),         intent(out)    :: ierr          ! error code
  character(*),         intent(out)    :: message       ! error message
  ! local variables
  integer(i4b)                         :: itdh          ! index loop for basin, time, respectively
  integer(i4b)                         :: ntdh          ! number of time step for future flow

  ierr=0; message='irf_conv/'

  ntdh = size(qfuture)

  ! place a fraction of runoff in future time steps and add to current state of q in the future
  do itdh=1,ntdh
   qfuture(itdh) = qfuture(itdh) + uh(itdh)*inq
  end do

  ! save the routed runoff
  delayq = qfuture(1)

  ! move array back
  do itdh=2,ntdh
   qfuture(itdh-1) = qfuture(itdh)
  end do
  qfuture(ntdh)    = 0._dp

  END SUBROUTINE irf_conv

END MODULE basinUH_module
