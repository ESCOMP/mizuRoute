module irf_route_module

!numeric type
USE nrtype
! data type
USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO        ! Network topology
! global parameters
USE public_var, ONLY: iulog          ! i/o logical unit number
USE public_var, ONLY: realMissing    ! missing value for real number
USE public_var, ONLY: integerMissing ! missing value for integer number
! subroutines: general
USE perf_mod,  ONLY: t_startf,t_stopf ! timing start/stop

! privary
implicit none
private

public::irf_route

contains

 ! *********************************************************************
 ! subroutine: perform network UH routing
 ! *********************************************************************
 subroutine irf_route(iEns,          &  ! input: index of runoff ensemble to be processed
                      river_basin,   &  ! input: river basin information (mainstem, tributary outlet etc.)
                      ixDesire,      &  ! input: reachID to be checked by on-screen pringing
                      NETOPO_in,     &  ! input: reach topology data structure
                      RCHFLX_out,    &  ! inout: reach flux data structure
                      ierr, message, &  ! output: error control
                      ixSubRch)         ! optional input: subset of reach indices to be processed

 ! global routing data
 USE dataTypes, ONLY: subbasin_omp   ! mainstem+tributary data structures
 USE model_utils, ONLY: handle_err

 implicit none
 ! Input
 integer(i4b),       intent(in)                  :: iEns                ! runoff ensemble to be routed
 type(subbasin_omp), intent(in),    allocatable  :: river_basin(:)      ! river basin information (mainstem, tributary outlet etc.)
 integer(i4b),       intent(in)                  :: ixDesire            ! index of the reach for verbose output ! Output
 type(RCHTOPO),      intent(in),    allocatable  :: NETOPO_in(:)        ! River Network topology
 ! inout
 TYPE(STRFLX),       intent(inout), allocatable  :: RCHFLX_out(:,:)     ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output variables
 integer(i4b),       intent(out)                 :: ierr                ! error code
 character(*),       intent(out)                 :: message             ! error message
 ! input (optional)
 integer(i4b),       intent(in), optional        :: ixSubRch(:)         ! subset of reach indices to be processed
 ! Local variables
 character(len=strLen)                           :: cmessage            ! error message from subroutine
 logical(lgt),                      allocatable  :: doRoute(:)          ! logical to indicate which reaches are processed
 integer(i4b)                                    :: nOrder              ! number of stream order
 integer(i4b)                                    :: nTrib               ! number of tributary basins
 integer(i4b)                                    :: nSeg                ! number of reaches in the network
 integer(i4b)                                    :: iSeg, jSeg          ! loop indices - reach
 integer(i4b)                                    :: iTrib               ! loop indices - branch
 integer(i4b)                                    :: ix                  ! loop indices stream order

 ! initialize error control
 ierr=0; message='irf_route/'

 ! number of reach check
 if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
  ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
 endif

 nSeg = size(RCHFLX_out(iens,:))

 allocate(doRoute(nSeg), stat=ierr)

 ! Initialize CHEC_IRF to False.
 RCHFLX_out(iEns,:)%isRoute=.False.

 if (present(ixSubRch))then
  doRoute(:)=.false.
  doRoute(ixSubRch) = .true. ! only subset of reaches are on
 else
  doRoute(:)=.true. ! every reach is on
 endif

 nOrder = size(river_basin)

 call t_startf('route/irf')

 do ix = 1,nOrder

   nTrib=size(river_basin(ix)%branch)

  ! 1. Route tributary reaches (parallel)
!$OMP PARALLEL DO schedule(dynamic,1)                   &
!$OMP          private(jSeg, iSeg)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(doRoute)                          & ! data array shared
!$OMP          shared(NETOPO_in)                        & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(ix, iEns, ixDesire)               & ! indices shared
!$OMP          firstprivate(nTrib)
   trib:do iTrib = 1,nTrib
     seg:do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
       jSeg = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
       if (.not. doRoute(jSeg)) cycle
       call segment_irf(iEns, jSeg, ixDesire, NETOPO_IN, RCHFLX_out, ierr, cmessage)
       if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
     end do seg
   end do trib
!$OMP END PARALLEL DO

 end do

 call t_stopf('route/irf')

 end subroutine irf_route


 ! *********************************************************************
 ! subroutine: perform one segment route UH routing
 ! *********************************************************************
 subroutine segment_irf(&
                        ! input
                        iEns,       &    ! input: index of runoff ensemble to be processed
                        segIndex,   &    ! input: index of runoff ensemble to be processed
                        ixDesire,   &    ! input: reachID to be checked by on-screen pringing
                        NETOPO_in,  &    ! input: reach topology data structure
                        ! inout
                        RCHFLX_out, &    ! inout: reach flux data structure
                        ! output
                        ierr, message)   ! output: error control

 implicit none
 ! Input
 INTEGER(I4B), intent(IN)                 :: iEns           ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)                 :: segIndex       ! segment where routing is performed
 INTEGER(I4B), intent(IN)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)   ! River Network topology
 ! inout
 TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! Output
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables to
 type(STRFLX), allocatable                :: uprflux(:)     ! upstream Reach fluxes
 INTEGER(I4B)                             :: nUps           ! number of upstream segment
 INTEGER(I4B)                             :: iUps           ! upstream reach index
 INTEGER(I4B)                             :: iRch_ups       ! index of upstream reach in NETOPO
 INTEGER(I4B)                             :: ntdh           ! number of time steps in IRF
 character(len=strLen)                    :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='segment_irf/'

 ! route streamflow through the river network
  if (.not.allocated(RCHFLX_out(iens,segIndex)%QFUTURE_IRF))then

   ntdh = size(NETOPO_in(segIndex)%UH)

   allocate(RCHFLX_out(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX_out(iens,segIndex)%QFUTURE_IRF'; return; endif

   RCHFLX_out(iens,segIndex)%QFUTURE_IRF(:) = 0._dp

  end if

  ! identify number of upstream segments of the reach being processed
  nUps = size(NETOPO_in(segIndex)%UREACHI)

  allocate(uprflux(nUps), stat=ierr, errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': uprflux'; return; endif

  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      uprflux(iUps) = RCHFLX_out(iens,iRch_ups)
    end do
  endif

  ! perform river network UH routing
  call conv_upsbas_qr(NETOPO_in(segIndex)%UH,    &    ! input: reach unit hydrograph
                      uprflux,                   &    ! input: upstream reach fluxes
                      RCHFLX_out(iens,segIndex), &    ! inout: updated fluxes at reach
                      ierr, message)                  ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! Check True since now this reach now routed
  RCHFLX_out(iEns,segIndex)%isRoute=.True.

  ! check
  if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   write(iulog,*) 'RCHFLX_out(iens,segIndex)%BASIN_QR(1),RCHFLX_out(iens,segIndex)%REACH_Q_IRF = ', &
                   RCHFLX_out(iens,segIndex)%BASIN_QR(1),RCHFLX_out(iens,segIndex)%REACH_Q_IRF
  endif

 end subroutine segment_irf


 ! *********************************************************************
 ! subroutine: Compute delayed runoff from the upstream segments
 ! *********************************************************************
 subroutine conv_upsbas_qr(reach_uh,   &    ! input: reach unit hydrograph
                           rflux_ups,  &    ! input: upstream reach fluxes
                           rflux,      &    ! input: input flux at reach
                           ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Details: Convolute runoff volume of upstream at one reach at one time step
 ! ----------------------------------------------------------------------------------------
 USE public_var, ONLY: dt
 implicit none
 ! Input
 real(dp),     intent(in)               :: reach_uh(:)  ! reach unit hydrograph
 type(STRFLX), intent(in)               :: rflux_ups(:) ! upstream Reach fluxes
 ! inout
 type(STRFLX), intent(inout)            :: rflux        ! current Reach fluxes
 ! Output
 integer(i4b), intent(out)              :: ierr         ! error code
 character(*), intent(out)              :: message      ! error message
 ! Local variables to
 real(dp)                               :: q_upstream   ! total discharge at top of the reach being processed
 INTEGER(I4B)                           :: nTDH         ! number of UH data
 INTEGER(I4B)                           :: iTDH         ! index of UH data (i.e.,future time step)
 INTEGER(I4B)                           :: nUps         ! number of all upstream segment
 INTEGER(I4B)                           :: iUps         ! loop indices for u/s reaches

 ! initialize error control
 ierr=0; message='conv_upsbas_qr/'

 ! identify number of upstream segments of the reach being processed
 nUps = size(rflux_ups)

 ! Find out total q at top of a segment
 q_upstream = 0.0_dp
 if(nUps>0)then
   do iUps = 1,nUps
     q_upstream = q_upstream + rflux_ups(iUps)%REACH_Q_IRF
   end do
 endif

 ! place a fraction of runoff in future time steps
 nTDH = size(reach_uh) ! identify the number of future time steps of UH for a given segment
 do iTDH=1,nTDH
   rflux%QFUTURE_IRF(iTDH) = rflux%QFUTURE_IRF(iTDH) &
                             + reach_uh(iTDH)*q_upstream
 enddo

 ! Volume stored in reach
 rflux%REACH_VOL(0) = rflux%REACH_VOL(1)
 rflux%REACH_VOL(1) = rflux%REACH_VOL(0) - (rflux%QFUTURE_IRF(1) - q_upstream)*dt

 ! Add local routed flow
 rflux%REACH_Q_IRF = rflux%QFUTURE_IRF(1) + rflux%BASIN_QR(1)

 ! move array back   use eoshift
 !rflux%QFUTURE_IRF=eoshift(rflux%QFUTURE_IRF,shift=1)

 do itdh=2,ntdh
  rflux%QFUTURE_IRF(itdh-1) = rflux%QFUTURE_IRF(itdh)
 enddo

 rflux%QFUTURE_IRF(ntdh) = 0._dp

 end subroutine conv_upsbas_qr

end module irf_route_module

