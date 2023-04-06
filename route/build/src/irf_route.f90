MODULE irf_route_module

USE nrtype
! data type
USE dataTypes, ONLY: STRFLX           ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO          ! Network topology
USE dataTypes, ONLY: RCHPRP           ! reach/lake property parameter
USE dataTypes, ONLY: irfRCH           ! irf specific state data structure
USE dataTypes, ONLY: subbasin_omp     ! mainstem+tributary data structures
! global parameters
USE public_var, ONLY: iulog           ! i/o logical unit number
USE public_var, ONLY: realMissing     ! missing value for real number
USE public_var, ONLY: integerMissing  ! missing value for integer number
USE public_var, ONLY: dt
USE globalData, ONLY: idxIRF          ! index of IRF method
! external subroutines
USE lake_route_module, ONLY: lake_route        ! lake route module
USE water_balance,     ONLY: comp_reach_wb     ! compute water balance error
USE perf_mod,          ONLY: t_startf,t_stopf  ! timing start/stop
USE model_utils,       ONLY: handle_err

implicit none

private
public::irf_route

CONTAINS

 ! *********************************************************************
 ! subroutine: perform network UH routing
 ! *********************************************************************
 SUBROUTINE irf_route(iEns,          &  ! input: index of runoff ensemble to be processed
                      river_basin,   &  ! input: river basin information (mainstem, tributary outlet etc.)
                      ixDesire,      &  ! input: reachID to be checked by on-screen pringing
                      NETOPO_in,     &  ! input: reach topology data structure
                      RPARAM_in,     &  ! input: reach parameter data structure
                      RCHFLX_out,    &  ! inout: reach flux data structure
                      ierr, message, &  ! output: error control
                      ixSubRch)         ! optional input: subset of reach indices to be processed

 USE public_var,  ONLY: is_lake_sim    ! logical whether or not lake should be simulated

 implicit none
 ! argument variables
 integer(i4b),       intent(in)                  :: iEns                ! runoff ensemble to be routed
 type(subbasin_omp), intent(in),    allocatable  :: river_basin(:)      ! river basin information (mainstem, tributary outlet etc.)
 integer(i4b),       intent(in)                  :: ixDesire            ! index of the reach for verbose output ! Output
 type(RCHTOPO),      intent(in),    allocatable  :: NETOPO_in(:)        ! River Network topology
 type(RCHPRP),       intent(inout), allocatable  :: RPARAM_in(:)        ! River Network parameters
 TYPE(STRFLX),       intent(inout)               :: RCHFLX_out(:,:)     ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b),       intent(out)                 :: ierr                ! error code
 character(*),       intent(out)                 :: message             ! error message
 integer(i4b),       intent(in),    optional     :: ixSubRch(:)         ! subset of reach indices to be processed
 ! Local variables
 character(len=strLen)                           :: cmessage            ! error message from subroutine
 logical(lgt),                      allocatable  :: doRoute(:)          ! logical to indicate which reaches are processed
 integer(i4b)                                    :: nOrder              ! number of stream order
 integer(i4b)                                    :: nTrib               ! number of tributary basins
 integer(i4b)                                    :: nSeg                ! number of reaches in the network
 integer(i4b)                                    :: iSeg, jSeg          ! loop indices - reach
 integer(i4b)                                    :: iTrib               ! loop indices - branch
 integer(i4b)                                    :: ix                  ! loop indices stream order

 ierr=0; message='irf_route/'

 ! number of reach check
 if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
  ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
 endif

 nSeg = size(RCHFLX_out(iens,:))

 allocate(doRoute(nSeg), stat=ierr)

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
!$OMP          shared(RPARAM_in)                        & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(ix, iEns, ixDesire)               & ! indices shared
!$OMP          firstprivate(nTrib)
   trib:do iTrib = 1,nTrib
     seg:do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
       jSeg = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
       if (.not. doRoute(jSeg)) cycle
       if ((NETOPO_in(jseg)%islake).and.(is_lake_sim)) then
        call lake_route (iEns, jSeg, ixDesire, NETOPO_in, RPARAM_in, RCHFLX_out, ierr, cmessage)
       else
        call irf_rch(iEns, jSeg, ixDesire, NETOPO_in,            RCHFLX_out, ierr, cmessage)
       endif
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
 SUBROUTINE irf_rch(iEns,         & ! input: index of runoff ensemble to be processed
                    segIndex,     & ! input: index of runoff ensemble to be processed
                    ixDesire,     & ! input: reachID to be checked by on-screen pringing
                    NETOPO_in,    & ! input: reach topology data structure
                    RCHFLX_out,   & ! inout: reach flux data structure
                    ierr, message)  ! output: error control

 USE public_var, ONLY:is_flux_wm   ! logical water management components fluxes should be read

 implicit none
 ! Argument variables
 integer(i4b), intent(in)                 :: iEns           ! runoff ensemble to be routed
 integer(i4b), intent(in)                 :: segIndex       ! segment where routing is performed
 integer(i4b), intent(in)                 :: ixDesire       ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:)   ! River Network topology
 type(STRFLX), intent(inout)              :: RCHFLX_out(:,:)! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables
 logical(lgt)                             :: verbose        ! check details of variables
 real(dp)                                 :: q_upstream     ! total discharge at top of the reach being processed
 real(dp)                                 :: Qabs           ! maximum allowable water abstraction rate [m3/s]
 real(dp)                                 :: Vmod           ! abstraction rate from storage [m3/s]
 integer(i4b)                             :: nUps           ! number of upstream segment
 integer(i4b)                             :: iUps           ! upstream reach index
 integer(i4b)                             :: iRch_ups       ! index of upstream reach in NETOPO
 integer(i4b)                             :: ntdh           ! number of time steps in IRF
 integer(i4b)                             :: itdh           ! loop index for unit hydrograph
 character(len=strLen)                    :: fmt1           ! format string
 character(len=strLen)                    :: cmessage       ! error message from subroutine

 ierr=0; message='irf_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   verbose = .true.
 end if

  ! initialize future discharge array at first time
  if (.not.allocated(RCHFLX_out(iens,segIndex)%QFUTURE_IRF))then
    ntdh = size(NETOPO_in(segIndex)%UH)
    allocate(RCHFLX_out(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX_out(iens,segIndex)%QFUTURE_IRF'; return; endif
    RCHFLX_out(iens,segIndex)%QFUTURE_IRF(:) = 0._dp
  end if

  ! get discharge coming from upstream
  nUps = size(NETOPO_in(segIndex)%UREACHI)
  q_upstream = 0.0_dp
  if (nUps>0) then
    do iUps = 1,nUps
      iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxIRF)%REACH_Q
    end do
  endif

  ! perform UH convolution
  call conv_upsbas_qr(NETOPO_in(segIndex)%UH,    &    ! input: reach unit hydrograph
                      q_upstream,                &    ! input: total discharge at top of the reach being processed
                      RCHFLX_out(iens,segIndex), &    ! inout: updated fluxes at reach
                      ierr, message)                  ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

!  ! take out the water from the reach if the wm flag is true and the value are not missing
!  ! here we should make sure the real missing is not injection (or negative abstration)
!  abstract_actual = 0._dp ! this can be removed later (after extensive testing)
!  init_STRQ = 0._dp ! this can be removed later (after extensive testing)
!  if((RCHFLX_out(iens,segIndex)%REACH_WM_FLUX /= realMissing).and.(is_flux_wm)) then
!    abstract_actual = RCHFLX_out(iens,segIndex)%REACH_Q ! get the reach streamflow as actual abstration
!    init_STRQ = RCHFLX_out(iens,segIndex)%REACH_Q ! TO BE DELETED
!    ! reach streamflow is updated based on abstration (positive) or injection (negative)
!    RCHFLX_out(iens,segIndex)%REACH_Q = RCHFLX_out(iens,segIndex)%REACH_Q - RCHFLX_out(iens,segIndex)%REACH_WM_FLUX
!    if (RCHFLX_out(iens,segIndex)%REACH_Q>0) then ! abstration was negative or smaller than reach streamflow
!      abstract_actual  =  RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! actual abstration will be equal to abstration value
!    else
!      RCHFLX_out(iens,segIndex)%REACH_Q = 0._dp ! all the water is taken and actual abstration is reach streamflow
!    endif
!  endif
!  WB_check = RCHFLX_out(iens,segIndex)%REACH_Q + abstract_actual - init_STRQ

  ! take out the water from the reach if the wm flag is true and the value are not missing
  ! here we should make sure the real missing is not injection (or negative abstration)
  if((RCHFLX_out(iens,segIndex)%REACH_WM_FLUX /= realMissing).and.(is_flux_wm)) then
    if (RCHFLX_out(iens,segIndex)%REACH_WM_FLUX <= 0) then ! negative/injection
      RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q = RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q - RCHFLX_out(iens,segIndex)%REACH_WM_FLUX
      RCHFLX_out(iens,segIndex)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX
    else ! positive/abstraction
      Qabs = min(RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1)/dt+RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q, &
                 RCHFLX_out(iens,segIndex)%REACH_WM_FLUX)
      Vmod = min(RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q-Qabs, 0._dp) ! is taken from storage and is <= 0

      RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) + Vmod*dt
      RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q      = RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q - (Qabs+Vmod)
      RCHFLX_out(iens,segIndex)%REACH_WM_FLUX_actual = Qabs

      !if (RCHFLX_out(iens,segIndex)%REACH_WM_FLUX <= RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q) then ! abstraction is smaller than water in the river
      !  RCHFLX_out(iens,segIndex)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX
      !else ! abstraction is larger than water in the river
      !  RCHFLX_out(iens,segIndex)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q
      !endif
    endif
    !RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q = RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q - RCHFLX_out(iens,segIndex)%REACH_WM_FLUX_actual
  endif

  if(verbose)then
    ntdh = size(NETOPO_in(segIndex)%UH)
    write(fmt1,'(A,I5,A)') '(A, 1X',ntdh,'(1X,F20.7))'
    write(*,'(2a)') new_line('a'),'** Check Impulse Response Function routing **'
    write(*,'(a,x,I10,x,I10)') ' Reach index & ID       =', segIndex, NETOPO_in(segIndex)%REACHID
    write(*,fmt1)              ' Unit-Hydrograph        =', (NETOPO_in(segIndex)%UH(itdh), itdh=1,ntdh)
    write(*,'(a)')             ' * total discharge from upstream(q_upstream) [m3/s], local area discharge [m3/s], and Final discharge [m3/s]:'
    write(*,'(a,x,F15.7)')     ' q_upstream             =', q_upstream
    write(*,'(a,x,F15.7)')     ' RCHFLX_out%BASIN_QR(1) =', RCHFLX_out(iens,segIndex)%BASIN_QR(1)
    write(*,'(a,x,F15.7)')     ' RCHFLX_out%REACH_Q =', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q
  endif

  if (RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) < 0) then
    write(iulog,'(A,X,G12.5,X,A,X,I9)') ' ---- NEGATIVE VOLUME [m3]= ', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1), 'at ', NETOPO_in(segIndex)%REACHID
!    RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) = 0._dp
  end if
  if (RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q < 0) then
    write(iulog,'(A,X,G12.5,X,A,X,I9)') ' ---- NEGATIVE FLOW [m3/s] = ', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q, 'at ', NETOPO_in(segIndex)%REACHID
!    RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q = 0._dp
  end if

  call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxIRF, q_upstream, RCHFLX_out(iens,segIndex), verbose, lakeFlag=.false.)

 END SUBROUTINE irf_rch


 ! *********************************************************************
 ! subroutine: Compute delayed runoff from the upstream segments
 ! *********************************************************************
 SUBROUTINE conv_upsbas_qr(reach_uh,   &    ! input: reach unit hydrograph
                           q_upstream, &    ! input:
                           rflux,      &    ! input: input flux at reach
                           ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Details: Convolute runoff volume of upstream at one reach at one time step
 ! ----------------------------------------------------------------------------------------

 implicit none
 ! Argument variables
 real(dp),     intent(in)               :: reach_uh(:)  ! reach unit hydrograph
 real(dp),     intent(in)               :: q_upstream   ! total discharge at top of the reach being processed
 type(STRFLX), intent(inout)            :: rflux        ! current Reach fluxes
 integer(i4b), intent(out)              :: ierr         ! error code
 character(*), intent(out)              :: message      ! error message
 ! Local variables
 integer(i4b)                           :: nTDH         ! number of UH data
 integer(i4b)                           :: iTDH         ! index of UH data (i.e.,future time step)

 ierr=0; message='conv_upsbas_qr/'


 ! place a fraction of runoff in future time steps
 nTDH = size(reach_uh) ! identify the number of future time steps of UH for a given segment
 do iTDH=1,nTDH
   rflux%QFUTURE_IRF(iTDH) = rflux%QFUTURE_IRF(iTDH) &
                             + reach_uh(iTDH)*q_upstream
 enddo

 ! compute volume in reach
 rflux%ROUTE(idxIRF)%REACH_VOL(0) = rflux%ROUTE(idxIRF)%REACH_VOL(1)
! ! For very low flow condition, outflow - inflow > current storage, so limit outflow and adjust rflux%QFUTURE_IRF(1)
! rflux%QFUTURE_IRF(1) = min(rflux%ROUTE(idxIRF)%REACH_VOL(0)/dt + q_upstream*0.999, rflux%QFUTURE_IRF(1))
 rflux%ROUTE(idxIRF)%REACH_VOL(1) = rflux%ROUTE(idxIRF)%REACH_VOL(0) - (rflux%QFUTURE_IRF(1) - q_upstream)*dt

 ! Add local routed flow at the bottom of reach
 rflux%ROUTE(idxIRF)%REACH_Q = rflux%QFUTURE_IRF(1) + rflux%BASIN_QR(1)

 ! move array back   use eoshift
 rflux%QFUTURE_IRF=eoshift(rflux%QFUTURE_IRF,shift=1)


 rflux%QFUTURE_IRF(nTDH) = 0._dp

 END SUBROUTINE conv_upsbas_qr

END MODULE irf_route_module
