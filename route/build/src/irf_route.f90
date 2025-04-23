MODULE irf_route_module

USE nrtype
USE dataTypes,         ONLY: STRFLX            ! fluxes in each reach
USE dataTypes,         ONLY: STRSTA            ! state in each reach
USE dataTypes,         ONLY: RCHTOPO           ! Network topology
USE dataTypes,         ONLY: RCHPRP            ! reach/lake property parameter
USE dataTypes,         ONLY: irfRCH            ! irf specific state data structure
USE public_var,        ONLY: iulog             ! i/o logical unit number
USE public_var,        ONLY: realMissing       ! missing value for real number
USE public_var,        ONLY: integerMissing    ! missing value for integer number
USE public_var,        ONLY: dt                ! simulation time step [sec]
USE public_var,        ONLY: qmodOption      ! qmod option (use 1==direct insertion)
USE public_var,        ONLY: hw_drain_point    ! headwater catchment pour point (top_reach==1 or bottom_reach==2)
USE public_var,        ONLY: min_length_route  ! minimum reach length for routing to be performed.
USE globalData,        ONLY: idxIRF            ! routing method index for IRF method
USE water_balance,     ONLY: comp_reach_wb     ! compute water balance error
USE base_route,        ONLY: base_route_rch    ! base (abstract) reach routing method class
USE data_assimilation, ONLY: direct_insertion  ! qmod option (use 1==direct insertion)

implicit none

private
public::irf_route_rch

integer(i4b), parameter :: top_reach=1
integer(i4b), parameter :: bottom_reach=2

type, extends(base_route_rch) :: irf_route_rch
 CONTAINS
   procedure, pass :: route => irf_rch
end type irf_route_rch

CONTAINS

 ! *********************************************************************
 ! subroutine: perform one segment route UH routing
 ! *********************************************************************
 SUBROUTINE irf_rch(this,         & ! irf_route_rch object to bound this procedure
                    iEns,         & ! input: index of runoff ensemble
                    segIndex,     & ! input: reach index
                    ixDesire,     & ! input: reachID to be checked by on-screen pringing
                    T0,T1,        & ! input: start and end of the time step
                    NETOPO_in,    & ! input: reach topology data structure
                    RPARAM_in,    & ! input: reach parameter data structure
                    RCHSTA_out,   & ! inout: reach state data structure
                    RCHFLX_out,   & ! inout: reach flux data structure
                    ierr, message)  ! output: error control

 USE public_var, ONLY:is_flux_wm   ! logical water management components fluxes should be read

 implicit none
 ! Argument variables
 class(irf_route_rch)                     :: this            ! irf_route_rch object to bound this procedure
 integer(i4b),  intent(in)                :: iEns            ! runoff ensemble to be routed
 integer(i4b),  intent(in)                :: segIndex        ! segment where routing is performed
 integer(i4b),  intent(in)                :: ixDesire        ! index of the reach for verbose output
 real(dp),      intent(in)                :: T0,T1           ! start and end of the time step (seconds)
 type(RCHTOPO), intent(in),   allocatable :: NETOPO_in(:)    ! River Network topology
 type(RCHPRP),  intent(inout),allocatable :: RPARAM_in(:)    ! River reach parameter
 type(STRSTA),  intent(inout)             :: RCHSTA_out(:,:) ! reach state data
 type(STRFLX),  intent(inout)             :: RCHFLX_out(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b),  intent(out)               :: ierr            ! error code
 character(*),  intent(out)               :: message         ! error message
 ! Local variables
 logical(lgt)                             :: verbose        ! check details of variables
 real(dp)                                 :: q_upstream     ! total discharge at top of the reach [m3/s]
 real(dp)                                 :: q_upstream_mod ! total discharge at top of the reach after water abstraction [m3/s]
 real(dp)                                 :: Qlat           ! lateral flow into channel [m3/s]
 real(dp)                                 :: Qabs           ! residual water abstraction [m3/s]
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

  ! get discharge coming from upstream
  nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
  q_upstream = 0.0_dp

  Qabs = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initial water abstraction (positive) or injection (negative)
  RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initialize actual water abstraction

  ! update volume at previous time step
  RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(0) = RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1)

  if (nUps>0) then ! this hru is not headwater
    do iUps = 1,nUps
      if (.not. NETOPO_in(segIndex)%goodBas(iUps)) cycle ! skip upstream reach which does not any flow due to zero total contributory areas
      iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
      q_upstream = q_upstream + RCHFLX_out(iens, iRch_ups)%ROUTE(idxIRF)%REACH_Q
    end do
    q_upstream_mod  = q_upstream
    Qlat = RCHFLX_out(iens,segIndex)%BASIN_QR(1)
  else ! head water
    if (verbose) then
      write(iulog,'(A)')            ' This is headwater '
    endif
    if (hw_drain_point==top_reach) then ! lateral flow is poured in a reach at the top
      q_upstream = q_upstream + RCHFLX_out(iens,segIndex)%BASIN_QR(1)
      q_upstream_mod = q_upstream
      Qlat = 0._dp
    else if (hw_drain_point==bottom_reach) then ! lateral flow is poured in a reach at the top
      q_upstream_mod = q_upstream
      Qlat = RCHFLX_out(iens,segIndex)%BASIN_QR(1)
    end if
  endif

  RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_INFLOW = q_upstream ! total inflow from the upstream reaches

  ! Water management - water injection or abstraction (irrigation or industrial/domestic water usage)
  ! For water abstraction, water is extracted from the following priorities:
  ! 1. existing storage(REACH_VOL(0), 2. upstream inflow , 3 lateral flow (BASIN_QR)
  if((RCHFLX_out(iens,segIndex)%REACH_WM_FLUX /= realMissing).and.(is_flux_wm)) then
    if (Qabs > 0) then ! positive == abstraction
      if (RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1)/dt > Qabs) then ! take out all abstraction from strorage
        RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) - Qabs*dt
      else ! if inital abstraction is greater than volume
        Qabs = Qabs - RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1)/dt ! get residual Qabs after extracting from strorage
        RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) = 0._dp ! volume gets 0
        if (q_upstream > Qabs) then ! then take out all residual abstraction from upstream inflow
          q_upstream_mod = q_upstream - Qabs
        else ! if residual abstraction is still greater than lateral flow
          Qabs = Qabs - q_upstream ! get residual abstraction after extracting upstream inflow and storage.
          q_upstream_mod = 0._dp ! upstream inflow gets 0 (all is gone to abstracted flow).
          if (Qlat > Qabs) then ! then take residual abstraction out from lateral flow
            Qlat = Qlat - Qabs
          else ! if residual abstraction is greater than upstream inflow
            Qabs = Qabs - Qlat ! take out residual abstraction from lateral flow
            Qlat = 0._dp ! lateral flow gets 0 (all are gone to abstraction)
            RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX - Qabs
          end if
        end if
      end if
    else ! negative == injection
      Qlat = Qlat - Qabs
    endif
  endif

  ! perform UH convolution
  call conv_upsbas_qr(RPARAM_in(segIndex),       &    ! input: parameter at segIndex reach
                      NETOPO_in(segIndex)%UH,    &    ! input: reach unit hydrograph
                      q_upstream_mod,            &    ! input: total discharge from the upstreams
                      Qlat,                      &    ! input: lateral flow [m3/s]
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

  if(verbose)then
    ntdh = size(NETOPO_in(segIndex)%UH)
    write(fmt1,'(A,I5,A)') '(A, 1X',ntdh,'(1X,F20.7))'
    write(*,'(2a)') new_line('a'),'** Check Impulse Response Function routing **'
    write(*,'(a,1x,I10,1x,I10)')' Reach index & ID       =', segIndex, NETOPO_in(segIndex)%REACHID
    write(*,fmt1)               ' Unit-Hydrograph        =', (NETOPO_in(segIndex)%UH(itdh), itdh=1,ntdh)
    write(*,'(a)')              ' * total discharge from upstream(q_upstream) [m3/s], local area discharge [m3/s]' // &
            ', and Final discharge [m3/s]:'
    write(*,'(a,1x,F15.7)')     ' q_upstream             =', q_upstream
    write(*,'(a,1x,F15.7)')     ' RCHFLX_out%BASIN_QR(1) =', RCHFLX_out(iens,segIndex)%BASIN_QR(1)
    write(*,'(a,1x,F15.7)')     ' RCHFLX_out%REACH_Q =', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q
  endif

  if (RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) < 0) then
    write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE VOLUME [m3]= ', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1), &
          'at ', NETOPO_in(segIndex)%REACHID
!    RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_VOL(1) = 0._dp
  end if
  if (RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q < 0) then
    write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE FLOW [m3/s] = ', RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q, &
           'at ', NETOPO_in(segIndex)%REACHID
!    RCHFLX_out(iens,segIndex)%ROUTE(idxIRF)%REACH_Q = 0._dp
  end if

  if (qmodOption==1) then
    call direct_insertion(iens, segIndex, & ! input: reach index
                          idxIRF,         & ! input: routing method id for Euler kinematic wave
                          ixDesire,       & ! input: verbose seg index
                          NETOPO_in,      & ! input: reach topology data structure
                          RCHSTA_out,     & ! inout: reach state data structure
                          RCHFLX_out,     & ! inout: reach fluxes datq structure
                          ierr, cmessage)   ! output: error control
    if(ierr/=0)then
      write(message,'(A,X,I12,X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage); return
    endif
  end if

  if (qmodOption==0) then ! check reach water balance only if data assimilation is off
    call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxIRF, q_upstream, Qlat, RCHFLX_out(iens,segIndex), verbose, lakeFlag=.false.)
  end if

 END SUBROUTINE irf_rch


 ! *********************************************************************
 ! subroutine: Compute delayed runoff from the upstream segments
 ! *********************************************************************
 SUBROUTINE conv_upsbas_qr(rch_param,  &    ! input: river parameter data structure
                           reach_uh,   &    ! input: reach unit hydrograph
                           q_upstream, &    ! input:
                           Qlat,       &    ! input:
                           rflux,      &    ! input: input flux at reach
                           ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Details: Convolute runoff volume of upstream at one reach at one time step
 ! ----------------------------------------------------------------------------------------

 implicit none
 ! Argument variables
 type(RCHPRP), intent(in)               :: rch_param    ! River reach parameter
 real(dp),     intent(in)               :: reach_uh(:)  ! reach unit hydrograph
 real(dp),     intent(in)               :: q_upstream   ! total discharge at top of the reach being processed
 real(dp),     intent(in)               :: Qlat         ! lataral flow
 type(STRFLX), intent(inout)            :: rflux        ! current Reach fluxes
 integer(i4b), intent(out)              :: ierr         ! error code
 character(*), intent(out)              :: message      ! error message
 ! Local variables
 integer(i4b)                           :: nTDH         ! number of UH data
 integer(i4b)                           :: iTDH         ! index of UH data (i.e.,future time step)

 ierr=0; message='conv_upsbas_qr/'

 if (rch_param%RLENGTH > min_length_route) then
 ! place a fraction of runoff in future time steps
 nTDH = size(reach_uh) ! identify the number of future time steps of UH for a given segment
 do iTDH=1,nTDH
   rflux%QFUTURE_IRF(iTDH) = rflux%QFUTURE_IRF(iTDH) &
                             + reach_uh(iTDH)*q_upstream
 enddo

 ! compute volume in reach
 ! For very low flow condition, outflow - inflow > current storage, so limit outflow and adjust rflux%QFUTURE_IRF(1)
! rflux%QFUTURE_IRF(1) = min(0.999*(rflux%ROUTE(idxIRF)%REACH_VOL(1)/dt + q_upstream), rflux%QFUTURE_IRF(1))
 rflux%QFUTURE_IRF(1) = min(rflux%ROUTE(idxIRF)%REACH_VOL(0)/dt + q_upstream*0.999, rflux%QFUTURE_IRF(1))
 rflux%ROUTE(idxIRF)%REACH_VOL(1) = rflux%ROUTE(idxIRF)%REACH_VOL(1) - (rflux%QFUTURE_IRF(1) - q_upstream)*dt

 ! Add local routed flow at the bottom of reach
 rflux%ROUTE(idxIRF)%REACH_Q = rflux%QFUTURE_IRF(1) + Qlat

 ! move array back   use eoshift
 rflux%QFUTURE_IRF=eoshift(rflux%QFUTURE_IRF,shift=1)

 rflux%QFUTURE_IRF(nTDH) = 0._dp
 else ! length < min_length_route: length is short enough to just pass upstream to downstream
   rflux%QFUTURE_IRF(:) = 0._dp
   rflux%QFUTURE_IRF(1) = q_upstream
   rflux%ROUTE(idxIRF)%REACH_Q = rflux%QFUTURE_IRF(1) + Qlat

   rflux%ROUTE(idxIRF)%REACH_VOL(0) = 0._dp
   rflux%ROUTE(idxIRF)%REACH_VOL(1) = 0._dp
 end if

 END SUBROUTINE conv_upsbas_qr

END MODULE irf_route_module
