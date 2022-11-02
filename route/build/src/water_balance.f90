MODULE water_balance

USE nrtype
! data type
USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
USE dataTypes, ONLY: RCHTOPO        ! Network topology
! global parameters
USE public_var, ONLY: iulog         ! i/o logical unit number
USE public_var, ONLY: dt            ! simulation time step

implicit none

private
public:: comp_reach_wb
public:: comp_global_wb

CONTAINS

  ! *********************************************************************
  ! public subroutine: compute reach water balance
  ! *********************************************************************
  SUBROUTINE comp_reach_wb(ixRoute,    &     ! input: index of routing method
                           Qupstream,  &     ! input: inflow from upstream
                           RCHFLX_in,  &     ! inout: reach flux data structure
                           verbose)

  ! Descriptions
  ! Compute water balance per simulation time step and reach/lake
  ! During a time step dt = t1 - t0 [sec]
  ! dVol     [m3]   = Vol(t1) - Vol(t0)
  ! flux_in  [m3/s] = inflow(dt) + lateral(dt) + Prec(dt)
  ! flux_out [m3/s] = outflow(dt) + abstract(dt) + Evap(dt)
  !
  ! Compare dVol vs (flux_in - flux_out)*dt

  implicit none
  ! Argument variables:
  integer(i4b), intent(in)                 :: ixRoute        ! input: routing method index
  real(dp),     intent(in)                 :: Qupstream      ! input: total inflow from upstream reaches
  type(STRFLX), intent(inout)              :: RCHFLX_in      ! inout: Reach fluxes data structure
  logical(lgt), intent(in)                 :: verbose        ! input: reach index to be examined
  ! Local variables:
  real(dp)                                 :: dVol           !
  real(dp)                                 :: Qin            !
  real(dp)                                 :: Qlateral       !
  real(dp)                                 :: Qout           !
  real(dp)                                 :: precip         !
  real(dp)                                 :: evapo          !
  real(dp)                                 :: Qtake          !

  ! volume change
  dVol     = RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(1)-RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(0)
  ! in flux
  Qin      = Qupstream *dt
  Qlateral = RCHFLX_in%BASIN_QR(1) *dt
  precip   = RCHFLX_in%basinprecip *dt
  ! out flux
  Qout     = -1._dp *RCHFLX_in%ROUTE(ixRoute)%REACH_Q *dt
  Qtake    = -1._dp *RCHFLX_in%REACH_WM_FLUX *dt
  evapo    = -1._dp *RCHFLX_in%basinevapo *dt

  RCHFLX_in%ROUTE(ixRoute)%WB = dVol - (Qin + Qlateral + precip + Qout + Qtake + evapo)

  if (verbose) then
    write(iulog,'(A,1PG15.7)') '  WBerr [m3]        = ', RCHFLX_in%ROUTE(ixRoute)%WB
    write(iulog,'(A,1PG15.7)') '  Vol at t0 [m3]    = ', RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(0)
    write(iulog,'(A,1PG15.7)') '  Vol at t1 [m3]    = ', RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(1)
    write(iulog,'(A,1PG15.7)') '  dVol [m3]         = ', dVol
    write(iulog,'(A,1PG15.7)') '  inflow [m3]       = ', Qin
    write(iulog,'(A,1PG15.7)') '  lateral flow [m3] = ', Qlateral
    write(iulog,'(A,1PG15.7)') '  precip [m3]       = ', precip
    write(iulog,'(A,1PG15.7)') '  outflow [m3]      = ', Qout
    write(iulog,'(A,1PG15.7)') '  abstraction [m3]  = ', Qtake
    write(iulog,'(A,1PG15.7)') '  evaporation [m3]  = ', evapo
  endif
  if (abs(RCHFLX_in%ROUTE(ixRoute)%WB) > 1.e-5_dp) then
    write(iulog,'(A,1PG15.7,1X,A)') ' WARNING: WB error [m3] = ', RCHFLX_in%ROUTE(ixRoute)%WB, '> 1.e-5 [m3]'
  end if

  END SUBROUTINE comp_reach_wb

  ! *********************************************************************
  ! public subroutine: compute global water balance
  ! *********************************************************************
  SUBROUTINE comp_global_wb(ixRoute, verbose, ierr, message)

    USE globalData, ONLY: multiProcs, masterproc
    USE globalData, ONLY: RCHFLX_trib
    USE globalData, ONLY: NETOPO_main
    USE globalData, ONLY: NETOPO_trib
    USE globalData, ONLY: nRch_mainstem        ! scalar data: number of mainstem reaches
    USE globalData, ONLY: nRch_trib            ! scalar data: number of tributary reaches
    USE globalData, ONLY: nTribOutlet
    USE mpi_utils,  ONLY: shr_mpi_reduce

    ! Descriptions
    ! Compute water balance per simulation time step over the whole domain
    ! During a time step dt = t1 - t0 [sec]
    ! dVol     [m3]   = Vol(t1) - Vol(t0)
    ! flux_in  [m3/s] = lateral(dt) + Prec(dt)
    ! flux_out [m3/s] = abstract(dt) + Evap(dt)
    !                   outflow(dt) where outlet
    !
    ! SUM(dVol) - SUM(flux_in - flux_out)*dt
    ! SUM is SUM over the whole domain

    implicit none
    ! Argument variables:
    integer(i4b),      intent(in)     :: ixRoute        ! input: routing method index
    logical(lgt),      intent(in)     :: verbose        ! input: reach index to be examined
    integer(i4b),      intent(out)    :: ierr
    character(strLen), intent(out)    :: message        ! error message
    ! Local variables:
    integer(i4b)                      :: lwr,upr        ! loop index
    real(dp)                          :: wb_local(6)
    real(dp)                          :: wb_global(6)
    real(dp)                          :: wb_mainstem(6)
    real(dp)                          :: wb_trib(6)
    real(dp)                          :: wb_error
    character(strLen)                 :: cmessage       ! error message from subroutine

    ierr=0; message='comp_global_wb/'

    if (multiProcs) then
      if (masterproc) then
        wb_trib     = 0._dp
        wb_mainstem = 0._dp
        if (nRch_mainstem > 0) then
          call accum_water_balance(NETOPO_main(1:nRch_mainstem), RCHFLX_trib(1,1:nRch_mainstem), wb_mainstem, ierr, cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        end if
        if (nRch_trib > 0) then
          lwr = nRch_mainstem + nTribOutlet + 1
          upr = nRch_mainstem + nTribOutlet + nRch_trib
          call accum_water_balance(NETOPO_trib, RCHFLX_trib(1,lwr:upr), wb_trib, ierr, cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        end if
        wb_local = wb_mainstem + wb_trib
      else ! non-main processors
        call accum_water_balance(NETOPO_trib, RCHFLX_trib(1,:), wb_local, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if
    else ! single core use
      call accum_water_balance(NETOPO_main, RCHFLX_trib(1,1:nRch_mainstem), wb_local, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end if

    ! sum water balance components across all the processors
    call shr_mpi_reduce(wb_local, 'sum', wb_global, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    wb_error = wb_global(1)-sum(wb_global(2:6))

    if (verbose) then
      if (masterproc) then
        write(iulog,'(A)') '  global water balance [m3] '
        write(iulog,'(A,1PG15.7)') '  dVol [m3]         = ', wb_global(1)
        write(iulog,'(A,1PG15.7)') '  lateral flow [m3] = ', wb_global(2)
        write(iulog,'(A,1PG15.7)') '  precip [m3]       = ', wb_global(3)
        write(iulog,'(A,1PG15.7)') '  abstraction [m3]  = ', wb_global(4)
        write(iulog,'(A,1PG15.7)') '  evaporation [m3]  = ', wb_global(5)
        write(iulog,'(A,1PG15.7)') '  outflow [m3]      = ', wb_global(6)
        write(iulog,'(A,1PG15.7)') '  WBerr [m3]        = ', wb_error
      end if
    endif

    if (abs(wb_error) > 1._dp) then ! tolerance is 1 [m3]
      write(iulog,'(A,1PG15.7,1X,A)') ' WARNING: global WB error [m3] = ', wb_error, '> 1.0 [m3]'
    end if

    CONTAINS

    SUBROUTINE accum_water_balance(NETOPO_in, RCHFLX_in, water_budget, ierr, message)

      USE dataTypes, ONLY: RCHTOPO   ! data structure - Network topology
      USE dataTypes, ONLY: STRFLX    ! data structure - fluxes in each reach

      implicit none
      ! Arguments:
      type(RCHTOPO),     intent(in)  :: NETOPO_in(:)
      type(STRFLX),      intent(in)  :: RCHFLX_in(:)
      real(dp),          intent(out) :: water_budget(6)
      integer(i4b),      intent(out) :: ierr
      character(strLen), intent(out) :: message           ! error message
      ! Local variables:
      integer(i4b)               :: ix          ! loop index
      integer(i4b)               :: nRch_local  ! number of reach

      ierr=0; message='accum_water_balance/'

      nRch_local = size(RCHFLX_in)
      if (nRch_local/=size(NETOPO_in)) then
        ierr=20; message=trim(message)//'RCHFLX and NETOPO has different sizes'; return
      end if

      water_budget = 0._dp
      do ix = 1, nRch_local
        water_budget(1) = water_budget(1) + (RCHFLX_in(ix)%ROUTE(ixRoute)%REACH_VOL(1)-RCHFLX_in(ix)%ROUTE(ixRoute)%REACH_VOL(0))
        ! in flux
        water_budget(2) = water_budget(2) + RCHFLX_in(ix)%BASIN_QR(1) *dt
        water_budget(3) = water_budget(3) + RCHFLX_in(ix)%basinprecip *dt
        ! out flux
        water_budget(4) = water_budget(4) - RCHFLX_in(ix)%basinevapo *dt
        water_budget(5) = water_budget(5) - RCHFLX_in(ix)%REACH_WM_FLUX *dt

        if (NETOPO_in(ix)%DREACHI==-1 .and. NETOPO_in(ix)%DREACHK<=0) then ! to-do: better way to detect outlet segments
          water_budget(6) = water_budget(6) - RCHFLX_in(ix)%ROUTE(ixRoute)%REACH_Q *dt
        end if
      end do
    END SUBROUTINE accum_water_balance

  END SUBROUTINE comp_global_wb

END MODULE water_balance
