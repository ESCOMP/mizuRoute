MODULE reach_flux
 USE nrtype
 IMPLICIT NONE
 SAVE
 REAL(DP), DIMENSION(:), ALLOCATABLE   :: FRAC_FUTURE       ! fraction of runoff in future time steps
 TYPE STRFLX
  REAL(DP), POINTER                    :: QFUTURE(:)        ! runoff volume in future time steps (m3/s)
  REAL(DP), POINTER                    :: QFUTURE_IRF(:)    ! runoff volume in future time steps for IRF routing (m3/s) added by NM
  REAL(DP)                             :: BASIN_QI          ! instantaneous runoff volume from the local basin (m3/s)
  REAL(DP)                             :: BASIN_QR(0:1)     ! routed runoff volume from the local basin (m3/s)
  REAL(DP)                             :: UPSBASIN_QR       ! routed runoff depth from the upstream basins (m/s) added by NM
  REAL(DP)                             :: BASIN_QR_IRF(0:1) ! routed runoff volume from all the upstream basin (m3/s) added by NM
  REAL(DP)                             :: REACH_Q           ! time-step average streamflow (m3/s)
  REAL(DP)                             :: REACH_Q_IRF       ! time-step average streamflow (m3/s) from IRF routing added by NM
  REAL(DP)                             :: UPSTREAM_QI       ! sum of upstream streamflow (m3/s)
  REAL(DP)                             :: TAKE              ! average take
 ENDTYPE STRFLX
 TYPE(STRFLX), DIMENSION(:,:), POINTER :: RCHFLX            ! Reach fluxes
END MODULE reach_flux

