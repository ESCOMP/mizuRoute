MODULE lake_route_module

USE nrtype
USE dataTypes,     ONLY: STRFLX          ! data struct: fluxes in each reach
USE dataTypes,     ONLY: RCHTOPO         ! data struct: Network topology
USE dataTypes,     ONLY: RCHPRP          ! data struct: Network parameter
USE public_var,    ONLY: iulog           ! parameter: i/o logical unit number
USE public_var,    ONLY: realMissing     ! parameter: missing value for real number
USE public_var,    ONLY: pi              ! parameter: pi value of 3.14159265359_dp
USE globalData,    ONLY: isColdStart     ! parameter: restart flag
USE water_balance, ONLY: comp_reach_wb   ! routine: compute water balance error
USE ascii_utils,   ONLY: lower           ! routine: convert string to lower case

implicit none
integer(i4b),parameter :: endorheic=0
integer(i4b),parameter :: doll03=1
integer(i4b),parameter :: hanasaki06=2
integer(i4b),parameter :: hype=3

private
public::lake_route

CONTAINS

  ! *********************************************************************
  ! subroutine: perform one segment lake routing
  ! *********************************************************************
  SUBROUTINE lake_route(iEns,       &    ! input: index of runoff ensemble to be processed
                        segIndex,   &    ! input: index of runoff ensemble to be processed
                        idxRoute,   &    ! input: routing method index
                        ixDesire,   &    ! input: reachID to be checked by on-screen pringing (here reachID can be lake)
                        NETOPO_in,  &    ! input: reach topology data structure
                        RPARAM_in,  &    ! input: reach parameter data strcuture
                        RCHFLX_out, &    ! inout: reach flux data structure
                        ierr, message)   ! output: error control

  USE globalData, ONLY: iTime               ! current model time step
  USE globalData, ONLY: simDatetime         ! previous and current model time
  USE public_var, ONLY: is_flux_wm          ! logical water management components fluxes should be read
  USE public_var, ONLY: dt                  ! simulation time step [sec]
  USE public_var, ONLY: lakeWBtol           ! lake water balance tolerance
  USE public_var, ONLY: is_vol_wm_jumpstart ! logical whether or not lake should be simulated
  USE public_var, ONLY: secprday            ! seconds per day = 86400
  USE public_var, ONLY: days_per_yr         ! days per a year = 365
  USE public_var, ONLY: months_per_yr       ! months per a year = 12
  USE public_var, ONLY: calendar            ! calendar name

  implicit none
  ! Argument variables:
  integer(i4b), intent(in)                 :: iEns           ! runoff ensemble to be routed
  integer(i4b), intent(in)                 :: segIndex       ! segment where routing is performed
  integer(i4b), intent(in)                 :: idxRoute       ! index of the reach for verbose output
  integer(i4b), intent(in)                 :: ixDesire       ! index of the reach for verbose output
  type(RCHTOPO), intent(in),   allocatable :: NETOPO_in(:)   ! River Network topology
  type(RCHPRP), intent(inout), allocatable :: RPARAM_in(:)   ! River Network topology
  TYPE(STRFLX), intent(inout)              :: RCHFLX_out(:,:)! Reach fluxes (ensembles, space [reaches]) for decomposed domains
  integer(i4b), intent(out)                :: ierr           ! error code
  character(*), intent(out)                :: message        ! error message
  ! Local variables:
  logical(lgt)                             :: verbose        ! check details of variables
  real(dp)                                 :: q_upstream     ! total discharge at top of the reach being processed
!  real(dp)                                 :: WB             ! water balance component in the lake
  integer(i4b)                             :: nUps           ! number of upstream segment
  integer(i4b)                             :: iUps           ! upstream reach index
  integer(i4b)                             :: iRch_ups       ! index of upstream reach in NETOPO
  integer(i4b)                             :: ntdh           ! number of time steps in IRF
  character(len=strLen)                    :: cmessage       ! error message from subroutine
  ! local variables for H06 routine
  real(dp)                                 :: c                   ! storage to yearly activity ratio
  real(dp)                                 :: I_yearly, D_yearly  ! mean annual inflow and demand
  real(dp), dimension(12)                  :: I_months, D_months  ! mean monthly inflow and demand
  integer(i4b), dimension(2)               :: array_size(2)       ! get the size of array_size
  integer(i4b)                             :: start_month=0       ! start month of the operational year
  integer(i4b)                             :: i                   ! index
  integer(i4b)                             :: past_length_I       ! pas length for inflow based on length in year and floor
  integer(i4b)                             :: past_length_D       ! pas length for demand based on length in year and floor
  real(dp)                                 :: target_r            ! target release
  ! local varibale for HYPE routine
  real(dp)                                 :: Day_of_year         ! the day number in a year
  integer(i4b)                             :: F_prim              ! factor for local flag is reservoir has primary spillway
  real(dp)                                 :: F_sin               ! factor for sin
  real(dp)                                 :: F_lin               ! factor for linear
  real(dp)                                 :: Q_prim              ! simulated outflow from main or primary spillway
  real(dp)                                 :: Q_spill             ! simulated outflow from emergency spillway
  real(dp)                                 :: Q_sim               ! simulated output from the reservoir

  ierr=0; message='lake_route/'

  verbose = .false.
  if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
    verbose = .true.
  end if

    ! identify number of upstream segments of the lake being processed
    nUps = size(NETOPO_in(segIndex)%UREACHI)

    q_upstream = 0.0_dp
    if(nUps>0)then
      do iUps = 1,nUps
        iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
        q_upstream = q_upstream + RCHFLX_out(iens,iRch_ups)%ROUTE(idxRoute)%REACH_Q
      end do
    endif

    ! no runoff input is added to the lake; the input are only precipitation and evaporation to the lake
    if (verbose) then
      write(iulog, '(A,1X,I15)') 'lake id=', NETOPO_in(segIndex)%REACHID ! to check the reach id of lake
      select case(NETOPO_in(segIndex)%LakeModelType)
        case(endorheic)
          write(iulog, '(A)') 'Lake type: Endorheic'
        case(doll03)
          write(iulog, '(A)') 'Lake type: Natural lake with Doll scheme'
          write(iulog, '(A,1X,G12.5)') 'D03_MaxStorage  = ', RPARAM_in(segIndex)%D03_MaxStorage
          write(iulog, '(A,1X,G12.5)') 'D03_Coefficient = ', RPARAM_in(segIndex)%D03_Coefficient
          write(iulog, '(A,1X,G12.5)') 'D03_Power       = ', RPARAM_in(segIndex)%D03_Power
        case(hanasaki06)
          write(iulog, '(A)') 'Lake type: Managed lake with Hanasaki scheme'
          write(iulog, '(A,1X,G12.5)') 'H06_Smax     =', RPARAM_in(segIndex)%H06_Smax
          write(iulog, '(A,1X,G12.5)') 'H06_alpha    =', RPARAM_in(segIndex)%H06_alpha
          write(iulog, '(A,1X,G12.5)') 'H06_envfact  =', RPARAM_in(segIndex)%H06_envfact
          write(iulog, '(A,1X,G12.5)') 'H06_S_ini    =', RPARAM_in(segIndex)%H06_S_ini
          write(iulog, '(A,1X,G12.5)') 'H06_c1       =', RPARAM_in(segIndex)%H06_c1
          write(iulog, '(A,1X,G12.5)') 'H06_c2       =', RPARAM_in(segIndex)%H06_c2
          write(iulog, '(A,1X,G12.5)') 'H06_exponent =', RPARAM_in(segIndex)%H06_exponent
        case(hype)
          write(iulog, '(A)') 'Lake type: Managed lake with Hype scheme'
          write(iulog, '(A,1X,G12.5)') 'HYP_E_emr .......= ', RPARAM_in(segIndex)%HYP_E_emr
          write(iulog, '(A,1X,G12.5)') 'HYP_E_lim .......= ', RPARAM_in(segIndex)%HYP_E_lim
          write(iulog, '(A,1X,G12.5)') 'HYP_E_min .......= ', RPARAM_in(segIndex)%HYP_E_min
          write(iulog, '(A,1X,G12.5)') 'HYP_E_zero ......= ', RPARAM_in(segIndex)%HYP_E_zero
          write(iulog, '(A,1X,G12.5)') 'HYP_Qrate_emr ...= ', RPARAM_in(segIndex)%HYP_Qrate_emr
          write(iulog, '(A,1X,G12.5)') 'HYP_Erate_emr ...= ', RPARAM_in(segIndex)%HYP_Erate_emr
          write(iulog, '(A,1X,G12.5)') 'HYP_Qrate_prim ..= ', RPARAM_in(segIndex)%HYP_Qrate_prim
          write(iulog, '(A,1X,G12.5)') 'HYP_Qrate_amp ...= ', RPARAM_in(segIndex)%HYP_Qrate_amp
          write(iulog, '(A,1X,G12.5)') 'HYP_Qrate_phs ...= ', RPARAM_in(segIndex)%HYP_Qrate_phs
          write(iulog, '(A,1X,G12.5)') 'HYP_prim_F ......= ', RPARAM_in(segIndex)%HYP_prim_F
          write(iulog, '(A,1X,G12.5)') 'HYP_A_avg .......= ', RPARAM_in(segIndex)%HYP_A_avg
        case default; ierr=20; message=trim(message)//'unable to identify the parametric lake model type'; return
      end select
      write(iulog,'(A,1X,G12.5)') 'total inflow [m3/s]  = ', q_upstream
      write(iulog,'(A,1X,G12.5)') 'precipitation [m3/s] = ', RCHFLX_out(iens,segIndex)%basinprecip
      write(iulog,'(A,1X,G12.5)') 'evaporation [m3/s]   = ', RCHFLX_out(iens,segIndex)%basinevapo
    end if

    ! jump start the lake volume to the target volume if provided for the first time step
    if (iTime==1) then
      if ((is_vol_wm_jumpstart).and.(NETOPO_in(segIndex)%LakeTargVol)) then
        RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_WM_VOL ! update the initial condition with first target volume value
      else ! the lake volume is not jump started based on lake target volume
        if (isColdStart) then    ! Cold start ....... initialize flux structures
          select case(NETOPO_in(segIndex)%LakeModelType)
            case(endorheic)
              RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RPARAM_in(segIndex)%D03_S0 ! currently assumes all endorheic max storage are provided under doll formulation
            case(doll03)
              RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RPARAM_in(segIndex)%D03_MaxStorage
            case(hanasaki06)
              RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RPARAM_in(segIndex)%H06_Smax
            case(hype)
              RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = (RPARAM_in(segIndex)%HYP_E_emr - RPARAM_in(segIndex)%HYP_E_zero) * RPARAM_in(segIndex)%HYP_A_avg
            case default; ierr=20; message=trim(message)//'unable to identify the parametric lake model type'; return
          end select
        endif
      endif
    endif

    ! add upstream, precipitation and subtract evaporation from the lake volume
    RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(0) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) ! updating storage at previous time step
    RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) + q_upstream * dt  ! input upstream discharge from m3/s to m3
    RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) + RCHFLX_out(iens,segIndex)%BASIN_QR(1) * dt ! add lateral flow
    RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) + RCHFLX_out(iens,segIndex)%basinprecip * dt ! input lake precipitation
    if (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) > RCHFLX_out(iens,segIndex)%basinevapo*dt) then ! enough water to evaporate
      RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%basinevapo * dt ! output lake evaporation
    else ! not enough water to evaporate (this cause water balance error in coupling mode)
      RCHFLX_out(iens,segIndex)%basinevapo = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1)/dt ! update basinevapo
      RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1)=0._dp
    endif

    RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! initialize actual water abstraction

    ! take out the water from the reach if the wm flag is true and the value are not missing
    ! here we should make sure the real missing is not injection (or negative abstration)
    if((RCHFLX_out(iens,segIndex)%REACH_WM_FLUX /= realMissing).and.(is_flux_wm)) then
      ! abstraction or injection
      if (RCHFLX_out(iens,segIndex)%REACH_WM_FLUX <= 0) then ! negative/injection
        RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_WM_FLUX*dt
        RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX
      else ! positive/abstraction
        if (RCHFLX_out(iens,segIndex)%REACH_WM_FLUX*dt <= RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1)) then ! abstraction is smaller than water in the lake
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_WM_FLUX*dt
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX
        else ! abstraction is larger than water in the river
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_WM_FLUX_actual = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1)/dt
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = 0._dp
        endif
      endif
    endif

    ! simulating lake/reservoir output
    if (NETOPO_in(segIndex)%LakeTargVol) then ! The lake should follow the given target volume

      if (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) < RCHFLX_out(iens,segIndex)%REACH_WM_VOL) then ! The storage is smaller, accumulate
        RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = 0 ! lake/reservoir output is set to zero (0)
      else ! The storage is largrt, release to get to the target level
        RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q      = (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_WM_VOL) /dt
        RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_WM_VOL
      endif

    else ! if the lake is parameteric
      select case(NETOPO_in(segIndex)%LakeModelType)
        case(endorheic)
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = 0._dp ! no outflow
        case(doll03)
          ! Check if the MaxStorage is larger than S0
          if (RPARAM_in(segIndex)%D03_MaxStorage < RPARAM_in(segIndex)%D03_S0) then !
            cmessage = 'Additional parameter of inactive storage is larger than MaxStorage of Doll formulation, please check and correct'
            ierr = 1; message=trim(message)//trim(cmessage);
          endif
          ! The D03_Coefficient is based on d**-1 meaning the result will be m**3 d**-1 and should be converter to m**3 s**-1
          if ((RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RPARAM_in(segIndex)%D03_S0) > 0) then
            RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = RPARAM_in(segIndex)%D03_Coefficient * &
            (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RPARAM_in(segIndex)%D03_S0) * &
            ((RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RPARAM_in(segIndex)%D03_S0) / &
            (RPARAM_in(segIndex)%D03_MaxStorage - RPARAM_in(segIndex)%D03_S0)) ** &
            RPARAM_in(segIndex)%D03_Power ! Q = AS(S/Smax)^B based on Eq. 1 Hanasaki et al., 2006 https://doi.org/10.1016/j.jhydrol.2005.11.011
          else
            RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = 0
          endif
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q / secprday ! conversion to m**3 s**-1
          ! in case is the output volume is more than lake volume
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = min(RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q, RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1)/dt)
          ! updating the storage
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q * dt

        case(hanasaki06)
          ! preserving the past upstrem discharge for lake models
          ! create memory of upstream inflow for Hanasaki formulation if memory flag is active for this lake
          if (RPARAM_in(segIndex)%H06_I_mem_F) then ! if memeory is acive for this case then allocate the past input
            if (.not.allocated(RCHFLX_out(iens,segIndex)%QPASTUP_IRF)) then
              past_length_I = floor(RPARAM_in(segIndex)%H06_I_mem_L * 31 * secprday / dt)
              allocate(RCHFLX_out(iens,segIndex)%QPASTUP_IRF(12,past_length_I),stat=ierr)
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 1,:) = RPARAM_in(segIndex)%H06_I_Jan
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 2,:) = RPARAM_in(segIndex)%H06_I_Feb
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 3,:) = RPARAM_in(segIndex)%H06_I_Mar
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 4,:) = RPARAM_in(segIndex)%H06_I_Apr
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 5,:) = RPARAM_in(segIndex)%H06_I_May
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 6,:) = RPARAM_in(segIndex)%H06_I_Jun
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 7,:) = RPARAM_in(segIndex)%H06_I_Jul
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 8,:) = RPARAM_in(segIndex)%H06_I_Aug
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 9,:) = RPARAM_in(segIndex)%H06_I_Sep
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF(10,:) = RPARAM_in(segIndex)%H06_I_Oct
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF(11,:) = RPARAM_in(segIndex)%H06_I_Nov
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF(12,:) = RPARAM_in(segIndex)%H06_I_Dec
            else
              array_size = shape(RCHFLX_out(iens,segIndex)%QPASTUP_IRF)
              past_length_I = array_size(2)
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF(simDatetime(1)%month(), 2:past_length_I) = RCHFLX_out(iens,segIndex)%QPASTUP_IRF(simDatetime(1)%month(), 1:past_length_I-1) ! shift the memory
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF(simDatetime(1)%month(), 1) = q_upstream ! allocate the current qupstream
            endif
            ! mean and updating the inflow parameters
            past_length_I = floor(RPARAM_in(segIndex)%H06_I_mem_L * 31    * secprday / dt)
            RPARAM_in(segIndex)%H06_I_Jan = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 1,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Mar = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 3,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_May = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 5,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Jul = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 7,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Aug = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 8,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Oct = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF(10,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Dec = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF(12,1:past_length_I))/past_length_I
            past_length_I = floor(RPARAM_in(segIndex)%H06_I_mem_L * 30    * secprday / dt)
            RPARAM_in(segIndex)%H06_I_Apr = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 4,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Jun = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 6,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Sep = sum(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 9,1:past_length_I))/past_length_I
            select case(trim(calendar))
              case('noleap','365_day')
                past_length_I = floor(RPARAM_in(segIndex)%H06_I_mem_L * 28    * secprday / dt)
              case ('standard','gregorian','proleptic_gregorian')
                past_length_I = floor(RPARAM_in(segIndex)%H06_I_mem_L * 28.25 * secprday / dt)
              case default;    ierr=20; message=trim(message)//'calendar name: '//trim(calendar)//' invalid'; return
            end select
            RPARAM_in(segIndex)%H06_I_Feb = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 2,1:past_length_I))/past_length_I
          endif

          ! create array with monthly inflow
          I_months = (/ RPARAM_in(segIndex)%H06_I_Jan, RPARAM_in(segIndex)%H06_I_Feb, &
                        RPARAM_in(segIndex)%H06_I_Mar, RPARAM_in(segIndex)%H06_I_Apr, &
                        RPARAM_in(segIndex)%H06_I_May, RPARAM_in(segIndex)%H06_I_Jun, &
                        RPARAM_in(segIndex)%H06_I_Jul, RPARAM_in(segIndex)%H06_I_Aug, &
                        RPARAM_in(segIndex)%H06_I_Sep, RPARAM_in(segIndex)%H06_I_Oct, &
                        RPARAM_in(segIndex)%H06_I_Nov, RPARAM_in(segIndex)%H06_I_Dec /)

          ! preserving the past demand discharge for lake models
          if ((RPARAM_in(segIndex)%H06_D_mem_F).and.(RCHFLX_out(iens,segIndex)%REACH_WM_FLUX/=realMissing).and.(is_flux_wm)) then ! if memeory is acive for this case then allocate the past input
            if (RCHFLX_out(iens,segIndex)%REACH_WM_FLUX<0) then ! demand cannot be negative for lake/reservoir
              RCHFLX_out(iens,segIndex)%REACH_WM_FLUX = 0._dp
            endif
            if (.not.allocated(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF)) then
              past_length_D = floor(RPARAM_in(segIndex)%H06_D_mem_L * 31 * secprday / dt)
              allocate(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(12,past_length_D),stat=ierr)
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 1,:) = RPARAM_in(segIndex)%H06_D_Jan
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 2,:) = RPARAM_in(segIndex)%H06_D_Feb
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 3,:) = RPARAM_in(segIndex)%H06_D_Mar
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 4,:) = RPARAM_in(segIndex)%H06_D_Apr
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 5,:) = RPARAM_in(segIndex)%H06_D_May
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 6,:) = RPARAM_in(segIndex)%H06_D_Jun
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 7,:) = RPARAM_in(segIndex)%H06_D_Jul
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 8,:) = RPARAM_in(segIndex)%H06_D_Aug
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 9,:) = RPARAM_in(segIndex)%H06_D_Sep
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(10,:) = RPARAM_in(segIndex)%H06_D_Oct
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(11,:) = RPARAM_in(segIndex)%H06_D_Nov
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(12,:) = RPARAM_in(segIndex)%H06_D_Dec
            else
              array_size = shape(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF)
              past_length_D = array_size(2)
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(simDatetime(1)%month(), 2:past_length_D) = RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(simDatetime(1)%month(), 1:past_length_D-1) ! shift the memory
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(simDatetime(1)%month(), 1) = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! allocate the current demand
            endif
            ! mean and updating the demand parameters
            past_length_D = floor(RPARAM_in(segIndex)%H06_D_mem_L * 31    * secprday / dt)
            RPARAM_in(segIndex)%H06_D_Jan = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 1,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Mar = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 3,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_May = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 5,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Jul = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 7,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Aug = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 8,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Oct = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(10,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Dec = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(12,1:past_length_D))/past_length_D
            past_length_D = floor(RPARAM_in(segIndex)%H06_D_mem_L * 30    * secprday / dt)
            RPARAM_in(segIndex)%H06_D_Apr = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 4,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Jun = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 6,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Sep = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 9,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Nov = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(11,1:past_length_D))/past_length_D
            select case(trim(calendar))
              case('noleap','365_day')
                past_length_D = floor(RPARAM_in(segIndex)%H06_D_mem_L * 28    * secprday / dt)
              case ('standard','gregorian','proleptic_gregorian')
                past_length_D = floor(RPARAM_in(segIndex)%H06_D_mem_L * 28.25 * secprday / dt)
              case default;    ierr=20; message=trim(message)//'calendar name: '//trim(calendar)//' invalid'; return
            end select
            RPARAM_in(segIndex)%H06_D_Feb = sum(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 2,1:past_length_D))/past_length_D
          endif
          ! create array with monthly demand
          D_months = (/ RPARAM_in(segIndex)%H06_D_Jan, RPARAM_in(segIndex)%H06_D_Feb, &
                        RPARAM_in(segIndex)%H06_D_Mar, RPARAM_in(segIndex)%H06_D_Apr, &
                        RPARAM_in(segIndex)%H06_D_May, RPARAM_in(segIndex)%H06_D_Jun, &
                        RPARAM_in(segIndex)%H06_D_Jul, RPARAM_in(segIndex)%H06_D_Aug, &
                        RPARAM_in(segIndex)%H06_D_Sep, RPARAM_in(segIndex)%H06_D_Oct, &
                        RPARAM_in(segIndex)%H06_D_Nov, RPARAM_in(segIndex)%H06_D_Dec /)

          ! calculate mean annual inflow and demand (to be integrated in condition not using of memory)
          I_yearly = sum(I_months)/months_per_yr
          D_yearly = sum(D_months)/months_per_yr

          ! calculate storage to yearly activity ratio
          c = RPARAM_in(segIndex)%H06_Smax/(I_yearly * days_per_yr * secprday)

          ! find start month of operational year
          do i=1,months_per_yr
            if (I_yearly <= I_months(i)) then
              start_month = i + 1
            endif
          enddo

          ! find start of operational year (add hour 1 when run hourly?) Once determined, this E_release should be communicated to the next timestep.
          if (simDatetime(1)%month() == start_month .AND. simDatetime(1)%day() == 1 ) then
             RPARAM_in(segIndex)%H06_E_rel_ini = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) / (RPARAM_in(segIndex)%H06_alpha * RPARAM_in(segIndex)%H06_Smax)
          endif

          ! Calculate target release
          if (RPARAM_in(segIndex)%H06_purpose == 1) then ! irrigation reservoir
            if (RPARAM_in(segIndex)%H06_envfact * I_yearly <= D_yearly) then ! larger demand than environmental flow requirement
              target_r = I_months(simDatetime(1)%month()) * RPARAM_in(segIndex)%H06_c1 + I_yearly * RPARAM_in(segIndex)%H06_c2 * (D_months(simDatetime(1)%month()) / D_yearly )
            else
              target_r = I_yearly + D_months(simDatetime(1)%month()) - D_yearly
            endif
          else ! non-irrigation reservoir
            target_r = I_yearly
          endif

          ! Calculate actual release
          if (c >= RPARAM_in(segIndex)%H06_c_compare) then ! multi-year reservoir
            RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = target_r * RPARAM_in(segIndex)%H06_E_rel_ini
          else if (0 <= c .AND. c < RPARAM_in(segIndex)%H06_c_compare) then
            RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = RPARAM_in(segIndex)%H06_E_rel_ini * target_r * (c / RPARAM_in(segIndex)%H06_denominator)**RPARAM_in(segIndex)%H06_exponent  + &
                                                   q_upstream * (1 - (c /RPARAM_in(segIndex)%H06_denominator)**RPARAM_in(segIndex)%H06_exponent)
          end if

          ! make sure reservoir volume does not drop below dead storage
          if (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) < (RPARAM_in(segIndex)%H06_Smax * RPARAM_in(segIndex)%H06_frac_Sdead)) then
            RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q - (RPARAM_in(segIndex)%H06_Smax * RPARAM_in(segIndex)%H06_frac_Sdead - RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) )/dt
            ! set negative outflow to zero
            if (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q<0) then
                RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q=0
            end if

          ! Account for spil overflow if reservoir is completely filled.
          else if (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) > RPARAM_in(segIndex)%H06_Smax) then
            RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q + (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RPARAM_in(segIndex)%H06_Smax)/ dt
          end if

          ! update the storage
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q * dt

        case(hype)
          ! update reach elevation
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_ELE = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) / RPARAM_in(segIndex)%HYP_A_avg &
                                                              + RPARAM_in(segIndex)%HYP_E_zero

          ! caclulate the day of calendar from 1st of January of current simulation year
          Day_of_year = simDatetime(1)%dayofyear()

          ! calculation of Fsin; sinusoidal aplication of flow
          F_sin = max(0._dp,(1+RPARAM_in(segIndex)%HYP_Qrate_amp*sin(2*pi*(Day_of_year+RPARAM_in(segIndex)%HYP_Qrate_phs)/365)))
          ! calculation of Flin; linear change in flow realted to the storage
          F_lin = min(max((RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_ELE-RPARAM_in(segIndex)%HYP_E_min)/(RPARAM_in(segIndex)%HYP_E_lim-RPARAM_in(segIndex)%HYP_E_min),0._dp),1._dp)
          ! flag for the model simulation in which is located in the area
          F_prim = 0
          if (RPARAM_in(segIndex)%HYP_prim_F) then
            F_prim = 1
          end if

          ! Q_main
          Q_prim = F_sin * F_lin * F_prim * RPARAM_in(segIndex)%HYP_Qrate_prim
          ! Q_spill
          Q_spill = 0._dp
          if (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_ELE > RPARAM_in(segIndex)%HYP_E_emr) then
            Q_spill = RPARAM_in(segIndex)%HYP_Qrate_emr* (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_ELE &
                      - RPARAM_in(segIndex)%HYP_E_emr)**RPARAM_in(segIndex)%HYP_Erate_emr
          end if

          if (RPARAM_in(segIndex)%HYP_Qsim_mode) then
            ! Q_sim
            Q_sim = Q_prim + Q_spill
          else
            ! Original implementation picks the maximum value of output
            ! from primary spillway and emergency spillway
            Q_sim = max(Q_prim, Q_spill)
          end if

          ! check if the output is not more than the existing stored water
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = min(Q_sim, max(0._dp,(RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_ELE-RPARAM_in(segIndex)%HYP_E_min)*RPARAM_in(segIndex)%HYP_A_avg)/dt)

          ! update the storage
          RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q * dt

        case default; ierr=20; message=trim(message)//'unable to identify the parametric lake model type'; return
      end select
    endif

!    TO-do: Make sure that WB check below is performed correctly in comp_reach_wb subroutine
!    ! calculate water balance (in this water balance we dont have the actual evaporation, assuming there is enough water for evaporation)
!    WB = q_upstream * dt + RCHFLX_out(iens,segIndex)%basinprecip * dt - RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q * dt &
!    - RCHFLX_out(iens,segIndex)%basinevapo * dt - RCHFLX_out(iens,segIndex)%REACH_WM_FLUX_actual * dt &
!    - (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(0))
!
!    if ((1._dp<WB).or.(segIndex==ixDesire)) then; ! larger than 1 cubic meter or desired lake
!      ! NOTE: The lake discharge and storage need to be solved iterative way to reduce water balance error
!      write(iulog,*) 'Water balance for lake ID = ', NETOPO_in(segIndex)%REACHID, ' excees the Tolerance'
!      write(iulog,'(A,1PG15.7)') 'WBerr [m3]       = ', WB
!      write(iulog,'(A,1PG15.7)') 'dS [m3]          = ', RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1)-RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(0)
!      write(iulog,'(A,1PG15.7)') 'S [m3]           = ', RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_VOL(1)
!      write(iulog,'(A,1PG15.7)') 'inflow [m3]      = ', q_upstream*dt
!      write(iulog,'(A,1PG15.7)') 'outflow [m3]     = ', RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q*dt
!      write(iulog,'(A,1PG15.7)') 'precip [m3]      = ', RCHFLX_out(iens,segIndex)%basinprecip*dt
!      write(iulog,'(A,1PG15.7)') 'evaporation [m3] = ', RCHFLX_out(iens,segIndex)%basinevapo*dt
!      if (is_flux_wm) then
!        write(iulog,'(A,1PG15.7)') 'actual abstraction [m3] = ', RCHFLX_out(iens,segIndex)%REACH_WM_FLUX_actual
!      endif
!      if (1<WB) then
!        cmessage = 'Water balance for lake ID = exceeds the Tolerance'
!        ierr = 1; message=trim(message)//trim(cmessage);
!      endif
!    endif

    call comp_reach_wb(NETOPO_in(segIndex)%REACHID, idxRoute, q_upstream, RCHFLX_out(iens,segIndex)%BASIN_QR(1), RCHFLX_out(iens,segIndex), &
                       verbose, lakeFlag=.true.,tolerance=lakeWBtol)

  END SUBROUTINE lake_route

END MODULE lake_route_module
