module lake_route_module

  !numeric type
 USE nrtype
 ! data type
 USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
 USE dataTypes, ONLY: RCHTOPO        ! Network topology
 USE dataTypes, ONLY: RCHPRP         ! Network parameter
 ! global parameters
 USE public_var, ONLY: iulog          ! i/o logical unit number
 USE public_var, ONLY: realMissing    ! missing value for real number
 USE public_var, ONLY: integerMissing ! missing value for integer number

  ! privary
 implicit none
 private

  public::lake_route

  contains

  ! *********************************************************************
  ! subroutine: perform one segment lake routing
  ! *********************************************************************
  subroutine lake_route(&
                         ! input
                         iEns,       &    ! input: index of runoff ensemble to be processed
                         segIndex,   &    ! input: index of runoff ensemble to be processed
                         ixDesire,   &    ! input: reachID to be checked by on-screen pringing (here reachID can be lake)
                         NETOPO_in,  &    ! input: reach topology data structure
                         RPARAM_in,  &    ! input: reach parameter data strcuture
                         ! inout
                         RCHFLX_out, &    ! inout: reach flux data structure
                         ! output
                         ierr, message)   ! output: error control

  USE globalData, ONLY: modTime           ! previous and current model time
  USE public_var, ONLY: is_flux_wm        ! logical water management components fluxes should be read
  USE public_var, ONLY: dt, lakeWBTol     ! lake water balance tolerance
  USE public_var, ONLY: lake_model_D03    ! logical whether or not lake should be simulated
  USE public_var, ONLY: lake_model_H06    ! logical whether or not lake should be simulated
  USE public_var, ONLY: secprday, days_per_yr, months_per_yr    ! time constants

  implicit none
  ! Input
  INTEGER(I4B), intent(in)                 :: iEns           ! runoff ensemble to be routed
  INTEGER(I4B), intent(in)                 :: segIndex       ! segment where routing is performed
  INTEGER(I4B), intent(in)                 :: ixDesire       ! index of the reach for verbose output
  type(RCHTOPO), intent(in),   allocatable :: NETOPO_in(:)   ! River Network topology
  type(RCHPRP), intent(inout), allocatable :: RPARAM_in(:)   ! River Network topology
  ! inout
  TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)! Reach fluxes (ensembles, space [reaches]) for decomposed domains
  ! Output
  integer(i4b), intent(out)                :: ierr           ! error code
  character(*), intent(out)                :: message        ! error message
  ! Local variables to
  real(dp)                                 :: q_upstream     ! total discharge at top of the reach being processed
  real(dp)                                 :: WB             ! water balance component in the lake
  type(STRFLX), allocatable                :: fluxstate(:)   ! upstream Reach fluxes
  INTEGER(I4B)                             :: nUps           ! number of upstream segment
  INTEGER(I4B)                             :: iUps           ! upstream reach index
  INTEGER(I4B)                             :: iRch_ups       ! index of upstream reach in NETOPO
  INTEGER(I4B)                             :: ntdh           ! number of time steps in IRF
  character(len=strLen)                    :: cmessage       ! error message from subroutine
  ! local variables for H06 routine
  real(dp)                                 :: c                   ! storage to yearly activity ratio
  real(dp)                                 :: I_yearly, D_yearly  ! mean annual inflow and demand
  real(dp), dimension(12)                  :: I_months, D_months  ! mean monthly inflow and demand
  INTEGER(I4B), allocatable                :: array_size (:)      ! get the size of array_size
  INTEGER(I4B)                             :: start_month=0       ! start month of the operational year
  INTEGER(I4B)                             :: i                   ! index
  INTEGER(I4B)                             :: past_length_I       ! pas length for inflow based on length in year and floor
  INTEGER(I4B)                             :: past_length_D       ! pas length for demand based on length in year and floor
  real(dp)                                 :: target_r            ! target release

  print*, 'inside lake, time at the model simulation',modTime(1)%iy,modTime(1)%im,modTime(1)%id,modTime(1)%ih,modTime(1)%imin,modTime(1)%dsec

    ! initialize error control
    ierr=0; message='lake_route/'

    ! identify number of upstream segments of the lake being processed
    nUps = size(NETOPO_in(segIndex)%UREACHI)

    allocate(fluxstate(nUps), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//': fluxstate'; return; endif

    ! loop to get the streamlow of the usptream discharges to the lake
    if (nUps>0) then
      do iUps = 1,nUps
        iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
        fluxstate(iUps) = RCHFLX_out(iens,iRch_ups)
      end do
    endif

    ! Find out total q draining into a lake
    q_upstream = 0.0_dp
    if(nUps>0)then
      do iUps = 1,nUps
        q_upstream = q_upstream + fluxstate(iUps)%REACH_Q_IRF
      end do
    endif

    ! perform lake routing based on a fixed storage discharge relationship Q=kS
    ! no runoff input is added to the lake; the input are only precipitation and evaporation to the lake

    ! if(NETOPO_in(segIndex)%REACHIX == ixDesire)then   ! uncommnet when the ixDesire is fixed and not -9999
    !print*, '------lake-simulation-------- '
    !print*, 'node id that is lake .......= ', NETOPO_in(segIndex)%REACHID ! to check the reach id of lake
    !print*, 'lake param D03_MaxStorage ..= ', RPARAM_in(segIndex)%D03_MaxStorage
    !print*, 'lake param D03_Coefficient .= ', RPARAM_in(segIndex)%D03_Coefficient
    !print*, 'lake param D03_Power .......= ', RPARAM_in(segIndex)%D03_Power
    !print*, 'lake param H06_Smax ........= ', RPARAM_in(segIndex)%H06_Smax
    !print*, 'lake param H06_alpha .......= ', RPARAM_in(segIndex)%H06_alpha
    !print*, 'lake param H06_S_ini .......= ', RPARAM_in(segIndex)%H06_S_ini
    !print*, 'lake target volum ..........= ', NETOPO_in(segIndex)%LakeTargVol
     print*, 'volume before simulation m3.= ', RCHFLX_out(iens,segIndex)%REACH_VOL(0)
    !print*, 'upstream streamflow m3/s ...= ', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
    !print*, 'upstream precipitation m3/s.= ', RCHFLX_out(iens,segIndex)%basinprecip
    !print*, 'upstream evaporation m3/s ..= ', RCHFLX_out(iens,segIndex)%basinevapo
    !print*, 'paraemters', RPARAM_in(segIndex)%D03_MaxStorage, RPARAM_in(segIndex)%D03_Coefficient, RPARAM_in(segIndex)%D03_Power, NETOPO_in(segIndex)%LakeTargVol, NETOPO_in(segIndex)%islake, NETOPO_in(segIndex)%LakeModelType


    ! add upstream, precipitation and subtract evaporation from the lake volume
    RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(0) ! updating storage for current time
    RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) + q_upstream * dt  ! input upstream discharge from m3/s to m3
    RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) + RCHFLX_out(iens,segIndex)%basinprecip * dt ! input lake precipitation
    RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%basinevapo * dt ! output lake evaporation
    if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) < 0) then; ! to avoid negative lake volume
       RCHFLX_out(iens,segIndex)%REACH_VOL(1)=0
    endif

    ! simulating lake/reservoir output
    if (NETOPO_in(segIndex)%LakeTargVol) then ! The lake should follow the given target volume

      if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) < RCHFLX_out(iens,segIndex)%REACH_WM_VOL) then ! The storage is smaller, accumulate
        RCHFLX_out(iens,segIndex)%REACH_Q_IRF  = 0 ! lake/reservoir output is set to zero (0)
      else ! The storage is largrt, release to get to the target level
        RCHFLX_out(iens,segIndex)%REACH_Q_IRF = (RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_WM_VOL) /dt
        RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_WM_VOL
      endif

      !print*, "inside the lake follow target", RCHFLX_out(iens,segIndex)%REACH_WM_VOL, RCHFLX_out(iens,segIndex)%REACH_VOL(1)

    else ! if the lake is paraemteric

      !print*, "lake model type", NETOPO_in(segIndex)%LakeModelType
      select case(NETOPO_in(segIndex)%LakeModelType)

        case (1)
          ! the model is Doll03
          ! print*, "lake model is Doll 2003"
          RCHFLX_out(iens,segIndex)%REACH_Q_IRF = RPARAM_in(segIndex)%D03_Coefficient * RCHFLX_out(iens,segIndex)%REACH_VOL(1) * &
                                                   (RCHFLX_out(iens,segIndex)%REACH_VOL(1) / RPARAM_in(segIndex)%D03_MaxStorage) ** &
                                                   RPARAM_in(segIndex)%D03_Power! Q = AS(S/Smax)^B based on Eq. 1 Hanasaki et al., 2006 https://doi.org/10.1016/j.jhydrol.2005.11.011
          ! in case is the output volume is more than lake volume
          RCHFLX_out(iens,segIndex)%REACH_Q_IRF = (min(RCHFLX_out(iens,segIndex)%REACH_Q_IRF * dt, RCHFLX_out(iens,segIndex)%REACH_VOL(1)) )/dt
          ! updating the storage
          RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_Q_IRF * dt
          if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) < 0) then; ! set the lake volume as 0 if it goes negative
            RCHFLX_out(iens,segIndex)%REACH_VOL(1)=0
          endif

        case (2)
          ! the model is Hanasaki06
          ! print*, "lake model is Hanasaki 2006"

          ! preserving the past upstrem discharge for lake models
          ! create memory of upstream inflow for Hanasaki formulation if memory flag is active for this lake
          if RCHFLX_out(iens,segIndex)%H06_I_mem_F then ! if memeory is acive for this case then allocate the past input
            if (.not.allocated(RCHFLX_out(iens,segIndex)%QPASTUP_IRF)) then
              past_length_I = floor(RCHFLX_out(iens,segIndex)%H06_I_mem_L * 31 * 3600 * 24 / dt)
              allocate(RCHFLX_out(iens,segIndex)%QPASTUP_IRF(12,past_length_I,stat=ierr)
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
              array_size = size(RCHFLX_out(iens,segIndex)%QPASTUP_IRF)
              past_length_I = array_size(2)
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF(modTime(1)%im, 2:past_length_I) = RCHFLX_out(iens,segIndex)%QPASTUP_IRF(modTime(1)%im, 1:past_length_I-1) ! shift the memory
              RCHFLX_out(iens,segIndex)%QPASTUP_IRF(modTime(1)%im, 1) = q_upstream ! allocate the current qupstream
            endif
            ! mean and updating the inflow parameters
            past_length_I = floor(RCHFLX_out(iens,segIndex)%H06_I_mem_L * 31 * 3600 * 24 / dt)
            RPARAM_in(segIndex)%H06_I_Jan = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 1,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Mar = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 3,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_May = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 5,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Jul = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 7,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Aug = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 8,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Oct = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF(10,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Dec = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF(12,1:past_length_I))/past_length_I
            past_length_I = floor(RCHFLX_out(iens,segIndex)%H06_I_mem_L * 28.25 * 3600 * 24 / dt)
            RPARAM_in(segIndex)%H06_I_Feb = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 2,1:past_length_I)/past_length_I
            past_length_I = floor(RCHFLX_out(iens,segIndex)%H06_I_mem_L * 30 * 3600 * 24 / dt)
            RPARAM_in(segIndex)%H06_I_Apr = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 4,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Jun = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 6,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Sep = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF( 9,1:past_length_I))/past_length_I
            RPARAM_in(segIndex)%H06_I_Nov = SUM(RCHFLX_out(iens,segIndex)%QPASTUP_IRF(11,1:past_length_I))/past_length_I
          endif
          !print*, RCHFLX_out(iens,segIndex)%QPASTUP_IRF
          ! create array with monthly inflow
          I_months = (/ RPARAM_in(segIndex)%H06_I_Jan, RPARAM_in(segIndex)%H06_I_Feb, &
                        RPARAM_in(segIndex)%H06_I_Mar, RPARAM_in(segIndex)%H06_I_Apr, &
                        RPARAM_in(segIndex)%H06_I_May, RPARAM_in(segIndex)%H06_I_Jun, &
                        RPARAM_in(segIndex)%H06_I_Jul, RPARAM_in(segIndex)%H06_I_Aug, &
                        RPARAM_in(segIndex)%H06_I_Sep, RPARAM_in(segIndex)%H06_I_Oct, &
                        RPARAM_in(segIndex)%H06_I_Nov, RPARAM_in(segIndex)%H06_I_Dec /)

          ! preserving the past demand discharge for lake models
          if ((RCHFLX_out(iens,segIndex)%H06_D_mem_F).and.(RCHFLX_out(iens,segIndex)%REACH_WM_FLUX/=realMissing).and.(is_flux_wm)) then ! if memeory is acive for this case then allocate the past input
            if (RCHFLX_out(iens,segIndex)%REACH_WM_FLUX<0) then ! demand cannot be negative for lake/reservoir
              RCHFLX_out(iens,segIndex)%REACH_WM_FLUX = 0._dp
            endif
            if (.not.allocated(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF)) then
              past_length_D = floor(RCHFLX_out(iens,segIndex)%H06_D_mem_L * 31 * 3600 * 24 / dt)
              allocate(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(12,past_length_D,stat=ierr)
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
              array_size = size(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF)
              past_length_D = array_size(2)
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(modTime(1)%im, 2:past_length_D) = RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(modTime(1)%im, 1:past_length_D-1) ! shift the memory
              RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(modTime(1)%im, 1) = RCHFLX_out(iens,segIndex)%REACH_WM_FLUX ! allocate the current demand
            endif
            ! mean and updating the demand parameters
            past_length_D = floor(RCHFLX_out(iens,segIndex)%H06_D_mem_L * 31 * 3600 * 24 / dt)
            RPARAM_in(segIndex)%H06_D_Jan = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 1,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Mar = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 3,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_May = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 5,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Jul = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 7,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Aug = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 8,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Oct = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(10,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Dec = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(12,1:past_length_D))/past_length_D
            past_length_D = floor(RCHFLX_out(iens,segIndex)%H06_D_mem_L * 28.25 * 3600 * 24 / dt)
            RPARAM_in(segIndex)%H06_D_Feb = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 2,1:past_length_D))/past_length_D
            past_length_D = floor(RCHFLX_out(iens,segIndex)%H06_D_mem_L * 30 * 3600 * 24 / dt)
            RPARAM_in(segIndex)%H06_D_Apr = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 4,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Jun = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 6,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Sep = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF( 9,1:past_length_D))/past_length_D
            RPARAM_in(segIndex)%H06_D_Nov = SUM(RCHFLX_out(iens,segIndex)%DEMANDPAST_IRF(11,1:past_length_D))/past_length_D
          endif
          ! create array with monthly demand
          D_months = (/ RPARAM_in(segIndex)%H06_D_Jan, RPARAM_in(segIndex)%H06_D_Feb, &
                        RPARAM_in(segIndex)%H06_D_Mar, RPARAM_in(segIndex)%H06_D_Apr, &
                        RPARAM_in(segIndex)%H06_D_May, RPARAM_in(segIndex)%H06_D_Jun, &
                        RPARAM_in(segIndex)%H06_D_Jul, RPARAM_in(segIndex)%H06_D_Aug, &
                        RPARAM_in(segIndex)%H06_D_Sep, RPARAM_in(segIndex)%H06_D_Oct, &
                        RPARAM_in(segIndex)%H06_D_Nov, RPARAM_in(segIndex)%H06_D_Dec /)

          ! calculate mean annual inflow and demand (to be integrated in condition not using of memory)
          I_yearly = SUM(I_months)/months_per_yr

          D_yearly = SUM(D_months)/months_per_yr

          ! calculate storage to yearly activity ratio
          c = RPARAM_in(segIndex)%H06_Smax/(I_yearly * days_per_yr * secprday)

          ! find start month of operational year
          do i=1,months_per_yr
            if (I_yearly <= I_months(i)) then
                start_month = i + 1
            endif
          enddo

          print*, 'start month', start_month

          ! find start of operational year (add hour 1 when run hourly?) Once determined, this E_release should be communicated to the next timestep.
          if (modTime(1)%im == start_month .AND. modTime(1)%id == 1 ) then
             RPARAM_in(segIndex)%H06_E_rel_ini = RCHFLX_out(iens,segIndex)%REACH_VOL(1) / (RPARAM_in(segIndex)%H06_alpha * RPARAM_in(segIndex)%H06_Smax)
          endif

          ! print*,'E_release ', RPARAM_in(segIndex)%H06_E_rel_ini

          ! Calculate target release
          if (RPARAM_in(segIndex)%H06_purpose == 1) then ! irrigation reservoir

            if (RPARAM_in(segIndex)%H06_envfact * I_yearly <= D_yearly) then ! larger demand than environmental flow requirement
                target_r = I_months(modTime(1)%im) * RPARAM_in(segIndex)%H06_c1 + I_yearly * RPARAM_in(segIndex)%H06_c2 * (D_months(modTime(1)%im) / D_yearly )
            else
                target_r = I_yearly + D_months(modTime(1)%im) - D_yearly
            endif

          else ! non-irrigation reservoir
            target_r = I_yearly

          endif

          ! Calculate actual release
          if (c >= RPARAM_in(segIndex)%H06_c_compare) then ! multi-year reservoir
            RCHFLX_out(iens,segIndex)%REACH_Q_IRF = target_r * RPARAM_in(segIndex)%H06_E_rel_ini
            print*,'multi-year reservoir'
          else if (0 <= c .AND. c < RPARAM_in(segIndex)%H06_c_compare) then
            RCHFLX_out(iens,segIndex)%REACH_Q_IRF = RPARAM_in(segIndex)%H06_E_rel_ini * target_r * (c / RPARAM_in(segIndex)%H06_denominator)**RPARAM_in(segIndex)%H06_exponent  + &
                                                   q_upstream * (1 - (c /RPARAM_in(segIndex)%H06_denominator)**RPARAM_in(segIndex)%H06_exponent)
            print*,'whithin-a-year reservoir'
          end if

          ! make sure reservoir volume does not drop below dead storage
          if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) < (RPARAM_in(segIndex)%H06_Smax * RPARAM_in(segIndex)%H06_frac_Sdead)) then
            RCHFLX_out(iens,segIndex)%REACH_Q_IRF = RCHFLX_out(iens,segIndex)%REACH_Q_IRF - (RPARAM_in(segIndex)%H06_Smax * RPARAM_in(segIndex)%H06_frac_Sdead - RCHFLX_out(iens,segIndex)%REACH_VOL(1) )/secprday
            print*, 'below dead storage'
            ! set negative outflow to zero
            if (RCHFLX_out(iens,segIndex)%REACH_Q_IRF<0) then
                RCHFLX_out(iens,segIndex)%REACH_Q_IRF=0
            end if

          ! Account for spil overflow if reservoir is completely filled.
          else if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) > RPARAM_in(segIndex)%H06_Smax) then
            print*, 'overflow evoked'
            RCHFLX_out(iens,segIndex)%REACH_Q_IRF = RCHFLX_out(iens,segIndex)%REACH_Q_IRF + (RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RPARAM_in(segIndex)%H06_Smax)/ secprday
          end if

          print*,'outflow ', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
          ! update the storage
          RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_Q_IRF * dt

          ! print*, modTime(1)%im ! month of the simulations
          ! print*, 'Hanasaki parameters'
          print*, RPARAM_in(segIndex)%H06_Smax, RPARAM_in(segIndex)%H06_alpha, RPARAM_in(segIndex)%H06_envfact, RPARAM_in(segIndex)%H06_S_ini, RPARAM_in(segIndex)%H06_c1, RPARAM_in(segIndex)%H06_c2, RPARAM_in(segIndex)%H06_exponent, RPARAM_in(segIndex)%H06_I_Feb, RPARAM_in(segIndex)%H06_D_Feb

        case default; ierr=20; message=trim(message)//'unable to identify the parametric lake model type'; return
      end select

    endif

    ! calculate water balance (in this water balance we dont have the actual evaporation, assuming there is enough water for evaporation)
    WB = q_upstream * dt + RCHFLX_out(iens,segIndex)%basinprecip * dt - RCHFLX_out(iens,segIndex)%REACH_Q_IRF * dt &
    - RCHFLX_out(iens,segIndex)%basinevapo * dt - (RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%REACH_VOL(0))


    !if(NETOPO_in(segIndex)%REACHIX == ixDesire)then    ! uncommnet when the ixDesire is fixed and not -9999
    !print*, 'lake simulated output m3/s .= ', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
    !print*, 'volume after simulation m3 .= ', RCHFLX_out(iens,segIndex)%REACH_VOL(1)
    !print*, 'read target volume     m3 .= ', RCHFLX_out(iens,segIndex)%REACH_WM_FLUX
    !print*, 'water balance error ........= ', WB
    !endif

    if(lakeWBTol<WB)then;
      ! NOTE: The lake discharge and storage need to be solved iterative way to reduce water balance error
      write(iulog,*) 'Water balance for lake ID = ', NETOPO_in(segIndex)%REACHID, ' excees the Tolerance'
      write(iulog,'(A,1PG15.7)') 'WBerr [m3]       = ', WB
      write(iulog,'(A,1PG15.7)') 'dS [m3]          = ', RCHFLX_out(iens,segIndex)%REACH_VOL(1)-RCHFLX_out(iens,segIndex)%REACH_VOL(0)
      write(iulog,'(A,1PG15.7)') 'inflow [m3]      = ', q_upstream*dt
      write(iulog,'(A,1PG15.7)') 'outflow [m3]     = ', RCHFLX_out(iens,segIndex)%REACH_Q_IRF*dt
      write(iulog,'(A,1PG15.7)') 'precip [m3]      = ', RCHFLX_out(iens,segIndex)%basinprecip*dt
      write(iulog,'(A,1PG15.7)') 'evaporation [m3] = ', RCHFLX_out(iens,segIndex)%basinevapo*dt
      cmessage = 'Water balance for lake ID = exceeds the Tolerance'
      ierr = 1; message=trim(message)//trim(cmessage);
    endif

    ! set the routed flag as .True.
    RCHFLX_out(iEns,segIndex)%isRoute=.True.

    ! pass the current storage for the past time step for the next time step simulation
    RCHFLX_out(iens,segIndex)%REACH_VOL(0) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) !shift on time step back

    ! assign the zero value as lake do not have a QFUTURE_IRF
    if (.not.allocated(RCHFLX_out(iens,segIndex)%QFUTURE_IRF))then
      ntdh = size(NETOPO_in(segIndex)%UH)
      allocate(RCHFLX_out(iens,segIndex)%QFUTURE_IRF(ntdh), stat=ierr, errmsg=cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage)//': RCHFLX_out(iens,segIndex)%QFUTURE_IRF'; return; endif
      RCHFLX_out(iens,segIndex)%QFUTURE_IRF(1:ntdh) = 0._dp
    end if

  end subroutine lake_route

end module lake_route_module
