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
  type(RCHPRP),  intent(in),   allocatable :: RPARAM_in(:)   ! River Network topology
  ! inout
  TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
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
  real(dp)                                 :: c              ! storage to yearly activity ratio
  real(dp)                                 :: I_yearly, D_yearly  ! mean annual inflow and demand
  real(dp), dimension(12)                  :: I_months, D_months  ! mean monthly inflow and demand  
  INTEGER(I4B)                             :: start_month=0     ! start month of the operational year
  INTEGER(I4B)                             :: i             ! index
  real(dp)                                 :: E_release     ! release coefficient
  real(dp)                                 :: target_r      ! target release


  print*, 'inside lake, time at the mdoel simulation',modTime(1)%iy,modTime(1)%im,modTime(1)%id,modTime(1)%ih,modTime(1)%imin,modTime(1)%dsec

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
    !print*, 'volume before simulation m3.= ', RCHFLX_out(iens,segIndex)%REACH_VOL(0)
    !print*, 'upstream streamflow m3/s ...= ', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
    !print*, 'upstream precipitation m3/s.= ', RCHFLX_out(iens,segIndex)%basinprecip
    !print*, 'upstream evaporation m3/s ..= ', RCHFLX_out(iens,segIndex)%basinevapo
    !print*, 'paraemters', RPARAM_in(segIndex)%D03_MaxStorage, RPARAM_in(segIndex)%D03_Coefficient, RPARAM_in(segIndex)%D03_Power, NETOPO_in(segIndex)%LakeTargVol, NETOPO_in(segIndex)%islake, NETOPO_in(segIndex)%LakeModelType


    ! add upstream, precipitation and subtract evaporation from the lake volume
    RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(0) ! updating storage for current time
    RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) + q_upstream * dt  ! input upstream discharge from m3/s to m3
    RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) + RCHFLX_out(iens,segIndex)%basinprecip * dt ! input lake precipitation
    RCHFLX_out(iens,segIndex)%REACH_VOL(1) = RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RCHFLX_out(iens,segIndex)%basinevapo * dt ! output lake evaporaiton
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

      !print*, "lake model is Doll 2003", NETOPO_in(segIndex)%LakeModelType
      select case(NETOPO_in(segIndex)%LakeModelType)

        case (1)
          ! the model is Doll03
          !print*, "lake model is Doll 2003"
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
          ! preserving the past upstrem discharge for lake models
          ! print*, "lake model is Hanasaki 2006"

          ! create memory of upstream inflow for Hanasaki formulation
          if (.not.allocated(RCHFLX_out(iens,segIndex)%QPASTUP_IRF)) then ! it is the first time step and should be allocated
            allocate(RCHFLX_out(iens,segIndex)%QPASTUP_IRF(10),stat=ierr)
          else
            RCHFLX_out(iens,segIndex)%QPASTUP_IRF(2:10) = RCHFLX_out(iens,segIndex)%QPASTUP_IRF(1:9) ! here the length is 10 as a randome varibale
            RCHFLX_out(iens,segIndex)%QPASTUP_IRF(1)    = q_upstream ! code needed to shift this as well.
          endif
          !print*, RCHFLX_out(iens,segIndex)%QPASTUP_IRF


          ! create array with monthly inflow
          I_months = (/ RPARAM_in(segIndex)%H06_I_Jan, RPARAM_in(segIndex)%H06_I_Feb, RPARAM_in(segIndex)%H06_I_Mar, RPARAM_in(segIndex)%H06_I_Apr, RPARAM_in(segIndex)%H06_I_May, RPARAM_in(segIndex)%H06_I_Jun, & 
                      RPARAM_in(segIndex)%H06_I_Jul, RPARAM_in(segIndex)%H06_I_Aug, RPARAM_in(segIndex)%H06_I_Sep, RPARAM_in(segIndex)%H06_I_Oct, RPARAM_in(segIndex)%H06_I_Nov, RPARAM_in(segIndex)%H06_I_Dec /)
          
          ! there is a problem with reading monthly inflow parameters, for testing puposes, replace by hardcoded values. 
          ! print*,'I_months',I_months 
          I_months = (/8.64614286,   8.58854822,   9.8415576 ,  18.64570952,  57.96402304, 132.58282381, 59.53373272, 41.27400922, 42.0379619, 27.56602765, 12.10262381, 9.80076959 /)
         
          ! Bhumiboi inflow
          ! I_months =  (/42.80004759, 23.10134666,  10.9559713,   15.51464843,  63.09036455, 121.77380353, 115.38119685, 291.89079013, 511.78662059, 423.13451822, 237.26578135,  98.9666038 /)
          

          ! create array with monthly inflow
          D_months = (/ RPARAM_in(segIndex)%H06_D_Jan, RPARAM_in(segIndex)%H06_D_Feb, RPARAM_in(segIndex)%H06_D_Mar, RPARAM_in(segIndex)%H06_D_Apr, RPARAM_in(segIndex)%H06_D_May, RPARAM_in(segIndex)%H06_D_Jun, & 
                      RPARAM_in(segIndex)%H06_D_Jul, RPARAM_in(segIndex)%H06_D_Aug, RPARAM_in(segIndex)%H06_D_Sep, RPARAM_in(segIndex)%H06_D_Oct, RPARAM_in(segIndex)%H06_D_Nov, RPARAM_in(segIndex)%H06_D_Dec /)
          ! print*,'D_months',D_months            
          ! D_months = (/0,150,200,250,200,75,50,50,20,0,30,0/)
          
          ! calculate mean annual inflow and demand (to be integrated in condition not using of memory)
          I_yearly = SUM(I_months)/months_per_yr
          print*,'I_yearly',I_yearly
 
          D_yearly = SUM(D_months)/months_per_yr
          
          ! calculate storage to yearly activity ratio
          c = RPARAM_in(segIndex)%H06_Smax/(I_yearly * days_per_yr * secprday)
          print*,'c', c

          ! find start month of operational year 
          do i=1,months_per_yr
            if (I_yearly <= I_months(i)) then 
                start_month = i + 1
            endif
          enddo
          
          print*, 'start month', start_month

          ! if operational year has not yet started (add condition!!), determine based on initial storage
          E_release = RPARAM_in(segIndex)%H06_S_ini/(RPARAM_in(segIndex)%H06_alpha * RPARAM_in(segIndex)%H06_Smax)
          

          ! find start of operational year (add hour 1 when run hourly?) Once determined, this E_release should be communicated to the next timestep. 
          if (modTime(1)%im == start_month .AND. modTime(1)%id == 1 ) then
             E_release = RCHFLX_out(iens,segIndex)%REACH_VOL(1) / (RPARAM_in(segIndex)%H06_alpha * RPARAM_in(segIndex)%H06_Smax)
          endif
          
          print*,'E_release', E_release
          
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
            RCHFLX_out(iens,segIndex)%REACH_Q_IRF = target_r * E_release
            print*,'multi-year reservoir'
            print*, 'target_r*E_release', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
          else if (0 <= c .AND. c < RPARAM_in(segIndex)%H06_c_compare) then
            RCHFLX_out(iens,segIndex)%REACH_Q_IRF = E_release * target_r * (c / RPARAM_in(segIndex)%H06_denominator)**RPARAM_in(segIndex)%H06_exponent  + & 
                                                   q_upstream * (1 - (c /RPARAM_in(segIndex)%H06_denominator)**RPARAM_in(segIndex)%H06_exponent) 
            print*,'whithin-a-year reservoir'
            print*, 'first part', E_release * target_r * (c / RPARAM_in(segIndex)%H06_denominator)**RPARAM_in(segIndex)%H06_exponent 
            print*, 'second part', q_upstream * (1 - (c /RPARAM_in(segIndex)%H06_denominator)**RPARAM_in(segIndex)%H06_exponent) 
            print*, 'q_upstream', q_upstream
            print*, 'RCHFLX_out(iens,segIndex)%REACH_Q_IRF', RCHFLX_out(iens,segIndex)%REACH_Q_IRF
          end if 
          
          ! make sure reservoir volume does not drop below dead storage
          if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) < (RPARAM_in(segIndex)%H06_Smax * RPARAM_in(segIndex)%H06_frac_Sdead)) then
            RCHFLX_out(iens,segIndex)%REACH_Q_IRF = RCHFLX_out(iens,segIndex)%REACH_Q_IRF - (RPARAM_in(segIndex)%H06_Smax * RPARAM_in(segIndex)%H06_frac_Sdead - RCHFLX_out(iens,segIndex)%REACH_VOL(1) )
          ! Account for spil overflow if reservoir is completely filled.
          else if (RCHFLX_out(iens,segIndex)%REACH_VOL(1) > RPARAM_in(segIndex)%H06_Smax) then
            print*, 'overflow evoked'
            RCHFLX_out(iens,segIndex)%REACH_Q_IRF = RCHFLX_out(iens,segIndex)%REACH_Q_IRF + (RCHFLX_out(iens,segIndex)%REACH_VOL(1) - RPARAM_in(segIndex)%H06_Smax)/ secprday
          end if
          
          !print*, modTime(1)%im ! month of the simulations
          print*, 'Hanasaki parameters'
          print*, RPARAM_in(segIndex)%H06_Smax, RPARAM_in(segIndex)%H06_alpha, RPARAM_in(segIndex)%H06_envfact, RPARAM_in(segIndex)%H06_S_ini, RPARAM_in(segIndex)%H06_c1, RPARAM_in(segIndex)%H06_c2, RPARAM_in(segIndex)%H06_exponent

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
      print*, 'Water balance for lake ID = ', NETOPO_in(segIndex)%REACHID, 'exceeds the Tolerance'
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
