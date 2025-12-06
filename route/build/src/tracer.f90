MODULE tracer_module

! Description:
! route a solute or constituent in water to downstream reach
!
! Note:
! As of Apr 17, 2025, in-stream processes - advection, dispersion and solute's reaction to environment
! are not implemented, though a draft subroutine - comp_advec_diffusion exsit.
! A solute moves with water movement (i.e., reach-to-reach advection),
! and solute concentration is assumed to be homogeneous within a reach

USE nrtype
USE dataTypes,     ONLY: STRFLX                ! fluxes in each reach
USE dataTypes,     ONLY: SUBRCH                ! states at sub-reach
USE dataTypes,     ONLY: STRSTA                ! state in each reach
USE dataTypes,     ONLY: RCHTOPO               ! Network topology
USE dataTypes,     ONLY: RCHPRP                ! Reach parameter
USE public_var,    ONLY: iulog                 ! i/o logical unit number
USE public_var,    ONLY: desireId              ! ID or reach where detailed reach state is print in log
USE public_var,    ONLY: dt                    ! simulation time step [sec]
USE public_var,    ONLY: hw_drain_point        ! headwater catchment pour point (top_reach==1 or bottom_reach==2)
USE public_var,    ONLY: min_length_route      ! minimum reach length for routing to be performed.
USE public_var,    ONLY: kinematicWave
USE public_var,    ONLY: muskingumCunge
USE public_var,    ONLY: diffusiveWave
USE water_balance, ONLY: comp_reach_mb         ! compute water balance error
USE hydraulic,     ONLY: flow_depth
USE hydraulic,     ONLY: flow_area

implicit none

private
public::constituent_rch

integer(i4b), parameter :: top_reach=1
integer(i4b), parameter :: bottom_reach=2

CONTAINS

 ! *********************************************************************
 ! subroutine: main subroutine for per-reach solute routing
 ! *********************************************************************
 SUBROUTINE constituent_rch(segIndex,         & ! input: index of reach to be processed
                            idxRoute,         & ! input: routing method index
                            NETOPO_in,        & ! input: reach topology data structure
                            RPARAM_in,        & ! input: reach parameter data structure
                            RCHSTA_out,       & ! inout: reach state data structure
                            RCHFLX_out,       & ! inout: reach flux data structure
                            ierr, message)      ! output: error control

 implicit none
 ! Argument variables
 integer(i4b),     intent(in)                 :: segIndex          ! segment where routing is performed
 integer(i4b),     intent(in)                 :: idxRoute          ! routing method index
 type(RCHTOPO),    intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
 type(RCHPRP),     intent(inout), allocatable :: RPARAM_in(:)      ! River reach parameter
 type(STRSTA),     intent(inout)              :: RCHSTA_out(:)     ! reach state data
 type(STRFLX),     intent(inout)              :: RCHFLX_out(:)     ! Reach fluxes for decomposed domains
 integer(i4b),     intent(out)                :: ierr              ! error code
 character(*),     intent(out)                :: message           ! error message
 ! Local variables
 logical(lgt)                                 :: verbose           ! check details of variables
 logical(lgt)                                 :: isHW              ! headwater basin?
 integer(i4b)                                 :: nUps              ! number of upstream segment
 integer(i4b)                                 :: iUps              ! upstream reach index
 integer(i4b)                                 :: iRch_ups          ! index of upstream reach in NETOPO
 real(dp)                                     :: Clat              ! lateral flow into channel [m3/s]
 real(dp)                                     :: Cupstream         ! total discharge at top of the reach [m3/s]
 character(len=strLen)                        :: cmessage          ! error message from subroutine

 ierr=0; message='constituent_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHID == desireId) verbose = .true.

 ! get mass flux from upstream
 nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
 isHW = .true.
 Cupstream = 0.0_dp  ! mass flux from upstream mg/s

 if (nUps>0) then ! this hru is not headwater
   isHW = .false.
   do iUps = 1,nUps
     if (.not. NETOPO_in(segIndex)%goodBas(iUps)) cycle ! skip upstream reach which does not any flow due to zero total contributory areas
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     Cupstream = Cupstream + RCHFLX_out(iRch_ups)%ROUTE(idxRoute)%reach_solute_flux   ! mass flux from upstream [mg/s]
   end do
   Clat = RCHFLX_out(segIndex)%basin_solute ! lateral mass flux [mg/s]
 else ! headwater
   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif
   if (hw_drain_point==top_reach) then ! lateral flow is poured in a reach at the top
     Cupstream = Cupstream + RCHFLX_out(segIndex)%basin_solute
     Clat = 0._dp
   else if (hw_drain_point==bottom_reach) then ! lateral flow is poured in a reach at the top
     Clat = RCHFLX_out(segIndex)%basin_solute
   end if
 end if

 if(verbose)then
   write(iulog,'(2A)') new_line('a'), '** CHECK tracer **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,1X,I12,1X,G15.4)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps), &
             RCHFLX_out(iRch_ups)%ROUTE(idxRoute)%reach_solute_flux
     enddo
   end if
   write(iulog,'(A,1X,G15.4)') ' RCHFLX_out(segIndex)%basin_solute=',RCHFLX_out(segIndex)%basin_solute
 endif

 ! compute mass flux out and average concentration
 call comp_mass_flux(idxRoute,                   &  ! input: routing method index
                     RPARAM_in(segIndex),        &  ! input: parameter at segIndex reach
                     Cupstream,                  &  ! input: total discharge at top of the reach being processed
                     Clat,                       &  ! input: lateral mass flow [mg/s]
                     isHW,                       &  ! input: is this headwater basin?
                     RCHFLX_out(segIndex),       &  ! inout: updated fluxes at reach
                     verbose,                    &  ! input: reach index to be examined
                     ierr, cmessage)                ! output: error control
 if(ierr/=0)then
    write(message, '(A,1X,I12,1X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage)
    return
 endif

 if(verbose)then
   write(iulog,'(A,1X,G15.4)') ' RCHFLX_out(segIndex)%reach_solute_flux=', RCHFLX_out(segIndex)%ROUTE(idxRoute)%reach_solute_flux
 endif

 if (RCHFLX_out(segIndex)%ROUTE(idxRoute)%reach_solute_flux < 0) then
   write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE MASS = ', RCHFLX_out(segIndex)%ROUTE(idxRoute)%reach_solute_flux, &
         'at ', NETOPO_in(segIndex)%REACHID
 end if

 call comp_reach_mb(NETOPO_in(segIndex)%REACHID, idxRoute, Cupstream, Clat, RCHFLX_out(segIndex), verbose, lakeFlag=.false.)

 END SUBROUTINE constituent_rch

 ! *********************************************************************
 ! subroutine: solve diffuisve wave equation
 ! *********************************************************************
 SUBROUTINE comp_mass_flux(idxRoute,       & ! input: routing method index
                           rch_param,      & ! input: river parameter data structure
                           Cupstream,      & ! input: discharge from upstream
                           Clat,           & ! input: lateral mass flux into chaneel [mg/s]
                           isHW,           & ! input: is this headwater basin?
                           rflux,          & ! inout: reach flux at a reach
                           verbose,        & ! input: reach index to be examined
                           ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Description:
 ! Compute a solute flux out of a reach and update solute mass in a reach
 ! ----------------------------------------------------------------------------------------
 implicit none
 ! Argument variables
 integer(i4b), intent(in)        :: idxRoute       ! routing method index
 type(RCHPRP), intent(in)        :: rch_param      ! River reach parameter
 real(dp),     intent(in)        :: Cupstream      ! total discharge at top of the reach being processed
 real(dp),     intent(in)        :: Clat           ! lateral discharge into chaneel [m3/s]
 logical(lgt), intent(in)        :: isHW           ! is this headwater basin?
 type(STRFLX), intent(inout)     :: rflux          ! current reach fluxes
 logical(lgt), intent(in)        :: verbose        ! reach index to be examined
 integer(i4b), intent(out)       :: ierr           ! error code
 character(*), intent(out)       :: message        ! error message
 ! Local variables
 real(dp)                        :: max_outMass    ! maximum possible constituent discharge [mg/s]
 real(dp)                        :: reach_vol      ! water volume [m3]
 real(dp)                        :: reach_mass     ! water volume [m3]
 real(dp)                        :: solute_out     ! mass flux out [mg/s]
 real(dp)                        :: solute_per_vol ! mass concentration in channel [mg/m3]
 character(len=strLen)           :: cmessage       ! error message from subroutine

 ierr=0; message='comp_mass_flux/'

 rflux%ROUTE(idxRoute)%reach_solute_mass(0) = rflux%ROUTE(idxRoute)%reach_solute_mass(1)

 if (.not. isHW .or. hw_drain_point==top_reach) then

   ! mass flux mg/s = discharge m3/s * concentration mg/m3
   reach_mass = Cupstream*dt + rflux%ROUTE(idxRoute)%reach_solute_mass(0)
   reach_vol  = rflux%ROUTE(idxRoute)%REACH_INFLOW*dt + rflux%ROUTE(idxRoute)%REACH_VOL(0)
   solute_per_vol = reach_mass/reach_vol

   solute_out = (rflux%ROUTE(idxRoute)%REACH_Q-rflux%BASIN_QR(1))*solute_per_vol

   max_outMass=rflux%ROUTE(idxRoute)%reach_solute_mass(1)/dt + Cupstream
   if (solute_out>max_outMass) then
     solute_out = max_outMass
     rflux%ROUTE(idxRoute)%reach_solute_mass(1) = 0
   else
     rflux%ROUTE(idxRoute)%reach_solute_mass(1) = rflux%ROUTE(idxRoute)%reach_solute_mass(1) + (Cupstream - solute_out)*dt
   end if

   rflux%ROUTE(idxRoute)%reach_solute_flux = solute_out + Clat

 else ! if head-water and pour runnof to the bottom of reach

   rflux%ROUTE(idxRoute)%reach_solute_flux = Clat
   rflux%ROUTE(idxRoute)%reach_solute_mass(1) = 0._dp

 endif

 END SUBROUTINE comp_mass_flux

END MODULE tracer_module
