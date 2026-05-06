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
USE globalData,    ONLY: nTracer
USE globalData,    ONLY: nMolecule             ! # of reach internal nodes for finite difference (including upstream and downstream boundaries)
USE water_balance, ONLY: comp_reach_mb         ! compute water balance error
USE hydraulic,     ONLY: flow_depth
USE hydraulic,     ONLY: flow_area

implicit none

private
public::constituent_rch

integer(i4b), parameter :: top_reach=1
integer(i4b), parameter :: bottom_reach=2

logical(lgt) :: instant_mix=.false.

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
 integer(i4b)                                 :: iTrace            ! loop index
 integer(i4b)                                 :: iRch_ups          ! index of upstream reach in NETOPO
 real(dp)                                     :: Clat(nTracer)     ! lateral flow into channel [m3/s]
 real(dp)                                     :: Cupstream(nTracer)! total discharge at top of the reach [m3/s]
 character(len=strLen)                        :: cmessage          ! error message from subroutine

 ierr=0; message='constituent_rch/'

 verbose = .false.
 if(NETOPO_in(segIndex)%REACHID == desireId) verbose = .true.

 ! get mass flux from upstream
 nUps = count(NETOPO_in(segIndex)%goodBas) ! reminder: goodBas is reach with >0 total contributory area
 isHW = .true.
 Cupstream = 0._dp
 Clat = 0._dp

 if (nUps>0) then ! this hru is not headwater
   isHW = .false.
   do iUps = 1,nUps
     if (.not. NETOPO_in(segIndex)%goodBas(iUps)) cycle ! skip upstream reach which does not any flow due to zero total contributory areas
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     Cupstream(1:nTracer) = Cupstream(1:nTracer) + &
                            RCHFLX_out(iRch_ups)%ROUTE(idxRoute)%reach_solute_flux(1:nTracer) ! mass flux from upstream [mg/s]
   end do
   Clat(1:nTracer) = RCHFLX_out(segIndex)%basin_solute(1:nTracer) ! lateral mass flux [mg/s]
 else ! headwater
   if (verbose) then
     write(iulog,'(A)')            ' This is headwater '
   endif
   if (hw_drain_point==top_reach) then ! lateral flow is poured in a reach at the top
     Cupstream(1:nTracer) = Cupstream(1:nTracer) + RCHFLX_out(segIndex)%basin_solute(1:nTracer)
   else if (hw_drain_point==bottom_reach) then ! lateral flow is poured in a reach at the top
     Clat(1:nTracer) = RCHFLX_out(segIndex)%basin_solute(1:nTracer)
   end if
 end if

 if(verbose)then
   write(iulog,'(2A)') new_line('a'), '** CHECK tracer **'
   if (nUps>0) then
     do iUps = 1,nUps
       iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
       write(iulog,'(A,1X,I12,1X,*(G15.4,:,": "))') ' UREACHK, uprflux=', NETOPO_in(segIndex)%UREACHK(iUps), &
             (RCHFLX_out(iRch_ups)%ROUTE(idxRoute)%reach_solute_flux(iTrace), iTrace=1,nTracer)
     enddo
   end if
   write(iulog,'(A,1X,*(G15.4,:,": "))') ' basin_solute=', (RCHFLX_out(segIndex)%basin_solute(iTrace), iTrace=1,nTracer)
 endif

 ! compute mass flux out and average concentration
 call comp_mass_flux(idxRoute,                   &  ! input: routing method index
                     RPARAM_in(segIndex),        &  ! input: parameter at segIndex reach
                     Cupstream,                  &  ! input: total discharge at top of the reach being processed
                     Clat,                       &  ! input: lateral mass flow [mg/s]
                     isHW,                       &  ! input: is this headwater basin?
                     RCHSTA_out(segIndex),       &  ! inout: updated reach states
                     RCHFLX_out(segIndex),       &  ! inout: updated fluxes at reach
                     verbose,                    &  ! input: reach index to be examined
                     ierr, cmessage)                ! output: error control
 if(ierr/=0)then
    write(message, '(A,1X,I12,1X,A)') trim(message)//'/segment=', NETOPO_in(segIndex)%REACHID, '/'//trim(cmessage)
    return
 endif

 if(verbose)then
   write(iulog,'(A,1X,*(G15.4,:,": "))') ' reach_solute_flux=', &
                                 (RCHFLX_out(segIndex)%ROUTE(idxRoute)%reach_solute_flux(iTrace), iTrace=1,nTracer)
 endif

 do iTrace=1,nTracer
   if (RCHFLX_out(segIndex)%ROUTE(idxRoute)%reach_solute_flux(iTrace) < 0) then
     write(iulog,'(A,1X,G12.5,1X,A,1X,I9)') ' ---- NEGATIVE MASS = ', &
           RCHFLX_out(segIndex)%ROUTE(idxRoute)%reach_solute_flux(iTrace), &
           'at ', NETOPO_in(segIndex)%REACHID
   end if
 end do

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
                           rstate,         & ! inout: reach states at a reach
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
 real(dp),     intent(in)        :: Cupstream(:)   ! total mass fluxes at top of the reach being processed [mg/s]
 real(dp),     intent(in)        :: Clat(:)        ! lateral mass fluxes into chanel [mg/s]
 logical(lgt), intent(in)        :: isHW           ! is this headwater basin?
 type(STRSTA), intent(inout)     :: rstate         ! current reach state
 type(STRFLX), intent(inout)     :: rflux          ! current reach fluxes
 logical(lgt), intent(in)        :: verbose        ! reach index to be examined
 integer(i4b), intent(out)       :: ierr           ! error code
 character(*), intent(out)       :: message        ! error message
 ! Local variables
 integer(i4b)                    :: iTrace         ! loop index
 real(dp),allocatable            :: Cinflux(:)     ! concentration influx [mg/m3/s]
 real(dp)                        :: max_outMass    ! maximum possible constituent discharge [mg/s]
 real(dp)                        :: reach_vol      ! water volume [m3]
 real(dp)                        :: reach_mass     ! water volume [m3]
 real(dp)                        :: solute_out     ! mass flux out [mg/s]
 real(dp)                        :: solute_per_vol ! mass concentration in channel [mg/m3]
 character(len=strLen)           :: cmessage       ! error message from subroutine

 ierr=0; message='comp_mass_flux/'

 rflux%ROUTE(idxRoute)%reach_solute_mass(0,1:nTracer) = rflux%ROUTE(idxRoute)%reach_solute_mass(1,1:nTracer)

 if (.not. isHW .or. hw_drain_point==top_reach) then

   allocate(Cinflux(ntracer))
   Cinflux=0._dp
   if (rflux%ROUTE(idxRoute)%REACH_INFLOW>0) then
      Cinflux(:) = Cupstream(:)/rflux%ROUTE(idxRoute)%REACH_INFLOW
   end if

   do iTrace=1,nTracer
     if (instant_mix) then ! To be removed
       ! mass flux mg/s = discharge m3/s * concentration mg/m3
       reach_mass = Cupstream(iTrace)*dt + rflux%ROUTE(idxRoute)%reach_solute_mass(1,iTrace)
       reach_vol  = rflux%ROUTE(idxRoute)%REACH_INFLOW*dt + rflux%ROUTE(idxRoute)%REACH_VOL(0)
       solute_per_vol=0._dp
       if (reach_vol>0._dp) then
         solute_per_vol = reach_mass/reach_vol
       end if
     else
       call comp_advec_diffusion(iTrace, &
                                 rch_param, &
                                 Cinflux, &
                                 Clat, &
                                 isHW, &
                                 rstate%DW_ROUTE%molecule, &
                                 verbose,   &
                                 ierr,cmessage)
       if(ierr/=0)then
          write(message,'(A,1X,I12,1X,A)') trim(message)//trim(cmessage);return
       endif
       solute_per_vol = rstate%DW_ROUTE%molecule%c_solute(nMolecule%DW_ROUTE-1,iTrace)
     end if
     solute_out = (rflux%ROUTE(idxRoute)%REACH_Q-rflux%BASIN_QR(1))*solute_per_vol

     ! limit maximum allowable mass flux out of the reach
     max_outMass=rflux%ROUTE(idxRoute)%reach_solute_mass(1,iTrace)/dt + Cupstream(iTrace)
     if (solute_out>max_outMass) then
       solute_out = max_outMass
       rflux%ROUTE(idxRoute)%reach_solute_mass(1,iTrace) = 0
     else
       rflux%ROUTE(idxRoute)%reach_solute_mass(1,iTrace) = rflux%ROUTE(idxRoute)%reach_solute_mass(1, iTrace) &
                                                         + (Cupstream(iTrace) - solute_out)*dt
     end if
     ! Update total mass flux from outlet of reach [mg/s]
     rflux%ROUTE(idxRoute)%reach_solute_flux(iTrace) = solute_out + Clat(iTrace)
   end do
 else ! if head-water and pour runnof to the bottom of reach
   rflux%ROUTE(idxRoute)%reach_solute_flux(:) = Clat(:)
   rflux%ROUTE(idxRoute)%reach_solute_mass(1,:) = 0._dp
 endif

 END SUBROUTINE comp_mass_flux

  ! *********************************************************************
  ! private subroutine: solve advection-diffusion equation within a reach
  ! *********************************************************************
  SUBROUTINE comp_advec_diffusion(iTrace,         & ! input: tracer index
                                  rch_param,      & ! input: river parameter data structure
                                  Cinflux,        & ! input: concentration influx [mg/m3/s]
                                  Clat,           & ! input: lateral discharge into chaneel [m3/s]
                                  isHW,           & ! input: is this headwater basin?
                                  subrch_state,   & ! inout: reach states at a reach
                                  verbose,        & ! input: reach index to be examined
                                  ierr,message)
  ! ----------------------------------------------------------------------------------------
  ! Description:
  !  Compute In-stream process - advection, dispersion, and reaction
  ! Note:
  !  Currently reaction is not implemented
  ! ----------------------------------------------------------------------------------------
  USE advection_diffusion, ONLY: solve_ade
  implicit none
  ! Argument variables
  integer(i4b), intent(in)        :: iTrace         ! tracer index
  type(RCHPRP), intent(in)        :: rch_param      ! River reach parameter
  real(dp),     intent(in)        :: Cinflux(:)     ! concentration influx [mg/m3/s]
  real(dp),     intent(in)        :: Clat(:)        ! lateral discharge into chaneel [m3/s]
  logical(lgt), intent(in)        :: isHW           ! is this headwater basin?
  type(SUBRCH), intent(inout)     :: subrch_state   ! curent states at sub-reaches
  logical(lgt), intent(in)        :: verbose        ! reach index to be examined
  integer(i4b), intent(out)       :: ierr           ! error code
  character(*), intent(out)       :: message        ! error message
  ! Local variables
  real(dp)                        :: Ybar           ! mean reach flow depth over sub-reach points [m]
  real(dp)                        :: Qbar           ! mean reach discharge over sub-reach points [m3/s]
  real(dp)                        :: Abar           ! mean reach flow area over sub-reach points [m2]
  real(dp)                        :: Vbar           ! cross-sectional mean water velocity [m/s]
  real(dp)                        :: dk             ! diffusivity [m2/s]
  real(dp), allocatable           :: Clocal(:,:)    ! sub-reach & sub-time step discharge at previous and current time step [m3/s]
  real(dp), allocatable           :: Cprev(:)       ! sub-reach discharge at previous time step [m3/s]
  real(dp)                        :: dTsub          ! time inteval for sub time-step [sec]
  integer(i4b)                    :: it             ! loop index
  integer(i4b)                    :: ntSub          ! number of sub time-step
  character(len=strLen)           :: cmessage       ! error message from subroutine

  ierr=0; message='comp_advec_diffusion/'

  ntSub = 1    ! number of sub-time step
  dk = 0._dp   ! dk=0 -> no diffusion, only advection

  associate(nMolecule  => nMolecule%DW_ROUTE,   &
            R_SLOPE    => rch_param%R_SLOPE,    & ! channel slope
            R_MAN_N    => rch_param%R_MAN_N,    & ! manning n
            R_WIDTH    => rch_param%R_WIDTH,    & ! channel bottom width
            R_DEPTH    => rch_param%R_DEPTH,    & ! bankfull depth
            SIDE_SLOPE => rch_param%SIDE_SLOPE, & ! channel side slope
            FLDP_SLOPE => rch_param%FLDP_SLOPE, & ! floodplain slope
            RLENGTH    => rch_param%RLENGTH)      ! channel length

  if (.not. isHW .or. hw_drain_point==top_reach) then

    allocate(Cprev(nMolecule), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! initialize previous time step flow
    Cprev = subrch_state%c_solute(:,iTrace)  ! soltute mass state at previous time step

    if (verbose) then
      write(iulog,'(A,1X,G12.5)') ' length [m] =',RLENGTH
    end if

    ! time-step adjustment so Courant number is less than 1
    dTsub = dt/ntSub

    if (verbose) then
      write(iulog,'(A,1X,I3,A,1X,G12.5)') ' No. sub timestep=',nTsub,' sub time-step [sec]=',dTsub
      write(iulog,'(A,1X,G12.5)') ' Mass concentration=',Cinflux(iTrace)
    end if

    allocate(Clocal(1:nMolecule, 0:1), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! use constant flow velocity
    Qbar = sum(subrch_state%Q)/nMolecule
    Vbar=0.0_dp
    if (abs(Qbar)>1.e-50_dp) then
      Ybar = flow_depth(abs(Qbar), R_WIDTH, SIDE_SLOPE, R_SLOPE, R_MAN_N, zf=FLDP_SLOPE, bankDepth=R_DEPTH)
      Abar = flow_area(abs(Ybar), R_WIDTH, SIDE_SLOPE, zf=FLDP_SLOPE, bankDepth=R_DEPTH)
      Vbar =  Qbar/Abar
    end if

    do it = 1, nTsub
      call solve_ade(RLENGTH,            & ! input: river parameter data structure
                     nMolecule,          & ! input: number of sub-segments
                     dTsub,              & ! input: time_step [sec]
                     Cinflux(iTrace),    & ! input: quantity from upstream [unit of quantity]
                     Vbar,               & ! input: reach mean flow velocity [m/s]
                     dk,                 & ! input: diffusivity [m2/s]
                     Clat(iTrace),       & ! input: lateral quantity into chaneel [unit of quantity]
                     Cprev,              & ! input: quantity at previous time step [unit of quantity]
                     Clocal,             & ! inout: quantity soloved at current time step [unit of quantity]
                     verbose,            & ! input: reach index to be examined
                     advec_scheme=1,     & ! optional input: reach index to be examined
                     downstreamBC=1)       ! optional input: reach index to be examined
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end do

    ! update state
    subrch_state%c_solute(1:nMolecule,iTrace) = Clocal(1:nMolecule,1)

  else ! if head-water and pour runnof to the bottom of reach

    subrch_state%c_solute(1:nMolecule,iTrace) = 0._dp
    if (verbose) then
      write(iulog,'(A)')            ' This is headwater '
    endif

  endif

  end associate

  END SUBROUTINE comp_advec_diffusion

END MODULE tracer_module
