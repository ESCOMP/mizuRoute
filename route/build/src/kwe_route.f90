MODULE kwe_route_module

! WORKING in PROGRESS

!numeric type
USE nrtype
! data types
USE dataTypes, ONLY: STRFLX            ! fluxes in each reach
USE dataTypes, ONLY: STRSTA            ! state in each reach
USE dataTypes, ONLY: RCHTOPO           ! Network topology
USE dataTypes, ONLY: RCHPRP            ! Reach parameter
USE dataTypes, ONLY: EKWRCH            ! Eulerian KW state
! global data
USE public_var, ONLY: iulog             ! i/o logical unit number
USE public_var, ONLY: verySmall         ! a very small value
USE public_var, ONLY: runoffMin         ! minimum runoff
USE public_var, ONLY: realMissing       ! missing value for real number
USE public_var, ONLY: integerMissing    ! missing value for integer number
! subroutines: general
USE perf_mod,  ONLY: t_startf,t_stopf   ! timing start/stop

! privary
implicit none
private

public::kwe_route

contains

 ! *********************************************************************
 ! subroutine: route kinematic waves with Euler solution through the river network
 ! *********************************************************************
 SUBROUTINE kwe_route(iens,                 & ! input: ensemble index
                      river_basin,          & ! input: river basin information (mainstem, tributary outlet etc.)
                      T0,T1,                & ! input: start and end of the time step
                      ixDesire,             & ! input: reachID to be checked by on-screen pringing
                      NETOPO_in,            & ! input: reach topology data structure
                      RPARAM_in,            & ! input: reach parameter data structure
                      RCHSTA_out,           & ! inout: reach state data structure
                      RCHFLX_out,           & ! inout: reach flux data structure
                      ierr,message,         & ! output: error control
                      ixSubRch)               ! optional input: subset of reach indices to be processed

   USE dataTypes, ONLY: subbasin_omp          ! mainstem+tributary data strucuture
   USE model_utils, ONLY: handle_err

   implicit none
   ! Input
   integer(i4b),       intent(in)                 :: iEns                 ! ensemble member
   type(subbasin_omp), intent(in),    allocatable :: river_basin(:)       ! river basin information (mainstem, tributary outlet etc.)
   real(dp),           intent(in)                 :: T0,T1                ! start and end of the time step (seconds)
   integer(i4b),       intent(in)                 :: ixDesire             ! index of the reach for verbose output
   type(RCHTOPO),      intent(in),    allocatable :: NETOPO_in(:)         ! River Network topology
   type(RCHPRP),       intent(in),    allocatable :: RPARAM_in(:)         ! River reach parameter
   ! inout
   type(STRSTA),       intent(inout), allocatable :: RCHSTA_out(:,:)      ! reach state data
   type(STRFLX),       intent(inout), allocatable :: RCHFLX_out(:,:)      ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   ! output variables
   integer(i4b),       intent(out)                :: ierr                 ! error code
   character(*),       intent(out)                :: message              ! error message
   ! input (optional)
   integer(i4b),       intent(in), optional       :: ixSubRch(:)          ! subset of reach indices to be processed
   ! local variables
   character(len=strLen)                          :: cmessage             ! error message for downwind routine
   logical(lgt),                      allocatable :: doRoute(:)           ! logical to indicate which reaches are processed
   integer(i4b)                                   :: LAKEFLAG=0           ! >0 if processing lakes
   integer(i4b)                                   :: nOrder               ! number of stream order
   integer(i4b)                                   :: nTrib                ! number of tributary basins
   integer(i4b)                                   :: nSeg                 ! number of reaches in the network
   integer(i4b)                                   :: iSeg, jSeg           ! loop indices - reach
   integer(i4b)                                   :: iTrib                ! loop indices - branch
   integer(i4b)                                   :: ix                   ! loop indices stream order

   ! initialize error control
   ierr=0; message='kwe_route/'

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

   call t_startf('route/kwe')

   ! route kinematic waves through the river network
   do ix = 1, nOrder

     nTrib=size(river_basin(ix)%branch)

!$OMP PARALLEL DO schedule(dynamic,1)                   & ! chunk size of 1
!$OMP          private(jSeg, iSeg)                      & ! private for a given thread
!$OMP          private(ierr, cmessage)                  & ! private for a given thread
!$OMP          shared(T0,T1)                            & ! private for a given thread
!$OMP          shared(LAKEFLAG)                         & ! private for a given thread
!$OMP          shared(river_basin)                      & ! data structure shared
!$OMP          shared(doRoute)                          & ! data array shared
!$OMP          shared(NETOPO_in)                        & ! data structure shared
!$OMP          shared(RPARAM_in)                        & ! data structure shared
!$OMP          shared(RCHSTA_out)                       & ! data structure shared
!$OMP          shared(RCHFLX_out)                       & ! data structure shared
!$OMP          shared(ix, iEns, ixDesire)               & ! indices shared
!$OMP          firstprivate(nTrib)
     trib:do iTrib = 1,nTrib
       seg:do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
         jSeg  = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
         if (.not. doRoute(jSeg)) cycle
         call kwe_rch(iEns,jSeg,           & ! input: array indices
                      ixDesire,            & ! input: index of the desired reach
                      T0,T1,               & ! input: start and end of the time step
                      LAKEFLAG,            & ! input: flag if lakes are to be processed
                      NETOPO_in,           & ! input: reach topology data structure
                      RPARAM_in,           & ! input: reach parameter data structure
                      RCHSTA_out,          & ! inout: reach state data structure
                      RCHFLX_out,          & ! inout: reach flux data structure
                      ierr,cmessage)         ! output: error control
         if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
       end do  seg
     end do trib
!$OMP END PARALLEL DO

   end do

   call t_stopf('route/kwe')

 END SUBROUTINE kwe_route

 ! *********************************************************************
 ! subroutine: perform one segment route KW routing
 ! *********************************************************************
 SUBROUTINE kwe_rch(iEns, segIndex, & ! input: index of runoff ensemble to be processed
                    ixDesire,       & ! input: reachID to be checked by on-screen pringing
                    T0,T1,          & ! input: start and end of the time step
                    LAKEFLAG,       & ! input: flag if lakes are to be processed
                    NETOPO_in,      & ! input: reach topology data structure
                    RPARAM_in,      & ! input: reach parameter data structure
                    RCHSTA_out,     & ! inout: reach state data structure
                    RCHFLX_out,     & ! inout: reach flux data structure
                    ierr, message)    ! output: error control

 implicit none

 ! Input
 integer(i4b),  intent(in)                 :: iEns              ! runoff ensemble to be routed
 integer(i4b),  intent(in)                 :: segIndex          ! segment where routing is performed
 integer(i4b),  intent(in)                 :: ixDesire          ! index of the reach for verbose output
 real(dp),      intent(in)                 :: T0,T1             ! start and end of the time step (seconds)
 integer(i4b),  intent(in)                 :: LAKEFLAG          ! >0 if processing lakes
 type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
 type(RCHPRP),  intent(in),    allocatable :: RPARAM_in(:)      ! River reach parameter
 ! inout
 type(STRSTA),  intent(inout), allocatable :: RCHSTA_out(:,:)   ! reach state data
 type(STRFLX),  intent(inout), allocatable :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! Output
 integer(i4b),  intent(out)                :: ierr              ! error code
 character(*),  intent(out)                :: message           ! error message
 ! Local variables to
 type(STRFLX), allocatable                 :: uprflux(:)        ! upstream Reach fluxes
 logical(lgt)                              :: doCheck           ! check details of variables
 integer(i4b)                              :: nUps              ! number of upstream segment
 integer(i4b)                              :: iUps              ! upstream reach index
 integer(i4b)                              :: iRch_ups          ! index of upstream reach in NETOPO
 character(len=strLen)                     :: cmessage          ! error message from subroutine

 ! initialize error control
 ierr=0; message='kwe_rch/'

 doCheck = .false.
 if(NETOPO_in(segIndex)%REACHIX == ixDesire)then
   doCheck = .true.
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

 if(doCheck)then
   write(iulog,'(A)') 'CHECK Euler Kinematic wave'
   if (nUps>0) then
     do iUps = 1,nUps
       write(iulog,'(A,X,I6,X,G12.5)') ' UREACHK, uprflux=',NETOPO_in(segIndex)%UREACHK(iUps),uprflux(iUps)%REACH_Q
     enddo
   end if
   write(iulog,'(A,X,G12.5)') ' RCHFLX_out(iEns,segIndex)%BASIN_QR(1)=',RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
 endif

 ! perform river network UH routing
 call kw_euler(RPARAM_in(segIndex),         &    ! input: parameter at segIndex reach
               T0,T1,                       &    ! input: start and end of the time step
               uprflux,                     &    ! input: upstream reach fluxes
               RCHSTA_out(iens,segIndex),   &    ! inout:
               RCHFLX_out(iens,segIndex),   &    ! inout: updated fluxes at reach
               doCheck,                    &    ! input: reach index to be examined
               ierr, message)                    ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! Check True since now this reach now routed
 RCHFLX_out(iEns,segIndex)%isRoute=.True.

 if(doCheck)then
  write(iulog,'(A,X,G12.5)') ' RCHFLX_out(iens,segIndex)%REACH_Q=', RCHFLX_out(iens,segIndex)%REACH_Q
 endif

 END SUBROUTINE kwe_rch


 ! *********************************************************************
 ! subroutine: route kinematic waves at one segment
 ! *********************************************************************
 SUBROUTINE kw_euler(rch_param,     & ! input: river parameter data structure
                     T0,T1,         & ! input: start and end of the time step
                     uprflux,       & ! input: upstream reach fluxes
                     rstate,        & ! inout: reach state at a reach
                     rflux,         & ! inout: reach flux at a reach
                     doCheck,       & ! input: reach index to be examined
                     ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Kinematic wave equation is solved based on conservative form the equation
 ! Reference: HEC-HMS technical reference mannual
 ! https://www.hec.usace.army.mil/software/hec-hms/documentation/HEC-HMS_Technical%20Reference%20Manual_(CPD-74B).pdf
 !
 ! Regarding state array:
 ! (time:0:1, loc:0:1) 0-previous time step/inlet, 1-current time step/outlet.
 ! Q or A(1,2,3,4): 1: (t=0,x=0), 2: (t=0,x=1), 3: (t=1,x=0), 4: (t=1,x=1)
 !
 ! TO-DO: 1. implement adaptive time step Euler

 IMPLICIT NONE
 ! Input
 type(RCHPRP), intent(in)                 :: rch_param    ! River reach parameter
 real(dp),     intent(in)                 :: T0,T1        ! start and end of the time step (seconds)
 type(STRFLX), intent(in), allocatable    :: uprflux(:)   ! upstream Reach fluxes
 logical(lgt), intent(in)                 :: doCheck      ! reach index to be examined
 ! Input/Output
 type(STRSTA), intent(inout)              :: rstate       ! curent reach states
 type(STRFLX), intent(inout)              :: rflux        ! current Reach fluxes
 ! Output
 integer(i4b), intent(out)                :: ierr         ! error code
 character(*), intent(out)                :: message      ! error message
 ! LOCAL VAIRABLES
 real(dP)                                 :: alpha        ! sqrt(slope)(/mannings N* width)
 real(dP)                                 :: beta         ! constant, 5/3
 real(dP)                                 :: alpha1       ! sqrt(slope)(/mannings N* width)
 real(dP)                                 :: beta1        ! constant, 5/3
 !real(dp)                                 :: theta        ! Courant number [-]
 real(dp)                                 :: dT           ! interval of time step [sec]
 real(dp)                                 :: dX           ! length of segment [m]
 real(dp)                                 :: A(0:1,0:1)   !
 real(dp)                                 :: Q(0:1,0:1)   !
 real(dp)                                 :: Abar         !
 real(dp)                                 :: Qbar         !
 integer(i4b)                             :: iUps,ix      ! loop index
 integer(i4b)                             :: nUps         ! number of all upstream segment
 character(len=strLen)                    :: cmessage     ! error message from subroutine

 ! initialize error control
 ierr=0; message='kw_euler/'

 nUps = size(uprflux)

 if (nUps>0) then

   !  current time and inlet  3 (1,0) -> previous time and inlet  1 (0,0)
   !  current time and outlet 4 (1,1) -> previous time and outlet 2 (0,1)
   do ix = 0,1
     A(0,ix) = rstate%EKW_ROUTE%A(ix+3)
     Q(0,ix) = rstate%EKW_ROUTE%Q(ix+3)
   end do

   A(1,1) = realMissing
   Q(1,1) = realMissing

   ! Get the reach parameters
   ! A = (Q/alpha)**(1/beta)
   ! Q = alpha*A**beta
   alpha = sqrt(rch_param%R_SLOPE)/(rch_param%R_MAN_N*rch_param%R_WIDTH**(2._dp/3._dp))
   beta  = 5._dp/3._dp
   beta1  = 1._dp/beta
   alpha1 = (1.0/alpha)**beta1
   dX = rch_param%RLENGTH
   dT = T1-T0

   ! compute flow rate and flow area at upstream end at current time step
   Q(1,0) = 0.0_dp
   do iUps = 1,nUps
     Q(1,0) = Q(1,0) + uprflux(iUps)%REACH_Q
   end do

   A(1,0) = (Q(1,0)/alpha)**(1/beta)

   if (doCheck) then
     write(iulog,'(3(A,X,G12.5))') ' R_SLOPE=',rch_param%R_SLOPE,' R_WIDTH=',rch_param%R_WIDTH,' R_MANN=',rch_param%R_MAN_N
     write(iulog,'(3(A,X,G12.5))') ' Q(0,0)=',Q(0,0),' Q(0,1)=',Q(0,1),' Q(1,0)=',Q(1,0)
     write(iulog,'(3(A,X,G12.5))') ' A(0,0)=',A(0,0),' A(0,1)=',A(0,1),' A(1,0)=',A(1,0)
   end if

! ----- Need adaptive time step forward Euler --------
!   ! Compute Courant number == ration of celarity to dx/dt
!   Qbar = 0.5*(Q(0,1)+Q(1,0))
!   Abar = 0.5*(A(0,1)+A(1,0))
!   if (Abar <= 0._dp) then
!     Qbar = runoffMin
!     Abar = (Qbar/alpha)**(1/beta)
!   endif
!   theta = beta*(dT/dX)*(Qbar)/(Abar)
!
!   if (doCheck) then
!     write(iulog,'(5(a,G12.5,x))') ' dT= ', dt, 'dX=', dX, 'Qbar= ', Qbar, 'Abar= ', ABar, 'theta= ', theta
!   end if

   ! ----------
   ! solve flow rate and flow area at downstream end at current time step
   ! ----------

! ----- Modify for adaptive time step forward Euler --------
!   if (theta >= 1.0) then
!     Q(1,1) = Q(1,0) - dX/dT*(A(1,0)-A(0,0))
!     A(1,1) = (Q(1,1)/alpha)**(1/beta)
!   else
!     A(1,1) = A(0,1) + dT/dX*(Q(0,0) - Q(0,1))
!     Q(1,1) = alpha*A(1,1)**beta
!   end if

   Qbar = (Q(0,1) + Q(1,0))/2
   Q(1,1) = dT/dX*Q(1,0) + alpha1*beta1*Qbar**(beta1-1)*Q(0,1)
   Q(1,1) = Q(1,1)/(dT/dX + alpha1*beta1*Qbar**(beta1-1))

   if (doCheck) then
     write(iulog,'(2(A,X,G12.5))') ' A(1,1)=',A(1,1),' Q(1,1)=',Q(1,1)
   end if

 else ! if head-water

   !  current time and inlet  3 (1,0) -> previous time and inlet (0,0)
   !  current time and outlet 4 (1,1) -> previous time and outlet (0,1)
   do ix = 0,1
     A(0,ix) = rstate%EKW_ROUTE%A(ix+3)
     Q(0,ix) = rstate%EKW_ROUTE%Q(ix+3)
   end do

   A(1,0) = 0._dp
   Q(1,0) = 0._dp

   Q(1,1) = 0._dp
   A(1,1) = 0._dp

   if (doCheck) then
     write(iulog,'(A)')            ' This is headwater '
     write(iulog,'(2(A,X,G12.5))') ' A(1,1)=',A(1,1),' Q(1,1)=',Q(1,1)
   endif

 endif

 ! add catchment flow
 rflux%REACH_Q = Q(1,1)+rflux%BASIN_QR(1)

 ! update state
 do ix = 0,1
   rstate%EKW_ROUTE%Q(ix+1) = Q(0,ix)
   rstate%EKW_ROUTE%Q(ix+3) = Q(1,ix)

   rstate%EKW_ROUTE%A(ix+1) = A(0,ix)
   rstate%EKW_ROUTE%A(ix+3) = A(1,ix)
 end do

 END SUBROUTINE kw_euler


END MODULE kwe_route_module
