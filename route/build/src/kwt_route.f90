module kwt_route

use nrtype
use nrutil, only : arth                                 ! Num. Recipies utilities

implicit none
private
public::reachorder
public::reach_list
public::qroute_rch

! common variables
integer(i4b),parameter  :: ixPrint = -9999  ! the desired reach (set to negative to avoid any printing)
real(dp),parameter      :: verySmall=tiny(1.0_dp)  ! a very small number

contains

 ! *********************************************************************
 ! subroutine: define processing order for the individual
 !                 stream segments in the river network
 ! *********************************************************************
 subroutine REACHORDER(NRCH, &           ! input
                       ierr, message)    ! error control
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   David Tarboton, 1997 (original code)
 !   Martyn Clark, 2007 (revised code, for use within TopNet)
 !   Martyn Clark, 2012 (stand-alone code)
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Defines the processing order for the individual stream segments in the river network
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O:
 !
 !   INPUTS:
 !    NRCH: Number of stream segments in the river network (reaches)
 !
 ! ----------------------------------------------------------------------------------------
 ! Structures modified:
 !
 !   Updates structure RHORDER in module reachparam
 !
 ! ----------------------------------------------------------------------------------------
 ! Source:
 !
 !   Subroutine MDDATA within TOPNET version 7
 !
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! ----------------------------------------------------------------------------------------
 USE reachparam
 IMPLICIT NONE
 ! input variables
 INTEGER(I4B), INTENT(IN)               :: NRCH            ! number of stream segments
 ! output variables
 integer(i4b), intent(out)              :: ierr            ! error code
 character(*), intent(out)              :: message         ! error message
 ! local variables
 INTEGER(I4B)                           :: IRCH,JRCH,KRCH  ! loop through reaches
 INTEGER(I4B)                           :: IUPS            ! loop through upstream reaches
 INTEGER(I4B)                           :: ICOUNT          ! counter for the gutters
 INTEGER(I4B)                           :: NASSIGN         ! # reaches currently assigned
 LOGICAL(LGT),DIMENSION(:),ALLOCATABLE  :: RCHFLAG         ! TRUE if reach is processed
 INTEGER(I4B)                           :: NUPS            ! number of upstream reaches
 INTEGER(I4B)                           :: UINDEX          ! upstream reach index
 ! initialize error control
 ierr=0; message='reachorder/'
 ! ----------------------------------------------------------------------------------------
 NASSIGN = 0
 ALLOCATE(RCHFLAG(NRCH),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for RCHFLAG'; return; endif
 RCHFLAG(1:NRCH) = .FALSE.
 ! ----------------------------------------------------------------------------------------
 ICOUNT=0
 DO  ! do until all reaches are assigned
  NASSIGN = 0
  DO IRCH=1,NRCH
   ! check if the reach is assigned yet
   IF(RCHFLAG(IRCH)) THEN
    NASSIGN = NASSIGN + 1
    CYCLE
   ENDIF
   ! climb upstream as far as possible
   JRCH = IRCH    ! the first reach under investigation 
   DO  ! do until get to a "most upstream" reach that is not assigned
    NUPS = SIZE(NETOPO(JRCH)%UREACHI)    ! get number of upstream reaches
    IF (NUPS.GE.1) THEN     ! (if NUPS = 0, then it is a first-order basin)
     KRCH = JRCH   ! the reach under investigation 
     ! loop through upstream reaches
     DO IUPS=1,NUPS
      UINDEX = NETOPO(JRCH)%UREACHI(IUPS)  ! POSITION of the upstream reach
      ! check if the reach is NOT assigned
      IF (.NOT.RCHFLAG(UINDEX)) THEN
       JRCH = UINDEX
       EXIT    ! exit IUPS loop
      END IF  ! if the reach is assigned
     END DO  ! (looping through upstream reaches)
     ! check if all upstream reaches are already assigned (only if KRCH=JRCH)
     IF (JRCH.EQ.KRCH) THEN
      ! assign JRCH
      ICOUNT=ICOUNT+1
      RCHFLAG(JRCH) = .TRUE.
      NETOPO(ICOUNT)%RHORDER = JRCH
      EXIT
     ENDIF  
     CYCLE   ! if jrch changes, keep looping (move upstream)
    ELSE    ! if the reach is a first-order basin
     ! assign JRCH
     ICOUNT=ICOUNT+1
     RCHFLAG(JRCH) = .TRUE.
     NETOPO(ICOUNT)%RHORDER = JRCH
     EXIT
    ENDIF
   END DO   !  climbing upstream (do-forever)
  END DO   ! looping through reaches
  IF (NASSIGN.EQ.NRCH) EXIT
 END DO  ! do forever (do until all reaches are assigned)
 DEALLOCATE(RCHFLAG,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for RCHFLAG'; return; endif
 ! ----------------------------------------------------------------------------------------
 end subroutine

 ! *********************************************************************
 ! subroutine: identify all reaches above the current reach
 ! *********************************************************************
 subroutine REACH_LIST(NRCH,NTOTAL,ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Martyn Clark, 2006
 !   Einar �~Vrn Hreinsson, 2009  -- adapt linked lists and selections of reaches
 !   Martyn Clark, 2014 -- modify to be used as a stand-alone module
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Generates a list of all reaches upstream of each reach (used to compute total runoff
 !     at each point in the river network)
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O:
 !
 !   INPUTS:
 !     NRCH: Number of stream segments in the river network (reaches)
 !
 !  OUTPUTS:
 !   NTOTAL: Total number of upstream reaches for all reaches
 !     ierr: Error code
 !  message: Error message
 !
 ! ----------------------------------------------------------------------------------------
 ! Structures modified:
 !
 !   Updates structure RCHLIST in module reachparam
 !
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! ----------------------------------------------------------------------------------------
 USE reachparam
 IMPLICIT NONE
 ! input variables
 INTEGER(I4B),INTENT(IN)                  :: NRCH            ! number of stream segments
 ! output variables
 integer(i4b), intent(out)                :: NTOTAL          ! total number of upstream reaches for all reaches
 integer(i4b), intent(out)                :: ierr            ! error code
 character(*), intent(out)                :: message         ! error message
 ! local variables
 INTEGER(I4B)                             :: IRCH,JRCH,KRCH  ! loop through reaches
 ! structure in a linked list
 TYPE NODE
  INTEGER(I4B)                            :: IUPS            ! index of reach u/s of station
  TYPE(NODE), POINTER                     :: NEXT            ! next node in linked list
 END TYPE NODE
 ! a structure object in a linked list
 TYPE(NODE), POINTER                      :: URCH            ! node in linked list
 ! structure in a array pointing to a linked list of upstream reaches
 TYPE UPSTR_RCH
  INTEGER(I4B)                            :: N_URCH          ! Number of upstream reaches
  TYPE(NODE), POINTER                     :: HPOINT          ! Head pointer in linked list
 END TYPE UPSTR_RCH
 TYPE(UPSTR_RCH),DIMENSION(:),ALLOCATABLE :: INTLIST         ! list of reaches u/s of each station
 INTEGER(I4B)                             :: NUMUPS          ! number of reaches upstream
 integer(i4b),parameter                   :: strLen=256      ! length of character string
 character(len=strLen)                    :: cmessage        ! error message of downwind routine
 ! ----------------------------------------------------------------------------------------
 message='REACH_LIST/'
 ! allocate space for intlist
 ALLOCATE(INTLIST(NRCH),STAT=IERR)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for INTLIST'; return; endif
 ! initialise the reach array
 DO IRCH=1,NRCH  ! loop through selected reaches
  INTLIST(IRCH)%N_URCH = 0       ! initialize the number of upstream reaches
  NULLIFY(INTLIST(IRCH)%HPOINT)  ! set pointer to a linked list to NULL
 END DO ! (irch)
 
 ! build the linked lists for all reaches
 DO KRCH=1,NRCH
  ! ensure take streamflow from surrounding basin (a reach is upstream of itself!)
  ! but only if reach is one of the selected
  CALL ADD2LIST(KRCH,KRCH,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! start with krch and loop reach by reach down through to the outlet
  IRCH = KRCH
  DO  ! (do-forever)
   JRCH = NETOPO(IRCH)%DREACHI
   IF (JRCH.GT.0) THEN !  (check that jrch is valid, negative means irch is the outlet)
    ! jrch is downstream of krch, which means that krch is upstream of jrch
    ! *** therefore, add the krch index to the jrch list of upstream reaches ***
    CALL ADD2LIST(JRCH,KRCH,ierr,cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    ! move the index downstream
    IRCH = JRCH
   ELSE
    EXIT   ! negative = missing, which means that irch is the outlet
   ENDIF
  END DO  ! do forever (go to the outlet)
 END DO  ! assess each reach

 NTOTAL=0 ! total number of upstream reaches for all reaches
 ! extract linked lists to allocated vectors
 DO JRCH=1,NRCH
  NUMUPS = INTLIST(JRCH)%N_URCH ! should be at least 1 (because reach is upstream of itself)
  NTOTAL = NTOTAL + NUMUPS
  ALLOCATE(NETOPO(JRCH)%RCHLIST(NUMUPS),STAT=IERR)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for RCHLIST'; return; endif
  ! copy list of upstream reaches to structures (and delete the list)
  CALL MOVE_LIST(JRCH,JRCH,NUMUPS,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  print*, 'jrch, numups, NETOPO(JRCH)%RCHLIST(:) = ', jrch, numups, NETOPO(JRCH)%RCHLIST(:)
 END DO  ! jrch
 
 ! free up memory
 DEALLOCATE(INTLIST,STAT=IERR)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem deallocating space for INTLIST'; return; endif
 ! ----------------------------------------------------------------------------------------
 ! ----------------------------------------------------------------------------------------
 contains
 
  ! For a down stream reach, add an upstream reach to its list of upstream reaches
  subroutine ADD2LIST(D_RCH,U_RCH,ierr,message)
    INTEGER(I4B),INTENT(IN)          :: U_RCH    ! upstream reach index
    INTEGER(I4B),INTENT(IN)          :: D_RCH    ! downstream reach index
    integer(i4b), intent(out)        :: ierr     ! error code
    character(*), intent(out)        :: message  ! error message
    message='ADD2LIST/'
    ! increment number of upstream reaches
    INTLIST(D_RCH)%N_URCH = INTLIST(D_RCH)%N_URCH + 1
    ! allocate node to be added to the list
    ALLOCATE(URCH,STAT=IERR)
    if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for linked list'; return; endif
    ! add upstream reach index
    URCH%IUPS = U_RCH
    ! insert at the beginning of existing list
    URCH%NEXT => INTLIST(D_RCH)%HPOINT
    INTLIST(D_RCH)%HPOINT => URCH
    ! nullify pointer
    NULLIFY(URCH)
  end subroutine

  ! Copy upstream reaches from linked list to structure and delete the list
  subroutine MOVE_LIST(IRCH,JRCH,NNODES,ierr,message)
    INTEGER(I4B),INTENT(IN)          :: IRCH   ! index in destination structure
    INTEGER(I4B),INTENT(IN)          :: JRCH   ! index in the intlist array
    INTEGER(I4B),INTENT(IN)          :: NNODES ! number of nodes in linked list
    TYPE(NODE), POINTER              :: URCH   ! node in linked list
    INTEGER(I4B)                     :: KRCH   ! index in structure
    integer(i4b), intent(out)        :: ierr     ! error code
    character(*), intent(out)        :: message  ! error message
    message='MOVE_LIST/'
    ! set urch to first node in the list of upstream reaches for reach jrch
    URCH => INTLIST(JRCH)%HPOINT
    ! set irch to number of upstream reaches
    KRCH = NNODES
    ! copy information in list to relevant structure
    DO WHILE (ASSOCIATED(URCH)) ! while URCH points to something
      NETOPO(IRCH)%RCHLIST(KRCH) = URCH%IUPS
      ! let urch point to next node in list
      URCH => URCH%NEXT
      ! deallocate the node (clean up)
      DEALLOCATE(INTLIST(JRCH)%HPOINT,STAT=IERR)
      if(ierr/=0)then; ierr=20; message=trim(message)//'problem deallocating space for linked list'; return; endif
      ! let hpoint follow urch along the list
      INTLIST(JRCH)%HPOINT => URCH
      ! decrement irch by one
      KRCH = KRCH - 1
    end do
  end subroutine

 end subroutine

 ! *********************************************************************
 ! subroutine: route kinematic waves through the river network
 ! *********************************************************************
 subroutine QROUTE_RCH(IENS,JRCH,    & ! input: array indices
                       ixOutlet,     & ! input: index of the outlet reach
                       T0,T1,        & ! input: start and end of the time step
                       MAXQPAR,      & ! input: maximum number of particle in a reach 
                       LAKEFLAG,     & ! input: flag if lakes are to be processed
                       ierr,message, & ! output: error control
                       RSTEP)          ! optional input: retrospective time step offset
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Ross Woods, 1997 (original code)
 !   Martyn Clark, 2006 (current version)
 !   Martyn Clark, 2012 (modified as a standalone module)
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Route kinematic waves through the river network
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O::
 !
 !   Input(s):
 !      IENS: ensemble member
 !      JRCH: index of stream segment
 !        T0: start of the time step (seconds)
 !        T1: end of the time step (seconds)
 !  LAKEFLAG: >0 if processing lakes 
 !     RSTEP: retrospective time step offset (optional)
 !
 !   Outputs (in addition to update of data structures):
 !      ierr: error code
 !   message: error message
 !
 ! ----------------------------------------------------------------------------------------
 ! Structures Used/Modified:
 !
 !   (1) MODEL_TIME
 !   --------------
 !   Uses the vector  MODTIM%TBOUNDS(:) to add time information to basin outflows
 !
 !   (2) BASIN_FLUX
 !   --------------
 !   Uses basin outflows [BASFLX(:)%INSTN_Q]
 !
 !   (3) REACHPARAM
 !   --------------
 !   Uses the network topology to process gutters [NETOPO(:)+] and the processing sequence
 !    to identify the reach to process and the last reach in the network [RHRODER(:)]
 !
 !   (4) REACHSTATE
 !   --------------
 !   Uses the data structure KROUTE to track flow particles through the river network.
 !    KROUTE is a data structure that contains the collection of flow points (in the
 !    structure KWAVE) for each stream segment.  Each flow point has attributes QF (flow),
 !    TI (time point entered a stream segment), TR (time point exited the stream segment,
 !    or is expected to exit), and RF (logical routing flag [.TRUE. if point has exited]).
 !    Hence [KROUTE(JRCH)%KWAVE(:)QF] defines the collection of flow points in the JRCH-th
 !    reach. KROUTE must be saved for model re-starts
 !
 !   (5) REACH_FLUX
 !   --------------
 !   Contains timestep-average flow for each stream segment (computed here)
 !
 !   (6) INTERBLOCK
 !   --------------
 !   Includes an explicit interface to the sub-programs
 !
 ! ----------------------------------------------------------------------------------------
 ! Method:
 !
 !   Flow routing is performed using a lagrangian one-dimensional Kinematic routing
 !   scheme.  Flow generated from each catchment is tracked through the river network.
 !   For each river segment, we compute the time each particle is expected to exit the
 !   river segment (travel time = length/celerity).  If the "exit time" occurs before
 !   the end of the time step, the particle is propagated to the downstream segment.
 !   The "exit time" then becomes the time the particle entered the downstream segment.
 !   This process is repeated for all river segments.  River segments are processed in
 !   order from upstream segments to downstream segments, meaning that it is possible
 !   for a particle to travel through several segments in a given time step.
 !
 !   Time step averages for each time step are computed by integrating over all the
 !   flow points that exit a river segment in a given time step.  To avoid boundary
 !   problems in the interpolation, the last routed wave from the previous time step
 !   and the next wave that is expected to exit the segment are used.
 !
 ! ----------------------------------------------------------------------------------------
 ! Source:
 !
 !   Most computations were originally performed within calcts in Topnet ver7, with calls
 !   to subroutines in kinwav_v7.f
 !      
 ! ----------------------------------------------------------------------------------------
 ! Modifications to Source (mclark@ucar.edu):
 !
 !   * code is more modular
 !
 !   * added additional comments
 !
 !   * all variables are defined (IMPLICIT NONE) and described (comments)
 !
 !   * use of a new data structure (KROUTE) to hold and update the flow particles
 !
 !   * upgrade to F90 (especially structured variables and dynamic memory allocation)
 !   
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! ----------------------------------------------------------------------------------------
   USE reachparam
   USE reachstate
   USE reach_flux
   implicit none
   ! Input
   INTEGER(I4B), INTENT(IN)                    :: IENS          ! ensemble member
   integer(i4b), intent(in)                    :: ixOutlet      ! index of the outlet reach
   INTEGER(I4B), INTENT(IN)                    :: JRCH          ! reach to process
   REAL(DP), INTENT(IN)                        :: T0,T1         ! start and end of the time step (seconds)
   INTEGER(I4B), INTENT(IN)                    :: MAXQPAR       ! maximum number of particles
   INTEGER(I4B), INTENT(IN)                    :: LAKEFLAG      ! >0 if processing lakes
   INTEGER(I4B), INTENT(IN), OPTIONAL          :: RSTEP         ! retrospective time step offset
   ! output variables
   integer(i4b), intent(out)                   :: ierr          ! error code
   character(*), intent(out)                   :: message       ! error message
   ! (1) extract flow from upstream reaches and append to the non-routed flow in JRCH
   INTEGER(I4B)                                :: NUPS          ! number of upstream reaches
   REAL(DP),DIMENSION(:),allocatable,SAVE      :: Q_JRCH        ! flow in downstream reach JRCH
   REAL(DP),DIMENSION(:),allocatable,SAVE      :: TENTRY        ! entry time to JRCH (exit time u/s)
   INTEGER(I4B)                                :: NQ1           ! # flow particles
   ! (2) route flow within the current [JRCH] river segment     
   INTEGER(I4B)                                :: ROFFSET       ! retrospective offset due to rstep
   REAL(DP)                                    :: T_START       ! start of time step
   REAL(DP)                                    :: T_END         ! end of time step
   REAL(DP),DIMENSION(:),allocatable,SAVE      :: T_EXIT        ! time particle expected exit JRCH
   LOGICAL(LGT),DIMENSION(:),allocatable,SAVE  :: FROUTE        ! routing flag .T. if particle exits
   INTEGER(I4B)                                :: NQ2           ! # flow particles (<=NQ1 b/c merge)
   ! (3) calculate time-step averages
   INTEGER(I4B)                                :: NR            ! # routed particles
   INTEGER(I4B)                                :: NN            ! # non-routed particles
   REAL(DP),DIMENSION(2)                       :: TNEW          ! start/end of time step
   REAL(DP),DIMENSION(1)                       :: QNEW          ! interpolated flow
   ! (4) housekeeping
   REAL(DP)                                    :: Q_END         ! flow at the end of the timestep
   REAL(DP)                                    :: TIMEI         ! entry time at the end of the timestep
   TYPE(FPOINT),allocatable,DIMENSION(:),SAVE  :: NEW_WAVE      ! temporary wave
   LOGICAL(LGT),SAVE                           :: INIT=.TRUE.   ! used to initialize pointers
   ! random stuff
   CHARACTER(LEN=256)                          :: CMESSAGE      ! error message for downwind routine

   ! initialize error control
   ierr=0; message='QROUTE_RCH/'
   ! ----------------------------------------------------------------------------------------
   ! (0) INITIALIZE POINTERS
   ! ----------------------------------------------------------------------------------------
   if(INIT) then
     INIT=.false.
     !NULLIFY(Q_JRCH,TENTRY,T_EXIT,FROUTE,NEW_WAVE)
     !deallocate(Q_JRCH,TENTRY,T_EXIT,FROUTE,NEW_WAVE)
   endif
   RCHFLX(IENS,JRCH)%TAKE=0.0_dp ! initialize take from this reach
    ! ----------------------------------------------------------------------------------------
    ! (1) EXTRACT FLOW FROM UPSTREAM REACHES & APPEND TO THE NON-ROUTED FLOW PARTICLES IN JRCH
    ! ----------------------------------------------------------------------------------------
    NUPS = count(NETOPO(JRCH)%goodBas)        ! number of desired upstream reaches
    !NUPS = size(NETOPO(JRCH)%UREACHI)        ! number of upstream reaches
    IF (NUPS.GT.0) THEN
      call GETUSQ_RCH(IENS,JRCH,LAKEFLAG,T0,T1, &                    ! input
                      Q_JRCH,TENTRY,T_EXIT,ierr,cmessage,RSTEP)      ! output
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      ! check for negative flow
      if (MINVAL(Q_JRCH).lt.0.0_dp) then
        ierr=20; message=trim(message)//'negative flow extracted from upstream reach'; return
      endif
      ! check
      if(JRCH==ixPrint)then
       print*, 'JRCH, Q_JRCH = ', JRCH, Q_JRCH
      endif 
    else
      ! set flow in headwater reaches to modelled streamflow from time delay histogram
      RCHFLX(IENS,JRCH)%REACH_Q = RCHFLX(IENS,JRCH)%BASIN_QR(1) 
      if (allocated(KROUTE(IENS,JRCH)%KWAVE)) THEN
        deallocate(KROUTE(IENS,JRCH)%KWAVE,STAT=IERR)
        if(ierr/=0)then; message=trim(message)//'problem deallocating space for KROUTE'; return; endif
      endif
      allocate(KROUTE(IENS,JRCH)%KWAVE(0:0),STAT=ierr)     
      if(ierr/=0)then; message=trim(message)//'problem allocating space for KROUTE(IENS,JRCH)%KWAVE(1)'; return; endif
      KROUTE(IENS,JRCH)%KWAVE(0)%QF=-9999 
      KROUTE(IENS,JRCH)%KWAVE(0)%TI=-9999 
      KROUTE(IENS,JRCH)%KWAVE(0)%TR=-9999 
      KROUTE(IENS,JRCH)%KWAVE(0)%RF=.False. 
      KROUTE(IENS,JRCH)%KWAVE(0)%QM=-9999
      return  ! no upstream reaches (routing for sub-basins done using time-delay histogram)
    endif
    ! ----------------------------------------------------------------------------------------
    ! (2) REMOVE FLOW PARTICLES (REDUCE MEMORY USAGE AND PROCESSING TIME)
    ! ----------------------------------------------------------------------------------------
    if (size(Q_JRCH).GT.MAXQPAR) then
      call REMOVE_RCH(MAXQPAR,Q_JRCH,TENTRY,T_EXIT,ierr,cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif
    NQ1 = SIZE(Q_JRCH)-1                                     ! -1 because of the zero element
    ! ----------------------------------------------------------------------------------------
    ! (3) ROUTE FLOW WITHIN THE CURRENT [JRCH] RIVER SEGMENT
    ! ----------------------------------------------------------------------------------------
    ! set the retrospective offset
    if (.not.present(RSTEP)) then
      ROFFSET = 0
    else
      ROFFSET = RSTEP
    endif
    ! set time boundaries
    T_START = T0 - (T1 - T0)*ROFFSET
    T_END   = T1 - (T1 - T0)*ROFFSET
    allocate(FROUTE(0:NQ1),STAT=IERR)
    if(ierr/=0)then; message=trim(message)//'problem allocating space for FROUTE'; return; endif
    FROUTE(0) = .TRUE.; FROUTE(1:NQ1)=.FALSE.  ! init. routing flags
    ! route flow through the current [JRCH] river segment (Q_JRCH in units of m2/s)
    call KINWAV_RCH(JRCH,T_START,T_END,Q_JRCH(1:NQ1),TENTRY(1:NQ1),&     ! (input)
                                       FROUTE(1:NQ1),T_EXIT(1:NQ1),NQ2,& ! (output)
                                       ierr,cmessage)                    ! (output)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    if(JRCH == ixPrint)then
      print*, 'Q_JRCH = ', Q_JRCH
      print*, 'FROUTE = ', FROUTE
      print*, 'TENTRY = ', TENTRY
      print*, 'T_EXIT = ', T_EXIT
    endif 
    ! ----------------------------------------------------------------------------------------
    ! (4) COMPUTE TIME-STEP AVERAGES
    ! ----------------------------------------------------------------------------------------
    NR = COUNT(FROUTE)-1   ! -1 because of the zero element (last routed)
    NN = NQ2-NR            ! number of non-routed points
    TNEW = (/T_START,T_END/)
    ! (zero position last routed; use of NR+1 instead of NR keeps next expected routed point) 
    call INTERP_RCH(T_EXIT(0:NR+1),Q_JRCH(0:NR+1),TNEW,QNEW,IERR,CMESSAGE)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    if(JRCH == ixPrint) print*, 'QNEW(1) = ', QNEW(1)
    ! m2/s --> m3/s + instantaneous runoff from basin
    RCHFLX(IENS,JRCH)%REACH_Q = QNEW(1)*RPARAM(JRCH)%R_WIDTH + RCHFLX(IENS,JRCH)%BASIN_QR(1)
    ! ----------------------------------------------------------------------------------------
    ! (5) HOUSEKEEPING
    ! ----------------------------------------------------------------------------------------
    ! compute the instantaneous flow at the end of the time step
    !   (last routed point)
    Q_END = Q_JRCH(NR) + &   !        (dQ/dT)                                 (dT)
             ( (Q_JRCH(NR+1)-Q_JRCH(NR)) / (T_EXIT(NR+1)-T_EXIT(NR)) ) * (T_END-T_EXIT(NR))
    ! compute an approximate entry time (needed for the remove routine later)
    TIMEI = TENTRY(NR) + &   !        (dT/dT)                                 (dT)
             ( (TENTRY(NR+1)-TENTRY(NR)) / (T_EXIT(NR+1)-T_EXIT(NR)) ) * (T_END-T_EXIT(NR))
    ! allocate space for the routed data (+1 to allocate space for the interpolated point)
    if (.not.allocated(KROUTE(IENS,JRCH)%KWAVE)) then
      ierr=20; message=trim(message)//'KROUTE is not associated'; return
    else
   ! if(allocated(KROUTE(IENS,JRCH)%KWAVE)) then
      deallocate(KROUTE(IENS,JRCH)%KWAVE, STAT=ierr)
   !   if (ierr.ne.0) then
   !     print *, '% CANNOT DEALLOCATE SPACE (FOR SOME UNKNOWN REASON)... TRY AND NULLIFY!'
   !     !NULLIFY(KROUTE(IENS,JRCH)%KWAVE)
   !     deallocate(KROUTE(IENS,JRCH)%KWAVE)
   !   endif
   ! endif
      allocate(KROUTE(IENS,JRCH)%KWAVE(0:NQ2+1),STAT=ierr)   ! NQ2 is number of points for kinematic routing
      if(ierr/=0)then; message=trim(message)//'problem allocating space for KROUTE(IENS,JRCH)%KWAVE(0:NQ2+1)'; return; endif
    endif
    ! insert the interpolated point (TI is irrelevant, as the point is "routed")
    KROUTE(IENS,JRCH)%KWAVE(NR+1)%QF=Q_END;   KROUTE(IENS,JRCH)%KWAVE(NR+1)%TI=TIMEI
    KROUTE(IENS,JRCH)%KWAVE(NR+1)%TR=T_END;   KROUTE(IENS,JRCH)%KWAVE(NR+1)%RF=.TRUE.
    ! add the output from kinwave...         - skip NR+1
    ! (when JRCH becomes IR routed points will be stripped out & the structures updated again)
    KROUTE(IENS,JRCH)%KWAVE(0:NR)%QF=Q_JRCH(0:NR); KROUTE(IENS,JRCH)%KWAVE(NR+2:NQ2+1)%QF=Q_JRCH(NR+1:NQ2)
    KROUTE(IENS,JRCH)%KWAVE(0:NR)%TI=TENTRY(0:NR); KROUTE(IENS,JRCH)%KWAVE(NR+2:NQ2+1)%TI=TENTRY(NR+1:NQ2)
    KROUTE(IENS,JRCH)%KWAVE(0:NR)%TR=T_EXIT(0:NR); KROUTE(IENS,JRCH)%KWAVE(NR+2:NQ2+1)%TR=T_EXIT(NR+1:NQ2)
    KROUTE(IENS,JRCH)%KWAVE(0:NR)%RF=FROUTE(0:NR); KROUTE(IENS,JRCH)%KWAVE(NR+2:NQ2+1)%RF=FROUTE(NR+1:NQ2)
    KROUTE(IENS,JRCH)%KWAVE(0:NQ2+1)%QM=-9999
    ! implement water use
    !IF (NUSER.GT.0.AND.UCFFLAG.GE.1) THEN
      !CALL EXTRACT_FROM_RCH(IENS,JRCH,NR,Q_JRCH,T_EXIT,T_END,TNEW)
    !ENDIF
    ! free up space for the next reach
    deallocate(Q_JRCH,TENTRY,T_EXIT,FROUTE,STAT=IERR)   ! FROUTE defined in this sub-routine
    if(ierr/=0)then; message=trim(message)//'problem deallocating space for [Q_JRCH, TENTRY, T_EXIT, FROUTE]'; return; endif
    ! ***
    ! remove flow particles from the most downstream reach
    if(jRch==ixOutlet)then
      if(jRch==ixPrint) print*, 'NETOPO(JRCH)%DREACHI = ', NETOPO(JRCH)%DREACHI
      if(NETOPO(JRCH)%DREACHI > 0)then
        message=trim(message)//'downstream reach index of outlet reach is greater than zero'
        ierr=20; return
      endif
    endif
    ! if the last reach or lake inlet (and lakes are enables), remove routed elements from memory
    IF (NETOPO(JRCH)%DREACHI.LT.0 .OR. &  ! if the last reach, then there is no downstream reach
        (LAKEFLAG.EQ.1.AND.NETOPO(JRCH)%LAKINLT)) THEN ! if lake inlet
      ! copy data to a temporary wave
      if (allocated(NEW_WAVE)) THEN
        DEALLOCATE(NEW_WAVE,STAT=IERR)
        if(ierr/=0)then; message=trim(message)//'problem deallocating space for NEW_WAVE'; return; endif
      endif
      ALLOCATE(NEW_WAVE(0:NN),STAT=IERR)  ! NN = number non-routed (the zero element is the last routed point)
      if(ierr/=0)then; message=trim(message)//'problem allocating space for NEW_WAVE'; return; endif
      NEW_WAVE(0:NN) = KROUTE(IENS,JRCH)%KWAVE(NR+1:NQ2+1)  ! +1 because of the interpolated point
      ! re-size wave structure
      if (allocated(KROUTE(IENS,JRCH)%KWAVE)) THEN
        deallocate(KROUTE(IENS,JRCH)%KWAVE,STAT=IERR)
        if(ierr/=0)then; message=trim(message)//'problem deallocating space for KROUTE'; return; endif
      endif
      allocate(KROUTE(IENS,JRCH)%KWAVE(0:NN),STAT=IERR)  ! again, the zero element for the last routed point
      if(ierr/=0)then; message=trim(message)//'problem allocating space for KROUTE'; return; endif
      ! copy data back to the wave structure and deallocate space for the temporary wave
      KROUTE(IENS,JRCH)%KWAVE(0:NN) = NEW_WAVE(0:NN)
      DEALLOCATE(NEW_WAVE,STAT=IERR)
      if(ierr/=0)then; message=trim(message)//'problem deallocating space for NEW_WAVE'; return; endif
    endif  ! (if JRCH is the last reach)
    ! get size of wave number for a reach
    return
  end subroutine

 ! *********************************************************************
 ! subroutine: extract flow from the reaches upstream of JRCH
 ! *********************************************************************
 subroutine GETUSQ_RCH(IENS,JRCH,LAKEFLAG,T0,T1,&                ! input
                       Q_JRCH,TENTRY,T_EXIT,ierr,message,&       ! output
                       RSTEP)                                    ! optional input
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Ross Woods, 2000; Martyn Clark, 2006
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Used to extract routed flow from the reaches upstream of JRCH, and append to the
 !    non-routed flow in JRCH
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O:
 !
 !   Input(s):
 !       IENS: Ensemble member
 !       JRCH: Index of the downstream reach
 !   LAKEFLAG: >0 means process lakes
 !         T0: Start of the time step
 !         T1: End of the time step
 !      RSTEP: Retrospective time step
 !
 !   Outputs:
 !  Q_JRCH(:): Vector of merged flow particles in reach JRCH
 !  TENTRY(:): Vector of times flow particles entered reach JRCH (exited upstream reaches)
 !  T_EXIT(:): Vector of times flow particles are expected to exit reach JRCH
 !
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! ----------------------------------------------------------------------------------------
 USE reachparam
 USE reachstate
 USE reach_flux
 USE lakes_flux
 USE lake_param
 IMPLICIT NONE
 ! Input
 INTEGER(I4B), INTENT(IN)                    :: IENS     ! ensemble member
 INTEGER(I4B), INTENT(IN)                    :: JRCH     ! reach to process
 INTEGER(I4B), INTENT(IN)                    :: LAKEFLAG ! >0 if processing lakes
 REAL(DP), INTENT(IN)                        :: T0,T1    ! start and end of the time step
 INTEGER(I4B), INTENT(IN), OPTIONAL          :: RSTEP    ! retrospective time step offset
 ! Local variables to hold the merged inputs to the downstream reach
 INTEGER(I4B)                                :: ROFFSET  ! retrospective offset due to rstep
 REAL(DP)                                    :: DT       ! model time step
 REAL(DP), DIMENSION(:),allocatable          :: QD       ! merged downstream flow
 REAL(DP), DIMENSION(:),allocatable          :: TD       ! merged downstream time
 INTEGER(I4B)                                :: ND       ! # points shifted downstream
 INTEGER(I4B)                                :: NJ       ! # points in the JRCH reach
 INTEGER(I4B)                                :: NK       ! # points for routing (NJ+ND)
 INTEGER(I4B)                                :: ILAK     ! lake index
 character(len=256)                          :: cmessage ! error message for downwind routine
 ! Output
 REAL(DP),DIMENSION(:),allocatable           :: Q_JRCH   ! merged (non-routed) flow in JRCH
 REAL(DP),DIMENSION(:),allocatable           :: TENTRY   ! time flow particles entered JRCH
 REAL(DP),DIMENSION(:),allocatable           :: T_EXIT   ! time flow is expected to exit JR
 integer(i4b), intent(out)                   :: ierr     ! error code
 character(*), intent(out)                   :: message  ! error message
 ! initialize error control
 ierr=0; message='GETUSQ_RCH/'
 ! ----------------------------------------------------------------------------------------
 ! (1) EXTRACT (AND MERGE) FLOW FROM UPSTREAM REACHES OR LAKE
 ! ----------------------------------------------------------------------------------------
 ! define dt
 DT = (T1 - T0)
 ! set the retrospective offset
 IF (.NOT.PRESENT(RSTEP)) THEN
   ROFFSET = 0
 ELSE
   ROFFSET = RSTEP
 END IF
 IF (LAKEFLAG.EQ.1) THEN                      ! lakes are enabled
  ! get lake outflow and only lake outflow if reach is a lake outlet reach, else do as normal 
  ILAK = NETOPO(JRCH)%LAKE_IX                 ! lake index
  IF (ILAK.GT.0) THEN                         ! part of reach is in lake
   IF (NETOPO(JRCH)%REACHIX.EQ.LKTOPO(ILAK)%DREACHI) THEN  ! we are in a lake outlet reach
    ND = 1
    ALLOCATE(QD(1),TD(1),STAT=IERR)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for QD and TD'; return; endif
    QD(1) = LAKFLX(IENS,ILAK)%LAKE_Q / RPARAM(JRCH)%R_WIDTH  ! lake outflow per unit reach width           
    TD(1) = T1 - DT*ROFFSET
   ELSE
    CALL QEXMUL_RCH(IENS,JRCH,T0,T1,ND,QD,TD,ierr,cmessage,RSTEP)        ! do as normal for unsubmerged part of inlet reach
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   ENDIF
  ELSE
   CALL QEXMUL_RCH(IENS,JRCH,T0,T1,ND,QD,TD,ierr,cmessage,RSTEP)         ! not in lake; do as normal
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  ENDIF
 ELSE                                         ! lakes disabled
  CALL QEXMUL_RCH(IENS,JRCH,T0,T1,ND,QD,TD,ierr,cmessage,RSTEP)         ! includes merging flow from different reaches
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  if(JRCH == ixPrint) print*, 'after QEXMUL_RCH: JRCH, ND, QD = ', JRCH, ND, QD
 ENDIF
 ! ----------------------------------------------------------------------------------------
 ! (2) EXTRACT NON-ROUTED FLOW FROM THE REACH JRCH & APPEND TO THE FLOW JUST ROUTED D/S
 ! ----------------------------------------------------------------------------------------
 ! check that the routing structure is associated
 if(allocated(KROUTE).eqv..FALSE.)THEN
  ierr=20; message='routing structure KROUTE is not associated'; return
 endif
 ! check that the wave has been initialized
 if (allocated(KROUTE(IENS,JRCH)%KWAVE).eqv..FALSE.) THEN
  ! if not initialized, then set initial flow to first flow
  ! (this will only occur for a cold start in the case of no streamflow observations)
  allocate(KROUTE(IENS,JRCH)%KWAVE(0:0),STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for KWAVE'; return; endif
  KROUTE(IENS,JRCH)%KWAVE(0)%QF = QD(1)
  KROUTE(IENS,JRCH)%KWAVE(0)%TI = T0 - DT - DT*ROFFSET
  KROUTE(IENS,JRCH)%KWAVE(0)%TR = T0      - DT*ROFFSET
  KROUTE(IENS,JRCH)%KWAVE(0)%RF = .TRUE.
 endif
 ! now extract the non-routed flow
 ! NB: routed flows were stripped out in the previous timestep when JRCH was index of u/s reach
 !  {only non-routed flows remain in the routing structure [ + zero element (last routed)]}
 NJ = SIZE(KROUTE(IENS,JRCH)%KWAVE) - 1           ! number of elements not routed (-1 for 0)
 NK = NJ + ND                                     ! pts still in reach + u/s pts just routed
 ALLOCATE(Q_JRCH(0:NK),TENTRY(0:NK),T_EXIT(0:NK),STAT=IERR) ! include zero element for INTERP later
 if(ierr/=0)then; message=trim(message)//'problem allocating array for [Q_JRCH, TENTRY, T_EXIT]'; return; endif
 Q_JRCH(0:NJ) = KROUTE(IENS,JRCH)%KWAVE(0:NJ)%QF  ! extract the non-routed flow from reach JR
 TENTRY(0:NJ) = KROUTE(IENS,JRCH)%KWAVE(0:NJ)%TI  ! extract the non-routed time from reach JR
 T_EXIT(0:NJ) = KROUTE(IENS,JRCH)%KWAVE(0:NJ)%TR  ! extract the expected exit time
 Q_JRCH(NJ+1:NJ+ND) = QD(1:ND)                    ! append u/s flow just routed downstream
 TENTRY(NJ+1:NJ+ND) = TD(1:ND)                    ! append u/s time just routed downstream
 T_EXIT(NJ+1:NJ+ND) = -9999.0D0                   ! set un-used T_EXIT to missing
 deallocate(QD,TD,STAT=IERR)                      ! routed flow appended, no longer needed
 if(ierr/=0)then; message=trim(message)//'problem deallocating array for QD and TD'; return; endif

 end subroutine

 ! *********************************************************************
 ! subroutine: extract flow from multiple reaches and merge into
 !                 a single series
 ! *********************************************************************
 subroutine QEXMUL_RCH(IENS,JRCH,T0,T1,&            ! input
                       ND,QD,TD,ierr,message,&      ! output
                       RSTEP)                       ! optional input
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Ross Woods, 2000 (two upstream reaches)
 !   Martyn Clark, 2006 (generalize to multiple upstream reaches, and major code overhaul)
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Used to extract routed flow from a multiple reaches upstream of JRCH.  Merges flow
 !    from multiple reaches into a single series.
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O:
 !
 !   Input(s):
 !       IENS: Ensemble member
 !       JRCH: Index of the downstream reach
 !         T0: Start of the time step
 !         T1: End of the time step
 !      RSTEP: Retrospective time step
 !
 !   Outputs:
 !      ND   : Number of routed particles
 !      QD(:): Vector of merged flow particles in reach JRCH
 !      TD(:): Vector of times flow particles entered reach JRCH (exited upstream reaches)
 !
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! ----------------------------------------------------------------------------------------
 USE reachparam
 USE reachstate
 USE reach_flux
 IMPLICIT NONE
 ! Input
 INTEGER(I4B), INTENT(IN)                    :: IENS      ! ensemble member
 INTEGER(I4B), INTENT(IN)                    :: JRCH      ! reach to process
 REAL(DP), INTENT(IN)                        :: T0,T1     ! start and end of the time step
 INTEGER(I4B), INTENT(IN), OPTIONAL          :: RSTEP     ! retrospective time step offset
 ! Local variables to hold flow/time from upstream reaches
 REAL(DP)                                    :: DT        ! model time step
 INTEGER(I4B)                                :: ROFFSET   ! retrospective offset due to rstep
 INTEGER(I4B)                                :: IUPS      ! loop through u/s reaches
 INTEGER(I4B)                                :: NUPB      ! number of upstream basins
 INTEGER(I4B)                                :: NUPR      ! number of upstream reaches
 INTEGER(I4B)                                :: INDX      ! index of the IUPS u/s reach
 INTEGER(I4B)                                :: MUPR      ! # reaches u/s of IUPS u/s reach
 INTEGER(I4B)                                :: NUPS      ! number of upstream elements
 TYPE(KREACH),DIMENSION(:),allocatable,SAVE  :: USFLOW    ! waves for all upstream segments
 REAL(DP), DIMENSION(:), ALLOCATABLE         :: UWIDTH    ! width of all upstream segments
 INTEGER(I4B)                                :: IMAX      ! max number of upstream particles
 INTEGER(I4B)                                :: IUPR      ! counter for reaches with particles
 INTEGER(I4B)                                :: IR        ! index of the upstream reach
 INTEGER(I4B)                                :: NS        ! size of  the wave
 INTEGER(I4B)                                :: NR        ! # routed particles in u/s reach
 INTEGER(I4B)                                :: NQ        ! NR+1, if non-routed particle exists 
 TYPE(FPOINT),DIMENSION(:),allocatable,SAVE  :: NEW_WAVE  ! temporary wave
 LOGICAL(LGT),SAVE                           :: INIT=.TRUE. ! used to initialize pointers
 ! Local variables to merge flow
 LOGICAL(LGT), DIMENSION(:), ALLOCATABLE     :: MFLG      ! T = all particles processed
 INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: ITIM      ! processing point for all u/s segments
 REAL(DP), DIMENSION(:), ALLOCATABLE         :: CTIME     ! central time for each u/s segment
 INTEGER(I4B)                                :: JUPS      ! index of reach with the earliest time
 REAL(DP)                                    :: Q_AGG     ! aggregarted flow at a given time
 INTEGER(I4B)                                :: IWAV      ! index of particle in the IUPS reach
 REAL(DP)                                    :: SCFAC     ! scale to conform to d/s reach width
 REAL(DP)                                    :: SFLOW     ! scaled flow at CTIME(JUPS)
 INTEGER(I4B)                                :: IBEG,IEND ! indices for particles that bracket time
 REAL(DP)                                    :: SLOPE     ! slope for the interpolation
 REAL(DP)                                    :: PREDV     ! value predicted by the interpolation
 INTEGER(I4B)                                :: IPRT      ! counter for flow particles
 INTEGER(I4B)                                :: JUPS_OLD  ! check that we don't get stuck in do-forever
 INTEGER(I4B)                                :: ITIM_OLD  ! check that we don't get stuck in do-forever
 REAL(DP)                                    :: TIME_OLD  ! previous time -- used to check for duplicates
 REAL(DP),DIMENSION(:),allocatable,SAVE      :: QD_TEMP   ! flow particles just enetered JRCH
 REAL(DP),DIMENSION(:),allocatable,SAVE      :: TD_TEMP   ! time flow particles entered JRCH
 ! Output
 INTEGER(I4B), INTENT(OUT)                   :: ND        ! number of routed particles
 REAL(DP), DIMENSION(:),allocatable          :: QD        ! flow particles just enetered JRCH
 REAL(DP), DIMENSION(:),allocatable          :: TD        ! time flow particles entered JRCH
 integer(i4b), intent(out)                   :: ierr     ! error code
 character(*), intent(out)                   :: message  ! error message
 ! initialize error control
 ierr=0; message='QEXMUL_RCH/'
 ! ----------------------------------------------------------------------------------------
 ! (0) INITIALIZE POINTERS
 ! ----------------------------------------------------------------------------------------
 IF(INIT) THEN
  INIT=.FALSE.
  !NULLIFY(USFLOW,NEW_WAVE,QD_TEMP,TD_TEMP)
  !deallocate(USFLOW,NEW_WAVE,QD_TEMP,TD_TEMP)
 ENDIF
 ! set the retrospective offset
 IF (.NOT.PRESENT(RSTEP)) THEN
   ROFFSET = 0
 ELSE
   ROFFSET = RSTEP
 END IF
 ! define dt
 DT = (T1 - T0)
 ! ----------------------------------------------------------------------------------------
 ! (1) DETERMINE THE NUMBER OF UPSTREAM REACHES
 ! ----------------------------------------------------------------------------------------
 ! ** SPECIAL CASE ** of no basins with any area
 if(count(NETOPO(JRCH)%goodBas)==0)return
 ! Need to extract and merge the runoff from all upstream BASINS as well as the streamflow
 !  from all upstream REACHES.  However, streamflow in headwater basins is undefined.  Thus
 !  the number of series merged from upstream reaches is the number of upstream basins +
 !  the number of upstream reaches that are not headwater basins.
 NUPR = 0                               ! number of upstream reaches
 NUPB = SIZE(NETOPO(JRCH)%UREACHI)      ! number of upstream basins
 !NUPB = count(NETOPO(JRCH)%goodBas)      ! number of upstream basins
 DO IUPS=1,NUPB
  INDX = NETOPO(JRCH)%UREACHI(IUPS)     ! index of the IUPS upstream reach
  !MUPR = SIZE(NETOPO(INDX)%UREACHI)     ! # reaches upstream of the IUPS upstream reach
  MUPR = count(NETOPO(INDX)%goodBas)     ! # reaches upstream of the IUPS upstream reach
  IF (MUPR.GT.0) NUPR = NUPR + 1        ! reach has streamflow in it, so add that as well
 END DO  ! iups
 NUPS = NUPB + NUPR                     ! number of upstream elements (basins + reaches)
 !print*, 'NUPB, NUPR, NUPS', NUPB, NUPR, NUPS
 !print*, 'NETOPO(JRCH)%UREACHK = ', NETOPO(JRCH)%UREACHK
 !print*, 'NETOPO(JRCH)%goodBas = ', NETOPO(JRCH)%goodBas
 ! if nups eq 1, then ** SPECIAL CASE ** of just one upstream basin that is a headwater
 IF (NUPS.EQ.1) THEN
  ND = 1
  ALLOCATE(QD(1),TD(1),STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem allocating array QD and TD'; return; endif
  ! get reach index
  IR = NETOPO(JRCH)%UREACHI(1)
  ! get flow in m2/s (scaled by with of downstream reach)
  QD(1) = RCHFLX(IENS,IR)%BASIN_QR(1)/RPARAM(JRCH)%R_WIDTH
  TD(1) = T1
  if(JRCH == ixPrint) print*, 'special case: JRCH, IR = ', JRCH, IR
  RETURN
 ENDIF
 ! allocate space for the upstream flow, time, and flags
 ALLOCATE(USFLOW(NUPS),UWIDTH(NUPS),CTIME(NUPS),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [USFLOW, UWIDTH, CTIME]'; return; endif
 ! define the minimum size of the routed data structure (number of flow particles)
 !  (IMAX is increased when looping through the reaches -- section 3 below)
 IMAX = NUPB                            ! flow from basins (one particle / timestep)
 ! ----------------------------------------------------------------------------------------
 ! (2) EXTRACT FLOW FROM UPSTREAM BASINS
 ! ----------------------------------------------------------------------------------------
 DO IUPS=1,NUPB
  ! identify the index for the IUPS upstream segment
  IR = NETOPO(JRCH)%UREACHI(IUPS)
  ! allocate space for the IUPS stream segment (flow, time, and flags)
  ALLOCATE(USFLOW(IUPS)%KWAVE(0:1),STAT=IERR)  ! basin, has flow @start and @end of the time step
  if(ierr>0)then; message=trim(message)//'problem allocating array USFLOW(IUPS)%KWAVE'; return; endif
  ! place flow and time in the KWAVE array (routing done with time-delay histogram in TIMDEL_BAS.F90)
  USFLOW(IUPS)%KWAVE(0:1)%QF = RCHFLX(IENS,IR)%BASIN_QR(0:1)          ! flow
  USFLOW(IUPS)%KWAVE(0:1)%TI = (/T0,T1/) - DT*ROFFSET                 ! entry time (not used)
  USFLOW(IUPS)%KWAVE(0:1)%TR = (/T0,T1/) - DT*ROFFSET                 ! exit time
  USFLOW(IUPS)%KWAVE(0:1)%RF = .TRUE.                                 ! routing flag
  !write(*,'(a,i4,1x,2(e20.10,1x))') 'IR, USFLOW(IUPS)%KWAVE(0:1)%QF = ', IR, USFLOW(IUPS)%KWAVE(0:1)%QF
  ! save the upstream width
  UWIDTH(IUPS) = 1.0D0                         ! basin = unit width
  ! save the the time for the first particle in each reach
  CTIME(IUPS) = USFLOW(IUPS)%KWAVE(1)%TR       ! central time
 END DO    ! (loop through upstream basins)
 ! ----------------------------------------------------------------------------------------
 ! (3) EXTRACT FLOW FROM UPSTREAM REACHES
 ! ----------------------------------------------------------------------------------------
 IUPR = 0
 DO IUPS=1,NUPB
  INDX = NETOPO(JRCH)%UREACHI(IUPS)     ! index of the IUPS upstream reach
  !MUPR = SIZE(NETOPO(INDX)%UREACHI)     ! # reaches upstream of the IUPS upstream reach
  MUPR = count(NETOPO(INDX)%goodBas)     ! # reaches upstream of the IUPS upstream reach
  IF (MUPR.GT.0) THEN                   ! reach has streamflow in it, so add that as well
   IUPR = IUPR + 1
   ! identify the index for the IUPS upstream segment
   IR = NETOPO(JRCH)%UREACHI(IUPS)
   ! identify the size of the wave
   NS = SIZE(KROUTE(IENS,IR)%KWAVE)
   ! identify number of routed flow elements in the IUPS upstream segment
   NR = COUNT(KROUTE(IENS,IR)%KWAVE(:)%RF)
   ! include a non-routed point, if it exists
   NQ = MIN(NR+1,NS)
   ! allocate space for the IUPS stream segment (flow, time, and flags)
   ALLOCATE(USFLOW(NUPB+IUPR)%KWAVE(0:NQ-1),STAT=IERR)  ! (zero position = last routed)
   if(ierr/=0)then; message=trim(message)//'problem allocating array USFLOW(NUPB+IUPR)%KWAVE(0:NQ-1)'; return; endif
   ! place data in the new arrays
   USFLOW(NUPB+IUPR)%KWAVE(0:NQ-1) = KROUTE(IENS,IR)%KWAVE(0:NQ-1)
   ! here a statement where we check for a modification in the upstream reach;
   ! if flow upstream is modified, then copy KROUTE(:,:)%KWAVE(:)%QM to USFLOW(..)%KWAVE%QF
   !IF (NUSER.GT.0.AND.SIMDAT%UCFFLAG.GE.1) THEN !if the irrigation module is active and there are users
   !  IF (RCHFLX(IENS,IR)%TAKE.GT.0._DP) THEN !if take from upstream reach is greater then zero
   !    ! replace QF with modified flow (as calculated in extract_from_rch)
   !    USFLOW(NUPB+IUPR)%KWAVE(0:NQ-1)%QF = KROUTE(IENS,IR)%KWAVE(0:NQ-1)%QM
   !  ENDIF
   !ENDIF
   ! ...and REMOVE the routed particles from the upstream reach
   ! (copy the wave to a temporary wave)
   IF (allocated(NEW_WAVE)) THEN
     DEALLOCATE(NEW_WAVE,STAT=IERR)    ! (so we can allocate)
     if(ierr/=0)then; message=trim(message)//'problem deallocating array NEW_WAVE'; return; endif
   END IF
   ALLOCATE(NEW_WAVE(0:NS-1),STAT=IERR)                 ! get new wave
   if(ierr/=0)then; message=trim(message)//'problem allocating array NEW_WAVE'; return; endif
   NEW_WAVE(0:NS-1) = KROUTE(IENS,IR)%KWAVE(0:NS-1)  ! copy
   ! (re-size wave structure)
   IF (.NOT.allocated(KROUTE(IENS,IR)%KWAVE))then; print*,' not allocated. in qex ';return; endif
   IF (allocated(KROUTE(IENS,IR)%KWAVE)) THEN
     deallocate(KROUTE(IENS,IR)%KWAVE,STAT=IERR)
     if(ierr/=0)then; message=trim(message)//'problem deallocating array KROUTE'; return; endif
   END IF
   ALLOCATE(KROUTE(IENS,IR)%KWAVE(0:NS-NR),STAT=IERR)   ! reduced size
   if(ierr/=0)then; message=trim(message)//'problem allocating array KROUTE'; return; endif
   ! (copy "last routed" and "non-routed" elements)
   KROUTE(IENS,IR)%KWAVE(0:NS-NR) = NEW_WAVE(NR-1:NS-1)
   ! (de-allocate temporary wave)
   DEALLOCATE(NEW_WAVE,STAT=IERR)
   if(ierr/=0)then; message=trim(message)//'problem deallocating array NEW_WAVE'; return; endif
  ! NULLIFY(NEW_WAVE)
   ! save the upstream width
   UWIDTH(NUPB+IUPR) = RPARAM(IR)%R_WIDTH            ! reach, width = parameter 
   ! save the time for the first particle in each reach
   CTIME(NUPB+IUPR) = USFLOW(NUPB+IUPR)%KWAVE(1)%TR  ! central time
   ! keep track of the total number of points that must be routed downstream
   IMAX = IMAX + (NR-1)     ! exclude zero point for the last routed
  ENDIF ! if reach has particles in it
 END DO  ! iups
 ! ----------------------------------------------------------------------------------------
 ! (4) MERGE FLOW FROM MULTIPLE UPSTREAM REACHES
 ! ----------------------------------------------------------------------------------------
 ! This is a bit tricky.  Consider a given upstream reach x.  For all upstream reaches
 !  *other than* x, we need to estimate (interpolate) flow for the *times* associted with
 !  each of the flow particles in reach x.  Then, at a given time, we can sum the flow
 !  (routed in reach x plus interpolated flow in all other reaches).  This needs to be done
 !  for all upstream reaches. 
 ! ----------------------------------------------------------------------------------------
 ! We accomplish this as follows.  We define a vector of indices (ITIM), where each
 !  element of ITIM points to a particle in a given upstream reach still to be processed.
 !  We also define a vector of times (CTIME), which is the time of the flow particles that
 !  relate to the vector of indices ITIM.  We identify upstream reach with the earliest
 !  time in CTIME, save the flow, and produce corresponding flow estimates for the same
 !  time in all other reaches.  We then scale the flow and flow estimates in all upstream
 !  reaches by the width of the downstream reach, and sum the flow over all upstream reaches.
 !  We then move the index forward in ITIM (for the upstream reach just processed), get a
 !  new vector CTIME, and process the next earliest particle.  We continue until all
 !  flow particles are processed in all upstream reaches.  
 ! ----------------------------------------------------------------------------------------
 IPRT = 0  ! initialize counter for flow particles in the output array
 ! allocate space for the merged flow at the downstream reach
 ALLOCATE(QD_TEMP(IMAX),TD_TEMP(IMAX),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [QD_TEMP, TD_TEMP]'; return; endif
 ! allocate positional arrays
 ALLOCATE(MFLG(NUPS),ITIM(NUPS),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [MFLG, ITIM]'; return; endif
 ! initalize the flag that defines whether all particles in a given reach are processed
 MFLG(1:NUPS)  = .FALSE.                     ! false until all particles are processed
 ! initialize the search vector
 ITIM(1:NUPS)  = 1                           ! start with the first element of the wave
 ! initialize jups_old and itim_old (used to check we don't get stuck in the do-forever loop)
 JUPS_OLD = HUGE(JUPS_OLD)
 ITIM_OLD = HUGE(ITIM_OLD)
 DO      ! loop through all the times in the upstream reaches until no more routed flows
  ! find the reach with the earliest time in all upstream reaches
  !  (NB: the time at the start of the timestep is the earliest possible time and
  !       the time at the end of the timestep is the latest possible time)
  JUPS  = MINLOC(CTIME,DIM=1)     ! JUPS = reach w/ earliest time
  ! check that we're not stuck in a continuous do loop
  IF (JUPS.EQ.JUPS_OLD .AND. ITIM(JUPS).EQ.ITIM_OLD) THEN
   ierr=20; message=trim(message)//'stuck in the continuous do-loop'; return
  ENDIF
  ! save jups and itim(jups) to check that we don't get stuck in a continuous do-loop
  JUPS_OLD = JUPS
  ITIM_OLD = ITIM(JUPS)
  ! check that there are still particles in the given reach that require processing
  IF (.NOT.MFLG(JUPS)) THEN
   ! check that the particle in question is a particle routed (if not, then don't process)
   IF (USFLOW(JUPS)%KWAVE(ITIM(JUPS))%RF.EQV..FALSE.) THEN
    MFLG(JUPS) = .TRUE. ! if routing flag is false, then have already processed all particles
    CTIME(JUPS) = HUGE(SFLOW)  ! largest possible number = ensure reach is not selected again
   ! the particle is in need of processing
   ELSE
    ! define previous time
    IF (IPRT.GE.1) THEN
      TIME_OLD = TD_TEMP(IPRT)
    ELSE ! (if no particles, set to largest possible negative number)
      TIME_OLD = -HUGE(SFLOW)
    END IF
    ! check that the particles are being processed in the correct order
    IF (CTIME(JUPS).LT.TIME_OLD) THEN
     ierr=30; message=trim(message)//'expect process in order of time'; return
    ENDIF
    ! don't process if time already exists
    IF (CTIME(JUPS).NE.TIME_OLD) THEN
     ! -------------------------------------------------------------------------------------
     ! compute sum of scaled flow for all reaches
     Q_AGG = 0.0D0
     DO IUPS=1,NUPS
      ! identify the element of the wave for the IUPS upstream reach
      IWAV = ITIM(IUPS)
      ! compute scale factor (scale upstream flow by width of downstream reach)
      SCFAC = UWIDTH(IUPS) / RPARAM(JRCH)%R_WIDTH
      ! case of the upstream reach with the minimum time (no interpolation required)
      IF (IUPS.EQ.JUPS) THEN
       SFLOW = USFLOW(IUPS)%KWAVE(IWAV)%QF * SCFAC  ! scaled flow
      ! case of all other upstream reaches (*** now, interpolate ***)
      ELSE
       ! identify the elements that bracket the flow particle in the reach JUPS
       ! why .GE.?  Why not .GT.??
       IBEG = IWAV; IF (USFLOW(IUPS)%KWAVE(IBEG)%TR.GE.CTIME(JUPS)) IBEG=IWAV-1
       IEND = IBEG+1  ! *** check the elements are ordered as we think ***
       ! test if we have bracketed properly
       IF (USFLOW(IUPS)%KWAVE(IEND)%TR.LT.CTIME(JUPS) .OR. &
           USFLOW(IUPS)%KWAVE(IBEG)%TR.GT.CTIME(JUPS)) THEN
            ierr=40; message=trim(message)//'the times are not ordered as we assume'; return
       ENDIF  ! test for bracketing
       ! estimate flow for the IUPS upstream reach at time CTIME(JUPS)
       SLOPE = (USFLOW(IUPS)%KWAVE(IEND)%QF - USFLOW(IUPS)%KWAVE(IBEG)%QF) / &
               (USFLOW(IUPS)%KWAVE(IEND)%TR - USFLOW(IUPS)%KWAVE(IBEG)%TR)
       PREDV =  USFLOW(IUPS)%KWAVE(IBEG)%QF + SLOPE*(CTIME(JUPS)-USFLOW(IUPS)%KWAVE(IBEG)%TR)
       SFLOW = PREDV * SCFAC  ! scaled flow
      ENDIF  ! (if interpolating)
      ! aggregate flow
      Q_AGG = Q_AGG + SFLOW
     END DO  ! looping through upstream elements
     ! -------------------------------------------------------------------------------------
     ! place Q_AGG and CTIME(JUPS) in the output arrays
     IPRT = IPRT + 1
     QD_TEMP(IPRT) = Q_AGG
     TD_TEMP(IPRT) = CTIME(JUPS)
    ENDIF  ! (check that time doesn't already exist)
    ! check if the particle just processed is the last element
    IF (ITIM(JUPS).EQ.SIZE(USFLOW(JUPS)%KWAVE)-1) THEN  ! -1 because of the zero element
     MFLG(JUPS) = .TRUE.            ! have processed all particles in a given u/s reach
     CTIME(JUPS) = HUGE(SFLOW)      ! largest possible number = ensure reach is not selected again
    ELSE
     ITIM(JUPS) = ITIM(JUPS) + 1                       ! move on to the next flow element
     CTIME(JUPS) = USFLOW(JUPS)%KWAVE(ITIM(JUPS))%TR   ! save the time
    ENDIF  ! (check if particle is the last element)
   ENDIF  ! (check if the particle is a routed element)
  ENDIF  ! (check that there are still particles to process)
  ! if processed all particles in all upstream reaches, then EXIT
  IF (COUNT(MFLG).EQ.NUPS) EXIT
 END DO   ! do-forever
 ! free up memory
 DO IUPS=1,NUPS  ! de-allocate each element of USFLOW
  DEALLOCATE(USFLOW(IUPS)%KWAVE,STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem deallocating array USFLOW(IUPS)%KWAVE'; return; endif
 END DO          ! looping thru elements of USFLOW
 DEALLOCATE(USFLOW,UWIDTH,CTIME,ITIM,MFLG,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [USFLOW, UWIDTH, CTIME, ITIM, MFLG]'; return; endif
 ! ...and, save reduced arrays in QD and TD
 ND = IPRT
 ALLOCATE(QD(ND),TD(ND),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [QD, TD]'; return; endif
 QD(1:ND) = QD_TEMP(1:ND)
 TD(1:ND) = TD_TEMP(1:ND)
 DEALLOCATE(QD_TEMP,TD_TEMP,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [QD_TEMP, TD_TEMP]'; return; endif
 ! ----------------------------------------------------------------------------------------
 end subroutine

 ! *********************************************************************
 !  subroutine: removes flow particles from the routing structure,
 !                 to reduce memory usage and processing time
 ! *********************************************************************
 subroutine REMOVE_RCH(MAXQPAR,&                           ! input
                       Q_JRCH,TENTRY,T_EXIT,ierr,message)  ! output
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Martyn Clark, 2006
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Used to remove flow particles from the routing structure, to decrease memory usage
 !    and processing time
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O:
 !
 !  Q_JRCH(:): Vector of merged flow particles in reach JRCH
 !  TENTRY(:): Vector of times flow particles entered reach JRCH (exited upstream reaches)
 !  T EXIT(:): Vector of times flow particles are EXPECTED to exit reach JRCH
 !
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! ----------------------------------------------------------------------------------------
 IMPLICIT NONE
 ! Input
 INTEGER(I4B), INTENT(IN)                    :: MAXQPAR  ! maximum number of flow particles allowed
 ! output
 REAL(DP), DIMENSION(:),allocatable          :: Q_JRCH   ! merged (non-routed) flow in JRCH
 REAL(DP), DIMENSION(:),allocatable          :: TENTRY   ! time flow particles entered JRCH
 REAL(DP), DIMENSION(:),allocatable          :: T_EXIT   ! time flow particles exited JRCH
 integer(i4b), intent(out)                   :: ierr     ! error code
 character(*), intent(out)                   :: message  ! error message
 ! Local variables
 INTEGER(I4B)                                :: NPRT     ! number of flow particles
 INTEGER(I4B)                                :: IPRT     ! loop through flow particles
 REAL(DP), DIMENSION(:), ALLOCATABLE         :: Q,T,Z    ! copies of Q_JRCH and T_JRCH
 LOGICAL(LGT), DIMENSION(:), ALLOCATABLE     :: PARFLG   ! .FALSE. if particle removed
 INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: INDEX0   ! indices of original vectors  
 REAL(DP), DIMENSION(:), ALLOCATABLE         :: ABSERR   ! absolute error btw interp and orig
 REAL(DP)                                    :: Q_INTP   ! interpolated particle
 INTEGER(I4B)                                :: MPRT     ! local number of flow particles
 INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: INDEX1   ! indices of particles retained
 REAL(DP), DIMENSION(:), ALLOCATABLE         :: E_TEMP   ! temp abs error btw interp and orig
 INTEGER(I4B), DIMENSION(1)                  :: ITMP     ! result of minloc function
 INTEGER(I4B)                                :: ISEL     ! index of local minimum value
 INTEGER(I4B)                                :: INEG     ! lower boundary for interpolation
 INTEGER(I4B)                                :: IMID     ! desired point for interpolation
 INTEGER(I4B)                                :: IPOS     ! upper boundary for interpolation
 ! initialize error control
 ierr=0; message='REMOVE_RCH/'
 ! ----------------------------------------------------------------------------------------
 ! (1) INITIALIZATION
 ! ----------------------------------------------------------------------------------------
 ! get the number of particles
 NPRT = SIZE(Q_JRCH)-1                       ! -1 because of zero element
 ! allocate and initialize arrays
 ALLOCATE(Q(0:NPRT),T(0:NPRT),Z(0:NPRT),PARFLG(0:NPRT),INDEX0(0:NPRT),ABSERR(0:NPRT),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [Q, T, Z, PARFLG, INDEX0, ABSERR]'; return; endif
 Q = Q_JRCH; T = TENTRY     ! get copies of Q_JRCH and TENTRY
 Z = T_EXIT                 ! (not used in the interp, but include for consistency)
 PARFLG = .TRUE.            ! particle flag = start with all points
 INDEX0 = arth(0,1,NPRT+1)  ! index = (0,1,2,...,NPRT)
 ABSERR = HUGE(Q)           ! largest possible double-precision number
 ! get the absolte difference between actual points and interpolated points
 DO IPRT=1,NPRT-1
  ! interpolate at point (iprt)
  Q_INTP = INTERP(T(IPRT),Q(IPRT-1),Q(IPRT+1),T(IPRT-1),T(IPRT+1))
  ! save the absolute difference between the actual value and the interpolated value
  ABSERR(IPRT) = ABS(Q_INTP-Q(IPRT))
 END DO
 ! ----------------------------------------------------------------------------------------
 ! (2) REMOVAL
 ! ----------------------------------------------------------------------------------------
 DO  ! continue looping until the number of particles is below the limit
  ! get the number of particles still in the structure
  MPRT = COUNT(PARFLG)-1       ! -1 because of the zero element
  ! get a copy of (1) indices of selected points, and (2) the interpolation errors
  ALLOCATE(INDEX1(0:MPRT),E_TEMP(0:MPRT),STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem allocating arrays [INDEX1, E_TEMP]'; return; endif
  INDEX1 = PACK(INDEX0,PARFLG) ! (restrict attention to the elements still present)
  E_TEMP = PACK(ABSERR,PARFLG)
  ! check for exit condition (exit after "pack" b.c. indices used to construct final vectors)
  IF (MPRT.LT.MAXQPAR) EXIT
  ! get the index of the minimum value
  ITMP = MINLOC(E_TEMP)
  ISEL = LBOUND(E_TEMP,DIM=1) + ITMP(1) - 1 ! MINLOC assumes count from 1, here (0,1,2,...NPRT)
  ! re-interpolate the point immediately before the point flagged for removal
  IF (INDEX1(ISEL-1).GT.0) THEN
   INEG=INDEX1(ISEL-2); IMID=INDEX1(ISEL-1); IPOS=INDEX1(ISEL+1)
   Q_INTP = INTERP(T(IMID),Q(INEG),Q(IPOS),T(INEG),T(IPOS))
   ABSERR(IMID) = ABS(Q_INTP-Q(IMID))
  ENDIF
  ! re-interpolate the point immediately after the point flagged for removal
  IF (INDEX1(ISEL+1).LT.NPRT) THEN
   INEG=INDEX1(ISEL-1); IMID=INDEX1(ISEL+1); IPOS=INDEX1(ISEL+2)
   Q_INTP = INTERP(T(IMID),Q(INEG),Q(IPOS),T(INEG),T(IPOS))
   ABSERR(IMID) = ABS(Q_INTP-Q(IMID))
  ENDIF
  ! flag the point as "removed"
  PARFLG(INDEX1(ISEL)) = .FALSE.
  ! de-allocate arrays
  DEALLOCATE(INDEX1,E_TEMP,STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [INDEX1, E_TEMP]'; return; endif
 END DO  ! keep looping until a sufficient number of points are removed
 ! ----------------------------------------------------------------------------------------
 ! (3) RE-SIZE DATA STRUCTURES
 ! ----------------------------------------------------------------------------------------
 DEALLOCATE(Q_JRCH,TENTRY,T_EXIT,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [Q_JRCH, TENTRY, T_EXIT]'; return; endif
 ALLOCATE(Q_JRCH(0:MPRT),TENTRY(0:MPRT),T_EXIT(0:MPRT),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [Q_JRCH, TENTRY, T_EXIT]'; return; endif
 Q_JRCH = Q(INDEX1)
 TENTRY = T(INDEX1)
 T_EXIT = Z(INDEX1)
 DEALLOCATE(INDEX1,E_TEMP,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [INDEX1, E_TEMP]'; return; endif
 DEALLOCATE(Q,T,Z,PARFLG,ABSERR,INDEX0,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [Q, T, Z, PARFLG, INDEX0, ABSERR]'; return; endif

 contains

  function INTERP(T0,Q1,Q2,T1,T2)
    REAL(DP),INTENT(IN)                        :: Q1,Q2  ! flow at neighbouring times 
    REAL(DP),INTENT(IN)                        :: T1,T2  ! neighbouring times 
    REAL(DP),INTENT(IN)                        :: T0     ! desired time
    REAL(DP)                                   :: INTERP ! function name
    INTERP = Q1 + ( (Q2-Q1) / (T2-T1) ) * (T0-T1)
  end function

 end subroutine

 ! *********************************************************************
 ! new subroutine: calculate the propagation of kinematic waves in a
 !                 single stream segment, including the formation and
 !                 propagation of a kinematic shock
 ! *********************************************************************
 subroutine KINWAV_RCH(JRCH,T_START,T_END,Q_JRCH,TENTRY,FROUTE,T_EXIT,NQ2,&
                       ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Derek Goring, 1994; modified by Ross Woods, 2002
 !   Martyn Clark, 2006: complete overhaul
 !   Martyn Clark, 2012: extract to stand-alone model
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Used to calculate the propagation of kinematic waves in an individual stream
 !   segment, including the formation and propagation of a kinematic shock.
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O:
 !
 !   Input(s):
 !   ---------
 !      JRCH: index of reach processed
 !   T_START: start of the time step
 !     T_END: end of the time step
 !    Q_JRCH: array of flow elements -- neighbouring times are merged if a shock forms, but
 !                                      then merged flows are dis-aggregated to save the
 !                                      flow either side of a shock -- thus we may have
 !                                      fewer elements on output if several particles are
 !                                      merged, INTENT(INOUT)
 !    TENTRY: array of time elements -- neighbouring times are merged if a shock forms, 
 !                                      then merged times are dis-aggregated, one second is
 !                                      added to the time corresponding to the higer merged
 !                                      flow (note also fewer elements), INTENT(INOUT)
 !
 !   Outputs:
 !   --------
 !    FROUTE: array of routing flags -- All inputs are .FALSE., but flags change to .TRUE.
 !                                      if element is routed INTENT(OUT)
 !    T_EXIT: array of time elements -- identify the time each element is EXPECTED to exit
 !                                      the stream segment, INTENT(OUT).  Used in INTERPTS
 !       NQ2: number of particles    -- <= input becuase multiple particles may merge
 !
 ! ----------------------------------------------------------------------------------------
 ! Structures Used:
 !
 !   REACHPARAM:  Use channel slope, width, length, etc.
 !
 ! ----------------------------------------------------------------------------------------
 ! Method:
 !
 !   Flow routing through an individual stream segment is performed using the method
 !   described by Goring (1994).  For each wave in the stream segment, Goring's
 !   method calculates the celerity (m/s) and travel time, length/celerity (s).  If
 !   the initial time of the wave plus the travel time is less than the time at the
 !   end of a time step, then the wave is routed to the end of the stream segment and
 !   the time stamp is updated.  Otherwise, the wave remains in the stream segment
 !   and the time stamp remains constant.  In this context "earlier" times imply that
 !   a kinematic shock is located nearer the downstream end of a stream segment.
 !
 !   A decrease in the time between kinematic waves deceases eventually produces a
 !   discontinuity known as a kinematic shock.  When this occurs the kinematic waves
 !   are merged and the celerity is modified.  See Goring (1994) for more details.
 !
 ! ----------------------------------------------------------------------------------------
 ! Source:
 ! 
 !   This routine is based on the subroutine kinwav, located in kinwav_v7.f
 !
 ! ----------------------------------------------------------------------------------------
 ! Modifications to source (mclark@ucar.edu):
 !
 !   * All variables are now defined (IMPLICIT NONE) and described (comments)
 !
 !   * Parameters are defined within the subroutine (for ease of readibility)
 !
 !   * Added many comments
 !
 !   * Replaced GOTO statements with DO-CYCLE loops and DO-FOREVER loops with EXIT clause
 !      (for ease of readability), and replaced some do-continue loops w/ array operations
 !      and use F90 dynamic memory features
 !
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! ----------------------------------------------------------------------------------------
 USE reachparam
 USE nrutil, ONLY : arth
 IMPLICIT NONE
 ! Input
 INTEGER(I4B), INTENT(IN)                    :: JRCH     ! Reach to process
 REAL(DP), INTENT(IN)                        :: T_START  ! start of the time step
 REAL(DP), INTENT(IN)                        :: T_END    ! end of the time step
 ! Input/Output
 REAL(DP), DIMENSION(:), INTENT(INOUT)       :: Q_JRCH   ! flow to be routed
 REAL(DP), DIMENSION(:), INTENT(INOUT)       :: TENTRY   ! time to be routed
 REAL(DP), DIMENSION(:), INTENT(INOUT)       :: T_EXIT   ! time pts expected exit segment
 LOGICAL(LGT), DIMENSION(:), INTENT(INOUT)   :: FROUTE   ! routing flag, T=routed
 ! Output
 INTEGER(I4B), INTENT(OUT)                   :: NQ2      ! # particles (<= input b/c merge)
 integer(i4b), intent(out)                   :: ierr     ! error code
 character(*), intent(out)                   :: message  ! error message
 ! Internal
 REAL(DP)                                    :: ALFA     ! constant, 5/3
 REAL(DP)                                    :: K        ! sqrt(slope)/mannings N
 REAL(DP)                                    :: XMX      ! length of the stream segment
 INTEGER(I4B)                                :: NN       ! number of input points
 INTEGER(I4B)                                :: NI       ! original size of the input
 INTEGER(I4B)                                :: NM       ! mumber of merged elements
 INTEGER(I4B), DIMENSION(SIZE(Q_JRCH))       :: IX       ! minimum index of each merged element
 INTEGER(I4B), DIMENSION(SIZE(Q_JRCH))       :: MF       ! index for input element merged
 REAL(DP), DIMENSION(SIZE(Q_JRCH))           :: T0,T1,T2 ! copy of input time
 REAL(DP), DIMENSION(SIZE(Q_JRCH))           :: Q0,Q1,Q2 ! flow series
 REAL(DP), DIMENSION(SIZE(Q_JRCH))           :: WC       ! wave celerity
 INTEGER(I4B)                                :: IW,JW    ! looping variables, break check
 REAL(DP)                                    :: X,XB     ! define smallest, biggest shock
 REAL(DP)                                    :: WDIFF    ! difference in wave celerity-1
 REAL(DP)                                    :: XXB      ! wave break
 INTEGER(I4B)                                :: IXB,JXB  ! define position of wave break
 REAL(DP)                                    :: A1,A2    ! stage - different sides of break
 REAL(DP)                                    :: CM       ! merged celerity 
 REAL(DP)                                    :: TEXIT    ! expected exit time of "current" particle
 REAL(DP)                                    :: TNEXT    ! expected exit time of "next" particle
 REAL(DP)                                    :: TEXIT2   ! exit time of "bottom" of merged element
 INTEGER(I4B)                                :: IROUTE   ! looping variable for routing
 INTEGER(I4B)                                :: JROUTE   ! looping variable for routing
 INTEGER(I4B)                                :: ICOUNT   ! used to account for merged pts
 character(len=256)                          :: cmessage ! error message of downwind routine
 ! ----------------------------------------------------------------------------------------
 ! NOTE: If merged particles DO NOT exit the reach in the current time step, they are 
 !       disaggregated into the original particles; if the merged particles DO exit the
 !       reach, then we save only the "slowest" and "fastest" particle.
 ! ----------------------------------------------------------------------------------------
 !       To disaggregate particles we need to keep track of the output element for which
 !       each input element is merged. This is done using the integer vector MF:  If, for
 !       example, MF = (1,2,2,2,3,4,5,5,5,5,6,7,8), this means the 2nd, 3rd, and 4th input
 !       particles have been merged into the 2nd particle of the output, and that the 7th,
 !       8th, 9th, and 10th input particles have been merged into the 5th particle of the
 !       output.  We only store the "slowest" and "fastest" particle within the merged set.
 ! ----------------------------------------------------------------------------------------
 !       Disaggregating the particles proceeds as follows:  If particles have been merged,
 !       then the flow in Q1 will be different from the flow in Q2, that is, continuing with
 !       the above example, Q1(IROUTE).NE.Q2(IROUTE), where IROUTE = 2 or 5.  In the case of
 !       a merged particle we identify all elements of MF that are equal to IROUTE (that is,
 !       the 2nd, 3rd, and 4th elements of MF = 2), populate the output vector with the
 !       selected elements (2,3,4) of the input vector.
 ! ----------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='KINWAV_RCH/'
 ! Get the reach parameters
 ALFA = 5._dp/3._dp        ! should this be initialized here or in a parameter file?
 K    = SQRT(RPARAM(JRCH)%R_SLOPE)/RPARAM(JRCH)%R_MAN_N
 XMX  = RPARAM(JRCH)%RLENGTH
 ! Identify the number of points to route
 NN = SIZE(Q1)                                ! modified when elements are merged
 NI = NN                                      ! original size of the input
 IF(NN.EQ.0) RETURN                           ! don't do anything if no points in the reach
 ! Initialize the vector that indicates which output element the input elements are merged
 MF = arth(1,1,NI)                            ! Num. Rec. intrinsic: see MODULE nrutil.f90
 ! Initialize the vector that indicates the minumum index of each merged element
 IX = arth(1,1,NI)                            ! Num. Rec. intrinsic: see MODULE nrutil.f90
 ! Get copies of the flow/time particles
 Q0=Q_JRCH; Q1=Q_JRCH; Q2=Q_JRCH
 T0=TENTRY; T1=TENTRY; T2=TENTRY
 ! compute wave celerity for all flow points (array operation)
 WC(1:NN) = ALFA*K**(1./ALFA)*Q1(1:NN)**((ALFA-1.)/ALFA)
 GT_ONE: IF(NN.GT.1) THEN                     ! no breaking if just one point
  X = 0.                                      ! altered later to describe "closest" shock
  GOTALL: DO                                  ! keep going until all shocks are merged 
   XB = XMX                                   ! initialized to length of the stream segment
   ! --------------------------------------------------------------------------------------
   ! check for breaking
   ! --------------------------------------------------------------------------------------
   WCHECK: DO IW=2,NN
    JW=IW-1
    IF(WC(IW).EQ.0. .OR. WC(JW).EQ.0.) CYCLE  ! waves not moving
    WDIFF = 1./WC(JW) - 1./WC(IW)             ! difference in wave celerity
    IF(WDIFF.EQ.0.) CYCLE                     ! waves moving at the same speed
    IF(WC(IW).EQ.WC(JW)) CYCLE                ! identical statement to the above?
    XXB = (T1(IW)-T1(JW)) / WDIFF             ! XXB is point of breaking in x direction
    IF(XXB.LT.X .OR. XXB.GT.XB) CYCLE         ! XB init at LENGTH, so > XB do in next reach
    ! if get to here, the wave is breaking
    XB  = XXB                                 ! identify break "closest to upstream" first
    IXB = IW
   END DO WCHECK
   ! --------------------------------------------------------------------------------------
   IF(XB.EQ.XMX) EXIT                         ! got all breaking waves, exit gotall
   ! --------------------------------------------------------------------------------------
   ! combine waves
   ! --------------------------------------------------------------------------------------
   NN  = NN-1
   JXB = IXB-1                                ! indices for the point of breaking
   NM  = NI-NN                                ! number of merged elements
   ! calculate merged shockwave celerity (CM) using finite-difference approximation
   Q2(JXB) =MAX(Q2(JXB),Q2(IXB))              ! flow of largest merged point
   Q1(JXB) =MIN(Q1(JXB),Q1(IXB))              ! flow of smallest merged point
   A2 = (Q2(JXB)/K)**(1./ALFA)                ! Q = (1./MAN_N) H**(ALFA) sqrt(SLOPE)
   A1 = (Q1(JXB)/K)**(1./ALFA)                ! H = (Q/K)**(1./ALFA) (K=sqrt(SLOPE)/MAN_N)
   CM = (Q2(JXB)-Q1(JXB))/(A2-A1)             ! NB:  A1,A2 are river stage
   ! update merged point
   T1(JXB) = T1(JXB) + XB/WC(JXB) - XB/CM     ! updated starting point
   WC(JXB) = CM
   ! if input elements are merged, then reduce index of merged element plus all remaining elements
   MF(IX(IXB):NI) = MF(IX(IXB):NI)-1         ! NI is the original size of the input
   ! re-number elements, ommitting the element just merged
   IX(IXB:NN) = IX(IXB+1:NN+1)               ! index (minimum index value of each merged particle)
   T1(IXB:NN) = T1(IXB+1:NN+1)               ! entry time
   WC(IXB:NN) = WC(IXB+1:NN+1)               ! wave celerity
   Q1(IXB:NN) = Q1(IXB+1:NN+1)               ! unmodified flows 
   Q2(IXB:NN) = Q2(IXB+1:NN+1)               ! unmodified flows
   ! update X - already got the "closest shock to start", see if there are any other shocks
   X = XB
   ! --------------------------------------------------------------------------------------
  END DO GOTALL
 ENDIF GT_ONE
 !
 ICOUNT=0
 ! ----------------------------------------------------------------------------------------
 ! perform the routing
 ! ----------------------------------------------------------------------------------------
 DO IROUTE = 1,NN    ! loop through the remaining particles (shocks,waves) (NM=NI-NN have been merged)
  ! check that we have non-zero flow
  if(WC(IROUTE) < verySmall)then
   ierr=20; message=trim(message)//'zero flow'; return
  endif
  ! compute the time the shock will exit the reach
  TEXIT = MIN(XMX/WC(IROUTE) + T1(IROUTE), HUGE(T1))
  ! compute the time the next shock will exit the reach
  IF (IROUTE.LT.NN) TNEXT = MIN(XMX/WC(IROUTE+1) + T1(IROUTE+1), HUGE(T1))
  IF (IROUTE.EQ.NN) TNEXT = HUGE(T1)
  ! check if element is merged
  MERGED: IF(Q1(IROUTE).NE.Q2(IROUTE)) THEN
   ! check if merged element has exited
   IF(TEXIT.LT.T_END) THEN
    ! when a merged element exits, save just the top and the bottom of the shock
    ! (identify the exit time for the "slower" particle)
    TEXIT2 = MIN(TEXIT+1.0D0, TEXIT + 0.5D0*(MIN(TNEXT,T_END)-TEXIT))
    ! unsure what will happen in the rare case if TEXIT and TEXIT2 are the same
    IF (TEXIT2.EQ.TEXIT) THEN
     ierr=30; message=trim(message)//'TEXIT equals TEXIT2 in kinwav'; return
    ENDIF
    ! fill output arrays
    CALL RUPDATE(Q1(IROUTE),T1(IROUTE),TEXIT,ierr,cmessage)    ! fill arrays w/ Q1, T1, + run checks
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    CALL RUPDATE(Q2(IROUTE),T1(IROUTE),TEXIT2,ierr,cmessage)   ! fill arrays w/ Q2, T1, + run checks
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   ELSE                                      ! merged elements have not exited
    ! when a merged element does not exit, need to disaggregate into original particles
    DO JROUTE=1,NI                           ! loop thru # original inputs    
     IF(MF(JROUTE).EQ.IROUTE) &
      CALL RUPDATE(Q0(JROUTE),T0(JROUTE),TEXIT,ierr,cmessage)  ! fill arrays w/ Q0, T0, + run checks
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    END DO  ! JROUTE
   ENDIF   ! TEXIT
  ! now process un-merged particles
  ELSE MERGED  ! (i.e., not merged)
   CALL RUPDATE(Q1(IROUTE),T1(IROUTE),TEXIT,ierr,cmessage)     ! fill arrays w/ Q1, T1, + run checks
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  ENDIF MERGED
 END DO
 ! update arrays
 NQ2 = ICOUNT
 RETURN
 ! ----------------------------------------------------------------------------------------
 contains

  subroutine RUPDATE(QNEW,TOLD,TNEW,ierr,message)
    REAL(DP),INTENT(IN)                        :: QNEW      ! Q0,Q1, or Q2
    REAL(DP),INTENT(IN)                        :: TOLD,TNEW ! entry/exit times
    integer(i4b), intent(out)                  :: ierr      ! error code
    character(*), intent(out)                  :: message   ! error message
    ! initialize error control
    ierr=0; message='RUPDATE/'
    ! ---------------------------------------------------------------------------------------
    ! Used to compute the time each element will exit stream segment & update routing flag
    ! NB: internal subroutine so all data from host is available
    ! ---------------------------------------------------------------------------------------
    ICOUNT=ICOUNT+1
    ! check for array bounds exceeded
    IF (ICOUNT.GT.SIZE(Q_JRCH)) THEN
     ierr=60; message=trim(message)//'array bounds exceeded'; return
    ENDIF
    ! fill output arrays
    Q_JRCH(ICOUNT) = QNEW                         ! flow (Q1 always smaller than Q2)
    TENTRY(ICOUNT) = TOLD                         ! time - note, T1 altered if element merged
    T_EXIT(ICOUNT) = TNEW
    ! time check -- occurs when disaggregating merged elements
    IF (ICOUNT.GT.1) THEN
     IF (T_EXIT(ICOUNT).LE.T_EXIT(ICOUNT-1)) T_EXIT(ICOUNT)=T_EXIT(ICOUNT-1)+1.
    ENDIF
    ! another time check -- rare problem when the shock can get the same time as tstart
    IF(ICOUNT.EQ.1.AND.T_EXIT(ICOUNT).LE.T_START) T_EXIT(ICOUNT)=T_START+1.
    ! update flag for routed elements
    IF(T_EXIT(ICOUNT).LT.T_END) FROUTE(ICOUNT) =.TRUE.
    return
  end subroutine

 end subroutine

 ! *********************************************************************
 ! new subroutine: calculate time-step averages from irregular values
 ! *********************************************************************
 subroutine INTERP_RCH(TOLD,QOLD,TNEW,QNEW,IERR,MESSAGE)
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Unknown (original Tideda routine?), fairly old
 !   Martyn Clark, 2006: major clean-up and re-write
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Used to calculate time step averages from a set of irregularly spaced values
 !    (provides one less output data value -- average BETWEEN time steps)
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O:
 !
 !   Input(s):
 !   ---------
 !     TOLD: vector of times
 !     QOLD: vector of flows
 !     TNEW: vector of desired output times
 !
 !   Outputs:
 !   --------
 !     QNEW: vector of interpolated output [dimension = TNEW-1 (average BETWEEN intervals)]
 !     IERR: error code (1 = bad bounds)
 !
 ! ----------------------------------------------------------------------------------------
 ! Structures Used:
 !
 !   NONE
 !
 ! ----------------------------------------------------------------------------------------
 ! Method:
 !
 !   Loop through all output times (can be just 2, start and end of the time step)
 !
 !   For each pair of output times, integrate over the irregularly-spaced data
 !   values.
 !
 !   First, estimate the area under the curve at the start and end of the time step.
 !   Identify index values IBEG and IEND that span the start and end of the time
 !   step: TOLD(IBEG-1) < T0 <= TOLD(IBEG), and TOLD(IEND-1) < T1 <= TOLD(IEND).
 !   Use linear interpolation to estimate the value at T0 (QEST0) and T1 (QEST1).
 !   The area under the curve at the start/end of the time step is then simply
 !   (TOLD(IBEG) - T0) and (T1 - TOLD(IEND-1)), and the corresponding data values
 !   are 0.5*(QEST0 + QOLD(IBEG)) and 0.5*(QOLD(IEND-1) + QEST1).
 !
 !   Second, integrate over the remaining datapoints between T0 and T1.  The area
 !   for the middle points is (TOLD(i) - TOLD(i-1)), and the corresponding data
 !   values are (QOLD(i) - QOLD(i-1)).
 !
 !   Finally, divide the sum of the data values * area by the sum of the areas
 !   (which in this case is simply T1-T0).
 !
 !   A special case occurs when both the start and the end points of the time step
 !   lie between two of the irregularly-spaced input data values.  Here, use linear
 !   interpolation to estimate values at T0 and T1, and simply average the
 !   interpolated values.
 !
 !   NB: values of TOLD with times <= T0 and >= T1 are required for reliable
 !   estimates of time step averages.
 !
 ! ----------------------------------------------------------------------------------------
 ! Source:
 !
 !   This routine is base on the subroutine interp, in file kinwav_v7.f
 !
 ! ----------------------------------------------------------------------------------------
 ! Modifications to source (mclark@ucar.edu):
 !
 !   * All variables are now defined (IMPLICIT NONE) and described (comments)
 !
 !   * Added extra comments
 !
 !   * Replaced GOTO statements with DO loops and IF statements
 !
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! --------------------------------------------------------------------------------------------
 IMPLICIT NONE
 ! Input
 REAL(DP), DIMENSION(:), INTENT(IN)          :: TOLD     ! input time array
 REAL(DP), DIMENSION(:), INTENT(IN)          :: QOLD     ! input flow array
 REAL(DP), DIMENSION(:), INTENT(IN)          :: TNEW     ! desired output times
 ! Output
 REAL(DP), DIMENSION(:), INTENT(OUT)         :: QNEW     ! flow averaged for desired times
 INTEGER(I4B), INTENT(OUT)                   :: IERR     ! error, 1= bad bounds
 character(*), intent(out)                   :: MESSAGE  ! error message
 ! Internal
 INTEGER(I4B)                                :: NOLD     ! number of elements in input array
 INTEGER(I4B)                                :: NNEW     ! number of desired new times
 INTEGER(I4B)                                :: IOLDLOOP ! loop through input times
 INTEGER(I4B)                                :: INEWLOOP ! loop through desired times
 REAL(DP)                                    :: T0,T1    ! time at start/end of the time step
 INTEGER(I4B)                                :: IBEG     ! identify input times spanning T0
 INTEGER(I4B)                                :: IEND     ! identify input times spanning T1
 INTEGER(I4B)                                :: IMID     ! input times in middle of the curve
 REAL(DP)                                    :: AREAB    ! area at the start of the time step
 REAL(DP)                                    :: AREAE    ! area at the end of the time step
 REAL(DP)                                    :: AREAM    ! area at the middle of the time step
 REAL(DP)                                    :: AREAS    ! sum of all areas
 REAL(DP)                                    :: SLOPE    ! slope between two input data values
 REAL(DP)                                    :: QEST0    ! flow estimate at point T0
 REAL(DP)                                    :: QEST1    ! flow estimate at point T1
 ! --------------------------------------------------------------------------------------------
 ! get array size
 NOLD = SIZE(TOLD); NNEW = SIZE(TNEW)

 !
 IERR=0
 message='INTERP_RCH/'
 ! check that the input time series starts before the first required output time
 ! and ends after the last required output time
 IF( (TOLD(1).GT.TNEW(1)) .OR. (TOLD(NOLD).LT.TNEW(NNEW)) ) THEN
  IERR=1; message=trim(message)//'bad bounds'; RETURN
 ENDIF
 !
 ! loop through the output times
 DO INEWLOOP=2,NNEW
  !
  T0 = TNEW(INEWLOOP-1)                      ! start of the time step
  T1 = TNEW(INEWLOOP)                        ! end of the time step
  !
  IBEG=1
  ! identify the index values that span the start of the time step
  BEG_ID: DO IOLDLOOP=2,NOLD
   IF(T0.LE.TOLD(IOLDLOOP)) THEN
    IBEG = IOLDLOOP
    EXIT
   ENDIF
  END DO BEG_ID
  !
  IEND=1
  ! identify the index values that span the end of the time step
  END_ID: DO IOLDLOOP=1,NOLD
   IF(T1.LE.TOLD(IOLDLOOP)) THEN
    IEND = IOLDLOOP
    EXIT
   ENDIF
  END DO END_ID
  !
  ! initialize the areas
  AREAB=0D0; AREAE=0D0; AREAM=0D0
  !
  ! special case: both TNEW(INEWLOOP-1) and TNEW(INEWLOOP) are within two original values
  ! (implies IBEG=IEND) -- estimate values at both end-points and average
  IF(T1.LT.TOLD(IBEG)) THEN
   SLOPE = (QOLD(IBEG)-QOLD(IBEG-1))/(TOLD(IBEG)-TOLD(IBEG-1))
   QEST0 = SLOPE*(T0-TOLD(IBEG-1)) + QOLD(IBEG-1)
   QEST1 = SLOPE*(T1-TOLD(IBEG-1)) + QOLD(IBEG-1)
   QNEW(INEWLOOP-1) = 0.5*(QEST0 + QEST1)
   CYCLE ! loop back to the next desired time
  ENDIF
  !
  ! estimate the area under the curve at the start of the time step
  IF(T0.LT.TOLD(IBEG)) THEN  ! if equal process as AREAM
   SLOPE = (QOLD(IBEG)-QOLD(IBEG-1))/(TOLD(IBEG)-TOLD(IBEG-1))
   QEST0 = SLOPE*(T0-TOLD(IBEG-1)) + QOLD(IBEG-1)
   AREAB = (TOLD(IBEG)-T0) * 0.5*(QEST0 + QOLD(IBEG))
  ENDIF
  !
  ! estimate the area under the curve at the end of the time step
  IF(T1.LT.TOLD(IEND)) THEN  ! if equal process as AREAM
   SLOPE = (QOLD(IEND)-QOLD(IEND-1))/(TOLD(IEND)-TOLD(IEND-1))
   QEST1 = SLOPE*(T1-TOLD(IEND-1)) + QOLD(IEND-1)
   AREAE = (T1-TOLD(IEND-1)) * 0.5*(QOLD(IEND-1) + QEST1)
  ENDIF
  !
  ! check if there are extra points to process
  IF(IBEG.LT.IEND) THEN
   ! loop through remaining points
   DO IMID=IBEG+1,IEND
    IF(IMID.LT.IEND .OR. &
      ! process the end slice as AREAM, but only if not already AREAB
      (IMID.EQ.IEND.AND.T1.EQ.TOLD(IEND).AND.T0.LT.TOLD(IEND-1)) ) THEN
       ! compute AREAM
       AREAM = AREAM + (TOLD(IMID) - TOLD(IMID-1)) * 0.5*(QOLD(IMID-1) + QOLD(IMID))
    ENDIF   ! if point is valid
   END DO  ! IMID
  ENDIF   ! If there is a possibility that middle points even exist
  !
  ! compute time step average
  AREAS = AREAB + AREAE + AREAM            ! sum of all areas
  QNEW(INEWLOOP-1) = AREAS / (T1-T0)       ! T1-T0 is the sum of all time slices
  !
 END DO
 return
 end subroutine

end module kwt_route
