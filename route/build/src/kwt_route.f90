MODULE kwt_route_module

!numeric type
USE nrtype
! data types
USE dataTypes, ONLY: FPOINT            ! particle
USE dataTypes, ONLY: STRFLX            ! fluxes in each reach
USE dataTypes, ONLY: STRSTA            ! states in each reach
USE dataTypes, ONLY: RCHTOPO           ! Network topology
USE dataTypes, ONLY: RCHPRP            ! Reach parameter
USE dataTypes, ONLY: kwtRCH            ! kwt specific state data structure 
! global data
USE public_var, ONLY: runoffMin         ! minimum runoff
USE public_var, ONLY: verySmall         ! a very small value
USE public_var, ONLY: realMissing       ! missing value for real number
USE public_var, ONLY: integerMissing    ! missing value for integer number
USE globalData, ONLY: idxKWT            ! index of KWT method 
! utilities
USE nr_utility_module, ONLY: arth       ! Num. Recipies utilities

implicit none

private

public::kwt_route

CONTAINS

 ! *********************************************************************
 ! subroutine: route kinematic waves through the river network
 ! *********************************************************************
 SUBROUTINE kwt_route(iens,                 & ! input: ensemble index
                      river_basin,          & ! input: river basin information (mainstem, tributary outlet etc.)
                      T0,T1,                & ! input: start and end of the time step
                      ixDesire,             & ! input: reachID to be checked by on-screen pringing
                      NETOPO_in,            & ! input: reach topology data structure
                      RPARAM_in,            & ! input: reach parameter data structure
                      RCHSTA_out,           & ! inout: reach state data structure
                      RCHFLX_out,           & ! inout: reach flux data structure
                      ierr,message,         & ! output: error control
                      ixSubRch)               ! optional input: subset of reach indices to be processed

   USE dataTypes,      ONLY : subbasin_omp        ! mainstem+tributary data strucuture
   USE model_finalize, ONLY : handle_err

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
   ! variables needed for timing
! integer(i4b)                             :: omp_get_thread_num
! integer(i4b), allocatable                :: ixThread(:)           ! thread id
! integer*8,    allocatable                :: openMPend(:)          ! time for the start of the parallelization section
! integer*8,    allocatable                :: timeTribStart(:)      ! time Tributaries start
! real(dp),     allocatable                :: timeTrib(:)           ! time spent on each Tributary

   ierr=0; message='kwt_route/'

   ! number of reach check
   if (size(NETOPO_in)/=size(RCHFLX_out(iens,:))) then
    ierr=20; message=trim(message)//'sizes of NETOPO and RCHFLX mismatch'; return
   endif

   nSeg = size(NETOPO_in)

   allocate(doRoute(nSeg), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating space for [doRoute]'; return; endif

   if (present(ixSubRch))then
    doRoute(:)=.false.
    doRoute(ixSubRch) = .true. ! only subset of reaches are on
   else
    doRoute(:)=.true. ! every reach is on
   endif

   nOrder = size(river_basin)

   do ix = 1, nOrder

     nTrib=size(river_basin(ix)%branch)

!   allocate(ixThread(nTrib), openMPend(nTrib), timeTrib(nTrib), timeTribStart(nTrib), stat=ierr)
!   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to allocate space for Trib timing'; return; endif
!   timeTrib(:) = realMissing
!   ixThread(:) = integerMissing

   ! 1. Route tributary reaches (parallel)
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
!!$OMP          shared(openMPend, nThreads)              & ! timing variables shared
!!$OMP          shared(timeTribStart)                    & ! timing variables shared
!!$OMP          shared(timeTrib)                         & ! timing variables shared
!!$OMP          shared(ixThread)                         & ! thread id array shared
     trib:do iTrib = 1,nTrib
!!$    ixThread(iTrib) = omp_get_thread_num()
!    call system_clock(timeTribStart(iTrib))
       seg:do iSeg=1,river_basin(ix)%branch(iTrib)%nRch
         jSeg  = river_basin(ix)%branch(iTrib)%segIndex(iSeg)
         if (.not. doRoute(jSeg)) cycle
         call qroute_rch(iEns,jSeg,           & ! input: array indices
                         ixDesire,            & ! input: index of verbose reach
                         T0,T1,               & ! input: start and end of the time step
                         LAKEFLAG,            & ! input: flag if lakes are to be processed
                         NETOPO_in,           & ! input: reach topology data structure
                         RPARAM_in,           & ! input: reach parameter data structure
                         RCHSTA_out,          & ! inout: reach state data structure
                         RCHFLX_out,          & ! inout: reach flux data structure
                         ierr,cmessage)         ! output: error control
         if(ierr/=0) call handle_err(ierr, trim(message)//trim(cmessage))
       end do seg
!     call system_clock(openMPend(iTrib))
!     timeTrib(iTrib) = real(openMPend(iTrib)-timeTribStart(iTrib), kind(dp))
     end do trib
!$OMP END PARALLEL DO

!   write(*,'(a)') 'iTrib nSeg ixThread nThreads StartTime EndTime'
!   do iTrib=1,nTrib
!     write(*,'(4(i5,1x),2(I20,1x))') iTrib, river_basin(ix)%branch(iTrib)%nRch, ixThread(iTrib), nThreads, timeTribStart(iTrib), openMPend(iTrib)
!   enddo
!   deallocate(ixThread, openMPend, timeTrib, timeTribStart, stat=ierr)
!   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': unable to deallocate space for Trib timing'; return; endif

   end do ! basin loop

 END SUBROUTINE kwt_route


 ! *********************************************************************
 ! subroutine: route kinematic waves at one segment
 ! *********************************************************************
 SUBROUTINE qroute_rch(IENS,JRCH,    & ! input: array indices
                       ixDesire,     & ! input: index of the reach for verbose output
                       T0,T1,        & ! input: start and end of the time step
                       LAKEFLAG,     & ! input: flag if lakes are to be processed
                       NETOPO_in,    & ! input: reach topology data structure
                       RPARAM_in,    & ! input: reach parameter data structure
                       RCHSTA_out,   & ! inout: reach state data structure
                       RCHFLX_out,   & ! inout: reach flux data structure
                       ierr,message, & ! output: error control
                       RSTEP)          ! optional input: retrospective time step offset
 USE public_var, ONLY : MAXQPAR        ! maximum number of waves per reach
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
 !   * all variables are defined (implicit none) and described (comments)
 !
 !   * use of a new data structure (RCHSTA_out) to hold and update the flow particles
 !
 !   * upgrade to F90 (especially structured variables and dynamic memory allocation)
 !
 ! ----------------------------------------------------------------------------------------
   implicit none
   ! Input
   integer(i4b), intent(in)                    :: IENS          ! ensemble member
   integer(i4b), intent(in)                    :: JRCH          ! reach to process
   integer(i4b), intent(in)                    :: ixDesire      ! index of the reach for verbose output
   real(dp),     intent(in)                    :: T0,T1         ! start and end of the time step (seconds)
   integer(i4b), intent(in)                    :: LAKEFLAG      ! >0 if processing lakes
   type(RCHTOPO),intent(in),    allocatable    :: NETOPO_in(:)  ! River Network topology
   type(RCHPRP), intent(in),    allocatable    :: RPARAM_in(:)  ! River reach parameter
   integer(i4b), intent(in), optional          :: RSTEP         ! retrospective time step offset
   ! inout
   type(STRSTA), intent(inout), allocatable    :: RCHSTA_out(:,:) ! reach state data
   type(STRFLX), intent(inout), allocatable    :: RCHFLX_out(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   ! output variables
   integer(i4b), intent(out)                   :: ierr          ! error code
   character(*), intent(out)                   :: message       ! error message
   ! (1) extract flow from upstream reaches and append to the non-routed flow in JRCH
   integer(i4b)                                :: NUPS          ! number of upstream reaches
   real(dp),dimension(:),allocatable           :: Q_JRCH        ! flow in downstream reach JRCH
   real(dp),dimension(:),allocatable           :: TENTRY        ! entry time to JRCH (exit time u/s)
   integer(i4b)                                :: NQ1           ! # flow particles
   ! (2) route flow within the current [JRCH] river segment
   integer(i4b)                                :: ROFFSET       ! retrospective offset due to rstep
   real(dp)                                    :: T_START       ! start of time step
   real(dp)                                    :: T_END         ! end of time step
   real(dp),dimension(:),allocatable           :: T_EXIT        ! time particle expected exit JRCH
   logical(LGT),dimension(:),allocatable       :: FROUTE        ! routing flag .T. if particle exits
   integer(i4b)                                :: NQ2           ! # flow particles (<=NQ1 b/c merge)
   ! (3) calculate time-step averages
   integer(i4b)                                :: NR            ! # routed particles
   integer(i4b)                                :: NN            ! # non-routed particles
   real(dp),dimension(2)                       :: TNEW          ! start/end of time step
   real(dp),dimension(1)                       :: QNEW          ! interpolated flow
   ! (4) housekeeping
   real(dp)                                    :: Q_END         ! flow at the end of the timestep
   real(dp)                                    :: TIMEI         ! entry time at the end of the timestep
   TYPE(FPOINT),allocatable,dimension(:)       :: NEW_WAVE      ! temporary wave
   ! random stuff
   integer(i4b)                                :: IWV           ! rech index
   character(len=strLen)                       :: fmt1,fmt2     ! format string
   character(len=strLen)                       :: CMESSAGE      ! error message for downwind routine

   ierr=0; message='qroute_rch/'

   if(JRCH==ixDesire) then
     write(*,'(2a)') new_line('a'),'** Check kinematic wave tracking routing **'
     write(*,"(a,x,I10,x,I10)")      ' Reach index & ID  =', JRCH, NETOPO_in(JRCH)%REACHID
     write(*,"(a,x,F20.7,1x,F20.7)") ' time step(T0,T1)  =', T0, T1
     write(*,'(a,x,F15.7)')          ' RPARAM_in%R_SLOPE =', RPARAM_in(JRCH)%R_SLOPE
     write(*,'(a,x,F15.7)')          ' RPARAM_in%R_MAN_N =', RPARAM_in(JRCH)%R_MAN_N
     write(*,'(a,x,F15.7)')          ' RPARAM_in%R_WIDTH =', RPARAM_in(JRCH)%R_WIDTH
   end if

   RCHFLX_out(IENS,JRCH)%TAKE=0.0_dp ! initialize take from this reach

    ! ----------------------------------------------------------------------------------------
    ! (1) EXTRACT FLOW FROM UPSTREAM REACHES & APPEND TO THE NON-ROUTED FLOW PARTICLES IN JRCH
    ! ----------------------------------------------------------------------------------------
    NUPS = count(NETOPO_in(JRCH)%goodBas)        ! number of desired upstream reaches
    !NUPS = size(NETOPO_in(JRCH)%UREACHI)        ! number of upstream reaches
    if (NUPS.GT.0) then
      call getusq_rch(IENS,JRCH,LAKEFLAG,T0,T1,ixDesire, & ! input
                      NETOPO_in,RPARAM_in,RCHFLX_out,    & ! input
                      RCHSTA_out,                        & ! inout
                      Q_JRCH,TENTRY,T_EXIT,ierr,cmessage,& ! output
                      RSTEP)                               ! optional input
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      if (minval(Q_JRCH).lt.0.0_dp) then
        ierr=20; message=trim(message)//'negative flow extracted from upstream reach'; return
      endif

      if(JRCH==ixDesire)then
        write(fmt1,'(A,I5,A)') '(A, 1X',size(Q_JRCH),'(1X,F20.7))'
        write(*,'(a)') ' * Wave discharge from upstream reaches (Q_JRCH) [m2/s]:'
        write(*,fmt1)  ' Q_JRCH=',(Q_JRCH(IWV), IWV=0,size(Q_JRCH)-1)
      endif
    else
      ! set flow in headwater reaches to modelled streamflow from time delay histogram
      RCHFLX_out(IENS,JRCH)%ROUTE(idxKWT)%REACH_Q = RCHFLX_out(IENS,JRCH)%BASIN_QR(1)
      if (allocated(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE)) then
        deallocate(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE,STAT=IERR)
        if(ierr/=0)then; message=trim(message)//'problem deallocating space for RCHSTA_out'; return; endif
      endif
      allocate(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:0),STAT=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating space for RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(1)'; return; endif
      RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%QF=-9999
      RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%TI=-9999
      RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%TR=-9999
      RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%RF=.False.
      RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%QM=-9999

      if(JRCH==ixDesire) then
        write(*,'(a)') ' * Final discharge (RCHFLX_out(IENS,JRCH)%REACH_Q) [m3/s]:'
        write(*,'(x,F20.7)') RCHFLX_out(IENS,JRCH)%ROUTE(idxKWT)%REACH_Q
      end if
      return  ! no upstream reaches (routing for sub-basins done using time-delay histogram)
    endif

    ! ----------------------------------------------------------------------------------------
    ! (2) REMOVE FLOW PARTICLES (REDUCE MEMORY USAGE AND PROCESSING TIME)
    ! ----------------------------------------------------------------------------------------
    if (size(Q_JRCH).GT.MAXQPAR) then
      call remove_rch(MAXQPAR,Q_JRCH,TENTRY,T_EXIT,ierr,cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif
    NQ1 = size(Q_JRCH)-1                                     ! -1 because of the zero element

    ! ----------------------------------------------------------------------------------------
    ! (x) Water use - take out (Qtake is negative)
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

    if (RPARAM_in(jrch)%QTAKE < 0._dp) then
      call extract_from_rch(iens, jrch,                       & ! input: ensemble and reach indices
                            T_START, T_END,                   & ! input: time [sec] of current time step bounds
                            RPARAM_in,                        & ! input: river reach parameters
                            RPARAM_in(jrch)%QTAKE,            & ! input: target Qtake (minus)
                            ixDesire,                         & ! input:
                            Q_JRCH, T_EXIT, TENTRY,           & ! inout: discharge and exit time for particle
                            ierr,cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end if

    ! ----------------------------------------------------------------------------------------
    ! (3) ROUTE FLOW WITHIN THE CURRENT [JRCH] RIVER SEGMENT
    ! ----------------------------------------------------------------------------------------
    allocate(FROUTE(0:NQ1),STAT=IERR)
    if(ierr/=0)then; message=trim(message)//'problem allocating space for FROUTE'; return; endif
    FROUTE(0) = .true.; FROUTE(1:NQ1)=.false.  ! init. routing flags

    ! route flow through the current [JRCH] river segment (Q_JRCH in units of m2/s)
    call kinwav_rch(JRCH,T_START,T_END,ixDesire,                             & ! input: location and time
                    NETOPO_in, RPARAM_in,                                    & ! input: river data structure
                    Q_JRCH(1:NQ1),TENTRY(1:NQ1),T_EXIT(1:NQ1),FROUTE(1:NQ1), & ! inout: kwt states
                    NQ2,ierr,cmessage)                                         ! output:
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    if(JRCH == ixDesire)then
      write(fmt1,'(A,I5,A)') '(A,1X',NQ1+1,'(1X,F20.7))'
      write(fmt2,'(A,I5,A)') '(A,1X',NQ1+1,'(1X,L))'
      write(*,'(a)') ' * After routed: wave discharge (Q_JRCH) [m2/s], isExit(FROUTE), entry time (TENTRY) [s], and exit time (T_EXIT) [s]:'
      write(*,fmt1)  ' Q_JRCH=',(Q_JRCH(IWV), IWV=0,NQ1)
      write(*,fmt1)  ' TENTRY=',(TENTRY(IWV), IWV=0,NQ1)
      write(*,fmt1)  ' T_EXIT=',(T_EXIT(IWV), IWV=0,NQ1)
      write(*,fmt2)  ' FROUTE=',(FROUTE(IWV), IWV=0,NQ1)
    end if

    ! ----------------------------------------------------------------------------------------
    ! (4) COMPUTE TIME-STEP AVERAGES
    ! ----------------------------------------------------------------------------------------
    NR = count(FROUTE)-1   ! -1 because of the zero element (last routed)
    NN = NQ2-NR            ! number of non-routed points
    TNEW = (/T_START,T_END/)

    ! (zero position last routed; use of NR+1 instead of NR keeps next expected routed point)
    call interp_rch(T_EXIT(0:NR+1),Q_JRCH(0:NR+1),TNEW,QNEW,IERR,CMESSAGE)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! m2/s --> m3/s + instantaneous runoff from basin
    RCHFLX_out(IENS,JRCH)%ROUTE(idxKWT)%REACH_Q = QNEW(1)*RPARAM_in(JRCH)%R_WIDTH + RCHFLX_out(IENS,JRCH)%BASIN_QR(1)

    if(JRCH == ixDesire)then
      write(*,'(a)')          ' * Time-ave. wave discharge that exit (QNEW(1)) [m2/s], local-area discharge (RCHFLX_out%BASIN_QR(1)) [m3/s] and Final discharge (RCHFLX_out%REACH_Q) [m3/s]:'
      write(*,"(A,1x,F15.7)") ' QNEW(1)                =', QNEW(1)
      write(*,"(A,1x,F15.7)") ' RCHFLX_out%BASIN_QR(1) =', RCHFLX_out(IENS,JRCH)%BASIN_QR(1)
      write(*,"(A,1x,F15.7)") ' RCHFLX_out%REACH_Q     =', RCHFLX_out(IENS,JRCH)%ROUTE(idxKWT)%REACH_Q
    endif

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
    if (.not.allocated(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE)) then
      ierr=20; message=trim(message)//'RCHSTA_out is not associated'; return
    else
      deallocate(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE, STAT=ierr)
      if(ierr/=0)then; message=trim(message)//'problem deallocating space for RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE'; return; endif
      allocate(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NQ2+1),STAT=ierr)   ! NQ2 is number of points for kinematic routing
      if(ierr/=0)then; message=trim(message)//'problem allocating space for RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NQ2+1)'; return; endif
    endif
    ! insert the interpolated point (TI is irrelevant, as the point is "routed")
    RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+1)%QF=Q_END;   RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+1)%TI=TIMEI
    RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+1)%TR=T_END;   RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+1)%RF=.true.
    ! add the output from kinwave...         - skip NR+1
    ! (when JRCH becomes IR routed points will be stripped out & the structures updated again)
    RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NR)%QF=Q_JRCH(0:NR); RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+2:NQ2+1)%QF=Q_JRCH(NR+1:NQ2)
    RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NR)%TI=TENTRY(0:NR); RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+2:NQ2+1)%TI=TENTRY(NR+1:NQ2)
    RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NR)%TR=T_EXIT(0:NR); RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+2:NQ2+1)%TR=T_EXIT(NR+1:NQ2)
    RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NR)%RF=FROUTE(0:NR); RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+2:NQ2+1)%RF=FROUTE(NR+1:NQ2)
    RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NQ2+1)%QM=-9999

    ! implement water use
    !IF (NUSER.GT.0.AND.UCFFLAG.GE.1) THEN
      !CALL EXTRACT_FROM_RCH(IENS,JRCH,NR,Q_JRCH,T_EXIT,T_END,TNEW)
    !ENDIF

    ! free up space for the next reach
    deallocate(Q_JRCH,TENTRY,T_EXIT,FROUTE,STAT=IERR)   ! FROUTE defined in this sub-routine
    if(ierr/=0)then; message=trim(message)//'problem deallocating space for [Q_JRCH, TENTRY, T_EXIT, FROUTE]'; return; endif

    ! ***
    ! remove flow particles from the most downstream reach
    ! if the last reach or lake inlet (and lakes are enabled), remove routed elements from memory
    if ((NETOPO_in(JRCH)%DREACHK<=0 ).OR. &  ! if the last reach (down reach ID:DREACHK is negative), then there is no downstream reach
        (LAKEFLAG.EQ.1.AND.NETOPO_in(JRCH)%LAKINLT)) then ! if lake inlet
      ! copy data to a temporary wave
      if (allocated(NEW_WAVE)) then
        deallocate(NEW_WAVE,STAT=IERR)
        if(ierr/=0)then; message=trim(message)//'problem deallocating space for NEW_WAVE'; return; endif
      endif
      allocate(NEW_WAVE(0:NN),STAT=IERR)  ! NN = number non-routed (the zero element is the last routed point)
      if(ierr/=0)then; message=trim(message)//'problem allocating space for NEW_WAVE'; return; endif
      NEW_WAVE(0:NN) = RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(NR+1:NQ2+1)  ! +1 because of the interpolated point
      ! re-size wave structure
      if (allocated(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE)) then
        deallocate(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE,STAT=IERR)
        if(ierr/=0)then; message=trim(message)//'problem deallocating space for RCHSTA_out'; return; endif
      endif
      allocate(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NN),STAT=IERR)  ! again, the zero element for the last routed point
      if(ierr/=0)then; message=trim(message)//'problem allocating space for RCHSTA_out'; return; endif
      ! copy data back to the wave structure and deallocate space for the temporary wave
      RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NN) = NEW_WAVE(0:NN)
    endif  ! (if JRCH is the last reach)

  END SUBROUTINE qroute_rch

  ! *********************************************************************
  ! subroutine: wave discharge mod to extract/infect water from the JRCH reach
  ! *********************************************************************
  SUBROUTINE extract_from_rch(iens, jrch,              & ! input: ensemble and reach indices
                              T_START, T_END,          & ! input: start and end time [sec] for this time step
                              RPARAM_in,               & ! input: river reach parameters
                              Qtake,                   & ! input: target Qtake (minus)
                              ixDesire,                & ! input:
                              Q_JRCH, T_EXIT, TENTRY,  & ! inout: discharge and exit time for particle
                              ierr,message)
  implicit none
  ! Input
  integer(i4b),              intent(in)    :: iens          ! ensemble member
  integer(i4b),              intent(in)    :: jRch          ! index of reach to process
  real(dp),                  intent(in)    :: T_START       ! start time [s]
  real(dp),                  intent(in)    :: T_END         ! end time [s]
  type(RCHPRP), allocatable, intent(in)    :: RPARAM_in(:)  ! River reach parameter
  real(dp),                  intent(in)    :: Qtake         ! target Q abstraction [m3/s]
  integer(i4b),              intent(in)    :: ixDesire      ! index of the reach for verbose output
  ! inout
  real(dp),     allocatable, intent(inout) :: Q_JRCH(:)     ! discharge of particle [m2/s] -- discharge for unit channel width
  real(dp),     allocatable, intent(inout) :: T_EXIT(:)     ! time flow is expected to exit JR
  real(dp),     allocatable, intent(inout) :: TENTRY(:)     ! time flow entered JR
  integer(i4b),              intent(out)   :: ierr          ! error code
  character(*),              intent(out)   :: message       ! error message
  ! Local variables
  real(dp)                                 :: totQ          ! total available flow
  real(dp)                                 :: Qabs          ! actual abstracted water
  real(dp)                                 :: Qfrac         ! fraction of target abstraction to total available flow
  real(dp)                                 :: alfa          ! constant = 5/3
  real(dp)                                 :: K             ! sqrt(slope)/mannings N
  real(dp)                                 :: TP(2)         ! start/end of time step
  real(dp)                                 :: Qavg(1)       ! time average flow [m2/s]
  real(dp), allocatable                    :: wc(:)         ! wave celerity [m/s]
  real(dp), allocatable                    :: q_jrch_mod(:) ! wave flow remaining after abstract [m2/s]
  real(dp), allocatable                    :: q_jrch_abs(:) ! wave abstracted flow [m2/s]
  real(dp), allocatable                    :: t_exit_mod(:) ! wave expected exit time [sec]
  integer(i4b)                             :: iw            ! loop index for wave
  integer(i4b)                             :: NR            ! number of particle (not including zero index)
  character(len=strLen)                    :: fmt1          ! format string
  character(len=strLen)                    :: cmessage      ! error message for downwind routine

  ierr=0; message='extract_from_rch/'

  ! uniform flow parameters
  alfa = 5._dp/3._dp
  K    = sqrt(RPARAM_in(JRCH)%R_SLOPE)/RPARAM_in(JRCH)%R_MAN_N

  ! number of waves
  NR = size(Q_JRCH)

  allocate(q_jrch_mod(0:NR-1), t_exit_mod(1:NR-1), wc(1:NR-1))

  ! total "available" discharge in current time step
  ! (zero position last routed; use of NR+1 instead of NR keeps next expected routed point)
  TP = [T_START,T_END]
  call interp_rch(TENTRY(0:NR-1),Q_JRCH(0:NR-1), TP, Qavg, ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  totQ = Qavg(1)*RPARAM_in(JRCH)%R_WIDTH

  if (abs(Qtake) < totQ) then

    ! modified wave Q [m3/s]
    Qfrac      = abs(Qtake)/TotQ    ! fraction of total wave Q to target abstracted Q
    Q_jrch_mod(0) = Q_JRCH(0)
    Q_jrch_mod(1:NR-1) = Q_JRCH(1:NR-1)*(1._dp-Qfrac) ! remaining wave Q after abstraction

  else ! everything taken....

    Q_jrch_mod = runoffMin  ! remaining wave Q after abstraction

  end if

  if(JRCH == ixDesire)then
    ! compute actual abstracted water
    allocate(Q_jrch_abs(0:NR-1))
    Q_jrch_abs = Q_JRCH - Q_jrch_mod
    call interp_rch(TENTRY(0:NR-1),Q_jrch_abs(0:NR-1), TP, Qavg, ierr,cmessage)
    Qabs = Qavg(1)*RPARAM_in(JRCH)%R_WIDTH
    write(*,'(a)')         ' * Target abstraction (Qtake) [m3/s], Available discharge (totQ) [m3/s], Actual abstraction (Qabs) [m3/s] '
    write(*,'(a,x,F15.7)') ' Qtake =', Qtake
    write(*,'(a,x,F15.7)') ' totQ  =', totQ
    write(*,'(a,x,F15.7)') ' Qabs  =', Qabs
  end if

  ! modify wave speed at modified wave discharge and re-compute exit time
  wc = alfa*K**(1._dp/alfa)*Q_jrch_mod**((alfa-1._dp)/alfa) ! modified wave celerity [m/s]
  do iw = 1,NR-1
    t_exit_mod(iw) = min(RPARAM_in(JRCH)%RLENGTH/wc(iw)+TENTRY(iw) , huge(TENTRY(iw)))
  end do

  ! modified final q_jrch and t_exit
  Q_JRCH(0:NR-1) = q_jrch_mod(0:NR-1)
  T_EXIT(1:NR-1) = t_exit_mod(1:NR-1)

  if(JRCH == ixDesire)then
    write(fmt1,'(A,I5,A)') '(A,1X',NR,'(1X,E15.7))'
    write(*,'(a)') ' * After abstracted: wave discharge (Q_JRCH) [m2/s], entry time (TENTRY) [s], and exit time (T_EXIT) [s]:'
    write(*,fmt1)  ' Q_JRCH=',(Q_JRCH(iw), iw=0,NR-1)
    write(*,fmt1)  ' TENTRY=',(TENTRY(iw), iw=0,NR-1)
    write(*,fmt1)  ' T_EXIT=',(T_EXIT(iw), iw=0,NR-1)
  endif

  END SUBROUTINE extract_from_rch


 ! *********************************************************************
 ! subroutine: extract flow from the reaches upstream of JRCH
 ! *********************************************************************
 subroutine getusq_rch(IENS,JRCH,LAKEFLAG,T0,T1,ixDesire, & ! input
                       NETOPO_in,RPARAM_in,RCHFLX_in,     & ! input
                       RCHSTA_out,                        & ! inout
                       Q_JRCH,TENTRY,T_EXIT,ierr,message, & ! output
                       RSTEP)                               ! optional input
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
 !  NETOPO_in: river network topology
 !  RPARAM_in: river reach parameter
 !      RSTEP: Retrospective time step
 !
 !      Inout:
 ! RCHSTA_out: reach wave data structures
 !
 !    Outputs:
 !  Q_JRCH(:): Vector of merged flow particles in reach JRCH
 !  TENTRY(:): Vector of times flow particles entered reach JRCH (exited upstream reaches)
 !  T_EXIT(:): Vector of times flow particles are expected to exit reach JRCH
 !
 ! ----------------------------------------------------------------------------------------
 USE globalData, ONLY: LKTOPO           ! Lake topology
 USE globalData, ONLY: LAKFLX           ! Lake fluxes
 implicit none
 ! Input
 integer(i4b), intent(in)                 :: IENS         ! ensemble member
 integer(i4b), intent(in)                 :: JRCH         ! reach to process
 integer(i4b), intent(in)                 :: LAKEFLAG     ! >0 if processing lakes
 real(dp),     intent(in)                 :: T0,T1        ! start and end of the time step
 integer(i4b), intent(in)                 :: ixDesire     ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable :: NETOPO_in(:) ! River Network topology
 type(RCHPRP), intent(in),    allocatable :: RPARAM_in(:) ! River reach parameter
 type(STRFLX), intent(in),    allocatable :: RCHFLX_in(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b), intent(in),    optional    :: RSTEP        ! retrospective time step offset
 ! inout
 type(STRSTA), intent(inout), allocatable :: RCHSTA_out(:,:) ! reach state data
 ! Output
 real(dp),allocatable, intent(out)        :: Q_JRCH(:)    ! merged (non-routed) flow in JRCH
 real(dp),allocatable, intent(out)        :: TENTRY(:)    ! time flow particles entered JRCH
 real(dp),allocatable, intent(out)        :: T_EXIT(:)    ! time flow is expected to exit JR
 integer(i4b),         intent(out)        :: ierr         ! error code
 character(*),         intent(out)        :: message      ! error message
 ! Local variables to hold the merged inputs to the downstream reach
 integer(i4b)                             :: ROFFSET      ! retrospective offset due to rstep
 real(dp)                                 :: DT           ! model time step
 real(dp), allocatable                    :: QD(:)        ! merged downstream flow
 real(dp), allocatable                    :: TD(:)        ! merged downstream time
 integer(i4b)                             :: iw           ! wave loop index
 integer(i4b)                             :: ND           ! # of waves routed from upstreams
 integer(i4b)                             :: NJ           ! # of waves in the JRCH reach
 integer(i4b)                             :: NK           ! # of waves for routing (NJ+ND)
 integer(i4b)                             :: ILAK         ! lake index
 character(len=strLen)                    :: fmt1         ! format string
 character(len=strLen)                    :: cmessage     ! error message for downwind routine

 ierr=0; message='getusq_rch/'
 ! ----------------------------------------------------------------------------------------
 ! (1) EXTRACT (AND MERGE) FLOW FROM UPSTREAM REACHES OR LAKE
 ! ----------------------------------------------------------------------------------------

 DT = (T1 - T0)
 ! set the retrospective offset
 if (.not.present(RSTEP)) then
   ROFFSET = 0
 else
   ROFFSET = RSTEP
 end if
 if (LAKEFLAG.EQ.1) then  ! lakes are enabled
  ! get lake outflow and only lake outflow if reach is a lake outlet reach, else do as normal
  ILAK = NETOPO_in(JRCH)%LAKE_IX                              ! lake index
  if (ILAK.GT.0) then                                         ! part of reach is in lake
   if (NETOPO_in(JRCH)%REACHIX.eq.LKTOPO(ILAK)%DREACHI) then  ! we are in a lake outlet reach
    ND = 1
    allocate(QD(1),TD(1),STAT=IERR)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for QD and TD'; return; endif
    QD(1) = LAKFLX(IENS,ILAK)%LAKE_Q / RPARAM_in(JRCH)%R_WIDTH  ! lake outflow per unit reach width
    TD(1) = T1 - DT*ROFFSET
   else
    call qexmul_rch(IENS,JRCH,T0,T1,ixDesire,     &   ! input
                   NETOPO_in,RPARAM_in,RCHFLX_in, &   ! input
                   RCHSTA_out,                    &   ! inout
                   ND,QD,TD,ierr,cmessage,        &   ! output
                   RSTEP)                             ! optional input
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   end if
  else
   call qexmul_rch(IENS,JRCH,T0,T1,ixDesire,      &   ! input
                   NETOPO_in,RPARAM_in,RCHFLX_in, &   ! input
                   RCHSTA_out,                    &   ! inout
                   ND,QD,TD,ierr,cmessage,        &   ! output
                   RSTEP)                             ! optional input
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif
 else  ! lakes disabled
  call qexmul_rch(IENS,JRCH,T0,T1,ixDesire,      &   ! input
                  NETOPO_in,RPARAM_in,RCHFLX_in, &   ! input
                  RCHSTA_out,                    &   ! inout
                  ND,QD,TD,ierr,cmessage,        &   ! output
                  RSTEP)                             ! optional input
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  if (JRCH == ixDesire) then
    write(fmt1,'(A,I5,A)') '(A,1X',ND,'(1X,F15.7))'
    write(*,'(a)')      ' * After qexmul_rch: # of routed wave from upstreams (ND) and wave discharge (QD) [m2/s]:'
    write(*,'(A,x,I5)') ' ND=', ND
    write(*,fmt1)       ' QD=', (QD(iw), iw=1,ND)
  end if
 end if

 ! ----------------------------------------------------------------------------------------
 ! (2) EXTRACT NON-ROUTED FLOW FROM THE REACH JRCH & APPEND TO THE FLOW JUST ROUTED D/S
 ! ----------------------------------------------------------------------------------------
 ! check that the routing structure is associated
 if (allocated(RCHSTA_out).eqv..false.) then
  ierr=20; message='routing structure RCHSTA_out is not associated'; return
 endif

 ! check that the wave has been initialized
 if (allocated(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE).eqv..false.) then
  ! if not initialized, then set initial flow to first flow
  ! (this will only occur for a cold start in the case of no streamflow observations)
  allocate(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:0),STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem allocating array for KWAVE'; return; endif
  RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%QF = QD(1)
  RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%TI = T0 - DT - DT*ROFFSET
  RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%TR = T0      - DT*ROFFSET
  RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0)%RF = .true.
 endif

 ! now extract the non-routed flow
 ! NB: routed flows were stripped out in the previous timestep when JRCH was index of u/s reach
 !  {only non-routed flows remain in the routing structure [ + zero element (last routed)]}
 NJ = size(RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE) - 1           ! number of elements not routed (-1 for 0)
 NK = NJ + ND                                     ! pts still in reach + u/s pts just routed

 allocate(Q_JRCH(0:NK),TENTRY(0:NK),T_EXIT(0:NK),STAT=IERR) ! include zero element for INTERP later
 if(ierr/=0)then; message=trim(message)//'problem allocating array for [Q_JRCH, TENTRY, T_EXIT]'; return; endif
 Q_JRCH(0:NJ) = RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NJ)%QF  ! extract the non-routed flow from reach JR
 TENTRY(0:NJ) = RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NJ)%TI  ! extract the non-routed time from reach JR
 T_EXIT(0:NJ) = RCHSTA_out(IENS,JRCH)%LKW_ROUTE%KWAVE(0:NJ)%TR  ! extract the expected exit time
 Q_JRCH(NJ+1:NJ+ND) = QD(1:ND)                                  ! append u/s flow just routed downstream
 TENTRY(NJ+1:NJ+ND) = TD(1:ND)                                  ! append u/s time just routed downstream
 T_EXIT(NJ+1:NJ+ND) = -9999.0D0                                 ! set un-used T_EXIT to missing

 END SUBROUTINE getusq_rch

 ! *********************************************************************
 ! subroutine: extract flow from multiple reaches and merge into
 !                 a single series
 ! *********************************************************************
 SUBROUTINE qexmul_rch(IENS,JRCH,T0,T1,ixDesire,      &   ! input
                       NETOPO_in,RPARAM_in,RCHFLX_in, &   ! input
                       RCHSTA_out,                    &   ! inout
                       ND,QD,TD,ierr,message,         &   ! output
                       RSTEP)                             ! optional input
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
 !  NETOPO_in: river network topology
 !  RPARAM_in: river reach parameter
 !  RCHFLX_in: reach flux
 !      RSTEP: Retrospective time step
 !
 !      Inout:
 ! RCHSTA_out: reach wave data structures
 !
 !   Outputs:
 !      ND   : Number of routed particles
 !      QD(:): Vector of merged flow particles in reach JRCH
 !      TD(:): Vector of times flow particles entered reach JRCH (exited upstream reaches)
 !
 ! ----------------------------------------------------------------------------------------
 implicit none
 ! Input
 integer(i4b), intent(in)                    :: IENS            ! ensemble member
 integer(i4b), intent(in)                    :: JRCH            ! reach to process
 real(dp),     intent(in)                    :: T0,T1           ! start and end of the time step
 integer(i4b), intent(in)                    :: ixDesire        ! index of the reach for verbose output
 type(RCHTOPO),intent(in), allocatable       :: NETOPO_in(:)    ! River Network topology
 type(RCHPRP), intent(in), allocatable       :: RPARAM_in(:)    ! River reach parameter
 type(STRFLX), intent(in), allocatable       :: RCHFLX_in(:,:)  ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b), intent(in), optional          :: RSTEP           ! retrospective time step offset
 ! Inout
 type(STRSTA), intent(inout), allocatable    :: RCHSTA_out(:,:) ! reach state data
 ! Output
 integer(i4b),          intent(out)          :: ND              ! number of routed particles
 real(dp), allocatable, intent(out)          :: QD(:)           ! flow particles just enetered JRCH
 real(dp), allocatable, intent(out)          :: TD(:)           ! time flow particles entered JRCH
 integer(i4b),          intent(out)          :: ierr            ! error code
 character(*),          intent(out)          :: message         ! error message
 ! Local variables to hold flow/time from upstream reaches
 real(dp)                                    :: DT        ! model time step
 integer(i4b)                                :: ROFFSET   ! retrospective offset due to rstep
 integer(i4b)                                :: IUPS      ! loop through u/s reaches
 integer(i4b)                                :: NUPB      ! number of upstream basins
 integer(i4b)                                :: NUPR      ! number of upstream reaches
 integer(i4b)                                :: INDX      ! index of the IUPS u/s reach
 integer(i4b)                                :: MUPR      ! # reaches u/s of IUPS u/s reach
 integer(i4b)                                :: NUPS      ! number of upstream elements
 TYPE(kwtRCH), allocatable                   :: USFLOW(:) ! waves for all upstream segments
 real(dp),     allocatable                   :: UWIDTH(:) ! width of all upstream segments
 integer(i4b)                                :: IMAX      ! max number of upstream particles
 integer(i4b)                                :: IUPR      ! counter for reaches with particles
 integer(i4b)                                :: IR        ! index of the upstream reach
 integer(i4b)                                :: NS        ! size of  the wave
 integer(i4b)                                :: NR        ! # routed particles in u/s reach
 integer(i4b)                                :: NQ        ! NR+1, if non-routed particle exists
 type(FPOINT), allocatable                   :: NEW_WAVE(:)  ! temporary wave
 ! Local variables to merge flow
 logical(lgt), dimension(:), allocatable     :: MFLG      ! T = all particles processed
 integer(i4b), dimension(:), allocatable     :: ITIM      ! processing point for all u/s segments
 real(dp),     dimension(:), allocatable     :: CTIME     ! central time for each u/s segment
 integer(i4b)                                :: JUPS      ! index of reach with the earliest time
 real(dp)                                    :: Q_AGG     ! aggregarted flow at a given time
 integer(i4b)                                :: IWAV      ! index of particle in the IUPS reach
 real(dp)                                    :: SCFAC     ! scale to conform to d/s reach width
 real(dp)                                    :: SFLOW     ! scaled flow at CTIME(JUPS)
 integer(i4b)                                :: IBEG,IEND ! indices for particles that bracket time
 real(dp)                                    :: SLOPE     ! slope for the interpolation
 real(dp)                                    :: PREDV     ! value predicted by the interpolation
 integer(i4b)                                :: IPRT      ! counter for flow particles
 integer(i4b)                                :: JUPS_OLD  ! check that we don't get stuck in do-forever
 integer(i4b)                                :: ITIM_OLD  ! check that we don't get stuck in do-forever
 real(dp)                                    :: TIME_OLD  ! previous time -- used to check for duplicates
 real(dp), allocatable                       :: QD_TEMP(:)! flow particles just enetered JRCH
 real(dp), allocatable                       :: TD_TEMP(:)! time flow particles entered JRCH

 ierr=0; message='qexmul_rch/'

 ! set the retrospective offset
 if (.NOT.present(RSTEP)) then
   ROFFSET = 0
 else
   ROFFSET = RSTEP
 end if
 DT = (T1 - T0)

 ! ----------------------------------------------------------------------------------------
 ! (1) DETERMINE THE NUMBER OF UPSTREAM REACHES
 ! ----------------------------------------------------------------------------------------
 ! ** SPECIAL CASE ** of no basins with any area
 if(count(NETOPO_in(JRCH)%goodBas)==0)return
 ! Need to extract and merge the runoff from all upstream BASINS as well as the streamflow
 !  from all upstream REACHES.  However, streamflow in headwater basins is undefined.  Thus
 !  the number of series merged from upstream reaches is the number of upstream basins +
 !  the number of upstream reaches that are not headwater basins.
 NUPR = 0                               ! number of upstream reaches
 NUPB = size(NETOPO_in(JRCH)%UREACHI)      ! number of upstream basins
 !NUPB = count(NETOPO_in(JRCH)%goodBas)      ! number of upstream basins
 do IUPS=1,NUPB
  INDX = NETOPO_in(JRCH)%UREACHI(IUPS)     ! index of the IUPS upstream reach
  !MUPR = SIZE(NETOPO_in(INDX)%UREACHI)     ! # reaches upstream of the IUPS upstream reach
  MUPR = count(NETOPO_in(INDX)%goodBas)     ! # reaches upstream of the IUPS upstream reach
  if (MUPR.GT.0) NUPR = NUPR + 1        ! reach has streamflow in it, so add that as well
 end do  ! iups
 NUPS = NUPB + NUPR                     ! number of upstream elements (basins + reaches)
 !print*, 'NUPB, NUPR, NUPS', NUPB, NUPR, NUPS
 !print*, 'NETOPO_in(JRCH)%UREACHK = ', NETOPO_in(JRCH)%UREACHK
 !print*, 'NETOPO_in(JRCH)%goodBas = ', NETOPO_in(JRCH)%goodBas

 ! ** SPECIAL CASE ** of just one upstream basin that is a headwater
 if (NUPS.EQ.1) then
  ND = 1
  allocate(QD(1),TD(1),STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem allocating array QD and TD'; return; endif
  ! get reach index
  IR = NETOPO_in(JRCH)%UREACHI(1)
  ! get flow in m2/s (scaled by with of downstream reach)
  QD(1) = RCHFLX_in(IENS,IR)%BASIN_QR(1)/RPARAM_in(JRCH)%R_WIDTH
  TD(1) = T1

  if(JRCH == ixDesire) then
    write(*,'(A,x,I8,x,I8)') ' * Special case - This reach has one headwater upstream: IR, NETOPO_in(IR)%REACHID = ', IR, NETOPO_in(IR)%REACHID
  end if

  return
 end if

 ! allocate space for the upstream flow, time, and flags
 allocate(USFLOW(NUPS),UWIDTH(NUPS),CTIME(NUPS),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [USFLOW, UWIDTH, CTIME]'; return; endif
 ! define the minimum size of the routed data structure (number of flow particles)
 !  (IMAX is increased when looping through the reaches -- section 3 below)
 IMAX = NUPB                            ! flow from basins (one particle / timestep)

 ! ----------------------------------------------------------------------------------------
 ! (2) EXTRACT FLOW FROM UPSTREAM BASINS
 ! ----------------------------------------------------------------------------------------
 do IUPS=1,NUPB
  ! identify the index for the IUPS upstream segment
  IR = NETOPO_in(JRCH)%UREACHI(IUPS)
  ! allocate space for the IUPS stream segment (flow, time, and flags)
  allocate(USFLOW(IUPS)%KWAVE(0:1),STAT=IERR)  ! basin, has flow @start and @end of the time step
  if(ierr>0)then; message=trim(message)//'problem allocating array USFLOW(IUPS)%KWAVE'; return; endif
  ! place flow and time in the KWAVE array (routing done with time-delay histogram in TIMDEL_BAS.F90)
  USFLOW(IUPS)%KWAVE(0:1)%QF = RCHFLX_in(IENS,IR)%BASIN_QR(0:1)      ! flow
  USFLOW(IUPS)%KWAVE(0:1)%TI = (/T0,T1/) - DT*ROFFSET                 ! entry time (not used)
  USFLOW(IUPS)%KWAVE(0:1)%TR = (/T0,T1/) - DT*ROFFSET                 ! exit time
  USFLOW(IUPS)%KWAVE(0:1)%RF = .true.                                 ! routing flag
  !write(*,'(a,i4,1x,2(e20.10,1x))') 'IR, USFLOW(IUPS)%KWAVE(0:1)%QF = ', IR, USFLOW(IUPS)%KWAVE(0:1)%QF
  ! save the upstream width
  UWIDTH(IUPS) = 1.0D0                         ! basin = unit width
  ! save the the time for the first particle in each reach
  CTIME(IUPS) = USFLOW(IUPS)%KWAVE(1)%TR       ! central time
 end do    ! (loop through upstream basins)

 ! ----------------------------------------------------------------------------------------
 ! (3) EXTRACT FLOW FROM UPSTREAM REACHES
 ! ----------------------------------------------------------------------------------------
 IUPR = 0
 do IUPS=1,NUPB
  INDX = NETOPO_in(JRCH)%UREACHI(IUPS)     ! index of the IUPS upstream reach
  !MUPR = SIZE(NETOPO_in(INDX)%UREACHI)     ! # reaches upstream of the IUPS upstream reach
  MUPR = count(NETOPO_in(INDX)%goodBas)     ! # reaches upstream of the IUPS upstream reach
  if (MUPR.GT.0) then                   ! reach has streamflow in it, so add that as well
   IUPR = IUPR + 1
   ! identify the index for the IUPS upstream segment
   IR = NETOPO_in(JRCH)%UREACHI(IUPS)
   ! identify the size of the wave
   NS = size(RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE)
   ! identify number of routed flow elements in the IUPS upstream segment
   NR = count(RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE(:)%RF)
   ! include a non-routed point, if it exists
   NQ = MIN(NR+1,NS)
   ! allocate space for the IUPS stream segment (flow, time, and flags)
   allocate(USFLOW(NUPB+IUPR)%KWAVE(0:NQ-1),STAT=IERR)  ! (zero position = last routed)
   if(ierr/=0)then; message=trim(message)//'problem allocating array USFLOW(NUPB+IUPR)%KWAVE(0:NQ-1)'; return; endif
   ! place data in the new arrays
   USFLOW(NUPB+IUPR)%KWAVE(0:NQ-1) = RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE(0:NQ-1)

   ! here a statement where we check for a modification in the upstream reach;
   ! if flow upstream is modified, then copy RCHSTA_out(:,:)%LKW_ROUTE%KWAVE(:)%QM to USFLOW(..)%KWAVE%QF
   !IF (NUSER.GT.0.AND.SIMDAT%UCFFLAG.GE.1) THEN !if the irrigation module is active and there are users
   !  IF (RCHFLX_out(IENS,IR)%TAKE.GT.0._DP) THEN !if take from upstream reach is greater then zero
   !    ! replace QF with modified flow (as calculated in extract_from_rch)
   !    USFLOW(NUPB+IUPR)%KWAVE(0:NQ-1)%QF = RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE(0:NQ-1)%QM
   !  ENDIF
   !ENDIF

   ! ...and REMOVE the routed particles from the upstream reach
   ! (copy the wave to a temporary wave)
   if (allocated(NEW_WAVE)) then
     deallocate(NEW_WAVE,STAT=IERR)    ! (so we can allocate)
     if(ierr/=0)then; message=trim(message)//'problem deallocating array NEW_WAVE'; return; endif
   end if
   allocate(NEW_WAVE(0:NS-1),STAT=IERR)                 ! get new wave
   if(ierr/=0)then; message=trim(message)//'problem allocating array NEW_WAVE'; return; endif
   NEW_WAVE(0:NS-1) = RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE(0:NS-1)  ! copy

   ! (re-size wave structure)
   if (.not.allocated(RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE))then
     ierr=20; message=trim(message)//'RCHSTA_out%LKW_ROUTE%KWAVE is not associated'; return
   end if
   if (allocated(RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE)) then
     deallocate(RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE,STAT=IERR)
     if(ierr/=0)then; message=trim(message)//'problem deallocating array RCHSTA_out'; return; endif
   end if
   allocate(RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE(0:NS-NR),STAT=IERR)   ! reduced size
   if(ierr/=0)then; message=trim(message)//'problem allocating array RCHSTA_out'; return; endif

   ! (copy "last routed" and "non-routed" elements)
   RCHSTA_out(IENS,IR)%LKW_ROUTE%KWAVE(0:NS-NR) = NEW_WAVE(NR-1:NS-1)

   ! (de-allocate temporary wave)
   deallocate(NEW_WAVE,STAT=IERR)
   if(ierr/=0)then; message=trim(message)//'problem deallocating array NEW_WAVE'; return; endif

   ! save the upstream width
   UWIDTH(NUPB+IUPR) = RPARAM_in(IR)%R_WIDTH            ! reach, width = parameter
   ! save the time for the first particle in each reach
   CTIME(NUPB+IUPR) = USFLOW(NUPB+IUPR)%KWAVE(1)%TR  ! central time

   ! keep track of the total number of points that must be routed downstream
   IMAX = IMAX + (NR-1)     ! exclude zero point for the last routed
  end if ! if reach has particles in it
 end do  ! iups

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
 allocate(QD_TEMP(IMAX),TD_TEMP(IMAX),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [QD_TEMP, TD_TEMP]'; return; endif
 ! allocate positional arrays
 allocate(MFLG(NUPS),ITIM(NUPS),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [MFLG, ITIM]'; return; endif

 ! initalize the flag that defines whether all particles in a given reach are processed
 MFLG(1:NUPS)  = .false.                     ! false until all particles are processed
 ! initialize the search vector
 ITIM(1:NUPS)  = 1                           ! start with the first element of the wave
 ! initialize jups_old and itim_old (used to check we don't get stuck in the do-forever loop)
 JUPS_OLD = HUGE(JUPS_OLD)
 ITIM_OLD = HUGE(ITIM_OLD)
 do      ! loop through all the times in the upstream reaches until no more routed flows
  ! find the reach with the earliest time in all upstream reaches
  !  (NB: the time at the start of the timestep is the earliest possible time and
  !       the time at the end of the timestep is the latest possible time)
  JUPS  = MINLOC(CTIME,DIM=1)     ! JUPS = reach w/ earliest time
  ! check that we're not stuck in a continuous do loop
  if (JUPS.EQ.JUPS_OLD .and. ITIM(JUPS).EQ.ITIM_OLD) then
   ierr=20; message=trim(message)//'stuck in the continuous do-loop'; return
  end if
  ! save jups and itim(jups) to check that we don't get stuck in a continuous do-loop
  JUPS_OLD = JUPS
  ITIM_OLD = ITIM(JUPS)
  ! check that there are still particles in the given reach that require processing
  if (.not.MFLG(JUPS)) then
   ! check that the particle in question is a particle routed (if not, then don't process)
   if (USFLOW(JUPS)%KWAVE(ITIM(JUPS))%RF.EQV..false.) then
    MFLG(JUPS) = .true. ! if routing flag is false, then have already processed all particles
    CTIME(JUPS) = HUGE(SFLOW)  ! largest possible number = ensure reach is not selected again
   ! the particle is in need of processing
   else
    ! define previous time
    if (IPRT.GE.1) then
      TIME_OLD = TD_TEMP(IPRT)
    else ! (if no particles, set to largest possible negative number)
      TIME_OLD = -HUGE(SFLOW)
    end if
    ! check that the particles are being processed in the correct order
    if (CTIME(JUPS).LT.TIME_OLD) then
     ierr=30; message=trim(message)//'expect process in order of time'; return
    end if
    ! don't process if time already exists
    if (CTIME(JUPS).NE.TIME_OLD) then
     ! -------------------------------------------------------------------------------------
     ! compute sum of scaled flow for all reaches
     Q_AGG = 0.0D0
     do IUPS=1,NUPS
      ! identify the element of the wave for the IUPS upstream reach
      IWAV = ITIM(IUPS)
      ! compute scale factor (scale upstream flow by width of downstream reach)
      SCFAC = UWIDTH(IUPS) / RPARAM_in(JRCH)%R_WIDTH
      ! case of the upstream reach with the minimum time (no interpolation required)
      if (IUPS.EQ.JUPS) then
       SFLOW = USFLOW(IUPS)%KWAVE(IWAV)%QF * SCFAC  ! scaled flow
      ! case of all other upstream reaches (*** now, interpolate ***)
      else
       ! identify the elements that bracket the flow particle in the reach JUPS
       ! why .GE.?  Why not .GT.??
       IBEG = IWAV; IF (USFLOW(IUPS)%KWAVE(IBEG)%TR.GE.CTIME(JUPS)) IBEG=IWAV-1
       IEND = IBEG+1  ! *** check the elements are ordered as we think ***
       ! test if we have bracketed properly
       if (USFLOW(IUPS)%KWAVE(IEND)%TR.LT.CTIME(JUPS) .or. &
           USFLOW(IUPS)%KWAVE(IBEG)%TR.GT.CTIME(JUPS)) then
            ierr=40; message=trim(message)//'the times are not ordered as we assume'; return
       end if  ! test for bracketing
       ! estimate flow for the IUPS upstream reach at time CTIME(JUPS)
       SLOPE = (USFLOW(IUPS)%KWAVE(IEND)%QF - USFLOW(IUPS)%KWAVE(IBEG)%QF) / &
               (USFLOW(IUPS)%KWAVE(IEND)%TR - USFLOW(IUPS)%KWAVE(IBEG)%TR)
       PREDV =  USFLOW(IUPS)%KWAVE(IBEG)%QF + SLOPE*(CTIME(JUPS)-USFLOW(IUPS)%KWAVE(IBEG)%TR)
       SFLOW = PREDV * SCFAC  ! scaled flow
      end if  ! (if interpolating)
      ! aggregate flow
      Q_AGG = Q_AGG + SFLOW
     end do  ! looping through upstream elements
     ! -------------------------------------------------------------------------------------
     ! place Q_AGG and CTIME(JUPS) in the output arrays
     IPRT = IPRT + 1
     QD_TEMP(IPRT) = Q_AGG
     TD_TEMP(IPRT) = CTIME(JUPS)
    end if  ! (check that time doesn't already exist)
    ! check if the particle just processed is the last element
    if (ITIM(JUPS).EQ.size(USFLOW(JUPS)%KWAVE)-1) then  ! -1 because of the zero element
     MFLG(JUPS) = .true.            ! have processed all particles in a given u/s reach
     CTIME(JUPS) = huge(SFLOW)      ! largest possible number = ensure reach is not selected again
    else
     ITIM(JUPS) = ITIM(JUPS) + 1                       ! move on to the next flow element
     CTIME(JUPS) = USFLOW(JUPS)%KWAVE(ITIM(JUPS))%TR   ! save the time
    end if  ! (check if particle is the last element)
   end if  ! (check if the particle is a routed element)
  end if  ! (check that there are still particles to process)
  ! if processed all particles in all upstream reaches, then EXIT
  IF (count(MFLG).EQ.NUPS) exit
 end do   ! do-forever

 ! free up memory
 do IUPS=1,NUPS  ! de-allocate each element of USFLOW
  deallocate(USFLOW(IUPS)%KWAVE,STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem deallocating array USFLOW(IUPS)%KWAVE'; return; endif
 end do          ! looping thru elements of USFLOW
 deallocate(USFLOW,UWIDTH,CTIME,ITIM,MFLG,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [USFLOW, UWIDTH, CTIME, ITIM, MFLG]'; return; endif

 ! ...and, save reduced arrays in QD and TD
 ND = IPRT
 allocate(QD(ND),TD(ND),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [QD, TD]'; return; endif
 QD(1:ND) = QD_TEMP(1:ND)
 TD(1:ND) = TD_TEMP(1:ND)

 END SUBROUTINE qexmul_rch

 ! *********************************************************************
 !  subroutine: removes flow particles from the routing structure,
 !                 to reduce memory usage and processing time
 ! *********************************************************************
 SUBROUTINE remove_rch(MAXQPAR,&                           ! input
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
 implicit none
 ! Input
 integer(i4b),          intent(in)           :: MAXQPAR  ! maximum number of flow particles allowed
 ! output
 real(dp), allocatable, intent(inout)        :: Q_JRCH(:)! merged (non-routed) flow in JRCH
 real(dp), allocatable, intent(inout)        :: TENTRY(:)! time flow particles entered JRCH
 real(dp), allocatable, intent(inout)        :: T_EXIT(:)! time flow particles exited JRCH
 integer(i4b),          intent(out)          :: ierr     ! error code
 character(*),          intent(out)          :: message  ! error message
 ! Local variables
 integer(i4b)                                :: NPRT     ! number of flow particles
 integer(i4b)                                :: IPRT     ! loop through flow particles
 real(dp),     dimension(:), allocatable     :: Q,T,Z    ! copies of Q_JRCH and T_JRCH
 logical(lgt), dimension(:), allocatable     :: PARFLG   ! .FALSE. if particle removed
 integer(i4b), dimension(:), allocatable     :: INDEX0   ! indices of original vectors
 real(dp),     dimension(:), allocatable     :: ABSERR   ! absolute error btw interp and orig
 real(dp)                                    :: Q_INTP   ! interpolated particle
 integer(i4b)                                :: MPRT     ! local number of flow particles
 integer(i4b), dimension(:), allocatable     :: INDEX1   ! indices of particles retained
 real(dp),     dimension(:), allocatable     :: E_TEMP   ! temp abs error btw interp and orig
 integer(i4b), dimension(1)                  :: ITMP     ! result of minloc function
 integer(i4b)                                :: ISEL     ! index of local minimum value
 integer(i4b)                                :: INEG     ! lower boundary for interpolation
 integer(i4b)                                :: IMID     ! desired point for interpolation
 integer(i4b)                                :: IPOS     ! upper boundary for interpolation

 ierr=0; message='remove_rch/'

 ! ----------------------------------------------------------------------------------------
 ! (1) INITIALIZATION
 ! ----------------------------------------------------------------------------------------
 ! get the number of particles
 NPRT = size(Q_JRCH)-1                       ! -1 because of zero element
 ! allocate and initialize arrays
 allocate(Q(0:NPRT),T(0:NPRT),Z(0:NPRT),PARFLG(0:NPRT),INDEX0(0:NPRT),ABSERR(0:NPRT),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [Q, T, Z, PARFLG, INDEX0, ABSERR]'; return; endif
 Q = Q_JRCH; T = TENTRY     ! get copies of Q_JRCH and TENTRY
 Z = T_EXIT                 ! (not used in the interp, but include for consistency)
 PARFLG = .true.            ! particle flag = start with all points
 INDEX0 = arth(0,1,NPRT+1)  ! index = (0,1,2,...,NPRT)
 ABSERR = HUGE(Q)           ! largest possible double-precision number
 ! get the absolte difference between actual points and interpolated points
 do IPRT=1,NPRT-1
  ! interpolate at point (iprt)
  Q_INTP = INTERP(T(IPRT),Q(IPRT-1),Q(IPRT+1),T(IPRT-1),T(IPRT+1))
  ! save the absolute difference between the actual value and the interpolated value
  ABSERR(IPRT) = abs(Q_INTP-Q(IPRT))
 end do

 ! ----------------------------------------------------------------------------------------
 ! (2) REMOVAL
 ! ----------------------------------------------------------------------------------------
 do  ! continue looping until the number of particles is below the limit
  ! get the number of particles still in the structure
  MPRT = count(PARFLG)-1       ! -1 because of the zero element
  ! get a copy of (1) indices of selected points, and (2) the interpolation errors
  allocate(INDEX1(0:MPRT),E_TEMP(0:MPRT),STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem allocating arrays [INDEX1, E_TEMP]'; return; endif
  INDEX1 = pack(INDEX0,PARFLG) ! (restrict attention to the elements still present)
  E_TEMP = pack(ABSERR,PARFLG)
  ! check for exit condition (exit after "pack" b.c. indices used to construct final vectors)
  if (MPRT.LT.MAXQPAR) exit
  ! get the index of the minimum value
  ITMP = minloc(E_TEMP)
  ISEL = lbound(E_TEMP,dim=1) + ITMP(1) - 1 ! MINLOC assumes count from 1, here (0,1,2,...NPRT)
  ! re-interpolate the point immediately before the point flagged for removal
  if (INDEX1(ISEL-1).GT.0) then
   INEG=INDEX1(ISEL-2); IMID=INDEX1(ISEL-1); IPOS=INDEX1(ISEL+1)
   Q_INTP = INTERP(T(IMID),Q(INEG),Q(IPOS),T(INEG),T(IPOS))
   ABSERR(IMID) = abs(Q_INTP-Q(IMID))
  end if
  ! re-interpolate the point immediately after the point flagged for removal
  if (INDEX1(ISEL+1).LT.NPRT) then
   INEG=INDEX1(ISEL-1); IMID=INDEX1(ISEL+1); IPOS=INDEX1(ISEL+2)
   Q_INTP = INTERP(T(IMID),Q(INEG),Q(IPOS),T(INEG),T(IPOS))
   ABSERR(IMID) = abs(Q_INTP-Q(IMID))
  end if
  ! flag the point as "removed"
  PARFLG(INDEX1(ISEL)) = .false.
  ! de-allocate arrays
  deallocate(INDEX1,E_TEMP,STAT=IERR)
  if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [INDEX1, E_TEMP]'; return; endif
 end do  ! keep looping until a sufficient number of points are removed

 ! ----------------------------------------------------------------------------------------
 ! (3) RE-SIZE DATA STRUCTURES
 ! ----------------------------------------------------------------------------------------
 deallocate(Q_JRCH,TENTRY,T_EXIT,STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem deallocating arrays [Q_JRCH, TENTRY, T_EXIT]'; return; endif
 allocate(Q_JRCH(0:MPRT),TENTRY(0:MPRT),T_EXIT(0:MPRT),STAT=IERR)
 if(ierr/=0)then; message=trim(message)//'problem allocating arrays [Q_JRCH, TENTRY, T_EXIT]'; return; endif
 Q_JRCH = Q(INDEX1)
 TENTRY = T(INDEX1)
 T_EXIT = Z(INDEX1)

 CONTAINS

  function INTERP(T0,Q1,Q2,T1,T2)
    real(dp),intent(in)                        :: Q1,Q2  ! flow at neighbouring times
    real(dp),intent(in)                        :: T1,T2  ! neighbouring times
    real(dp),intent(in)                        :: T0     ! desired time
    real(dp)                                   :: INTERP ! function name
    INTERP = Q1 + ( (Q2-Q1) / (T2-T1) ) * (T0-T1)
  end function INTERP

 END SUBROUTINE

 ! *********************************************************************
 ! new subroutine: calculate the propagation of kinematic waves in a
 !                 single stream segment, including the formation and
 !                 propagation of a kinematic shock
 ! *********************************************************************
 SUBROUTINE kinwav_rch(JRCH,T_START,T_END,ixDesire,               & ! input: location and time
                       NETOPO_in, RPARAM_in,                      & ! input: river data structure
                       Q_JRCH,TENTRY,T_EXIT,FROUTE,               & ! inout: kwt states
                       NQ2,                                       & ! output:
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
 !                                      if element is routed intent(out)
 !    T_EXIT: array of time elements -- identify the time each element is EXPECTED to exit
 !                                      the stream segment, intent(out).  Used in INTERPTS
 !       NQ2: number of particles    -- <= input becuase multiple particles may merge
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
 !   * All variables are now defined (implicit none) and described (comments)
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
 implicit none
 ! Input
 integer(i4b), intent(in)                    :: JRCH     ! Reach to process
 real(dp),     intent(in)                    :: T_START  ! start of the time step
 real(dp),     intent(in)                    :: T_END    ! end of the time step
 integer(i4b), intent(in)                    :: ixDesire ! index of the reach for verbose output
 type(RCHTOPO),intent(in),    allocatable    :: NETOPO_in(:)    ! River Network topology
 type(RCHPRP), intent(in),    allocatable    :: RPARAM_in(:)    ! River reach parameter
 ! Input/Output
 real(dp),     intent(inout)                 :: Q_JRCH(:)! flow to be routed
 real(dp),     intent(inout)                 :: TENTRY(:)! time to be routed
 real(dp),     intent(inout)                 :: T_EXIT(:)! time pts expected exit segment
 logical(lgt), intent(inout)                 :: FROUTE(:)! routing flag, T=routed
 ! Output
 integer(i4b), intent(out)                   :: NQ2      ! # particles (<= input b/c merge)
 integer(i4b), intent(out)                   :: ierr     ! error code
 character(*), intent(out)                   :: message  ! error message
 ! Internal
 real(dp)                                    :: ALFA     ! constant, 5/3
 real(dp)                                    :: K        ! sqrt(slope)/mannings N
 real(dp)                                    :: XMX      ! length of the stream segment
 integer(i4b)                                :: NN       ! number of input points
 integer(i4b)                                :: NI       ! original size of the input
 integer(i4b)                                :: NM       ! mumber of merged elements
 integer(i4b), dimension(size(Q_JRCH))       :: IX       ! minimum index of each merged element
 integer(i4b), dimension(size(Q_JRCH))       :: MF       ! index for input element merged
 real(dp),     dimension(size(Q_JRCH))       :: T0,T1,T2 ! copy of input time
 real(dp),     dimension(size(Q_JRCH))       :: Q0,Q1,Q2 ! flow series
 real(dp),     dimension(size(Q_JRCH))       :: WC       ! wave celerity
 integer(i4b)                                :: IW,JW    ! looping variables, break check
 real(dp)                                    :: X,XB     ! define smallest, biggest shock
 real(dp)                                    :: WDIFF    ! difference in wave celerity-1
 real(dp)                                    :: XXB      ! wave break
 integer(i4b)                                :: IXB,JXB  ! define position of wave break
 real(dp)                                    :: A1,A2    ! stage - different sides of break
 real(dp)                                    :: CM       ! merged celerity
 real(dp)                                    :: TEXIT    ! expected exit time of "current" particle
 real(dp)                                    :: TNEXT    ! expected exit time of "next" particle
 real(dp)                                    :: TEXIT2   ! exit time of "bottom" of merged element
 integer(i4b)                                :: IROUTE   ! looping variable for routing
 integer(i4b)                                :: JROUTE   ! looping variable for routing
 integer(i4b)                                :: ICOUNT   ! used to account for merged pts
 character(len=strLen)                       :: fmt1     ! format string
 character(len=strLen)                       :: cmessage ! error message of downwind routine
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

 ierr=0; message='kinwav_rch/'

 ! Get the reach parameters
 ALFA = 5._dp/3._dp        ! should this be initialized here or in a parameter file?
 K    = sqrt(RPARAM_in(JRCH)%R_SLOPE)/RPARAM_in(JRCH)%R_MAN_N
 XMX  = RPARAM_in(JRCH)%RLENGTH

 ! Identify the number of points to route
 NN = size(Q1)                                ! modified when elements are merged
 NI = NN                                      ! original size of the input
 if (NN.EQ.0) return                           ! don't do anything if no points in the reach

 ! Initialize the vector that indicates which output element the input elements are merged
 MF = arth(1,1,NI)                            ! Num. Rec. intrinsic: see MODULE nrutil.f90
 ! Initialize the vector that indicates the minumum index of each merged element
 IX = arth(1,1,NI)                            ! Num. Rec. intrinsic: see MODULE nrutil.f90
 ! Get copies of the flow/time particles
 Q0=Q_JRCH; Q1=Q_JRCH; Q2=Q_JRCH
 T0=TENTRY; T1=TENTRY; T2=TENTRY
 ! compute wave celerity for all flow points (array operation)
 WC(1:NN) = ALFA*K**(1./ALFA)*Q1(1:NN)**((ALFA-1.)/ALFA)

 if(jRch==ixDesire) then
   write(fmt1,'(A,I5,A)') '(A,1X',NN,'(1X,F15.7))'
   write(*,'(a)')      ' * Wave discharge (q1) [m2/s] and wave celertiy (wc) [m/s]:'
   write(*,'(a,x,I3)') ' Number of wave =', NN
   write(*,fmt1)       ' q1=', (q1(iw), iw=1,NN)
   write(*,fmt1)       ' wc=', (wc(iw), iw=1,NN)
 end if

 ! handle breaking waves
 GT_ONE: if (NN.GT.1) then                     ! no breaking if just one point
  X = 0.                                      ! altered later to describe "closest" shock
  GOTALL: do                                  ! keep going until all shocks are merged
   XB = XMX                                   ! initialized to length of the stream segment
   ! --------------------------------------------------------------------------------------
   ! check for breaking
   ! --------------------------------------------------------------------------------------
   WCHECK: do IW=2,NN
    JW=IW-1
    if (WC(IW).EQ.0. .or. WC(JW).EQ.0.) cycle  ! waves not moving
    WDIFF = 1./WC(JW) - 1./WC(IW)             ! difference in wave celerity
    if (WDIFF.EQ.0.) cycle                     ! waves moving at the same speed
    if (WC(IW).EQ.WC(JW)) cycle                ! identical statement to the above?
    XXB = (T1(IW)-T1(JW)) / WDIFF             ! XXB is point of breaking in x direction
    if (XXB.LT.X .or. XXB.GT.XB) cycle         ! XB init at LENGTH, so > XB do in next reach
    ! if get to here, the wave is breaking
    XB  = XXB                                 ! identify break "closest to upstream" first
    IXB = IW
   end do WCHECK
   ! --------------------------------------------------------------------------------------
   if (XB.EQ.XMX) exit                         ! got all breaking waves, exit gotall
   ! --------------------------------------------------------------------------------------
   ! combine waves
   ! --------------------------------------------------------------------------------------
   NN  = NN-1
   JXB = IXB-1                                ! indices for the point of breaking
   NM  = NI-NN                                ! number of merged elements
   ! calculate merged shockwave celerity (CM) using finite-difference approximation
   Q2(JXB) =max(Q2(JXB),Q2(IXB))              ! flow of largest merged point
   Q1(JXB) =min(Q1(JXB),Q1(IXB))              ! flow of smallest merged point
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
  end do GOTALL
 end if GT_ONE

 ! check
 if(jRch==ixDesire) then
   write(fmt1,'(A,I5,A)') '(A,1X',NN,'(1X,F15.7))'
   write(*,'(a)')      ' * After wave merge: wave celertiy (wc) [m/s]:'
   write(*,'(a,x,I3)') ' Number of wave =', NN
   write(*,fmt1)       ' wc=', (wc(iw), iw=1,NN)
 end if

 ICOUNT=0
 ! ----------------------------------------------------------------------------------------
 ! perform the routing
 ! ----------------------------------------------------------------------------------------
 do IROUTE = 1,NN    ! loop through the remaining particles (shocks,waves) (NM=NI-NN have been merged)
  ! check that we have non-zero flow
  if(WC(IROUTE) < verySmall)then
   write(message,'(a,i0)') trim(message)//'zero flow for reach id ', NETOPO_in(jRch)%REACHID
   ierr=20; return
  endif
  ! compute the time the shock will exit the reach
  TEXIT = min(XMX/WC(IROUTE) + T1(IROUTE), huge(T1))
  ! compute the time the next shock will exit the reach
  if (IROUTE.LT.NN) TNEXT = min(XMX/WC(IROUTE+1) + T1(IROUTE+1), huge(T1))
  if (IROUTE.EQ.NN) TNEXT = huge(T1)
  ! check if element is merged
  MERGED: if (Q1(IROUTE).NE.Q2(IROUTE)) then
   ! check if merged element has exited
   if (TEXIT.LT.T_END) then
    ! when a merged element exits, save just the top and the bottom of the shock
    ! (identify the exit time for the "slower" particle)
    TEXIT2 = min(TEXIT+1.0D0, TEXIT + 0.5D0*(min(TNEXT,T_END)-TEXIT))
    ! unsure what will happen in the rare case if TEXIT and TEXIT2 are the same
    if (TEXIT2.EQ.TEXIT) then
     ierr=30; message=trim(message)//'TEXIT equals TEXIT2 in kinwav'; return
    end if
    ! fill output arrays
    call rUpdate(Q1(IROUTE),T1(IROUTE),TEXIT,ierr,cmessage)    ! fill arrays w/ Q1, T1, + run checks
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call rUpdate(Q2(IROUTE),T1(IROUTE),TEXIT2,ierr,cmessage)   ! fill arrays w/ Q2, T1, + run checks
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   else                                      ! merged elements have not exited
    ! when a merged element does not exit, need to disaggregate into original particles
    do JROUTE=1,NI                           ! loop thru # original inputs
     if (MF(JROUTE).EQ.IROUTE) &
      call rUpdate(Q0(JROUTE),T0(JROUTE),TEXIT,ierr,cmessage)  ! fill arrays w/ Q0, T0, + run checks
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    end do  ! JROUTE
   end if   ! TEXIT
  ! now process un-merged particles
  else MERGED  ! (i.e., not merged)
   call rUpdate(Q1(IROUTE),T1(IROUTE),TEXIT,ierr,cmessage)     ! fill arrays w/ Q1, T1, + run checks
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if MERGED
 end do
 ! update arrays
 NQ2 = ICOUNT

 CONTAINS

  SUBROUTINE rUpdate(QNEW,TOLD,TNEW,ierr,message)
    real(dp),intent(in)                        :: QNEW      ! Q0,Q1, or Q2
    real(dp),intent(in)                        :: TOLD,TNEW ! entry/exit times
    integer(i4b), intent(out)                  :: ierr      ! error code
    character(*), intent(out)                  :: message   ! error message

    ierr=0; message='RUPDATE/'
    ! ---------------------------------------------------------------------------------------
    ! Used to compute the time each element will exit stream segment & update routing flag
    ! NB: internal subroutine so all data from host is available
    ! ---------------------------------------------------------------------------------------
    ICOUNT=ICOUNT+1
    ! check for array bounds exceeded
    if (ICOUNT.GT.size(Q_JRCH)) then
     ierr=60; message=trim(message)//'array bounds exceeded'; return
    end if
    ! fill output arrays
    Q_JRCH(ICOUNT) = QNEW                         ! flow (Q1 always smaller than Q2)
    TENTRY(ICOUNT) = TOLD                         ! time - note, T1 altered if element merged
    T_EXIT(ICOUNT) = TNEW
    ! time check -- occurs when disaggregating merged elements
    if (ICOUNT.GT.1) then
     if (T_EXIT(ICOUNT).LE.T_EXIT(ICOUNT-1)) T_EXIT(ICOUNT)=T_EXIT(ICOUNT-1)+1.
    end if
    ! another time check -- rare problem when the shock can get the same time as tstart
    if (ICOUNT.EQ.1.and.T_EXIT(ICOUNT).LE.T_START) T_EXIT(ICOUNT)=T_START+1.
    ! update flag for routed elements
    if (T_EXIT(ICOUNT).LT.T_END) FROUTE(ICOUNT) =.true.
  END SUBROUTINE rUpdate

 END SUBROUTINE kinwav_rch

 ! *********************************************************************
 ! new subroutine: calculate time-step averages from irregular values
 ! *********************************************************************
 SUBROUTINE interp_rch(TOLD,QOLD,TNEW,QNEW,IERR,MESSAGE)
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
 !   * All variables are now defined (implicit none) and described (comments)
 !
 !   * Added extra comments
 !
 !   * Replaced GOTO statements with DO loops and IF statements
 !
 ! --------------------------------------------------------------------------------------------
 implicit none
 ! Input
 real(dp), dimension(:), intent(in)          :: TOLD     ! input time array
 real(dp), dimension(:), intent(in)          :: QOLD     ! input flow array
 real(dp), dimension(:), intent(in)          :: TNEW     ! desired output times
 ! Output
 real(dp), dimension(:), intent(out)         :: QNEW     ! flow averaged for desired times
 integer(i4b), intent(out)                   :: IERR     ! error, 1= bad bounds
 character(*), intent(out)                   :: MESSAGE  ! error message
 ! Internal
 integer(i4b)                                :: NOLD     ! number of elements in input array
 integer(i4b)                                :: NNEW     ! number of desired new times
 integer(i4b)                                :: IOLDLOOP ! loop through input times
 integer(i4b)                                :: INEWLOOP ! loop through desired times
 real(dp)                                    :: T0,T1    ! time at start/end of the time step
 integer(i4b)                                :: IBEG     ! identify input times spanning T0
 integer(i4b)                                :: IEND     ! identify input times spanning T1
 integer(i4b)                                :: IMID     ! input times in middle of the curve
 real(dp)                                    :: AREAB    ! area at the start of the time step
 real(dp)                                    :: AREAE    ! area at the end of the time step
 real(dp)                                    :: AREAM    ! area at the middle of the time step
 real(dp)                                    :: AREAS    ! sum of all areas
 real(dp)                                    :: SLOPE    ! slope between two input data values
 real(dp)                                    :: QEST0    ! flow estimate at point T0
 real(dp)                                    :: QEST1    ! flow estimate at point T1

 IERR=0; message='interp_rch/'

 ! get array size
 NOLD = size(TOLD); NNEW = size(TNEW)

 ! check that the input time series starts before the first required output time
 ! and ends after the last required output time
 if( (TOLD(1).GT.TNEW(1)) .or. (TOLD(NOLD).LT.TNEW(NNEW)) ) then
  IERR=1; message=trim(message)//'bad bounds'; return
 end if

 ! loop through the output times
 do INEWLOOP=2,NNEW

  T0 = TNEW(INEWLOOP-1)                      ! start of the time step
  T1 = TNEW(INEWLOOP)                        ! end of the time step

  IBEG=1
  ! identify the index values that span the start of the time step
  BEG_ID: do IOLDLOOP=2,NOLD
   if (T0.LE.TOLD(IOLDLOOP)) then
    IBEG = IOLDLOOP
    exit
   end if
  end do BEG_ID

  IEND=1
  ! identify the index values that span the end of the time step
  END_ID: do IOLDLOOP=1,NOLD
   if (T1.LE.TOLD(IOLDLOOP)) then
    IEND = IOLDLOOP
    exit
   end if
  end do END_ID

  ! initialize the areas
  AREAB=0D0; AREAE=0D0; AREAM=0D0

  ! special case: both TNEW(INEWLOOP-1) and TNEW(INEWLOOP) are within two original values
  ! (implies IBEG=IEND) -- estimate values at both end-points and average
  if (T1.LT.TOLD(IBEG)) then
   SLOPE = (QOLD(IBEG)-QOLD(IBEG-1))/(TOLD(IBEG)-TOLD(IBEG-1))
   QEST0 = SLOPE*(T0-TOLD(IBEG-1)) + QOLD(IBEG-1)
   QEST1 = SLOPE*(T1-TOLD(IBEG-1)) + QOLD(IBEG-1)
   QNEW(INEWLOOP-1) = 0.5*(QEST0 + QEST1)
   cycle ! loop back to the next desired time
  end if

  ! estimate the area under the curve at the start of the time step
  if (T0.LT.TOLD(IBEG)) then  ! if equal process as AREAM
   SLOPE = (QOLD(IBEG)-QOLD(IBEG-1))/(TOLD(IBEG)-TOLD(IBEG-1))
   QEST0 = SLOPE*(T0-TOLD(IBEG-1)) + QOLD(IBEG-1)
   AREAB = (TOLD(IBEG)-T0) * 0.5*(QEST0 + QOLD(IBEG))
  end if

  ! estimate the area under the curve at the end of the time step
  if (T1.LT.TOLD(IEND)) then  ! if equal process as AREAM
   SLOPE = (QOLD(IEND)-QOLD(IEND-1))/(TOLD(IEND)-TOLD(IEND-1))
   QEST1 = SLOPE*(T1-TOLD(IEND-1)) + QOLD(IEND-1)
   AREAE = (T1-TOLD(IEND-1)) * 0.5*(QOLD(IEND-1) + QEST1)
  end if

  ! check if there are extra points to process
  if (IBEG.LT.IEND) then
   ! loop through remaining points
   do IMID=IBEG+1,IEND
    if (IMID.LT.IEND .or. &
      ! process the end slice as AREAM, but only if not already AREAB
      (IMID.EQ.IEND.and.T1.EQ.TOLD(IEND).and.T0.LT.TOLD(IEND-1)) ) then
       ! compute AREAM
       AREAM = AREAM + (TOLD(IMID) - TOLD(IMID-1)) * 0.5*(QOLD(IMID-1) + QOLD(IMID))
    end if   ! if point is valid
   end do  ! IMID
  end if   ! If there is a possibility that middle points even exist

  ! compute time step average
  AREAS = AREAB + AREAE + AREAM            ! sum of all areas
  QNEW(INEWLOOP-1) = AREAS / (T1-T0)       ! T1-T0 is the sum of all time slices

 end do

 END SUBROUTINE interp_rch

END MODULE kwt_route_module
