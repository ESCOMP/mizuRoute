MODULE network_topo
 ! *********************************************************************
 ! NOTE
 ! *********************************************************************
 ! This module is copied/edited based on network_route.f90 and vic_route.f90
 ! Removed all subroutines in network_rout.f90 EXCEPT REACH_LIST
 ! Added upstrm_length from vic_route.f90

USE nrtype
implicit none
private
public::reach_list
public::upstrm_length

! define the desired reach (set to negative to avoid any printing)
integer(i4b),parameter  :: ixDesire = -9999
contains

 ! *********************************************************************
 ! new subroutine: identify all reaches above the current reach
 ! *********************************************************************
 SUBROUTINE REACH_LIST(NRCH,NTOTAL,ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Martyn Clark, 2006
 !   Einar Ã~Vrn Hreinsson, 2009  -- adapt linked lists and selections of reaches
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
 USE nrtype
 USE reachparam
 USE nrutil, ONLY : arth                                     ! Num. Recipies utilities
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
    !print*, 'irch, jrch, krch = ', irch, jrch, krch
    if (jrch.eq.irch) THEN !  (check that donwstream reach index is the same as current reach index, which means basin w/o reach)
      exit
    endif
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
  !print*, 'jrch, numups, NETOPO(JRCH)%RCHLIST(:) = ', jrch, numups, NETOPO(JRCH)%RCHLIST(:)
 END DO  ! jrch

 ! free up memory
 DEALLOCATE(INTLIST,STAT=IERR)
 if(ierr/=0)then; ierr=20; message=trim(message)//'problem deallocating space for INTLIST'; return; endif
 ! ----------------------------------------------------------------------------------------
 ! ----------------------------------------------------------------------------------------
 CONTAINS

  ! For a down stream reach, add an upstream reach to its list of upstream reaches
  SUBROUTINE ADD2LIST(D_RCH,U_RCH,ierr,message)
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
  END SUBROUTINE ADD2LIST

  ! Copy upstream reaches from linked list to structure and delete the list
  SUBROUTINE MOVE_LIST(IRCH,JRCH,NNODES,ierr,message)
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
  END DO
  END SUBROUTINE MOVE_LIST
 ! ----------------------------------------------------------------------------------------
 END SUBROUTINE REACH_LIST

 ! *********************************************************************
 ! subroutine: Update total stream network length from a each esgment
 ! *********************************************************************
  subroutine upstrm_length(nSeg, &          ! input
                           ierr, message)   ! error control
  ! ----------------------------------------------------------------------------------------
  ! Creator(s):
  !
  !   Naoki Mizukami
  !
  ! ----------------------------------------------------------------------------------------
  ! Purpose:
  !
  !   Calculate total length of upstream reach network from each segment
  !
  ! ----------------------------------------------------------------------------------------
  ! I/O:
  !
  !   INPUTS:
  !    nSeg: Number of stream segments in the river network (reaches)
  !
  ! ----------------------------------------------------------------------------------------
  ! Structures modified:
  !
  !   Updates structure NETOPO%UPSLENG (module reachparam)
  !
  ! ----------------------------------------------------------------------------------------
  USE reachparam

  implicit none
  ! input variables
  integer(I4B), intent(in)               :: nSeg          ! number of stream segments
  ! output variables
  integer(i4b), intent(out)              :: ierr          ! error code
  character(*), intent(out)              :: message       ! error message
  ! local variables
  integer(I4B)                           :: iSeg          ! index for segment loop
  integer(I4B)                           :: iUps          ! index for upstream segment loop
  integer(I4B)                           :: jUps          ! index for upstream segment loop
  INTEGER(I4B)                           :: NASSIGN       ! # reaches currently assigned
  logical(LGT),dimension(:),allocatable  :: RCHFLAG       ! TRUE if reach is processed
  integer(I4B)                           :: nUps          ! number of upstream reaches
  real(DP)                               :: xLocal        ! length of one segement
  real(DP)                               :: xTotal        ! total length of upstream segment

  ! initialize error control
  ierr=0; message='strmlength/'
  ! ----------------------------------------------------------------------------------------
  NASSIGN = 0
  allocate(RCHFLAG(nSeg),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for RCHFLAG'; return; endif
  RCHFLAG(1:nSeg) = .FALSE.
  ! ----------------------------------------------------------------------------------------

  seg_loop: do iSeg=1,nSeg !Loop through each segment

    nUps = size(NETOPO(iSeg)%RCHLIST) ! size of upstream segment
    allocate(NETOPO(iSeg)%UPSLENG(nUps),stat=ierr)
    !print *,'--------------------------------------------'
    !print *,'Seg ID, Num of upstream', iSeg, nUps
    !print *,'list of upstream index',NETOPO(iSeg)%RCHLIST
    upstrms_loop: do iUps=1,nUps !Loop through upstream segments of current segment
      jUps=NETOPO(iSeg)%RCHLIST(iUps) !index of upstreamf segment
      xTotal = 0.0 !Initialize total length of upstream segments
      do
        xLocal=RPARAM(jUps)%RLENGTH ! Get a segment length
        xTotal=xTotal+xLocal
        if (jUps.eq.NETOPO(iSeg)%REACHIX) exit
        jUps = NETOPO(jUps)%DREACHI ! Get index of immediate downstream segment
      enddo
      NETOPO(iSeg)%UPSLENG(iUps)=xTotal
    enddo upstrms_loop
    !print*, 'iSeg,  NETOPO(iSeg)%UPSLENG(:) = ', iSeg, NETOPO(iSeg)%UPSLENG(:)
  enddo seg_loop

  end subroutine upstrm_length

end module network_topo
