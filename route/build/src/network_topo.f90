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

! define the desired reach (set to negative to avoid any printing)
integer(i4b),parameter  :: ixDesire = -9999
contains

 ! *********************************************************************
 ! new subroutine: identify all reaches above the current reach
 ! *********************************************************************
 SUBROUTINE REACH_LIST(ixDesire,NETOPO,ixUpstream,ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Creator(s):
 !   Martyn Clark, 2006
 !   Einar Ã~Vrn Hreinsson, 2009  -- adapt linked lists and selections of reaches
 !   Martyn Clark, 2014 -- modify to be used as a stand-alone module
 !
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Generates a list of all reaches upstream of a given reach
 !
 ! ----------------------------------------------------------------------------------------
 ! I/O:
 !
 !   INPUTS:
 !     ixDesire:   Index of the desired stream segment
 !     NETOPO:     Network topology structure
 !
 !  OUTPUTS:
 !     ixUpstream: Indices of upstream reaches
 !     ierr:       Error code
 !     message:    Error message
 !
 ! ----------------------------------------------------------------------------------------
 ! Structures modified:
 !
 !   None
 !
 ! ----------------------------------------------------------------------------------------
 ! Future revisions:
 !
 !   (none planned)
 !
 ! ----------------------------------------------------------------------------------------
 USE nrtype
 USE reachparam, only : RCHTOPO                              ! Network topology structure
 USE nr_utility_module, ONLY : arth                          ! Num. Recipies utilities
 IMPLICIT NONE
 ! input variables
 integer(i4b), intent(in)                 :: ixDesire        ! index of the desired reach
 type(RCHTOPO), intent(in)                :: NETOPO(:)       ! River Network topology
 ! output variables
 integer(i4b), intent(out)                :: ixUpstream(:)   ! indices of upstream reaches for reach ixDesire
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
 TYPE(UPSTR_RCH),ALLOCATABLE              :: INTLIST         ! list of reaches u/s of each station
 INTEGER(I4B)                             :: NUMUPS          ! number of reaches upstream
 integer(i4b),parameter                   :: strLen=256      ! length of character string
 character(len=strLen)                    :: cmessage        ! error message of downwind routine
 ! ----------------------------------------------------------------------------------------
 message='REACH_LIST/'
 ! initialise the reach array
 INTLIST%N_URCH = 0       ! initialize the number of upstream reaches
 NULLIFY(INTLIST%HPOINT)  ! set pointer to a linked list to NULL

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

  ! allocate space for the upstream reach indices
  NUMUPS = INTLIST%N_URCH ! should be at least 1 (because reach is upstream of itself)
  ALLOCATE(ixUpstream(NUMUPS),STAT=IERR)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for RCHLIST'; return; endif

  ! copy list of upstream reaches to structures (and delete the list)
  CALL MOVE_LIST(JRCH,JRCH,NUMUPS,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  !print*, 'jrch, numups, ixUpstream(:) = ', jrch, numups, ixUpstream(:)

 END DO  ! jrch

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
   ixUpstream(KRCH) = URCH%IUPS
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

end module network_topo
