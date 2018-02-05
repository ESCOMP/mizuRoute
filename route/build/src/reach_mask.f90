MODULE reach_mask_module
 ! *********************************************************************
 ! NOTE
 ! *********************************************************************
 ! This module is copied/edited based on network_route.f90 and vic_route.f90
 ! Removed all subroutines in network_rout.f90 EXCEPT REACH_LIST
 ! Added upstrm_length from vic_route.f90

! data types
USE nrtype
USE nrtype, only : strLen          ! string length
USE nrtype, only : integerMissing  ! missing values for integers
USE dataTypes, only : var_ilength  ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength  ! double precision type: var(:)%dat

! named variables
USE var_lookup,only:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology

implicit none
private
public::reach_mask

contains

 ! *********************************************************************
 ! new subroutine: identify all reaches above the current reach
 ! *********************************************************************
 SUBROUTINE REACH_MASK(desireId,structNTOPO,nRch,isDesired,ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Generates a list of all reaches upstream of a given reach
 !
 ! ----------------------------------------------------------------------------------------
 USE nrtype
 USE nr_utility_module, ONLY : arth                            ! Num. Recipies utilities
 IMPLICIT NONE
 ! input variables
 integer(i4b)      , intent(in)           :: desireId          ! id of the desired reach
 type(var_ilength) , intent(in)           :: structNTOPO(:)    ! network topology structure
 integer(i4b)      , intent(in)           :: nRch              ! number of reaches
 ! output variables
 logical(lgt)      , intent(out)          :: isDesired(:)      ! flags to define that the reach is desired
 integer(i4b)      , intent(out)          :: ierr              ! error code
 character(*)      , intent(out)          :: message           ! error message
 ! ----------------------------------------------------------------------------------------
 ! general local variables
 integer(i4b)                             :: nDesire           ! num desired reaches
 INTEGER(I4B)                             :: IRCH,JRCH,KRCH    ! loop through reaches
 logical(lgt)                             :: isTested(nRch)    ! flags to define that the reach has been tested
 character(len=strLen)                    :: cmessage          ! error message of downwind routine
 ! structure in a linked list
 TYPE NODE_STRUCTURE
  INTEGER(I4B)                            :: IX                ! index of downstream reach
  TYPE(NODE_STRUCTURE), POINTER           :: NEXT              ! next node in linked list
 END TYPE NODE_STRUCTURE
 ! a structure object in a linked list
 TYPE(NODE_STRUCTURE), POINTER            :: NODE              ! node in linked list
 TYPE(NODE_STRUCTURE), POINTER            :: HPOINT            ! Head pointer in linked list
 ! vectors of downstream reaches
 integer(i4b)                             :: nDown             ! number d/s reaches
 integer(i4b)                             :: ixDesire          ! index of desired reach
 integer(i4b),allocatable                 :: ixDownstream(:)   ! indices of downstream reaches
 ! ----------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='REACH_MASK/'

 ! check if we actually want the mask
 if(desireId<0)then
  isDesired(:) = .true.
  return
 endif

 ! initialize logical vectors
 isTested(:)  = .false. ! .true. if we have processed a reach already
 isDesired(:) = .false. ! .true. if reach is desired

 ! loop through all reaches
 DO KRCH=1,nRch

  ! print progress
  !if(mod(kRch,100000)==0) print*, 'Getting reach subset: kRch, nRch = ', kRch, nRch

  ! skip if we have already tested krch
  if(isTested(kRch)) cycle

  ! ---------- get a vector of reaches downstream of a given reach -----------------------------------------------------------------------

  ! initialise the reach array
  nDown = 0        ! initialize the number of upstream reaches
  NULLIFY(HPOINT)  ! set pointer to a linked list to NULL

  ! ensure take streamflow from surrounding basin (a reach is upstream of itself!)
  ! but only if reach is one of the selected
  CALL ADD2LIST(KRCH,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! climb down the network
  IRCH = KRCH
  DO  ! (do-forever)

   ! index of downstream reach
   JRCH = structNTOPO(iRch)%var(ixNTOPO%downSegIndex)%dat(1)

   ! check if complete
   if(jRch<=0) exit         ! negative = missing, which means that the reach is the outlet
   if(isTested(jRch)) exit  ! if we have tested jRch then we have also tested all reaches downstream

   ! add downstream reach to the list
   CALL ADD2LIST(JRCH,ierr,cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! move the index downstream
   IRCH = JRCH

  END DO  ! do forever (go to the outlet)

  ! ---------- extract the vector from the linked list -----------------------------------------------------------------------------------

  ! initialize desired index
  ixDesire = integerMissing

  ! allocate space
  allocate(ixDownstream(nDown),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for downstream reaches'; return; endif

  ! set node to first node in the list of reaches
  node => HPOINT

  jRch = nDown
  ! keep going while node still points to something
  do while (associated(node))

   ! extract the reach index from the linked list
   ixDownstream(jRch) = node%ix

   ! if desired reach is contained in downstream vector, then set the index in downstream vector
   if(structNTOPO(node%ix)%var(ixNTOPO%segId)%dat(1)==desireId) ixDesire=jRch

   ! point to the next node
   node=>node%next

   ! remove the previous node
   deallocate(HPOINT, STAT=IERR)
   if(ierr/=0)then; ierr=20; message=trim(message)//'problem deallocating space for linked list'; return; endif

   ! reset head pointer
   HPOINT=>node
   jRch = jRch-1

  end do  ! while

  ! set the flags to define that these reaches have been tested
  isTested(ixDownstream) = .true.

  ! set the flags to denote the reach is desired
  if(ixDesire/=integerMissing)then
   isDesired(ixDownstream(1:ixDesire)) = .true.
  endif

  ! deallocate space
  deallocate(ixDownstream,stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for downstream reaches'; return; endif

 END DO  ! assess each reach

 ! ---------- check for errors ----------------------------------------------------------------------------------------------------------

 ! check that we processed all reaches
 if(count(isTested)/=nRch)then
  message=trim(message)//'did not process all reaches'
  ierr=20; return
 endif

 ! check that the desired reach exists
 nDesire = count(isDesired)
 if(nDesire==0)then
  message=trim(message)//'desired reach does not exist: NOTE: reach Ids must be >0'
  ierr=20; return
 endif

 ! ----------------------------------------------------------------------------------------
 ! ----------------------------------------------------------------------------------------
 CONTAINS

  ! for a given upstream reach, add a downstream reach index to the list
  subroutine add2list(ixDown,ierr,message)
  integer(i4b) ,intent(in)         :: ixDown   ! downstream reach index
  integer(i4b), intent(out)        :: ierr     ! error code
  character(*), intent(out)        :: message  ! error message
  message='ADD2LIST/'

  ! increment number of downstream reaches
  nDown = nDown+1

  ! allocate node to be added to the list
  allocate(node, stat=ierr)
  if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for linked list'; return; endif

  ! add downstream reach index
  node%ix = ixDown

  ! insert downstream reach at the beginning of existing list
  node%next => HPOINT
  HPOINT    => node

  ! nullify pointer
  NULLIFY(node)

  END SUBROUTINE ADD2LIST

 ! ----------------------------------------------------------------------------------------
 END SUBROUTINE REACH_MASK

end module reach_mask_module
