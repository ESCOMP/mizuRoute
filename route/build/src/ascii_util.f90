MODULE ascii_util_module

USE nrtype

implicit none

private

public::file_open
public::split_line
public::get_vlines
public::lower
public::upper

CONTAINS

 ! **********************************************************************************************
 ! new subroutine: get unused file unit (modified from DMSL)
 ! **********************************************************************************************
 subroutine getSpareUnit(unt,err,message)
 ! Purpose: returns un-used file unit
 ! Unit will be in the range 7->2.2billion (see comment below).
 ! Comment:
 ! * Can not be pure as it contains the inquire function,tough life...
 ! * In Fortran-95,available units range from 0 to 2,147,483,640.
 !   Preconnected units 5 (keyboard) and 6 (screen) can be re-connected
 !   to a different logical device but to avoid silly errors this is avoided
 !   in this procedure.
 implicit none
 ! dummies
 integer(i4b),intent(out)::unt
 integer(i4b),intent(out)::err
 character(*),intent(out)::message
 ! locals
 integer(i4b)::i
 logical(lgt)::opened,xist
 integer(i4b),parameter::minUnits=7,maxUnits=2147483639
 ! Start procedure here
 do i=minUnits,maxUnits
  inquire(unit=i,opened=opened,exist=xist) ! check unit status
  if(.not.opened.and.xist)then ! un-opened existing unit found
   unt=i; err=0; message="getSpareUnit/ok"
   exit
  endif
  if(i==maxUnits)then  ! all units in use
   unt=-1; err=-10; message="getSpareUnit/allUnitsInUse&"//&
       "&(all 2.2billion-u've goda b jokin')"
  endif
 enddo
 endsubroutine getSpareUnit


 ! **********************************************************************************************
 ! new subroutine: open file
 ! **********************************************************************************************
 subroutine file_open(infile,unt,err,message)
 implicit none
 ! declare dummy variables
 character(*),intent(in)              :: infile      ! filename
 integer(i4b),intent(out)             :: unt         ! file unit
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 logical(lgt)                         :: xist        ! .TRUE. if the file exists
 logical(lgt)                         :: xopn        ! .TRUE. if the file is already open
 character(len=256)                   :: cmessage    ! error message of downwind routine
 ! initialize errors
 err=0; message="f-file_open/"
 ! check if the file exists
 inquire(file=trim(infile),exist=xist)
 if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   err=10; return
 endif
 ! check if the file is already open
 inquire(file=trim(infile),opened=xopn)
 if(xopn)then
   message=trim(message)//"FileIsAlreadyOpen[file='"//trim(infile)//"']"
   err=20; return
 endif
 ! get spare file unit
 call getSpareUnit(unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! open file
 open(unt,file=trim(infile),status="old",action="read",iostat=err)
 if(err/=0)then
   message=trim(message)//"OpenError['"//trim(infile)//"']"
   err=30; return
 endif
 end subroutine file_open


 ! **********************************************************************************************
 ! new subroutine: split a line of characters into an vector of "words"
 ! **********************************************************************************************
 subroutine split_line(inline,words,err,message)
 ! do not know how many "words", so use linked lists
 implicit none
 ! declare dummy arguments
 character(*),intent(in)              :: inline     ! line of characters
 character(*),intent(out),allocatable :: words(:) ! vector of "words"
 integer(i4b),intent(out)             :: err      ! error code
 character(*),intent(out)             :: message  ! error message
 ! declare local variables
 character(len=256)      :: temp                  ! temporary line of characters
 integer(i4b)            :: iword                 ! loop through words
 integer(i4b),parameter  :: maxWords=100          ! maximum number of words in a line
 integer(i4b)            :: i1                    ! index at the start of a given word
 character(len=256)      :: cword                 ! the current word
 integer(i4b)            :: nWords                ! number of words in the character string
 ! define pointers for linked list
 type node
  character(len=256)     :: chardat
  integer(i4b)           :: ix
  type(node),pointer     :: next=>null()
 end type node
 type(node),pointer      :: list=>null()
 type(node),pointer      :: current=>null()
 type(node),pointer      :: previous=>null()
 ! start procedure here
 err=0; message='split_line/'
 temp=inline  ! initialize string of characters
 i1=1         ! initialize the index at the start of the first word
 ! ***** loop through the character string
 do iword=1,maxWords
  ! extract a given "word"
  temp=adjustl(temp(i1:len_trim(temp))); if(len_trim(temp)==0) exit
  read(temp,*) cword
  i1  =len_trim(cword)+1
  ! add the variable to the linked list
  if(iword==1)then
   allocate(list); list=node(cword,iword,null())
   current=>list
  else
   allocate(current%next); current%next=node(cword,iword,null())
   current=>current%next
  endif
  ! check that the line has fewer words than maxWords
  if (iword==maxWords)then; err=20; message=trim(message)//"exceedMaxWords"; return; endif
 end do
 ! ***** allocate space for the list of words
 nWords = current%ix
 allocate(words(nWords),stat=err)
 if(err/=0)then; err=30; message=trim(message)//"problemAllocateWords"; return; endif
 ! ***** save the list in a vector, and deallocate space as we go...
 current=>list
 do while(associated(current))
  words(current%ix) = current%chardat
  previous=>current; current=>current%next
  deallocate(previous)
 end do
 end subroutine split_line


 ! **********************************************************************************************
 ! new subroutine: get valid lines of data from file and store as a vector of charater strings
 ! **********************************************************************************************
 subroutine get_vlines(unt,vlines,err,message)
 ! do not know how many valid lines, so use linked lists
 implicit none
 ! declare dummy arguments
 integer(i4b),intent(in)              :: unt         ! file unit
 character(*),intent(out),allocatable :: vlines(:)   ! vector of character strings
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 integer(i4b)            :: iline                    ! loop through lines in the file
 integer(i4b),parameter  :: maxLines=1000            ! maximum number of valid lines in a file
 character(len=256)      :: temp                     ! character data or a given line
 integer(i4b)            :: icount                   ! counter for the valid lines
 integer(i4b)            :: iend                     ! index to indicate end of the file
 ! define pointers for linked list
 type node
  character(len=256)     :: chardat
  integer(i4b)           :: ix
  type(node),pointer     :: next=>null()
 end type node
 type(node),pointer      :: list=>null()
 type(node),pointer      :: current=>null()
 type(node),pointer      :: previous=>null()
 ! start procedure here
 err=0; message='get_vlines/'
 ! ***** get the valid lines of data from the file and store in linked lists *****
 icount=0  ! initialize the counter for the valid lines
 do iline=1,maxLines
  read(unt,'(a)',iostat=iend)temp; if(iend/=0)exit    ! read line of data
  if (temp(1:1)=='!')cycle
  icount = icount+1
  ! add the variable to the linked list
  if(.not.associated(list))then
   allocate(list,previous,current); list=node(temp,icount,null())
   current=>list
  else
   allocate(current%next)
   current%next=node(temp,icount,null())
   current=>current%next
  endif
  if (iline==maxLines)then; err=20; message=trim(message)//"exceedMaxLines"; return; endif
 end do  ! looping through the lines in the file (exit clause above will kick in)
 ! ***** allocate space for the valid lines *****
 allocate(vlines(icount),stat=err)
 if(err/=0)then; err=30; message=trim(message)//"problemAllocateVlines"; return; endif
 ! ***** save the list in a vector, and deallocate space as we go... *****
 current=>list
 do while(associated(current))
  vlines(current%ix) = current%chardat
  previous=>current; current=>current%next
  deallocate(previous)
 end do
 if(associated(list)) nullify(list)
 end subroutine get_vlines


  FUNCTION upper(strIn) RESULT(strOut)
    ! convert string to upper-case
    ! only ASCII character code works
    implicit none

    character(*), intent(in) :: strIn
    character(len(strIn))    :: strOut
    integer, parameter       :: DUC = ichar('A') - ichar('a')
    character                :: ch
    integer                  :: i

    do i = 1, len(strIn)
      ch = strIn(i:i)
      if (ch>='a' .and. ch<='z') ch = char(ichar(ch)+DUC)
      strOut(i:i) = ch
    end do

  END FUNCTION upper


  pure FUNCTION lower(strIn) RESULT(strOut)
    ! convert string to lower-case
    ! only ASCII character code works
     implicit none

    character(*), intent(in)  :: strIn
    character(len(strIn))     :: strOut
    integer, parameter        :: DUC = ichar('A') - ichar('a')
    character                 :: ch
    integer                   :: i

    do i = 1,len(strIn)
      ch = strIn(i:i)
      if (ch>='A' .and. ch<='Z') ch = char(ichar(ch)-DUC)
      strOut(i:i) = ch
    end do

  END FUNCTION lower


END MODULE ascii_util_module
