MODULE ascii_utils

USE nrtype

implicit none

private
public::file_open
public::split_line
public::get_vlines
public::lower

CONTAINS

 ! **********************************************************************************************
 ! public subroutine: open file
 ! **********************************************************************************************
 SUBROUTINE file_open(infile,unt,err,message)
 implicit none
 ! declare dummy variables
 character(*),intent(in)              :: infile      ! filename
 integer(i4b),intent(out)             :: unt         ! file unit
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 logical(lgt)                         :: xist        ! .TRUE. if the file exists
 logical(lgt)                         :: xopn        ! .TRUE. if the file is already open
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
 ! open file
 open(newunit=unt,file=trim(infile),status="old",action="read",iostat=err)
 if(err/=0)then
   message=trim(message)//"OpenError['"//trim(infile)//"']"
   err=30; return
 endif
 end SUBROUTINE file_open


 ! **********************************************************************************************
 ! public function: split a line of characters into an vector of "words"
 ! **********************************************************************************************
 function split_line(inline, delim) result(words)
 ! do not know how many "words", so use linked lists
 implicit none
 ! declare dummy arguments
 character(*),intent(in)              :: inline   ! line of characters
 character(*),optional,intent(in)     :: delim    ! one character delimiter (optional)
 ! declare local variables
 character(len=256),allocatable :: words(:)             ! returned object: vector of "words"
 character(len=256)      :: temp                  ! temporary line of characters
 character(len=256)      :: delimiter             ! delimiter actually used
 integer(i4b)            :: p1,p2                 !
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
 if (present(delim)) then
   delimiter=delim
 else
   delimiter=' '
 endif

 temp=inline  ! initialize string of characters
 i1=1         ! initialize the index at the start of the first word
 ! ***** loop through the character string
 do iword=1,maxWords
  temp = temp(i1:); if(len_trim(temp)==0) exit
  ! skip leading delimiters
  p1 = verify(temp, trim(delimiter)) ! find first index of character not maching delimiter
  if (p1 == 0) exit ! if all the delimiter character....

  temp = temp(p1:)
  p2 = scan(temp,trim(delimiter))
  if (p2==0) then
    cword=temp
    i1=len_trim(inline)+1
  else
    cword=temp(:p2-1)
    i1=p1+p2
  endif

  ! add the variable to the linked list
  if(iword==1)then
   allocate(list); list=node(cword,iword,null())
   current=>list
  else
   allocate(current%next); current%next=node(cword,iword,null())
   current=>current%next
  endif
  ! check that the line has fewer words than maxWords
  if (iword==maxWords)then; print*, "WARNING: split_line exceeding 100 words"; endif
 end do
 ! ***** allocate space for the list of words
 nWords = current%ix
 allocate(words(nWords))
 ! ***** save the list in a vector, and deallocate space as we go...
 current=>list
 do while(associated(current))
  words(current%ix) = current%chardat
  previous=>current; current=>current%next
  deallocate(previous)
 end do
 end function split_line


 ! **********************************************************************************************
 ! public subroutine: get valid lines of data from file and store as a vector of charater strings
 ! **********************************************************************************************
 SUBROUTINE get_vlines(unt,vlines,err,message)
 ! do not know how many valid lines, so use linked lists
 implicit none
 ! declare dummy arguments
 integer(i4b),intent(in)              :: unt         ! file unit
 character(*),intent(out),allocatable :: vlines(:)   ! vector of character strings
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 integer(i4b)            :: iline                    ! loop through lines in the file
 integer(i4b),parameter  :: maxLines=50000           ! maximum number of valid lines in a file
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
  nullify(previous)
 end do
 nullify(list)
 END SUBROUTINE get_vlines

 ! **********************************************************************************************
 ! public function: convert string to lower case
 ! **********************************************************************************************
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

END MODULE ascii_utils
