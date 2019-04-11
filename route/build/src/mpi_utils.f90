Module mpi_mod

  !-------------------------------------------------------------------------------
  ! PURPOSE: general layer on MPI functions
  !-------------------------------------------------------------------------------
  USE mpi

  USE nrtype

  implicit none
  private

  public :: shr_mpi_chkerr
  public :: shr_mpi_commsize
  public :: shr_mpi_commrank
  public :: shr_mpi_initialized
  public :: shr_mpi_abort
  public :: shr_mpi_init
  public :: shr_mpi_finalize

CONTAINS

  SUBROUTINE shr_mpi_chkerr(rcode,string)

    IMPLICIT none

    !----- arguments ---
    integer(i4b),  intent(in)   :: rcode  ! input MPI error code
    character(*),  intent(in)   :: string ! message

    !----- local ---
    character(strLen),parameter :: subName = 'shr_mpi_chkerr/'
    character(strLen)           :: lstring
    integer(i4b)                :: slen
    integer(i4b)                :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: layer on MPI error checking
    !-------------------------------------------------------------------------------

    if (rcode /= MPI_SUCCESS) then
       call MPI_ERROR_STRING(rcode,lstring,slen,ierr)
       write(*,*) trim(subName),lstring(1:slen)
       call shr_mpi_abort(string,rcode)
    endif

  END SUBROUTINE shr_mpi_chkerr

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_commsize(comm,size,string)

    IMPLICIT none

    !----- arguments ---
    integer,intent(in)                 :: comm
    integer,intent(out)                :: size
    character(*),optional,intent(in)   :: string   ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_commsize/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI commsize
    !-------------------------------------------------------------------------------

    call MPI_COMM_SIZE(comm,size,ierr)
    if (present(string)) then
       call shr_mpi_chkerr(ierr,subName//trim(string))
    else
       call shr_mpi_chkerr(ierr,subName)
    endif

  END SUBROUTINE shr_mpi_commsize

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_commrank(comm,rank,string)

    IMPLICIT none

    !----- arguments ---
    integer(i4b),intent(in)            :: comm
    integer(i4b),intent(out)           :: rank
    character(*),optional,intent(in)   :: string   ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_commrank/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI commrank
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(comm,rank,ierr)
    if (present(string)) then
       call shr_mpi_chkerr(ierr,subName//trim(string))
    else
       call shr_mpi_chkerr(ierr,subName)
    endif

  END SUBROUTINE shr_mpi_commrank

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_initialized(flag,string)

    IMPLICIT none

    !----- arguments ---
    logical,intent(out)                :: flag
    character(*),optional,intent(in)   :: string   ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_initialized/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI initialized
    !-------------------------------------------------------------------------------

    call MPI_INITIALIZED(flag,ierr)
    if (present(string)) then
       call shr_mpi_chkerr(ierr,subName//trim(string))
    else
       call shr_mpi_chkerr(ierr,subName)
    endif

  END SUBROUTINE shr_mpi_initialized

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_abort(string,rcode)

    IMPLICIT none

    !----- arguments ---
    character(*),optional,intent(in)   :: string   ! message
    integer(i4b),optional,intent(in)   :: rcode    ! optional code

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_abort/'
    integer(i4b)                       :: ierr
    integer                            :: rc       ! return code

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI abort
    !-------------------------------------------------------------------------------

    if ( present(string) .and. present(rcode) ) then
       write(*,*) trim(subName),trim(string),rcode
    endif
    if ( present(rcode) )then
       rc = rcode
    else
       rc = 1001
    end if
    call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
    stop

  END SUBROUTINE shr_mpi_abort

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_init(string)

    IMPLICIT none

    !----- arguments ---
    character(*),optional,intent(in)   :: string   ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_init/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI init
    !-------------------------------------------------------------------------------

    call MPI_INIT(ierr)
    if (present(string)) then
       call shr_mpi_chkerr(ierr,subName//trim(string))
    else
       call shr_mpi_chkerr(ierr,subName)
    endif

  END SUBROUTINE shr_mpi_init

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_finalize(string)

    IMPLICIT none

    !----- arguments ---
    character(*),optional,intent(in)   :: string   ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_finalize/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI finalize
    !-------------------------------------------------------------------------------

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    if (present(string)) then
       call shr_mpi_chkerr(ierr,subName//trim(string))
    else
       call shr_mpi_chkerr(ierr,subName)
    endif

  END SUBROUTINE shr_mpi_finalize

  !===============================================================================
  !===============================================================================

END MODULE mpi_mod
