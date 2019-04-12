Module mpi_mod

  !-------------------------------------------------------------------------------
  ! PURPOSE: general layer on MPI functions
  !-------------------------------------------------------------------------------
  USE mpi

  USE nrtype

  implicit none
  private

  public :: shr_mpi_gatherV
  public :: shr_mpi_scatterV
  public :: shr_mpi_chkerr
  public :: shr_mpi_commsize
  public :: shr_mpi_commrank
  public :: shr_mpi_initialized
  public :: shr_mpi_abort
  public :: shr_mpi_init
  public :: shr_mpi_finalize

  interface shr_mpi_scatterV ; module procedure &
    shr_mpi_scatterIntV,    &
    shr_mpi_scatterRealV,   &
    shr_mpi_scatterLogicalV
  end interface

  interface shr_mpi_gatherV ; module procedure &
    shr_mpi_gatherIntV,    &
    shr_mpi_gatherRealV,   &
    shr_mpi_gatherLogicalV
  end interface

CONTAINS

  ! ----------------------------------
  ! 1D integer array scatter
  ! ----------------------------------
  SUBROUTINE shr_mpi_scatterIntV(globalArray, num_per_proc, & ! input
                                localArray, ierr, message)   ! output
    USE globalData,  ONLY: pid, nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    integer(i4b),              intent(in)  :: globalArray(:)            ! input: global array at root proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! input: number of elements per proc
    integer(i4b), allocatable, intent(out) :: localArray(:)             ! output: scattered array for each proc
    ! Output error handling variables
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                           :: displs(0:nNodes-1)
    integer(i4b)                           :: myid

    ierr=0; message='shr_mpi_scatterIntV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(localArray(num_per_proc(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_SCATTERV(globalArray, num_per_proc(0:nNodes-1), displs, MPI_INT,       & ! flows from proc
                      localArray,  num_per_proc(pid),                MPI_INT, root, & ! scattered flows at root node
                      MPI_COMM_WORLD, ierr)

  END SUBROUTINE shr_mpi_scatterIntV

  ! ----------------------------------
  ! 1D double precision array scatter
  ! ----------------------------------
  SUBROUTINE shr_mpi_scatterRealV(globalArray, num_per_proc, & ! input
                                 localArray, ierr, message)   ! output
    USE globalData,  ONLY: pid, nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    real(dp),              intent(in)  :: globalArray(:)            ! input: global array at root proc
    integer(i4b),          intent(in)  :: num_per_proc(0:nNodes-1)  ! input: number of elements per proc
    real(dp), allocatable, intent(out) :: localArray(:)             ! output: scattered array for each proc
    ! Output error handling variables
    integer(i4b),          intent(out) :: ierr
    character(strLen),     intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                       :: displs(0:nNodes-1)
    integer(i4b)                       :: myid

    ierr=0; message='shr_mpi_scatterRealV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(localArray(num_per_proc(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_SCATTERV(globalArray, num_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION,       & ! flows from proc
                      localArray,  num_per_proc(pid),                MPI_DOUBLE_PRECISION, root, & ! scattered flows at root node
                      MPI_COMM_WORLD, ierr)

  END SUBROUTINE shr_mpi_scatterRealV

  ! ----------------------------------
  ! SCATTER - 1D logical array
  ! ----------------------------------
  SUBROUTINE shr_mpi_scatterLogicalV(globalArray, num_per_proc, & ! input
                                     localArray, ierr, message)   ! output
    USE globalData,  ONLY: pid, nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    logical(lgt),              intent(in)  :: globalArray(:)            ! input: global array at root proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! input: number of elements per proc
    logical(lgt), allocatable, intent(out) :: localArray(:)             ! output: scattered array for each proc
    ! Output error handling variables
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                           :: displs(0:nNodes-1)
    integer(i4b)                           :: myid

    ierr=0; message='shr_mpi_scatterLogicalV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(localArray(num_per_proc(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_SCATTERV(globalArray, num_per_proc(0:nNodes-1), displs, MPI_LOGICAL,       & ! flows from proc
                      localArray,  num_per_proc(pid),                MPI_LOGICAL, root, & ! scattered flows at root node
                      MPI_COMM_WORLD, ierr)

  END SUBROUTINE shr_mpi_scatterLogicalV


  ! ----------------------------------
  ! GATHER - 1D integer array
  ! ----------------------------------
  SUBROUTINE shr_mpi_gatherIntV(localArray, num_per_proc, & ! input
                                globalArray, ierr, message) ! output
    USE globalData,  ONLY: pid, nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    integer(i4b),              intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
    ! Output
    integer(i4b), allocatable, intent(out) :: globalArray(:)            ! gathered  array at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                           :: displs(0:nNodes-1)
    integer(i4b)                           :: myid

    ierr=0; message='shr_mpi_gatherIntV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_GATHERV(localArray,  num_per_proc(pid),                MPI_INT,       & ! local array stuff
                     globalArray, num_per_proc(0:nNodes-1), displs, MPI_INT, root, & ! global array stuff
                     MPI_COMM_WORLD, ierr)

  END SUBROUTINE shr_mpi_gatherIntV

  ! ----------------------------------
  ! GATHER - 1D real8 array
  ! ----------------------------------
  SUBROUTINE shr_mpi_gatherRealV(localArray, num_per_proc,  & ! input
                                 globalArray, ierr, message)  ! output
    USE globalData,  ONLY: pid, nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    real(dp),              intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),          intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
    ! Output
    real(dp), allocatable, intent(out) :: globalArray(:)            ! gathered  array at root proc
    integer(i4b),          intent(out) :: ierr
    character(strLen),     intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                       :: displs(0:nNodes-1)
    integer(i4b)                       :: myid

    ierr=0; message='shr_mpi_gatherRealV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_GATHERV(localArray,  num_per_proc(pid),                MPI_DOUBLE_PRECISION,       & ! local array stuff
                     globalArray, num_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, root, & ! global array stuff
                     MPI_COMM_WORLD, ierr)

  END SUBROUTINE shr_mpi_gatherRealV

  ! ----------------------------------
  ! GATHER - 1D logical array
  ! ----------------------------------
  SUBROUTINE shr_mpi_gatherLogicalV(localArray, num_per_proc,  & ! input
                                    globalArray, ierr, message)  ! output
    USE globalData,  ONLY: pid, nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    logical(lgt),              intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
    ! Output
    logical(lgt), allocatable, intent(out) :: globalArray(:)            ! gathered  array at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                           :: displs(0:nNodes-1)
    integer(i4b)                           :: myid

    ierr=0; message='shr_mpi_gatherLogicalV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_GATHERV(localArray,  num_per_proc(pid),                MPI_LOGICAL,       & ! local array stuff
                     globalArray, num_per_proc(0:nNodes-1), displs, MPI_LOGICAL, root, & ! global array stuff
                     MPI_COMM_WORLD, ierr)

  END SUBROUTINE shr_mpi_gatherLogicalV


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
