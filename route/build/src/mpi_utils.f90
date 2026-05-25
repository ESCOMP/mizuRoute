MODULE mpi_utils
  !-------------------------------------------------------------------------------
  ! PURPOSE: general layer on MPI functions
  !-------------------------------------------------------------------------------
  USE mpi
  USE nrtype
  USE globalData, ONLY: pid, nNodes
  USE globalData, ONLY: mpicom_route
  USE globalData, ONLY: masterproc
  USE public_var, ONLY: root
  USE public_var, ONLY: iulog

  implicit none

  private
  public :: shr_mpi_send
  public :: shr_mpi_bcast
  public :: shr_mpi_gatherV
  public :: shr_mpi_scatterV
  public :: shr_mpi_allgather
  public :: shr_mpi_reduce
  public :: shr_mpi_sparse_distribute
  public :: shr_mpi_chkerr
  public :: shr_mpi_commsize
  public :: shr_mpi_commrank
  public :: shr_mpi_initialized
  public :: shr_mpi_barrier
  public :: shr_mpi_abort
  public :: shr_mpi_init
  public :: shr_mpi_finalize

  INTERFACE shr_mpi_send; module procedure &
    shr_mpi_sendReal
  END INTERFACE

  INTERFACE shr_mpi_bcast; module procedure &
    shr_mpi_bcastInt_scalar,  &
    shr_mpi_bcastChar_scalar, &
    shr_mpi_bcastInt,         &
    shr_mpi_bcastLong,        &
    shr_mpi_bcastReal,        &
    shr_mpi_bcastChar,        &
    shr_mpi_bcastLogical
  END INTERFACE

  INTERFACE shr_mpi_scatterV ; module procedure &
    shr_mpi_scatterIntV,    &
    shr_mpi_scatterRealV,   &
    shr_mpi_scatterLogicalV
  END INTERFACE

  INTERFACE shr_mpi_gatherV ; module procedure &
    shr_mpi_gatherIntV,    &
    shr_mpi_gatherRealV,   &
    shr_mpi_gatherLogicalV
  END INTERFACE

  INTERFACE shr_mpi_allgather ; module procedure &
    shr_mpi_allgatherIntV,    &
    shr_mpi_allgatherRealV,   &
    shr_mpi_allgatherLogicalV,&
    shr_mpi_allgatherInt,     &
    shr_mpi_allgatherReal,    &
    shr_mpi_allgatherLogical
  END INTERFACE

  INTERFACE shr_mpi_reduce ; module procedure &
    shr_mpi_reduceReal,   &
    shr_mpi_reduceInt,    &
    shr_mpi_reduceRealV
  END INTERFACE

  integer(i4b), parameter :: send_data_tag=2001
  integer(i4b), parameter :: return_data_tag=2002
  integer(i4b)            :: status(MPI_STATUS_SIZE)

CONTAINS

  ! ----------------------------------
  ! SEND-RECV - real array element
  ! ----------------------------------
  SUBROUTINE shr_mpi_sendReal(sendArray,        & ! in:    an array to be sent
                              pidSend, ixSend,  & ! in:    source task and index of the array to be sent
                              recArray,         & ! inout: array to be modified after communication
                              pidRec, ixRec,    & ! in:    destination task and index of the array to be changed
                              ierr, message)

    ! Descriptions:
    !  sending one real element of array in a particular task to the other task
    !  and add the element value to the values in destination array

    implicit none
    ! Argument variables:
    real(dp),         intent(inout) :: sendArray(:)     ! in:    array to be sent
    integer(i4b),     intent(in)    :: pidSend, ixSend  ! in:    task holding an array to be sent and index of the array
    real(dp),         intent(inout) :: recArray(:)      ! inout: array to be modified after communication
    integer(i4b),     intent(in)    :: pidRec, ixRec    ! in:    task holding modifying array and index of the array
    integer(i4b),     intent(out)   :: ierr             ! out:   error code
    character(strLen),intent(out)   :: message          ! out:   error message
    ! local variable
    real(dp)                        :: dummyArray(1)! temporary reciving element

    ierr=0; message='shr_mpi_sendReal/'

    if (pid==pidSend) then
      call MPI_SEND(sendArray(ixSend), 1, MPI_DOUBLE_PRECISION, pidRec, send_data_tag, mpicom_route, ierr)
    else if (pid==pidRec) then
      call MPI_RECV(dummyArray,  1, MPI_DOUBLE_PRECISION, pidSend, send_data_tag, mpicom_route, status, ierr)
      recArray(ixRec) = recArray(ixRec) + dummyArray(1)
    endif

  END SUBROUTINE shr_mpi_sendReal

  ! ----------------------------------
  ! BROADCAST - integer scalar
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastInt_scalar(scalar,       & ! inout:  array to be broadcasted to each proc
                                     ierr, message)  ! output: error handling

    implicit none
    ! Argument variables:
    integer(i4b),              intent(inout) :: scalar        ! inout:  array to be sent to proc
    integer(i4b),              intent(out)   :: ierr
    character(strLen),         intent(out)   :: message       ! error message

    ierr=0; message='shr_mpi_bcastInt/'

    call MPI_BCAST(scalar, 1, MPI_INTEGER, root, mpicom_route, ierr)

  END SUBROUTINE shr_mpi_bcastInt_scalar

  ! ----------------------------------
  ! BROADCAST - character scalar
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastChar_scalar(scalar,       & ! inout:  array to be broadcasted to each proc
                                      ierr, message)  ! output: error handling

    implicit none
    ! Argument variables:
    character(strLen),         intent(inout) :: scalar        ! inout:  array to be sent to proc
    integer(i4b),              intent(out)   :: ierr
    character(strLen),         intent(out)   :: message       ! error message

    ierr=0; message='shr_mpi_bcastInt/'

    call MPI_BCAST(scalar, strLen, MPI_CHARACTER, root, mpicom_route, ierr)

  END SUBROUTINE shr_mpi_bcastChar_scalar

  ! ----------------------------------
  ! BROADCAST - integer allocatable array
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastInt(allocArray,   & ! inout:  array to be broadcasted to each proc
                              ierr, message)  ! output: error handling
    implicit none
    ! Argument variables:
    integer(i4b), allocatable, intent(inout) :: allocArray(:) ! inout:  array to be sent to proc
    integer(i4b),              intent(out)   :: ierr
    character(strLen),         intent(out)   :: message       ! error message
    ! local variable
    integer(i4b)                             :: myid
    integer(i4b)                             :: array_size
    integer(i4b)                             :: bound(2)

    ierr=0; message='shr_mpi_bcastInt/'

    if (masterproc) then

      bound = (/lbound(allocArray), ubound(allocArray)/)

      array_size = bound(2)-bound(1)+1

      do myid = 1, nNodes-1
       call MPI_SEND(bound(1),             2,          MPI_INTEGER, myid, send_data_tag, mpicom_route, ierr)
       call MPI_SEND(allocArray(bound(1)), array_size, MPI_INTEGER, myid, send_data_tag, mpicom_route, ierr)
      end do

    else

       call MPI_RECV(bound, 2, MPI_INTEGER, root, send_data_tag, mpicom_route, status, ierr)

       array_size = bound(2)-bound(1)+1

       if (allocated(allocArray)) then
         deallocate(allocArray, stat=ierr)
         if(ierr/=0)then; message=trim(message)//'probleml de-allocating array for [allocArray]'; return; endif
       endif
       allocate(allocArray(bound(1):bound(2)), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating array for [allocArray]'; return; endif

       call MPI_RECV(allocArray, array_size, MPI_INTEGER, root, send_data_tag, mpicom_route, status, ierr)

    endif

  END SUBROUTINE shr_mpi_bcastInt

  ! ----------------------------------
  ! BROADCAST - integer allocatable array
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastLong(allocArray,   & ! inout:  array to be broadcasted to each proc
                              ierr, message)  ! output: error handling
    implicit none
    ! Argument variables:
    integer(i8b), allocatable, intent(inout) :: allocArray(:) ! inout:  array to be sent to proc
    integer(i4b),              intent(out)   :: ierr
    character(strLen),         intent(out)   :: message       ! error message
    ! local variable
    integer(i4b)                             :: myid
    integer(i4b)                             :: array_size
    integer(i4b)                             :: bound(2)

    ierr=0; message='shr_mpi_bcastLong/'

    if (masterproc) then

      bound = (/lbound(allocArray), ubound(allocArray)/)

      array_size = bound(2)-bound(1)+1

      do myid = 1, nNodes-1
       call MPI_SEND(bound(1),             2,          MPI_INTEGER, myid, send_data_tag, mpicom_route, ierr)
       call MPI_SEND(allocArray(bound(1)), array_size, MPI_INTEGER8, myid, send_data_tag, mpicom_route, ierr)
      end do

    else

       call MPI_RECV(bound, 2, MPI_INTEGER, root, send_data_tag, mpicom_route, status, ierr)

       array_size = bound(2)-bound(1)+1

       if (allocated(allocArray)) then
         deallocate(allocArray, stat=ierr)
         if(ierr/=0)then; message=trim(message)//'probleml de-allocating array for [allocArray]'; return; endif
       endif
       allocate(allocArray(bound(1):bound(2)), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating array for [allocArray]'; return; endif

       call MPI_RECV(allocArray, array_size, MPI_INTEGER8, root, send_data_tag, mpicom_route, status, ierr)

    endif

  END SUBROUTINE shr_mpi_bcastLong

  ! ----------------------------------
  ! BROADCAST - real allocatable array
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastReal(allocArray,   & ! input:  array to be broadcasted to each proc
                               ierr, message)  ! output: error handling
    implicit none
    ! Argument variables:
    real(dp), allocatable, intent(inout) :: allocArray(:) ! inout:  array to be sent to proc
    integer(i4b),          intent(out)   :: ierr
    character(strLen),     intent(out)   :: message       ! error message
    ! local variable
    integer(i4b)                         :: myid
    integer(i4b)                         :: array_size
    integer(i4b)                         :: bound(2)

    ierr=0; message='shr_mpi_bcastReal/'

    if (masterproc) then

      bound = (/lbound(allocArray), ubound(allocArray)/)

      array_size = bound(2)-bound(1)+1

      do myid = 1, nNodes-1
       call MPI_SEND(bound(1),             2,          MPI_INTEGER,              myid, send_data_tag, mpicom_route, ierr)
       call MPI_SEND(allocArray(bound(1)), array_size, MPI_DOUBLE_PRECISION, myid, send_data_tag, mpicom_route, ierr)
      end do

    else

       call MPI_RECV(bound, 2, MPI_INTEGER, root, send_data_tag, mpicom_route, status, ierr)

       array_size = bound(2)-bound(1)+1

       if (allocated(allocArray)) then
         deallocate(allocArray, stat=ierr)
         if(ierr/=0)then; message=trim(message)//'probleml de-allocating array for [allocArray]'; return; endif
       endif
       allocate(allocArray(bound(1):bound(2)), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating array for [allocArray]'; return; endif

       call MPI_RECV(allocArray, array_size, MPI_DOUBLE_PRECISION, root, send_data_tag, mpicom_route, status, ierr)

    endif

  END SUBROUTINE shr_mpi_bcastReal

  ! ----------------------------------
  ! BROADCAST - logical allocatable array
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastLogical(allocArray,    & ! inout: array to be broadcasted to each proc
                                  ierr, message)   ! output: error handling
    implicit none
    ! Argument variables:
    logical(lgt), allocatable, intent(inout) :: allocArray(:)   ! input:  array to be sent to proc
    integer(i4b),              intent(out)   :: ierr
    character(strLen),         intent(out)   :: message         ! error message
    ! local variable
    integer(i4b)                             :: myid
    integer(i4b)                             :: array_size
    integer(i4b)                             :: bound(2)

    ierr=0; message='shr_mpi_bcastLogical/'

    ! send allocatable arrays
    if (masterproc) then

      bound = (/lbound(allocArray), ubound(allocArray)/)

      array_size = bound(2)-bound(1)+1

      do myid = 1, nNodes-1
       call MPI_SEND(bound(1),             2,          MPI_INTEGER,     myid, send_data_tag, mpicom_route, ierr)
       call MPI_SEND(allocArray(bound(1)), array_size, MPI_LOGICAL, myid, send_data_tag, mpicom_route, ierr)
      end do

    else

       call MPI_RECV(bound, 2, MPI_INTEGER, root, send_data_tag, mpicom_route, status, ierr)

       array_size = bound(2)-bound(1)+1

       if (allocated(allocArray)) then
         deallocate(allocArray, stat=ierr)
         if(ierr/=0)then; message=trim(message)//'probleml de-allocating array for [allocArray]'; return; endif
       endif
       allocate(allocArray(bound(1):bound(2)), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating array for [allocArray]'; return; endif

       call MPI_RECV(allocArray, array_size, MPI_LOGICAL, root, send_data_tag, mpicom_route, status, ierr)

    endif

  END SUBROUTINE shr_mpi_bcastLogical

  ! ----------------------------------
  ! BROADCAST - character allocatable array (fixed character length)
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastChar(allocArray,    & ! inout: array to be broadcasted to each proc
                               ierr, message)   ! output: error handling
    implicit none
    ! Argument variables:
    character(strLen), allocatable, intent(inout) :: allocArray(:) ! inout: array to be sent to proc
    integer(i4b),                   intent(out)   :: ierr
    character(strLen),              intent(out)   :: message       ! error message
    ! local variable
    integer(i4b)                                  :: flat_size
    integer(i4b)                                  :: num_elements

    ierr=0; message='shr_mpi_bcastChar/'

    if (masterproc) then
      num_elements = size(allocArray)
    end if
    call MPI_BCAST(num_elements, 1, MPI_INTEGER, root, mpicom_route, ierr)
    flat_size = num_elements*strLen

    ! Allocate allocArray in other ranks after receiving num_elements
    if (.not. masterproc) then
      if (allocated(allocArray)) then
        deallocate(allocArray, stat=ierr)
        if(ierr/=0)then; message=trim(message)//'probleml de-allocating array for [allocArray]'; return; endif
      endif
      allocate(allocArray(num_elements), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating array for [allocArray]'; return; endif
    end if

    ! Broadcast the flattened character array
    call MPI_BCAST(allocArray, flat_size, MPI_CHARACTER, root, mpicom_route, ierr)

  END SUBROUTINE shr_mpi_bcastChar

  ! ----------------------------------
  ! SCATTERV - 1D integer array
  ! ----------------------------------
  SUBROUTINE shr_mpi_scatterIntV(globalArray, num_per_proc, & ! input
                                localArray, ierr, message)   ! output
    implicit none
    ! Argument variables:
    integer(i4b),              intent(in)  :: globalArray(:)            ! input: global array at root proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! input: number of elements per proc
    integer(i4b), allocatable, intent(out) :: localArray(:)             ! output: scattered array for each proc
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

    if (allocated(localArray)) then
     deallocate(localArray, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [localArray]'; return; endif
    endif
    allocate(localArray(num_per_proc(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_SCATTERV(globalArray, num_per_proc(0:nNodes-1), displs, MPI_INTEGER,       & ! flows from proc
                      localArray,  num_per_proc(pid),                MPI_INTEGER, root, & ! scattered flows at root node
                      mpicom_route, ierr)

  END SUBROUTINE shr_mpi_scatterIntV

  ! ----------------------------------
  ! SCATTERV - 1D double precision array
  ! ----------------------------------
  SUBROUTINE shr_mpi_scatterRealV(globalArray, num_per_proc, & ! input
                                 localArray, ierr, message)   ! output
    implicit none
    ! Argument variables:
    real(dp),              intent(in)  :: globalArray(:)            ! input: global array at root proc
    integer(i4b),          intent(in)  :: num_per_proc(0:nNodes-1)  ! input: number of elements per proc
    real(dp), allocatable, intent(out) :: localArray(:)             ! output: scattered array for each proc
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

    if (allocated(localArray)) then
     deallocate(localArray, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [localArray]'; return; endif
    endif
    allocate(localArray(num_per_proc(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_SCATTERV(globalArray, num_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION,       & ! flows from proc
                      localArray,  num_per_proc(pid),                MPI_DOUBLE_PRECISION, root, & ! scattered flows at root node
                      mpicom_route, ierr)

  END SUBROUTINE shr_mpi_scatterRealV

  ! ----------------------------------
  ! SCATTERV - 1D logical array
  ! ----------------------------------
  SUBROUTINE shr_mpi_scatterLogicalV(globalArray, num_per_proc, & ! input
                                     localArray, ierr, message)   ! output
    implicit none
    ! Argument variables:
    logical(lgt),              intent(in)  :: globalArray(:)            ! input: global array at root proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! input: number of elements per proc
    logical(lgt), allocatable, intent(out) :: localArray(:)             ! output: scattered array for each proc
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

    if (allocated(localArray)) then
     deallocate(localArray, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [localArray]'; return; endif
    endif
    allocate(localArray(num_per_proc(pid)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_SCATTERV(globalArray, num_per_proc(0:nNodes-1), displs, MPI_LOGICAL,       & ! flows from proc
                      localArray,  num_per_proc(pid),                MPI_LOGICAL, root, & ! scattered flows at root node
                      mpicom_route, ierr)

  END SUBROUTINE shr_mpi_scatterLogicalV

  ! ----------------------------------
  ! GATHERV - 1D integer array
  ! ----------------------------------
  SUBROUTINE shr_mpi_gatherIntV(localArray, num_per_proc, & ! input
                                globalArray, ierr, message) ! output
    implicit none
    ! Argument variables:
    integer(i4b),              intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
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

    if (allocated(globalArray)) then
     deallocate(globalArray, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [globalArray]'; return; endif
    endif
    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [globalArray]'; return; endif

    call MPI_GATHERV(localArray,  num_per_proc(pid),                MPI_INTEGER,       & ! local array stuff
                     globalArray, num_per_proc(0:nNodes-1), displs, MPI_INTEGER, root, & ! global array stuff
                     mpicom_route, ierr)

  END SUBROUTINE shr_mpi_gatherIntV

  ! ----------------------------------
  ! GATHERV - 1D real8 array
  ! ----------------------------------
  SUBROUTINE shr_mpi_gatherRealV(localArray, num_per_proc,  & ! input
                                 globalArray, ierr, message)  ! output
    implicit none
    ! Argument variables:
    real(dp),              intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),          intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
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

    if (allocated(globalArray)) then
     deallocate(globalArray, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [globalArray]'; return; endif
    endif
    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_GATHERV(localArray,  num_per_proc(pid),                MPI_DOUBLE_PRECISION,       & ! local array stuff
                     globalArray, num_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, root, & ! global array stuff
                     mpicom_route, ierr)

  END SUBROUTINE shr_mpi_gatherRealV

  ! ----------------------------------
  ! GATHERV - 1D logical array
  ! ----------------------------------
  SUBROUTINE shr_mpi_gatherLogicalV(localArray, num_per_proc,  & ! input
                                    globalArray, ierr, message)  ! output
    implicit none
    ! Argument variables:
    logical(lgt),              intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
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

    if (allocated(globalArray)) then
     deallocate(globalArray, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem de-allocating array for [globalArray]'; return; endif
    endif
    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_GATHERV(localArray,  num_per_proc(pid),                MPI_LOGICAL,       & ! local array stuff
                     globalArray, num_per_proc(0:nNodes-1), displs, MPI_LOGICAL, root, & ! global array stuff
                     mpicom_route, ierr)

  END SUBROUTINE shr_mpi_gatherLogicalV

  ! ALLGATHER
  ! ----------------------------------
  ! ALLGATHER - 1D integer array
  ! ----------------------------------
  SUBROUTINE shr_mpi_allgatherIntV(localArray,  num_per_proc, & ! input
                                   globalArray, ierr, message) ! output
    implicit none
    ! Argument variables:
    integer(i4b),              intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
    integer(i4b), allocatable, intent(out) :: globalArray(:)            ! gathered  array at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                           :: displs(0:nNodes-1)
    integer(i4b)                           :: myid

    ierr=0; message='shr_mpi_allgatherIntV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_ALLGATHERV(localArray,  num_per_proc(pid),                MPI_INTEGER, & ! local array stuff
                        globalArray, num_per_proc(0:nNodes-1), displs, MPI_INTEGER, & ! global array stuff
                        mpicom_route, ierr)

  END SUBROUTINE shr_mpi_allgatherIntV

  ! ----------------------------------
  ! ALLGATHER - real8 array
  ! ----------------------------------
  SUBROUTINE shr_mpi_allgatherRealV(localArray,  num_per_proc, & ! input
                                    globalArray, ierr, message) ! output
    implicit none
    ! Argument variables:
    real(dp),                  intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
    real(dp),     allocatable, intent(out) :: globalArray(:)            ! gathered  array at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                           :: displs(0:nNodes-1)
    integer(i4b)                           :: myid

    ierr=0; message='shr_mpi_allgatherRealV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_ALLGATHERV(localArray,  num_per_proc(pid),                MPI_DOUBLE_PRECISION, & ! local array stuff
                        globalArray, num_per_proc(0:nNodes-1), displs, MPI_DOUBLE_PRECISION, & ! global array stuff
                        mpicom_route, ierr)

  END SUBROUTINE shr_mpi_allgatherRealV

  ! ----------------------------------
  ! ALLGATHER - logical array
  ! ----------------------------------
  SUBROUTINE shr_mpi_allgatherLogicalV(localArray,  num_per_proc, & ! input
                                       globalArray, ierr, message) ! output
    implicit none
    ! Argument variables:
    logical(lgt),              intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
    logical(lgt), allocatable, intent(out) :: globalArray(:)            ! gathered  array at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message                   ! error message
    ! local variable
    integer(i4b)                           :: displs(0:nNodes-1)
    integer(i4b)                           :: myid

    ierr=0; message='shr_mpi_allgatherLogicalV/'

    displs(0) = 0
    do myid = 1, nNodes-1
     displs(myid) = sum(num_per_proc(0:myid-1))
    end do

    allocate(globalArray(sum(num_per_proc)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_ALLGATHERV(localArray,  num_per_proc(pid),                MPI_LOGICAL, & ! local array stuff
                        globalArray, num_per_proc(0:nNodes-1), displs, MPI_LOGICAL, & ! global array stuff
                        mpicom_route, ierr)

  END SUBROUTINE shr_mpi_allgatherLogicalV

  ! ----------------------------------
  ! ALLGATHER - integer scalar
  ! ----------------------------------
  SUBROUTINE shr_mpi_allgatherInt(localScalar,  num,        & ! input
                                  globalArray, ierr, message) ! output

    implicit none
    ! Argument variables:
    integer(i4b),              intent(in)  :: localScalar        ! local array at each proc
    integer(i4b),              intent(in)  :: num                ! number of elements per proc (i.e., size of localArray)
    integer(i4b), allocatable, intent(out) :: globalArray(:)     ! gathered  array at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message            ! error message
    ! local variable
    ! None

    ierr=0; message='shr_mpi_allgatherInt/'

    allocate(globalArray(num*nNodes), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_ALLGATHER(localScalar, num, MPI_INTEGER, & ! local array stuff
                       globalArray, num, MPI_INTEGER, & ! global array stuff
                       mpicom_route, ierr)

  END SUBROUTINE shr_mpi_allgatherInt

  ! ----------------------------------
  ! ALLGATHER - integer scalar
  ! ----------------------------------
  SUBROUTINE shr_mpi_allgatherReal(localScalar,  num,        & ! input
                                   globalArray, ierr, message) ! output

    implicit none
    ! Argument variables:
    real(dp),                  intent(in)  :: localScalar        ! local array at each proc
    integer(i4b),              intent(in)  :: num                ! number of elements per proc (i.e., size of localArray)
    real(dp), allocatable,     intent(out) :: globalArray(:)     ! gathered  array at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message            ! error message
    ! local variable
    ! None

    ierr=0; message='shr_mpi_allgatherReal/'

    allocate(globalArray(num*nNodes), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_ALLGATHER(localScalar,  num, MPI_DOUBLE_PRECISION,  & ! local array stuff
                        globalArray, num, MPI_DOUBLE_PRECISION,  & ! global array stuff
                        mpicom_route, ierr)

  END SUBROUTINE shr_mpi_allgatherReal

  ! ----------------------------------
  ! ALLGATHER - logical scalar
  ! ----------------------------------
  SUBROUTINE shr_mpi_allgatherLogical(localScalar,  num,        & ! input
                                      globalArray, ierr, message) ! output

    implicit none
    ! Argument variables:
    logical(lgt),              intent(in)  :: localScalar        ! local array at each proc
    integer(i4b),              intent(in)  :: num                ! number of elements per proc (i.e., size of localArray)
    logical(lgt), allocatable, intent(out) :: globalArray(:)     ! gathered  array at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message            ! error message
    ! local variable
    ! None

    ierr=0; message='shr_mpi_allgatherLogical/'

    allocate(globalArray(num*nNodes), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating array for [localArray]'; return; endif

    call MPI_ALLGATHER(localScalar,  num, MPI_LOGICAL,  & ! local array stuff
                        globalArray, num, MPI_LOGICAL,  & ! global array stuff
                        mpicom_route, ierr)

  END SUBROUTINE shr_mpi_allgatherLogical

  ! ----------------------------------
  ! REDUCE - real scalar
  ! ----------------------------------
  SUBROUTINE shr_mpi_reduceReal(localScalar,  method,       & ! input
                                reducedScalar, ierr, message) ! output

    implicit none
    ! Argument variables:
    real(dp),                  intent(in)  :: localScalar        ! local scalar at each proc
    character(*),              intent(in)  :: method             ! reduction operation
    real(dp),                  intent(out) :: reducedScalar      ! reduced scalar at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message            ! error message
    ! local variable
    ! None

    ierr=0; message='shr_mpi_reduceReal/'

    select case(trim(method))
      case('sum')
        call MPI_REDUCE(localScalar, reducedScalar, 1, MPI_DOUBLE_PRECISION,  &
                        MPI_SUM, root, mpicom_route, ierr)
      case('max')
        call MPI_REDUCE(localScalar, reducedScalar, 1, MPI_DOUBLE_PRECISION,  &
                        MPI_MAX, root, mpicom_route, ierr)
      case('min')
        call MPI_REDUCE(localScalar, reducedScalar, 1, MPI_DOUBLE_PRECISION,  &
                        MPI_MIN, root, mpicom_route, ierr)
      case default
        ierr=20; message=trim(message)//'reduction operation: '//trim(method)//' not available'; return
    end select

  END SUBROUTINE shr_mpi_reduceReal

  ! ----------------------------------
  ! REDUCE - integer scalar
  ! ----------------------------------
  SUBROUTINE shr_mpi_reduceInt(localScalar,  method,       & ! input
                                reducedScalar, ierr, message) ! output

    implicit none
    ! Argument variables:
    integer(i4b),              intent(in)  :: localScalar        ! local array at each proc
    character(*),              intent(in)  :: method             ! reduction operation
    integer(i4b),              intent(out) :: reducedScalar      ! reduced scalar at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message            ! error message
    ! local variable
    ! None

    ierr=0; message='shr_mpi_reduceInt/'

    select case(trim(method))
      case('sum')
        call MPI_REDUCE(localScalar, reducedScalar, 1, MPI_INTEGER,  &
                        MPI_SUM, root, mpicom_route, ierr)
      case('max')
        call MPI_REDUCE(localScalar, reducedScalar, 1, MPI_INTEGER,  &
                        MPI_MAX, root, mpicom_route, ierr)
      case('min')
        call MPI_REDUCE(localScalar, reducedScalar, 1, MPI_INTEGER,  &
                        MPI_MIN, root, mpicom_route, ierr)
      case default
        ierr=20; message=trim(message)//'reduction operation: '//trim(method)//' not available'; return
    end select

  END SUBROUTINE shr_mpi_reduceInt

  ! ----------------------------------
  ! REDUCE - real array
  ! ----------------------------------
  SUBROUTINE shr_mpi_reduceRealV(localArray,  method,       & ! input
                                 reducedArray, ierr, message) ! output

    implicit none
    ! Argument variables:
    real(dp),                  intent(in)  :: localArray(:)     ! local scalar at each proc
    character(*),              intent(in)  :: method            ! reduction operation
    real(dp),                  intent(out) :: reducedArray(:)   ! reduced scalar at root proc
    integer(i4b),              intent(out) :: ierr
    character(strLen),         intent(out) :: message           ! error message
    ! local variable
    integer(i4b)                           :: n

    ierr=0; message='shr_mpi_reduceRealV/'

    n=size(localArray)
    if (size(reducedArray)/=n)then
      ierr=20; message=trim(message)//'two array input has different sizes'; return
    end if

    select case(trim(method))
      case('sum')
        call MPI_REDUCE(localArray, reducedArray, n, MPI_DOUBLE_PRECISION,  &
                        MPI_SUM, root, mpicom_route, ierr)
      case('max')
        call MPI_REDUCE(localArray, reducedArray, n, MPI_DOUBLE_PRECISION,  &
                        MPI_MAX, root, mpicom_route, ierr)
      case('min')
        call MPI_REDUCE(localArray, reducedArray, n, MPI_DOUBLE_PRECISION,  &
                        MPI_MIN, root, mpicom_route, ierr)
      case default
        ierr=20; message=trim(message)//'reduction operation: '//trim(method)//' not available'; return
    end select

  END SUBROUTINE shr_mpi_reduceRealV

  ! -----------------------------------------------
  ! Sparse real array element transfer across tasks
  ! ------------------------------------------------
  subroutine shr_mpi_sparse_distribute(src_array, dest_task, dest_index, & ! input
                                       dest_array,                       & ! output
                                       agg, fillvalue)                     ! optional arguments
    ! Description:
    !   Every element in real src_array is sent to a specified task and index of dest_array
    !   if there are multiple source elements from different src_array(s), aggregate them with max/min/sum/mean.
    !   If there is no source elements sent to any indices in dest_array(s), put fillvalue
    implicit none
    ! Argument variables:
    real(dp),    intent(in)           :: src_array(:)   ! source array to be distributed
    integer,     intent(in)           :: dest_task(:)   ! destination task for each element
    integer,     intent(in)           :: dest_index(:)  ! destination index for for each element
    real(dp),    intent(out)          :: dest_array(:)  ! final array after transferred
    character(*),intent(in),optional  :: agg            ! aggregation methods at overlapped transferred elements
    real(dp),    intent(in),optional  :: fillvalue      ! fill value (missing value) in dest_array
    ! Local variables:
    integer                       :: i,r,k
    integer,  allocatable         :: send_counts(:)         ! total number of destination tasks in a source array per task
    integer,  allocatable         :: recv_counts(:)         ! total number of source tasks in a dest array per task
    integer,  allocatable         :: sdispls(:), rdispls(:) ! sent/recive displacemnet in MPR_ALLtoallV function
    integer                       :: total_send, total_recv ! total number of elements sents/recieved per task
    integer,  allocatable         :: send_indices(:)        ! index of packed source array
    integer,  allocatable         :: recv_indices(:)        ! index of parked recived array
    integer,  allocatable         :: cursor(:)              !
    integer,  allocatable         :: elem_counts(:)         ! number of elements sent at each dest array index
    real(dp), allocatable         :: recv_values(:)         ! packed recived array after trasfer
    real(dp), allocatable         :: send_values(:)         ! source array packed
    real(dp)                      :: fillvalue_local        ! fill values to be filled in missing dest array elements
    character(4)                  :: agg_method             ! aggregation method: sum, max, min, or mean
    integer                       :: ierr

    agg_method='sum'
    if (present(agg)) then
      agg_method=trim(agg)
    end if
    fillvalue_local=-9999._dp
    if (present(fillvalue)) then
      fillvalue_local=fillvalue
    end if

    ! Count how many to send to each task
    allocate(send_counts(nNodes), recv_counts(nNodes))
    send_counts(:) = 0
    do i = 1, size(dest_task)
      send_counts(dest_task(i)+1) = send_counts(dest_task(i)+1) + 1 ! task starts with 0, but array index starts with 1
    end do

    ! Exchange counts
    call MPI_Alltoall(send_counts, 1, MPI_INTEGER, &
                      recv_counts, 1, MPI_INTEGER, mpicom_route, ierr)
    total_send = sum(send_counts)
    total_recv = sum(recv_counts)

    ! Build displacements
    allocate(sdispls(nNodes), rdispls(nNodes))

    sdispls(1) = 0
    rdispls(1) = 0
    do r = 2, nNodes
      sdispls(r) = sdispls(r-1) + send_counts(r-1)
      rdispls(r) = rdispls(r-1) + recv_counts(r-1)
    end do

    ! Pack send buffers (value + index)
    allocate(send_values(total_send), send_indices(total_send))
    allocate(cursor(nNodes))

    cursor = sdispls
    do i = 1, size(dest_task)
      r = dest_task(i) + 1
      k = cursor(r) + 1

      send_values(k)  = src_array(i)
      send_indices(k) = dest_index(i)

      cursor(r) = cursor(r) + 1
    end do

    ! All-to-all exchange
    allocate(recv_values(total_recv))
    allocate(recv_indices(total_recv))

    call MPI_Alltoallv(send_values, send_counts, sdispls, MPI_DOUBLE_PRECISION, &
                       recv_values, recv_counts, rdispls, MPI_DOUBLE_PRECISION, &
                       mpicom_route, ierr)

    call MPI_Alltoallv(send_indices, send_counts, sdispls, MPI_INTEGER, &
                       recv_indices, recv_counts, rdispls, MPI_INTEGER, &
                       mpicom_route, ierr)

    ! aggregate contributions
    allocate(elem_counts(size(dest_array)), source=0)
    if (trim(agg_method)=="sum" .or. trim(agg_method)=="mean") then
      dest_array = 0._dp
      do r = 1, total_recv
        elem_counts(recv_indices(r)) = elem_counts(recv_indices(r)) + 1
        dest_array(recv_indices(r)) = dest_array(recv_indices(r)) + recv_values(r)
      end do
      if (trim(agg_method)=="mean") then
        where (elem_counts>0)
          dest_array = dest_array / elem_counts
        end where
      end if
    else if (trim(agg_method)=="min") then
      dest_array = huge(0._dp)
      do r = 1, total_recv
        elem_counts(recv_indices(r)) = elem_counts(recv_indices(r)) + 1
        dest_array(recv_indices(r)) = min(dest_array(recv_indices(r)), recv_values(r))
      end do
    else if (trim(agg_method)=="max") then
      dest_array = -huge(0._dp)
      do r = 1, total_recv
        elem_counts(recv_indices(r)) = elem_counts(recv_indices(r)) + 1
        dest_array(recv_indices(r)) = max(dest_array(recv_indices(r)), recv_values(r))
      end do
    end if

    ! finalize array by filling missing element
    where (elem_counts==0)
      dest_array = fillvalue_local
    end where

  end subroutine shr_mpi_sparse_distribute

  !-------------------------------------------------------------------------------
  ! PURPOSE: MPI number of tasks
  !-------------------------------------------------------------------------------
  SUBROUTINE shr_mpi_commsize(comm, ntasks, message)

    implicit none
    ! Argument variables:
    integer,              intent(in)  :: comm
    integer,              intent(out) :: ntasks
    character(*),optional,intent(in)  :: message   ! message
    ! local variable
    character(strLen),parameter       :: subName = 'shr_mpi_commsize/'
    integer(i4b)                      :: ierr

    call MPI_COMM_SIZE(comm, ntasks, ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, comm=comm, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, comm=comm, message=subName)
    endif

  END SUBROUTINE shr_mpi_commsize

  !-------------------------------------------------------------------------------
  ! PURPOSE: MPI rank
  !-------------------------------------------------------------------------------
  SUBROUTINE shr_mpi_commrank(comm, rank, message)

    implicit none
    ! Argument variables:
    integer(i4b),intent(in)            :: comm
    integer(i4b),intent(out)           :: rank
    character(*),optional,intent(in)   :: message   ! message
    ! local variable
    character(strLen),parameter        :: subName = 'shr_mpi_commrank/'
    integer(i4b)                       :: ierr

    call MPI_COMM_RANK(comm,rank,ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, comm=comm, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, comm=comm, message=subName)
    endif

  END SUBROUTINE shr_mpi_commrank

  !-------------------------------------------------------------------------------
  ! PURPOSE: MPI initialized
  !-------------------------------------------------------------------------------
  SUBROUTINE shr_mpi_initialized(flag, message)

    implicit none
    ! Argument variables:
    logical,              intent(out)  :: flag
    character(*),optional,intent(in)   :: message   ! message
    ! local variable
    character(strLen),parameter        :: subName = 'shr_mpi_initialized/'
    integer(i4b)                       :: ierr

    call MPI_INITIALIZED(flag,ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, message=subName)
    endif

  END SUBROUTINE shr_mpi_initialized

  !-------------------------------------------------------------------------------
  ! PURPOSE: MPI init
  !-------------------------------------------------------------------------------
  SUBROUTINE shr_mpi_init(comm, message)

    implicit none
    ! Argument variables:
    integer(i4b),         intent(out)  :: comm      ! communicator
    character(*),optional,intent(in)   :: message   ! message
    ! local variable
    character(strLen),parameter        :: subName = 'shr_mpi_init/'
    integer(i4b)                       :: ierr

    call MPI_INIT(ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, message=subName)
    endif

    comm = MPI_COMM_WORLD

  END SUBROUTINE shr_mpi_init

  !-------------------------------------------------------------------------------
  ! PURPOSE: mpi barrier (wait for the other procs to catch up)
  !-------------------------------------------------------------------------------
  SUBROUTINE shr_mpi_barrier(comm, message)

    implicit none
    ! Argument variables:
    integer(i4b),           intent(in) :: comm         ! communicator
    character(*), optional, intent(in) :: message      ! error message
    ! local variable
    character(strLen),parameter        :: subName = 'shr_mpi_barrier/'
    integer(i4b)                       :: ierr         ! error code

    call MPI_BARRIER(comm, ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, comm=comm, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, comm=comm, message=subName)
    endif

  END SUBROUTINE shr_mpi_barrier

  !-------------------------------------------------------------------------------
  ! PURPOSE: MPI finalize
  !-------------------------------------------------------------------------------
  SUBROUTINE shr_mpi_finalize(comm, message)

    implicit none
    ! Argument variables:
    integer(i4b),           intent(in) :: comm      ! communicator
    character(*), optional, intent(in) :: message   ! message
    ! local variable
    character(strLen)                  :: cmessage
    character(strLen),parameter        :: subName = 'shr_mpi_finalize/'
    integer(i4b)                       :: ierr

    cmessage=trim(subName)//'shr_mpi_barrier'
    call shr_mpi_barrier(comm, cmessage)

    call MPI_FINALIZE(ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, comm=comm, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, comm=comm, message=subName)
    endif

  END SUBROUTINE shr_mpi_finalize

  !-------------------------------------------------------------------------------
  ! PURPOSE: layer on MPI error checking
  !-------------------------------------------------------------------------------
  SUBROUTINE shr_mpi_chkerr(ierr, comm, message)

    implicit none
    ! Argument variables:
    integer(i4b),           intent(in) :: ierr    ! input MPI error code
    integer(i4b), optional, intent(in) :: comm    ! communicator
    character(*), optional, intent(in) :: message ! message
    ! local variable
    character(strLen),parameter        :: subName = 'shr_mpi_chkerr/'
    character(strLen)                  :: cmessage
    character(strLen)                  :: errMsg
    integer(i4b)                       :: errLen
    integer(i4b)                       :: jerr

    if (ierr /= MPI_SUCCESS) then
      call MPI_ERROR_STRING(ierr, errMsg, errLen, jerr)
      if (present(message)) then
        cmessage = trim(subName)//trim(message)//errMsg(1:errLen)
      else
        cmessage = trim(subName)//errMsg(1:errLen)
      endif
      if (present(comm)) then
        call shr_mpi_abort(cmessage, ierr, comm)
      else
        call shr_mpi_abort(cmessage, ierr)
      endif
    endif

  END SUBROUTINE shr_mpi_chkerr

  !-------------------------------------------------------------------------------
  ! PURPOSE: MPI abort
  !-------------------------------------------------------------------------------
  SUBROUTINE shr_mpi_abort(message, ierr, comm)

    implicit none
    ! Argument variables:
    character(*),           intent(in) :: message  ! message
    integer(i4b),           intent(in) :: ierr     ! error code
    integer(i4b), optional, intent(in) :: comm     ! communicator
    ! local variables
    character(strLen),parameter :: subName = 'shr_mpi_abort/'
    integer(i4b)                :: jerr

    write(iulog,*) trim(subName),trim(message)
    flush(iulog)

    if (present(comm)) then
      call MPI_ABORT(comm, ierr, jerr)
    else
      call MPI_ABORT(MPI_COMM_WORLD, ierr, jerr)
    endif

    stop

  END SUBROUTINE shr_mpi_abort

  !-------------------------------------------------------------------------------
  ! handle non-MPI error codes
  !-------------------------------------------------------------------------------
  SUBROUTINE mpi_handle_err(ierr,pid)

    implicit none
    ! Argument variables:
    integer(i4b),intent(in):: ierr   ! error code
    integer(i4b),intent(in):: pid    ! process ID
    ! local variables
    integer(i4b)           :: jerr   ! error code with message string call
    integer(i4b)           :: errLen ! length of error message
    character(len=strLen)  :: errMsg ! error message

    ! check errors
    if(ierr/=0)then
      ! get error string
      call MPI_Error_String(ierr, errMsg, errLen, jerr)
      if(jerr==0)errMsg='problemIdentifyingErrorMessage'
      if(errLen>strLen)errMsg='errorMessageLengthTooLong'

      ! include process ID
      write(iulog,'(a,1x,i4)') 'FATAL ERROR (MPI): '//trim(errMsg)//' for process ID ', pid

      ! finalize MPI
      call MPI_FINALIZE(jerr)
      flush(iulog)
      stop
    endif

  END SUBROUTINE mpi_handle_err

END MODULE mpi_utils
