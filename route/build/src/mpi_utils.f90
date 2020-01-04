Module mpi_mod
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

  public :: shr_mpi_bcast
  public :: shr_mpi_gatherV
  public :: shr_mpi_scatterV
  public :: shr_mpi_allgather

  public :: shr_mpi_chkerr
  public :: shr_mpi_commsize
  public :: shr_mpi_commrank
  public :: shr_mpi_initialized
  public :: shr_mpi_barrier
  public :: shr_mpi_abort
  public :: shr_mpi_init
  public :: shr_mpi_finalize

  interface shr_mpi_bcast; module procedure &
    shr_mpi_bcastInt,    &
    shr_mpi_bcastReal,   &
    shr_mpi_bcastLogical
  end interface

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

  interface shr_mpi_allgather ; module procedure &
    shr_mpi_allgatherIntV,    &
    shr_mpi_allgatherRealV,   &
    shr_mpi_allgatherLogicalV,&
    shr_mpi_allgatherInt,     &
    shr_mpi_allgatherReal,    &
    shr_mpi_allgatherLogical
  end interface

  integer(i4b), parameter :: send_data_tag=2001
  integer(i4b), parameter :: return_data_tag=2002
  integer(i4b)            :: status(MPI_STATUS_SIZE)

CONTAINS

  ! ----------------------------------
  ! BROADCAST - integer allocatable array
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastInt(allocArray,   & ! inout:  array to be broadcasted to each proc
                              ierr, message)  ! output: error handling
    implicit none
    ! Input
    integer(i4b), allocatable, intent(inout) :: allocArray(:) ! inout:  array to be sent to proc
    ! Output error handling variables
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
  ! BROADCAST - real allocatable array
  ! ----------------------------------
  SUBROUTINE shr_mpi_bcastReal(allocArray,   & ! input:  array to be broadcasted to each proc
                               ierr, message)  ! output: error handling
    implicit none
    ! Input
    real(dp), allocatable, intent(inout) :: allocArray(:) ! inout:  array to be sent to proc
    ! Output error handling variables
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
    ! Input
    logical(lgt), allocatable, intent(inout) :: allocArray(:)   ! input:  array to be sent to proc
    ! Output error handling variables
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
  ! SCATTERV - 1D integer array
  ! ----------------------------------
  SUBROUTINE shr_mpi_scatterIntV(globalArray, num_per_proc, & ! input
                                localArray, ierr, message)   ! output
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
    ! Input
    real(dp),                  intent(in)  :: localArray(:)             ! local array at each proc
    integer(i4b),              intent(in)  :: num_per_proc(0:nNodes-1)  ! number of elements per proc (i.e., size of localArray)
    ! Output
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
    USE globalData,  ONLY: nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    integer(i4b),              intent(in)  :: localScalar        ! local array at each proc
    integer(i4b),              intent(in)  :: num                ! number of elements per proc (i.e., size of localArray)
    ! Output
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
    USE globalData,  ONLY: nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    real(dp),                  intent(in)  :: localScalar        ! local array at each proc
    integer(i4b),              intent(in)  :: num                ! number of elements per proc (i.e., size of localArray)
    ! Output
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
    USE globalData,  ONLY: nNodes
    USE public_var,  ONLY: root
    implicit none
    ! Input
    logical(lgt),              intent(in)  :: localScalar        ! local array at each proc
    integer(i4b),              intent(in)  :: num                ! number of elements per proc (i.e., size of localArray)
    ! Output
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


  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_commsize(comm, ntasks, message)

    IMPLICIT none

    !----- arguments ---
    integer,              intent(in)  :: comm
    integer,              intent(out) :: ntasks
    character(*),optional,intent(in)  :: message   ! message

    !----- local ---
    character(strLen),parameter       :: subName = 'shr_mpi_commsize/'
    integer(i4b)                      :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI number of tasks
    !-------------------------------------------------------------------------------

    call MPI_COMM_SIZE(comm, ntasks, ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, comm=comm, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, comm=comm, message=subName)
    endif

  END SUBROUTINE shr_mpi_commsize

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_commrank(comm, rank, message)

    IMPLICIT none

    !----- arguments ---
    integer(i4b),intent(in)            :: comm
    integer(i4b),intent(out)           :: rank
    character(*),optional,intent(in)   :: message   ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_commrank/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI rank
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(comm,rank,ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, comm=comm, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, comm=comm, message=subName)
    endif

  END SUBROUTINE shr_mpi_commrank

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_initialized(flag, message)

    IMPLICIT none

    !----- arguments ---
    logical,              intent(out)  :: flag
    character(*),optional,intent(in)   :: message   ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_initialized/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI initialized
    !-------------------------------------------------------------------------------

    call MPI_INITIALIZED(flag,ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, message=subName)
    endif

  END SUBROUTINE shr_mpi_initialized

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_init(comm, message)

    IMPLICIT none

    !----- arguments ---
    integer(i4b),         intent(out)  :: comm      ! communicator
    character(*),optional,intent(in)   :: message   ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_init/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI init
    !-------------------------------------------------------------------------------

    call MPI_INIT(ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, message=subName)
    endif

    comm = MPI_COMM_WORLD

  END SUBROUTINE shr_mpi_init

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_barrier(comm, message)

    IMPLICIT none

    !----- argument ---
    integer(i4b),           intent(in) :: comm         ! communicator
    character(*), optional, intent(in) :: message      ! error message
    !----- local variables ---
    character(strLen),parameter        :: subName = 'shr_mpi_barrier/'
    integer(i4b)                       :: ierr         ! error code

    call MPI_BARRIER(comm, ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, comm=comm, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, comm=comm, message=subName)
    endif

  END SUBROUTINE shr_mpi_barrier

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_finalize(comm, message)

    IMPLICIT none

    !----- arguments ---
    integer(i4b),           intent(in) :: comm      ! communicator
    character(*), optional, intent(in) :: message   ! message

    !----- local ---
    character(strLen)                  :: cmessage
    character(strLen),parameter        :: subName = 'shr_mpi_finalize/'
    integer(i4b)                       :: ierr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI finalize
    !-------------------------------------------------------------------------------
    cmessage=trim(subName)//'shr_mpi_barrier'
    call shr_mpi_barrier(comm, cmessage)

    call MPI_FINALIZE(ierr)
    if (present(message)) then
       call shr_mpi_chkerr(ierr, comm=comm, message=subName//trim(message))
    else
       call shr_mpi_chkerr(ierr, comm=comm, message=subName)
    endif

  END SUBROUTINE shr_mpi_finalize


  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_chkerr(ierr, comm, message)

    IMPLICIT none

    !----- arguments ---
    integer(i4b),           intent(in) :: ierr    ! input MPI error code
    integer(i4b), optional, intent(in) :: comm    ! communicator
    character(*), optional, intent(in) :: message ! message

    !----- local ---
    character(strLen),parameter        :: subName = 'shr_mpi_chkerr/'
    character(strLen)                  :: cmessage
    character(strLen)                  :: errMsg
    integer(i4b)                       :: errLen
    integer(i4b)                       :: jerr

    !-------------------------------------------------------------------------------
    ! PURPOSE: layer on MPI error checking
    !-------------------------------------------------------------------------------

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

  !===============================================================================
  !===============================================================================

  SUBROUTINE shr_mpi_abort(message, ierr, comm)

    IMPLICIT none

    !----- arguments ---
    character(*),           intent(in) :: message  ! message
    integer(i4b),           intent(in) :: ierr     ! error code
    integer(i4b), optional, intent(in) :: comm     ! communicator

    !----- local ---
    character(strLen),parameter :: subName = 'shr_mpi_abort/'
    integer(i4b)                :: jerr

    !-------------------------------------------------------------------------------
    ! PURPOSE: MPI abort
    !-------------------------------------------------------------------------------

    write(iulog,*) trim(subName),trim(message),ierr
    call flush(6)

    if (present(comm)) then
      call MPI_ABORT(comm, ierr, jerr)
    else
      call MPI_ABORT(MPI_COMM_WORLD, ierr, jerr)
    endif

    stop

  END SUBROUTINE shr_mpi_abort

  !===============================================================================
  !===============================================================================

  SUBROUTINE mpi_handle_err(ierr,pid)
  ! handle non-MPI error codes
  implicit none
  ! arguments
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
   write(*,'(a,1x,i4)') 'FATAL ERROR (MPI): '//trim(errMsg)//' for process ID ', pid

   ! finalize MPI
   call MPI_FINALIZE(jerr)
   call flush(6)
   stop

  endif

  END SUBROUTINE mpi_handle_err

  !===============================================================================
  !===============================================================================

END MODULE mpi_mod
