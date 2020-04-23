MODULE model_utils

! data types
USE nrtype,    ONLY: i4b        ! variable types, etc.

implicit none

! privacy -- everything private unless declared explicitly
private

public :: model_finalize
public :: handle_err

CONTAINS

 ! *********************************************************************
 ! public subroutine: finalize model
 ! *********************************************************************
 SUBROUTINE model_finalize(comm)

  USE perf_mod,   ONLY: t_prf            ! timing output
  USE perf_mod,   ONLY: t_finalizef      ! finalize timing routines
  USE globalData, ONLY: masterproc       ! root proc logical
  USE public_var, ONLY: iulog            ! i/o logical unit number
  USE mpi_mod,    ONLY: shr_mpi_finalize ! mpi utilities: shut down mpi

  implicit none

  integer(i4b), intent(in) :: comm   ! communicator

  if (masterproc) then
    write(iulog,*) '----------------------'
    write(iulog,*) ' SUCCESSFUL EXECUTION '
    write(iulog,*) '----------------------'
  end if

  call t_prf(mpicom=comm)
  call t_finalizef()
  call shr_mpi_finalize(comm)

  stop

 END SUBROUTINE model_finalize


 ! *********************************************************************
 ! public subroutine: error handling (maybe place different module)
 ! *********************************************************************
 SUBROUTINE handle_err(err,message)

  USE mpi_mod,    ONLY: shr_mpi_abort     ! mpi utilities: abort mpi

  implicit none

  integer(i4b),intent(in)::err             ! error code
  character(*),intent(in)::message         ! error message

  if(err/=0)then
   call shr_mpi_abort('FATAL ERROR: '//trim(message), err)
  endif

 END SUBROUTINE handle_err


END MODULE model_utils
