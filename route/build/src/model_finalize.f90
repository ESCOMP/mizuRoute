MODULE model_finalize

USE nrtype,     ONLY: i4b
USE public_var, ONLY: iulog            ! i/o logical unit number

implicit none

private

public :: finalize
public :: handle_err

CONTAINS

 ! *********************************************************************
 ! public subroutine: finalize model
 ! *********************************************************************
 SUBROUTINE finalize()
  implicit none
  write(iulog,'(a)') new_line('a'), '--------------------'
  write(iulog,'(a)')                'Finished simulation'
  write(iulog,'(a)')                '--------------------'
  stop
 END SUBROUTINE finalize


 ! *********************************************************************
 ! public subroutine: error handling
 ! *********************************************************************
 SUBROUTINE handle_err(err,message)
 implicit none
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 if(err/=0)then
   write(iulog,*) 'FATAL ERROR: '//trim(message)
   call flush(6)
   stop
 endif
 END SUBROUTINE handle_err


END MODULE model_finalize
