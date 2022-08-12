MODULE model_finalize

USE nrtype
USE public_var, ONLY: iulog            ! i/o logical unit number
USE public_var, ONLY: qmodOption       !
USE globalData, ONLY: gage_obs_data
USE globalData, ONLY: rch_qtake_data

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
  integer(i4b)      :: ierr             ! error code
  character(strLen) :: cmessage         ! error message

  if (qmodOption==1) then
    call gage_obs_data%closeNC(ierr, cmessage)
  else if (qmodOption==2) then
    call rch_qtake_data%closeNC(ierr, cmessage)
  end if

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
