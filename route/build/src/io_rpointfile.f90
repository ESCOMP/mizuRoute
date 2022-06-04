MODULE io_rpointfile

USE nrtype
USE public_var,  ONLY: gageOutput           ! ascii containing last restart and history files
USE public_var,  ONLY: restart_dir              ! restart directory
USE public_var,  ONLY: rpntfil                  !
USE globalData,  ONLY: masterproc
USE globalData,  ONLY: hfileout, hfileout_gage
USE globalData,  ONLY: rfileout

implicit none

private
public::io_rpfile

CONTAINS

  ! *********************************************************************
  ! public subroutine: main routine to define new output file
  ! *********************************************************************
  SUBROUTINE io_rpfile(mode, ierr, message)

    implicit none
    ! argument variables
    character(1),   intent(in)     :: mode        ! io mode: 'r' or 'w'
    integer(i4b),   intent(out)    :: ierr        ! error code
    character(*),   intent(out)    :: message     ! error message

    ierr=0; message='io_rpfile/'

    select case(mode)
      case("w")
        if (masterproc) then
          open(1, file = trim(restart_dir)//trim(rpntfil), status='replace', action='write')
          write(1,'(a)') trim(rfileout)
          write(1,'(a)') trim(hfileout)
          if (gageOutput) then
            write(1,'(a)') trim(hfileout_gage)
          end if
          close(1)
        end if
      case("r")
        open(1, file = trim(restart_dir)//trim(rpntfil), status='old', action='read')
        read(1, '(A)') rfileout
        read(1, '(A)') hfileout
        if (gageOutput) then
          read(1, '(A)') hfileout_gage
        end if
        close(1)
      case default
        ierr=20; message=trim(message)//"Mode argument must be in 'w','r'!"
    end select

  END SUBROUTINE io_rpfile

END MODULE io_rpointfile
