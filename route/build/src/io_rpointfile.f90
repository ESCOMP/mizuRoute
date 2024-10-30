MODULE io_rpointfile

USE nrtype
USE public_var,  ONLY: outputAtGage             ! ascii containing last restart and history files
USE public_var,  ONLY: restart_dir              ! restart directory
USE public_var,  ONLY: rpntfil                  !
USE public_var,  ONLY: secprmin, secprhour
USE globalData,  ONLY: masterproc
USE globalData,  ONLY: hfileout, hfileout_gage
USE globalData,  ONLY: rfileout
USE datetime_data, ONLY: datetime               ! datetime data class

implicit none

private
public::io_rpfile

CONTAINS

  ! *********************************************************************
  ! public subroutine: main routine to define new output file
  ! *********************************************************************
  SUBROUTINE io_rpfile(mode, ierr, message, curDatetime)

    implicit none
    ! argument variables
    character(1),   intent(in)           :: mode          ! io mode: 'r' or 'w'
    integer(i4b),   intent(out)          :: ierr          ! error code
    character(*),   intent(out)          :: message       ! error message
    type(datetime), intent(in), optional :: curDatetime   ! current datetime
    integer                              :: sec_in_day    ! time in second
    character(len=17)                    :: timestamp     ! datetime string in yyyy-mm-dd-sssss
    character(len=strLen)                :: rpntfil_path  ! rpointer file path

    ierr=0; message='io_rpfile/'

    if (present(curDatetime)) then
      sec_in_day = curDatetime%hour()*nint(secprhour)+curDatetime%minute()*nint(secprmin)+nint(curDatetime%sec())
      write(timestamp,'(".",I4.4,"-",I2.2,"-",I2.2,"-",I5.5)') &
            curDatetime%year(),curDatetime%month(),curDatetime%day(),sec_in_day
      rpntfil_path = trim(restart_dir)//trim(rpntfil)//timestamp
    else
      rpntfil_path = trim(restart_dir)//trim(rpntfil)
    end if

    select case(mode)
      case("w")
        if (masterproc) then
          open(1, file = trim(rpntfil_path), status='replace', action='write')
          write(1,'(a)') trim(rfileout)
          write(1,'(a)') trim(hfileout)
          if (outputAtGage) then
            write(1,'(a)') trim(hfileout_gage)
          end if
          close(1)
        end if
      case("r")
        open(1, file = trim(rpntfil_path), status='old', action='read')
        read(1, '(A)') rfileout
        read(1, '(A)') hfileout
        if (outputAtGage) then
          read(1, '(A)') hfileout_gage
        end if
        close(1)
      case default
        ierr=20; message=trim(message)//"Mode argument must be in 'w','r'!"
    end select

  END SUBROUTINE io_rpfile

END MODULE io_rpointfile
