module model_setup

! variable types
USE nrtype                                    ! variable types, etc.
! global data
USE public_var
! provide access to desired subroutines...
USE mpi                                       ! MPI

implicit none

private
public :: model_setup

contains

 subroutine model_setup(pid, ierr, cmessage)

  ! subroutines: populate metadata
  USE popMetadat_module,   only : popMetadat       ! populate metadata
  ! subroutines: model control
  USE read_control_module, only : read_control     ! read the control file
  USE read_param_module,   only : read_param       ! read the routing parameters

  implicit none

  integer(i4b), intent(in)    :: pid                 ! process id
  ! output: error control
  integer(i4b), intent(out)   :: ierr             ! error code
  character(*), intent(out)   :: message          ! error message
  ! local variables
  character(len=strLen)       :: cfile_name          ! name of the control file
  integer(i4b)                :: iens                ! ensemble member
  integer(i4b)                :: nHRU                ! number of HRUs
  integer(i4b)                :: nRch                ! number of desired reaches
  character(len=strLen)       :: cmessage            ! error message of downwind routine

  ! initialize error control
  ierr=0; message='model_setup/'

  ! if the master node
  if (pid==0) then

   ! populate the metadata files
   call popMetadat(ierr,cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! get command-line argument defining the full path to the control file
   call getarg(1,cfile_name)
   if(len_trim(cfile_name)==0) ierr=50;message=trim(message)//'need to supply name of the control file as a command-line argument'; return

   ! read the control file
   call read_control(trim(cfile_name), ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! read the routing parameter namelist
   call read_param(trim(ancil_dir)//trim(param_nml),ierr,cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif  ! if the master node
 end subroutine model_setup

end module model_setup
