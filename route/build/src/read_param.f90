MODULE read_param_module

 USE nrtype, ONLY : i4b,   &
                    strLen

 private
 public :: read_param

 CONTAINS

  SUBROUTINE read_param(fname, ierr, message)
   USE ascii_util_module, ONLY : file_open          ! open file (performs a few checks as well)
   USE globalData,        ONLY : fshape, tscale, &  ! basin IRF routing parameters
                                 velo, diff,     &  ! IRF routing parameters
                                 mann_n, wscale     ! KWT routing parameters
   implicit none
   ! input variables
   character(*),intent(in)    :: fname              ! parameter namelist name
   ! output variables
   integer(i4b),intent(out)   :: ierr               ! error code
   character(*),intent(out)   :: message            ! error message
   ! local variables
   character(len=strLen)      :: cmessage           ! error message from subroutine
   integer(i4b)               :: iunit              ! file unit
   namelist /HSLOPE/fshape,tscale  ! route simulated runoff through the local basin
   namelist /IRF_UH/velo,diff      ! route delayed runoff through river network with St.Venant UH
   namelist /KWT/mann_n,wscale     ! route kinematic waves through the river network

   ! error control initialization
   ierr=0; message='read_param'

   ! read the name list
   call file_open(trim(fname),iunit, ierr, cmessage)
   if(ierr/=0)then; ierr=30; message=trim(cmessage)//': '//trim(fname); return;endif
   read(iunit, nml=HSLOPE)              ! basin IRF routing parameters
   read(iunit, nml=IRF_UH)              ! IRF routing parameters
   read(iunit, nml=KWT)                 ! kinematic wave routing parameters
   close(iunit)

  END SUBROUTINE read_param

END MODULE read_param_module
