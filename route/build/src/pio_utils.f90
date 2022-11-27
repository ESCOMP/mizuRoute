MODULE pio_utils

  USE mpi
  USE nrtype
  USE pio

  implicit none

  private
  public::pio_sys_init             ! Define pio system descriptor (iosystem_desc_t)
  public::pio_decomp               ! Define decomposition (io_desc_t)
  public::createFile               ! Create netcdf
  public::def_dim                  ! Define dimension
  public::def_var                  ! Define variable
  public::end_def                  ! End definition mode
  public::inq_dim_len              ! inquire dimension size
  public::openFile                 ! Open netcdf
  public::freeDecomp               ! free decomposition (io_desc_t)
  public::closeFile                ! close netcdf (if it's open) and clean file_desc_t
  public::finalizeSystem           ! free pio system descriptor (iosystem_desc_t)
  public::sync_file                ! Write out data into disk
  public::write_scalar_netcdf      ! write non-distributed data
  public::write_netcdf             ! write non-distributed data
  public::write_pnetcdf            ! write distributed data without record dimension
  public::write_pnetcdf_recdim     ! write distributed data at a specified index of record dimension

  ! public PIO parameter
  integer(i4b),parameter,public :: ncd_int       = pio_int
  integer(i4b),parameter,public :: ncd_float     = pio_real
  integer(i4b),parameter,public :: ncd_double    = pio_double
  integer(i4b),parameter,public :: ncd_char      = pio_char
  integer(i4b),parameter,public :: ncd_global    = pio_global
  integer(i4b),parameter,public :: ncd_write     = pio_write
  integer(i4b),parameter,public :: ncd_nowrite   = pio_nowrite
  integer(i4b),parameter,public :: ncd_clobber   = pio_clobber
  integer(i4b),parameter,public :: ncd_noclobber = pio_noclobber
  integer(i4b),parameter,public :: ncd_nofill    = pio_nofill
  integer(i4b),parameter,public :: ncd_unlimited = pio_unlimited

  ! PIO types needed for public interface calls
  public iosystem_desc_t
  public file_desc_t
  public var_desc_t
  public io_desc_t

  INTERFACE write_scalar_netcdf
    module procedure write_int_scalar
    module procedure write_double_scalar
    module procedure write_char_scalar
  END INTERFACE

  INTERFACE write_netcdf
    module procedure write_int_array1D
    module procedure write_int_array2D
    module procedure write_double_array1D
    module procedure write_double_array2D
  END INTERFACE

  INTERFACE write_pnetcdf
    module procedure write_int_darray1D
    module procedure write_int_darray2D
    module procedure write_int_darray3D
    module procedure write_double_darray1D
    module procedure write_double_darray2D
    module procedure write_double_darray3D
  END INTERFACE

  INTERFACE write_pnetcdf_recdim
    module procedure write_int_darray2D_recdim
    module procedure write_int_darray3D_recdim
    module procedure write_int_darray4D_recdim
    module procedure write_float_darray2D_recdim
    module procedure write_float_darray3D_recdim
    module procedure write_float_darray4D_recdim
    module procedure write_double_darray2D_recdim
    module procedure write_double_darray3D_recdim
    module procedure write_double_darray4D_recdim
  END INTERFACE

CONTAINS

  ! *********************************************************************
  ! subroutine: initialize ParallelIO system
  ! *********************************************************************
  SUBROUTINE pio_sys_init(pid,       &  !input
                          comm,      &  !input
                          stride,    &  !input
                          nIOtasks,  &  !input
                          rearranger,&  !input
                          idxBase,   &  !input
                          pioIOsystem)  !output
    implicit none
    ! ARGUMENTS:
    integer(i4b),          intent(in)  :: pid
    integer(i4b),          intent(in)  :: comm
    integer(i4b),          intent(in)  :: stride       ! stride in MPI rank between IO tasks.
    integer(i4b),          intent(in)  :: nIOtasks     ! number of IO tasks
    integer(i4b),          intent(in)  :: rearranger   !
    integer(i4b),          intent(in)  :: idxBase      ! Start index of IO tasks.
    type(iosystem_desc_t), intent(out) :: pioIoSystem  ! pio system descriptor
    ! Local variables
    integer(i4b)                       :: nAggregator  ! MPI aggregator count

    ! set up PIO for rest of parameters
    nAggregator = 0
    call pio_init(pid,              & ! input: MPI rank
                  comm,             & ! input: MPI communicator
                  nIOtasks,         & ! input: Number of iotasks
                  nAggregator,      & ! input: number of aggregators to use
                  stride,           & ! input: MPI rank stride between IO tasks
                  rearranger,       & ! input: do not use any form of rearrangement
                  pioIoSystem,      & ! output: pio system descriptor
                  base=idxBase)       ! base (optional argument)

  END SUBROUTINE pio_sys_init

  ! *********************************************************************
  ! subroutine: PIO domain decomposition data
  ! *********************************************************************
  SUBROUTINE pio_decomp(pioIoSystem,  & ! input: pio system descriptor
                        piotype,      & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                        dimLen,       & ! input: dimension length for global array
                        gidx_local,   & ! input: local global-index array at one dimension
                        iodesc)         ! output:
    ! Details:
    ! Initialize domain descomposition for PIO (gidx_local is defined only in the 1st dimension.
    !
    ! Note:
    ! gidx_local doesn't need to be sequential but it will work a lot better if it's monotonically increasing on each task

    implicit none
    ! ARGUMENTS:
    type(iosystem_desc_t),intent(inout) :: pioIoSystem    ! pio system descriptor
    integer(i4b),         intent(in)    :: piotype        ! pio data type
    integer(i4b),         intent(in)    :: dimLen(:)      ! length of dimension for global array
    integer(i4b),         intent(in)    :: gidx_local(:)  ! indices of global array describing the decomposition of the data
    type(io_desc_t),      intent(inout) :: iodesc         ! io descriptor handle that is generated in PIO_initdecomp
    ! Local variables
    integer(i4b),         allocatable   :: compdof(:)     ! indices of global array describing the decomposition of the data
    integer(i4b),         allocatable   :: compdof2d(:,:) ! indices of global array describing the decomposition of the data
    integer(i4b)                        :: ndims          ! number of dimensions
    integer(i4b)                        :: totnum         ! total number of global array elements
    integer(i4b)                        :: gsize          ! global array 1st dimension size
    integer(i4b)                        :: lsize          ! local array 1st dimension size
    integer(i4b)                        :: ix             ! counter
    integer(i4b)                        :: nn             !

    ndims = size(dimLen)
    gsize = dimLen(1)        ! 1st dimension size for global array
    lsize = size(gidx_local) ! local size = 1st dimension size for local array

    totnum = lsize
    nn     = 1

    ! if array dimension is more than 1
    if (ndims>1) then
      do ix = 2, ndims
        totnum = totnum*dimLen(ix)
      end do

      do ix = 2, ndims
        nn = nn*dimLen(ix)
      enddo

      ! fill global indices in higher dimension
      allocate(compdof2d(lsize,nn))
      do ix = 1,nn
        compdof2d(1:lsize, ix) = (ix-1)*gsize + gidx_local(1:lsize)
      end do

      ! take care of ghost indices (0 index value)
      do ix = 1, lsize
       if (gidx_local(ix) == 0) then
         compdof2d(ix,:) = 0
       endif
      enddo

!   do ix = 1,lsize
!    write(*,*) (compdof2d(ix,jx),jx=1,nn)
!   enddo

      allocate(compdof(totnum))
      compdof = reshape(compdof2d,[totnum])

    else

      allocate(compdof(totnum))
      compdof(1:lsize) = gidx_local(1:lsize)

    endif

!    write(*,*) (compdof(ix),ix=1,totnum)

    call pio_initdecomp(pioIoSystem,      & ! input: pio system descriptor
                        piotype,          & ! input: data type
                        dimLen(1:ndims),  & ! input: dimension length in global array
                        compdof,          & ! input: decomposition for local array
                        iodesc)             ! output:

  END SUBROUTINE pio_decomp

  !-----------------------------------------------------------------------
  FUNCTION iotype_id(iotype,        &  ! input:  netcdf type name
                     ierr, message)    ! output: error handling

    ! Valid netcdf type: "netcdf", "netcdf4c", "netcdf4p", "pnetcdf"
    !
    !  - pnetcdf  ==> pio_iotype_pnetcdf  = 1   Parallel Netcdf  (parallel)
    !  - netcdf   ==> pio_iotype_netcdf   = 2   Netcdf3 Classic format (serial)
    !  - netcdf4c ==> pio_iotype_netcdf4c = 3   NetCDF4 (HDF5) compressed format (serial)
    !  - netcdf4p ==> pio_iotype_NETCDF4p = 4   NetCDF4 (HDF5) parallel

    implicit none
    ! ARGUMENTS:
    character(*), intent(in)  :: iotype    ! Input: netcdf type
    ! local variable
    integer(i4b), intent(out) :: ierr      ! error code
    character(*), intent(out) :: message   ! error message
    integer(i4b)              :: iotype_id !

    ierr=0; message='iotype_id'

    select case(trim(iotype))
      case('netcdf');   iotype_id = pio_iotype_netcdf
      case('pnetcdf');  iotype_id = pio_iotype_pnetcdf
      case('netcdf4c'); iotype_id = pio_iotype_netcdf4c
      case('netcdf4p'); iotype_id = pio_iotype_NETCDF4p
      case default
        message=trim(message)//'unexpected netcdf type name '//trim(iotype)
        ierr=20; return
    end select

  END FUNCTION iotype_id

  !-----------------------------------------------------------------------
  FUNCTION ioformat_id(ioformat,     &  ! input:  netcdf format name, default="64bit_offset"
                       ierr, message)   ! output: error handling
    ! Valid netcdf format: "64bit_offset"
    !
    !  - 64bit_offset,

    implicit none
    ! ARGUMENTS:
    character(*), intent(in)  :: ioformat    ! Input: netcdf format name
    integer(i4b), intent(out) :: ierr        ! error code
    character(*), intent(out) :: message     ! error message
    integer(i4b)              :: ioformat_id ! netcdf format id

    ierr=0; message='ioformat_id'

    select case(trim(ioformat))
      case('64bit_data');   ioformat_id = PIO_64BIT_DATA
      case('64bit_offset'); ioformat_id = PIO_64BIT_OFFSET
      case default
        message=trim(message)//'unexpected netcdf format name '//trim(ioformat)
        ierr=20; return
    end select

  END FUNCTION ioformat_id

  !-----------------------------------------------------------------------
  SUBROUTINE createFile(pioIoSystem,   &  ! inout:  pio system descriptor (initialized by pio_init)
                        fileName,      &  ! input:  output netcdf name
                        netcdf_type,   &  ! input:  netcdf type name,   default="netcdf"
                        netcdf_format, &  ! input:  netcdf format name, default="64bit_offset"
                        pioFileDesc,   &  ! output: pio file descriptor (use in writing function)
                        ierr, message)    ! output: error handling

    implicit none
    ! ARGUMENTS:
    type(iosystem_desc_t), intent(inout):: pioIoSystem   ! input: pio system descriptor
    character(*),          intent(in)   :: fileName      ! input: netcdf name
    character(*),          intent(in)   :: netcdf_type   ! input: netcdf type name
    character(*),          intent(in)   :: netcdf_format ! input: netcdf format name
    type(file_desc_t),     intent(out)  :: pioFileDesc   ! contains data identifying the file.
    integer(i4b),          intent(out)  :: ierr          ! error code
    character(*),          intent(out)  :: message       ! error message
    ! local variable
    integer(i4b)                        :: iotype        ! netcdf type ID
    integer(i4b)                        :: ioformat      ! netcdf format ID
    integer(i4b)                        :: mode
    character(len=strLen)               :: cmessage      ! error message from subroutine

   iotype = iotype_id(netcdf_type, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ioformat = ioformat_id(netcdf_format, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   mode = ior(ioformat,PIO_CLOBBER)

   ierr = pio_createfile(pioIOsystem,    & ! input:
                         pioFileDesc,    & ! output:
                         iotype,         & ! input:
                         trim(fileName), & ! input: input file name
                         mode)             ! append
   if(ierr/=pio_noerr)then; message=trim(message)//'cannot create netCDF'; return; endif

  END SUBROUTINE createFile

  !-----------------------------------------------------------------------
  SUBROUTINE openFile(pioIoSystem, pioFileDesc, fname, netcdf_type, mode, fileOpen, ierr, message)
    !
    ! DESCRIPTION:
    ! Open a NetCDF PIO file
    !
    implicit none
    ! ARGUMENTS:
    type(iosystem_desc_t), intent(inout) :: pioIoSystem
    type(file_desc_t),     intent(inout) :: pioFileDesc  ! contains data identifying the file.
    character(*),          intent(in)    :: fname        ! filename
    character(*),          intent(in)    :: netcdf_type  ! input: netcdf type name
    integer(i4b),          intent(in)    :: mode         ! file mode: pio_nowrite or pio_write
    logical(lgt),          intent(inout) :: fileOpen     ! logical to indicate file open or not
    integer(i4b),          intent(out)   :: ierr         ! error status
    character(*),          intent(out)   :: message      ! error message
    ! local variable
    integer(i4b)                         :: iotype       ! netcdf type ID
    character(len=strLen)                :: cmessage     ! error message from subroutine

    ierr=0; message='openFile/'

    iotype = iotype_id(netcdf_type, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ierr = pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode)
    if(ierr/=pio_noerr)then; message=trim(message)//'Could not open netCDF'; return; endif

    fileOpen = .true.

  END SUBROUTINE openFile

  !-----------------------------------------------------------------------
  SUBROUTINE closeFile(pioFileDesc, fileOpen)
    ! !DESCRIPTION:
    ! Close a NetCDF PIO file
    !
    implicit none
    ! ARGUMENTS:
    type(file_desc_t), intent(inout) :: pioFileDesc   ! PIO file handle to close
    logical(lgt),      intent(inout) :: fileOpen      ! logical to indicate file open or not

    call pio_closefile(pioFileDesc)

    fileOpen = .false.

  END SUBROUTINE closeFile

  !-----------------------------------------------------------------------
  SUBROUTINE freeDecomp(pioFileDesc, iodesc)
    ! !DESCRIPTION:
    ! Free decomposition
    !
    implicit none
    ! ARGUMENTS:
    type(file_desc_t), intent(inout) :: pioFileDesc   ! PIO file handle to close
    type(io_desc_t),   intent(inout) :: iodesc

    call pio_freedecomp(pioFileDesc, ioDesc)

  END SUBROUTINE freeDecomp

  !-----------------------------------------------------------------------
  SUBROUTINE finalizeSystem(pioIOsystem)
    ! !DESCRIPTION:
    ! Free IO system and releasing all resources
    !
    implicit none
    ! ARGUMENTS:
    type(iosystem_desc_t),intent(inout) :: pioIOsystem   !
    integer(i4b)                        :: ierr

    call pio_finalize(pioIOsystem, ierr)

  END SUBROUTINE finalizeSystem

  !-----------------------------------------------------------------------
  SUBROUTINE sync_file(pioFileDesc, ierr, message)
    ! !DESCRIPTION:
    ! sync a file to disc (write out all)
    !
    implicit none
    ! ARGUMENTS:
    type(file_desc_t), intent(inout) :: pioFileDesc  ! netcdf file id
    integer(i4b),      intent(out)   :: ierr         ! error status
    character(*),      intent(out)   :: message      ! error message

    ierr=0; message='sync_file/'

    call PIO_syncfile(pioFileDesc)

  END SUBROUTINE sync_file

  !-----------------------------------------------------------------------
  SUBROUTINE end_def(pioFileDesc, ierr, message)
    ! !DESCRIPTION:
    ! end definition of netcdf file
    !
    implicit none
    ! ARGUMENTS:
    type(file_desc_t), intent(inout) :: pioFileDesc  ! netcdf file id
    integer(i4b),      intent(out)   :: ierr         ! error status
    character(*),      intent(out)   :: message      ! error message

    ierr=0; message='end_def/'

    ierr = pio_enddef(pioFileDesc)
    if(ierr/=pio_noerr)then; message=trim(message)//'Could not end define mode'; return; endif

  END SUBROUTINE end_def

  !-----------------------------------------------------------------------
  SUBROUTINE def_dim(pioFileDesc,  & ! input: file descriptor
                    dimname,       & ! input: dimension name
                    dimlen,        & ! input: dimension size negative => record dimension
                    dimid)           ! output: dimension id
    implicit none
    ! ARGUMENTS:
    type(file_desc_t),      intent(in) :: pioFileDesc  ! netcdf file id
    character(len=*),       intent(in) :: dimname      ! netcdf attrib
    integer(i4b),           intent(in) :: dimlen       ! netcdf attrib value
    integer(i4b),           intent(out):: dimid        ! netcdf dimension id
    integer(i4b)                       :: ierr

    if (dimlen>0) then
      ierr = pio_def_dim(pioFileDesc, trim(dimname), dimlen, dimid)
    else
      ierr = pio_def_dim(pioFileDesc, trim(dimname), pio_unlimited, dimid)
    endif

  END SUBROUTINE def_dim

  !-----------------------------------------------------------------------
  SUBROUTINE inq_dim_len(pioFileDesc,  & ! input: file descriptor
                         dimname,      & ! input: dimension name
                         dimlen)         ! output: dimension id
    implicit none
    ! ARGUMENTS:
    type(file_desc_t),      intent(in)  :: pioFileDesc  ! netcdf file id
    character(len=*),       intent(in)  :: dimname      ! name of a dimension
    integer(i4b),           intent(out) :: dimlen       ! length of a dimension
    ! local variables
    integer(i4b)                        :: dimid        ! dimension IDs
    integer(i4b)                        :: ierr

    ierr = pio_inq_dimid(pioFileDesc, trim(dimname), dimid)
    ierr = pio_inq_dimlen(pioFileDesc, dimid, dimlen)

  END SUBROUTINE inq_dim_len

  !-----------------------------------------------------------------------
  SUBROUTINE def_var(pioFileDesc,   & ! input: file descriptor
                     vname,         & ! input: variable name
                     ivtype,        & ! input: variable type. pio_int, pio_real, pio_double, pio_char
                     ierr, message, & ! output: error code and message
                     pioDimId,      & ! input: optional. dimension ID(s)
                     vdesc,         & ! input: optional. long_name
                     vunit,         & ! input: optional. unit
                     vcal)            ! input: optional. calendar type
  implicit none
  ! ARGUMENTS:
  type(file_desc_t),intent(inout)         :: pioFileDesc            ! contains data identifying the file.
  character(*),     intent(in)            :: vname                  ! Input: variable name
  integer(i4b),     intent(in)            :: ivtype                 ! Input: variable type. pio_int, pio_real, pio_double, pio_char
  integer(i4b),     intent(in), optional  :: pioDimId(:)            ! Input: variable dimension IDs
  character(*),     intent(in), optional  :: vdesc                  ! Input: variable description
  character(*),     intent(in), optional  :: vunit                  ! Input: variable units
  character(*),     intent(in), optional  :: vcal                   ! Input: calendar (if time variable)
  integer(i4b),     intent(out)           :: ierr
  character(*),     intent(out)           :: message      ! error message
  ! local variables
  integer(i4b)                            :: dimid0(1)           ! dimid for no dimension
  type(var_desc_t)                        :: pioVarId

  dimid0 = 0

  if (present(pioDimId)) then
    ierr = pio_def_var(pioFileDesc, trim(vname), ivtype, pioDimId, pioVarId)
  else
    ierr = pio_def_var(pioFileDesc, trim(vname), ivtype, pioVarId)
  end if
  if(ierr/=pio_noerr)then; message=trim(message)//'Could not define variable'; return; endif

  if (present(vdesc)) then ! add long_name
    ierr = pio_put_att(pioFileDesc, pioVarId, 'long_name', trim(vdesc))
    if(ierr/=0)then; message=trim(message)//'ERROR: adding long_name'; return; endif
  end if

  if (present(vunit)) then ! add variable unit
    ierr = pio_put_att(pioFileDesc, pioVarId, 'units', trim(vunit))
    if(ierr/=0)then; message=trim(message)//'ERROR: adding unit'; return; endif
  end if

  if (present(vcal)) then ! add time calendar
    ierr = pio_put_att(pioFileDesc, pioVarId, 'calendar', trim(vcal))
    if(ierr/=0)then; message=trim(message)//'ERROR: adding calendar'; return; endif
  end if

  END SUBROUTINE def_var

  ! -----------------------------
  ! Writing routine
  ! -----------------------------

  ! ---------------------------------------------------------------
  ! write global integer scalar
  SUBROUTINE write_int_scalar(pioFileDesc,     &
                              vname,           &  ! input: variable name
                              scalar,          &  ! input: variable data
                              ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),  intent(inout) :: pioFileDesc  ! pio file handle
  character(*),       intent(in)    :: vname        ! variable name
  integer(i4b),       intent(in)    :: scalar       ! variable data
  integer(i4b),       intent(out)   :: ierr
  character(*),       intent(out)   :: message      ! error message
  ! local variables
  type(var_desc_t)                  :: pioVarId

  ierr=0; message='write_int_scalar/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, [scalar])
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_scalar

  ! ---------------------------------------------------------------
  ! write global real scalar
  SUBROUTINE write_double_scalar(pioFileDesc,     &
                               vname,           &  ! input: variable name
                               scalar,          &  ! input: variable data
                               ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),  intent(inout) :: pioFileDesc  ! pio file handle
  character(*),       intent(in)    :: vname        ! variable name
  real(dp),           intent(in)    :: scalar       ! variable data
  integer(i4b),       intent(out)   :: ierr
  character(*),       intent(out)   :: message      ! error message
  ! local variables
  type(var_desc_t)                  :: pioVarId

  ierr=0; message='write_double_scalar/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, [scalar])
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_scalar

  ! ---------------------------------------------------------------
  ! write global character scalar
  SUBROUTINE write_char_scalar(pioFileDesc,    &
                              vname,           &  ! input: variable name
                              scalar,          &  ! input: variable data
                              ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),  intent(inout) :: pioFileDesc  ! pio file handle
  character(*),       intent(in)    :: vname        ! variable name
  character(*),       intent(in)    :: scalar       ! variable data
  integer(i4b),       intent(out)   :: ierr
  character(*),       intent(out)   :: message      ! error message
  ! local variables
  type(var_desc_t)                  :: pioVarId

  ierr=0; message='write_char_scalar/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, [scalar])
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_char_scalar

  ! ---------------------------------------------------------------
  ! write global integer vector into 1D variable
  SUBROUTINE write_int_array1D(pioFileDesc,     &
                               vname,           &  ! input: variable name
                               array,           &  ! input: variable data
                               iStart,          &  ! input: start index
                               iCount,          &  ! input: length of vector
                               ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout) :: pioFileDesc  ! pio file handle
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:)     ! variable data
  integer(i4b),          intent(in)    :: iStart(:)    ! start index
  integer(i4b),          intent(in)    :: iCount(:)    ! length of vector
  integer(i4b),     intent(out)        :: ierr
  character(*),     intent(out)        :: message      ! error message
  ! local variables
  type(var_desc_t)                     :: pioVarId

  ierr=0; message='write_int_array1D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, iStart, iCount, array)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_array1D

  ! ---------------------------------------------------------------
  ! write global real vector into 1D variable
  SUBROUTINE write_double_array1D(pioFileDesc,    &
                                 vname,           &  ! input: variable name
                                 array,           &  ! input: variable data
                                 iStart,          &  ! input: start index
                                 iCount,          &  ! input: length of vector
                                 ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout) :: pioFileDesc  ! pio file handle
  character(*),          intent(in)    :: vname        ! variable name
  real(dp),              intent(in)    :: array(:)     ! variable data
  integer(i4b),          intent(in)    :: iStart(:)    ! start index
  integer(i4b),          intent(in)    :: iCount(:)    ! length of vector
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message      ! error message
  ! local variables
  type(var_desc_t)                     :: pioVarId

  ierr=0; message='write_double_array1D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, iStart, iCount, array)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_array1D

  ! ---------------------------------------------------------------
  ! write global integer 2D into 2D variable
  SUBROUTINE write_int_array2D(pioFileDesc,     &
                               vname,           &  ! input: variable name
                               array,           &  ! input: variable data
                               iStart,          &  ! input: start index
                               iCount,          &  ! input: length of vector
                               ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout) :: pioFileDesc  ! pio file handle
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:,:)   ! variable data
  integer(i4b),          intent(in)    :: iStart(:)    ! start index
  integer(i4b),          intent(in)    :: iCount(:)    ! length of vector
  integer(i4b),     intent(out)        :: ierr
  character(*),     intent(out)        :: message      ! error message
  ! local variables
  type(var_desc_t)                     :: pioVarId

  ierr=0; message='write_int_array2D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, iStart, iCount, array)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_array2D

  ! ---------------------------------------------------------------
  ! write global real 2D into 2D variable
  SUBROUTINE write_double_array2D(pioFileDesc,     &
                                  vname,           &  ! input: variable name
                                  array,           &  ! input: variable data
                                  iStart,          &  ! input: start index
                                  iCount,          &  ! input: length of vector
                                  ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout) :: pioFileDesc  ! pio file handle
  character(*),          intent(in)    :: vname        ! variable name
  real(dp),              intent(in)    :: array(:,:)   ! variable data
  integer(i4b),          intent(in)    :: iStart(:)    ! start index
  integer(i4b),          intent(in)    :: iCount(:)    ! length of vector
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message      ! error message
  ! local variables
  type(var_desc_t)                     :: pioVarId

  ierr=0; message='write_double_array2D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, iStart, iCount, array)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_array2D

  ! ---------------------------------------------------------------
  ! write distributed integer vector into 1D variable
  SUBROUTINE write_int_darray1D(pioFileDesc,     &
                                vname,           &  ! input: variable name
                                array,           &  ! input: variable data
                                iodesc,          &  ! input: ??? it is from initdecomp routine
                                ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout) :: pioFileDesc  ! pio file handle
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:)     ! variable data
  type(io_desc_t),       intent(inout) :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b), intent(out)            :: ierr         ! error code
  character(*), intent(out)            :: message      ! error message
  ! local variables
  type(var_desc_t)                     :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_int_darray1D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_darray1D

  ! ---------------------------------------------------------------
  ! write distributed 2D integer array into 2D variable
  SUBROUTINE write_int_darray2D(pioFileDesc,     &
                                vname,           &  ! input: variable name
                                array,           &  ! input: variable data
                                iodesc,          &  ! input: ??? it is from initdecomp routine
                                ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout) :: pioFileDesc  ! pio file handle
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout) :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(out)   :: ierr         ! error code
  character(*),          intent(out)   :: message      ! error message
  ! local variables
  type(var_desc_t)                     :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_int_darray2D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_darray2D

  ! ---------------------------------------------------------------
  ! write distributed 3D integer array into 2D variable
  SUBROUTINE write_int_darray3D(pioFileDesc,     &
                                vname,           &  ! input: variable name
                                array,           &  ! input: variable data
                                iodesc,          &  ! input: ??? it is from initdecomp routine
                                ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout) :: pioFileDesc  ! pio file handle
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:,:,:) ! variable data
  type(io_desc_t),       intent(inout) :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(out)   :: ierr         ! error code
  character(*),          intent(out)   :: message      ! error message
  ! local variables
  type(var_desc_t)                     :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_int_darray3D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_darray3D

  ! ---------------------------------------------------------------
  ! write distributed double vector into 1D data
  SUBROUTINE write_double_darray1D(pioFileDesc,   &
                                 vname,           &  ! input: variable name
                                 array,           &  ! input: variable data
                                 iodesc,          &  ! input: ??? it is from initdecomp routine
                                 ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout) :: pioFileDesc  ! pio file handle
  character(*),          intent(in)    :: vname        ! variable name
  real(dp),              intent(in)    :: array(:)     ! variable data
  type(io_desc_t),       intent(inout) :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b), intent(out)            :: ierr         ! error code
  character(*), intent(out)            :: message      ! error message
  ! local variables
  type(var_desc_t)                     :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_double_darray1D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_darray1D

  ! ---------------------------------------------------------------
  ! write distributed double 2D integer array into 2D variable
  SUBROUTINE write_double_darray2D(pioFileDesc,   &
                                 vname,           &  ! input: variable name
                                 array,           &  ! input: variable data
                                 iodesc,          &  ! input: ??? it is from initdecomp routine
                                 ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  real(dp),              intent(in)     :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_double_darray2D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_darray2D

  ! ---------------------------------------------------------------
  ! write distributed double 3D integer array into 2D variable
  SUBROUTINE write_double_darray3D(pioFileDesc,   &
                                 vname,           &  ! input: variable name
                                 array,           &  ! input: variable data
                                 iodesc,          &  ! input: ??? it is from initdecomp routine
                                 ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  real(dp),              intent(in)     :: array(:,:,:) ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_double_darray3D/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_darray3D

  !
  ! Writing distributed data in nc variable with record dimension
  !
  ! ---------------------------------------------------------------
  ! write distributed integer vector into 2D variable with record dimension
  SUBROUTINE write_int_darray2D_recdim(pioFileDesc,     &
                                       vname,           &  ! input: variable name
                                       array,           &  ! input: variable data
                                       iodesc,          &  ! input: ??? it is from initdecomp routine
                                       nr,              &  ! input: index of record dimension
                                       ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  integer(i4b),          intent(in)     :: array(:)     ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_int_darray2D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr,kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_darray2D_recdim

  ! ---------------------------------------------------------------
  ! write distributed float vector into 2D variable with record dimension
  ! ---------------------------------------------------------------
  SUBROUTINE write_float_darray2D_recdim(pioFileDesc,     &
                                         vname,           &  ! input: variable name
                                         array,           &  ! input: variable data
                                         iodesc,          &  ! input: it is from initdecomp routine
                                         nr,              &  ! input: index of record dimension
                                         ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  real(sp),              intent(in)     :: array(:)     ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_float_darray2D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr, kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_float_darray2D_recdim

  ! ---------------------------------------------------------------
  ! write distributed real vector into 2D variable with record dimension
  ! ---------------------------------------------------------------
  SUBROUTINE write_double_darray2D_recdim(pioFileDesc,     &
                                          vname,           &  ! input: variable name
                                          array,           &  ! input: variable data
                                          iodesc,          &  ! input: it is from initdecomp routine
                                          nr,              &  ! input: index of record dimension
                                          ierr, message)      ! output: error control
  implicit none
  ! Argument variables:
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  real(dp),              intent(in)     :: array(:)     ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_double_darray2D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr, kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_darray2D_recdim

  ! ---------------------------------------------------------------
  ! write distributed integer 2D array into 3D variable with record dimension
  SUBROUTINE write_int_darray3D_recdim(pioFileDesc,     &
                                       vname,           &  ! input: variable name
                                       array,           &  ! input: variable data
                                       iodesc,          &  ! input: ??? it is from initdecomp routine
                                       nr,              &  ! input: index of record dimension
                                       ierr, message)      ! output: error control
  implicit none
  ! Argument variables
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  integer(i4b),          intent(in)     :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_int_darray3D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr,kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_darray3D_recdim

  ! ---------------------------------------------------------------
  ! write distributed float 2D array into 3D variable with record dimension
  ! ---------------------------------------------------------------
  SUBROUTINE write_float_darray3D_recdim(pioFileDesc,     &
                                         vname,           &  ! input: variable name
                                         array,           &  ! input: variable data
                                         iodesc,          &  ! input: ??? it is from initdecomp routine
                                         nr,              &  ! input: index of record dimension
                                         ierr, message)      ! output: error control
  implicit none
  ! Argument variables
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  real(sp),              intent(in)     :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_float_darray3D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr, kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_float_darray3D_recdim

  ! ---------------------------------------------------------------
  ! write distributed double 2D array into 3D variable with record dimension
  ! ---------------------------------------------------------------
  SUBROUTINE write_double_darray3D_recdim(pioFileDesc,     &
                                          vname,           &  ! input: variable name
                                          array,           &  ! input: variable data
                                          iodesc,          &  ! input: ??? it is from initdecomp routine
                                          nr,              &  ! input: index of record dimension
                                          ierr, message)      ! output: error control
  implicit none
  ! Argument variables
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  real(dp),              intent(in)     :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_double_darray3D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr, kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_darray3D_recdim

  ! ---------------------------------------------------------------
  ! write distributed integer 3D array into 4D variable with record dimension
  SUBROUTINE write_int_darray4D_recdim(pioFileDesc,     &
                                       vname,           &  ! input: variable name
                                       array,           &  ! input: variable data
                                       iodesc,          &  ! input: ??? it is from initdecomp routine
                                       nr,              &  ! input: index of record dimension
                                       ierr, message)      ! output: error control
  implicit none
  ! Argument variables
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  integer(i4b),          intent(in)     :: array(:,:,:) ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_int_darray4D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr,kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_int_darray4D_recdim

  ! ---------------------------------------------------------------
  ! write distributed float 3D array into 4D variable with record dimension
  ! ---------------------------------------------------------------
  SUBROUTINE write_float_darray4D_recdim(pioFileDesc,     &
                                         vname,           &  ! input: variable name
                                         array,           &  ! input: variable data
                                         iodesc,          &  ! input: ??? it is from initdecomp routine
                                         nr,              &  ! input: index of record dimension
                                         ierr, message)      ! output: error control
  implicit none
  ! Argument variables
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  real(sp),              intent(in)     :: array(:,:,:) ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_float_darray4D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr, kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_float_darray4D_recdim

  ! ---------------------------------------------------------------
  ! write distributed double 3D array into 4D variable with record dimension
  ! ---------------------------------------------------------------
  SUBROUTINE write_double_darray4D_recdim(pioFileDesc,     &
                                          vname,           &  ! input: variable name
                                          array,           &  ! input: variable data
                                          iodesc,          &  ! input: ??? it is from initdecomp routine
                                          nr,              &  ! input: index of record dimension
                                          ierr, message)      ! output: error control
  implicit none
  ! Argument variables
  type(file_desc_t),     intent(inout)  :: pioFileDesc  ! pio file handle
  character(*),          intent(in)     :: vname        ! variable name
  real(dp),              intent(in)     :: array(:,:,:) ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_double_darray4D_recdim/'

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr, kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  END SUBROUTINE write_double_darray4D_recdim

END MODULE pio_utils
