module pio_utils

  USE mpi
  USE nrtype
  USE pio

  implicit none

  private

  public::pio_sys_init
  public::pio_decomp
  public::createFile
  public::defDim
  public::defVar
  public::endDef
  public::closeFile
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

  interface write_netcdf
    module procedure write_int_array1D
    module procedure write_real_array1D
    module procedure write_int_array2D
    module procedure write_real_array2D
  end interface

  interface write_pnetcdf
    module procedure write_int_darray1D
    module procedure write_int_darray2D
    module procedure write_real_darray1D
    module procedure write_real_darray2D
  end interface

  interface write_pnetcdf_recdim
    module procedure write_int_darray2D_recdim
    module procedure write_int_darray3D_recdim
    module procedure write_real_darray2D_recdim
    module procedure write_real_darray3D_recdim
  end interface

  !! pio_iotype_pnetcdf = 1     Parallel Netcdf  (parallel)
  !! pio_iotype_netcdf = 2      Netcdf3 Classic format (serial)
  !! pio_iotype_netcdf4c = 3    NetCDF4 (HDF5) compressed format (serial)
  !! pio_iotype_NETCDF4p = 4    NetCDF4 (HDF5) parallel
  integer, parameter :: iotype = pio_iotype_netcdf

contains

  ! *********************************************************************
  ! subroutine: initialize ParallelIO system
  ! *********************************************************************
  subroutine pio_sys_init(pid,       &  !input
                          nNodes,    &  !input
                          pioIOsystem)  !output
    implicit none
    ! input variables
    integer(i4b),          intent(in)  :: pid
    integer(i4b),          intent(in)  :: nNodes
    type(iosystem_desc_t), intent(out) :: pioIoSystem  ! pio system descriptor
    ! Local variables
    integer(i4b)                       :: stride       ! stride in MPI rank between IO tasks.
    integer(i4b)                       :: nIOtasks     ! number of IO tasks
    integer(i4b)                       :: nAggregator  ! MPI aggregator count
    integer(i4b),parameter             :: idxBase=1    ! Start index of IO tasks.

    ! set up PIO for rest of parameters
    stride      = 1
    nIOtasks    = nNodes/stride
    nAggregator = 0
    call pio_init(pid,              & ! input: MPI rank
                  MPI_COMM_WORLD,   & ! input: MPI communicator
                  nIOtasks,         & ! input: Number of iotasks
                  nAggregator,      & ! input: number of aggregators to use
                  stride,           & ! input: MPI rank stride between IO tasks
                  pio_rearr_subset, & ! input: do not use any form of rearrangement
                  pioIoSystem,      & ! output: pio system descriptor
                  base=idxBase)       ! base (optional argument)

  end subroutine pio_sys_init

  ! *********************************************************************
  ! subroutine: Initialize decompoistion
  ! *********************************************************************
  subroutine pio_decomp(pioIoSystem,  & ! input: pio system descriptor
                        piotype,      & ! input: data type (pio_int, pio_real, pio_double, pio_char)
                        dimLen,       & ! input: dimension length for global array
                        gidx_local,   & ! input: local global-index array at one dimension
                        iodesc)         ! output:
    ! Details:
    !
    !

    implicit none
    ! input variables
    type(iosystem_desc_t),intent(inout) :: pioIoSystem    ! pio system descriptor
    integer,              intent(in)    :: piotype        ! pio data type
    integer(i4b),         intent(in)    :: dimLen(:)      ! length of dimension for global array
    integer(i4b),         intent(in)    :: gidx_local(:)  ! indices of global array describing the decomposition of the data
    ! output variables
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
    endif

    allocate(compdof2d(lsize,nn))
    do ix = 1,nn
      compdof2d(1:lsize, ix) = (ix-1)*gsize + gidx_local(1:lsize)
    end do

    allocate(compdof(totnum))
    compdof = reshape(compdof2d,[totnum])

    call pio_initdecomp(pioIoSystem,      & ! input: pio system descriptor
                        piotype,          & ! input: data type
                        dimLen(1:ndims),  & ! input: dimension length in global array
                        compdof,          & ! input: decomposition for local array
                        iodesc)             ! output:

  end subroutine pio_decomp


  ! *********************************************************************
  ! subroutine: create an netCDF file
  ! *********************************************************************
  subroutine createFile(pioIoSystem, fileName, pioFileDesc, ierr, message)

    implicit none
    ! input variables
    type(iosystem_desc_t), intent(inout):: pioIoSystem
    character(*),          intent(in)   :: fileName
    type(file_desc_t),     intent(out)  :: pioFileDesc            ! contains data identifying the file.
    integer(i4b),          intent(out)  :: ierr
    character(*),          intent(out)  :: message      ! error message

    ierr = pio_createfile(pioIOsystem,    & ! input:
                          pioFileDesc,    & ! output:
                          iotype,         & ! input: ??
                          trim(fileName), & ! input: input file name
                          pio_clobber)      ! append
    if(ierr/=pio_noerr)then; message=trim(message)//'cannot create netCDF'; return; endif

  end subroutine createFile

  !-----------------------------------------------------------------------
  subroutine closeFile(file)
    ! !DESCRIPTION:
    ! Close a NetCDF PIO file
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: file   ! PIO file handle to close
    !-----------------------------------------------------------------------

    call pio_closefile(file)

  end subroutine closeFile

  !-----------------------------------------------------------------------
  subroutine endDef(pioFileDesc, ierr, message)
    ! !DESCRIPTION:
    ! end definition of netcdf file
    !
    ! input
    type(file_desc_t), intent(inout) :: pioFileDesc  ! netcdf file id
    ! output
    integer(i4b),      intent(out)   :: ierr         ! error status
    character(*),      intent(out)   :: message      ! error message
    ! Local variable
    integer(i4b)                     :: status

    ierr=0; message='endDef/'

    status = pio_enddef(pioFileDesc)
    if(status/=pio_noerr)then; ierr=10;message=trim(message)//'Could not end define mode'; return; endif

  end subroutine endDef

  ! *********************************************************************
  ! subroutine: define dimension
  ! *********************************************************************
  subroutine defDim(pioFileDesc,   & ! input: file descriptor
                    dimname,       & ! input: dimension name
                    dimlen,        & ! input: dimension size negative => record dimension
                    dimid)           ! output: dimension id
    implicit none
    type(file_desc_t),      intent(in) :: pioFileDesc  ! netcdf file id
    character(len=*),       intent(in) :: dimname      ! netcdf attrib
    integer(i4b),           intent(in) :: dimlen       ! netcdf attrib value
    integer(i4b),           intent(out):: dimid        ! netcdf dimension id
    integer(i4b)                       :: ierr
    !-----------------------------------------------------------------------
    if (dimlen>0) then
      ierr = pio_def_dim(pioFileDesc, dimname, dimlen, dimid)
    else
      ierr = pio_def_dim(pioFileDesc, dimname, pio_unlimited, dimid)
    endif

  end subroutine defDim

  ! *********************************************************************
  ! subroutine: define variable
  ! *********************************************************************
  subroutine defVar(pioFileDesc,   & ! input: file descriptor
                    vname,         & ! input: variable name
                    pioDimId,      & ! input: dimension ID(s)
                    ivtype,        & ! input: variable type. pio_int, pio_real, pio_double, pio_char
                    ierr, message, & ! output: error code and message
                    vdesc,         & ! input: optional. long_name
                    vunit,         & ! input: optional. unit
                    vcal)            ! input: optional. calendar type
  implicit none
  ! input
  type(file_desc_t),intent(inout)         :: pioFileDesc            ! contains data identifying the file.
  character(*),     intent(in)            :: vname                  ! Input: variable name
  integer(i4b),     intent(in)            :: pioDimId(:)            ! Input: variable dimension IDs
  integer(i4b),     intent(in)            :: ivtype                 ! Input: variable type. pio_int, pio_real, pio_double, pio_char
  character(*),     intent(in), optional  :: vdesc                  ! Input: variable description
  character(*),     intent(in), optional  :: vunit                  ! Input: variable units
  character(*),     intent(in), optional  :: vcal                   ! Input: calendar (if time variable)
  ! output
  integer(i4b),     intent(out)           :: ierr
  character(*),     intent(out)           :: message      ! error message
  ! local variables
  type(var_desc_t)                        :: pioVarId

  ierr = pio_def_var(pioFileDesc, trim(vname), ivtype, pioDimId, pioVarId)
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

  end subroutine defVar

  ! -----------------------------
  ! Writing routine
  ! -----------------------------

  ! ---------------------------------------------------------------
  ! write global integer vector into 1D variable
  subroutine write_int_array1D(pioIoSystem,     &
                               fname,           &  ! input: filename
                               vname,           &  ! input: variable name
                               array,           &  ! input: variable data
                               iStart,          &  ! input: start index
                               iCount,          &  ! input: length of vector
                               ierr, message)      ! output: error control
  implicit none
  ! input
  type(iosystem_desc_t), intent(inout) :: pioIoSystem
  character(*),          intent(in)    :: fname        ! filename
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:)     ! variable data
  integer(i4b),          intent(in)    :: iStart(:)    ! start index
  integer(i4b),          intent(in)    :: iCount(:)    ! length of vector
  ! output
  integer(i4b),     intent(out)        :: ierr
  character(*),     intent(out)        :: message      ! error message
  ! local variables
  type(file_desc_t)                    :: pioFileDesc  ! pio file handle
  type(var_desc_t)                     :: pioVarId

  ierr=0; message='write_int_array1D/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, iStart, iCount, array)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_closefile(pioFileDesc)

  end subroutine write_int_array1D

  ! ---------------------------------------------------------------
  ! write global real vector into 1D variable
  subroutine write_real_array1D(pioIoSystem,     &
                                fname,           &  ! input: filename
                                vname,           &  ! input: variable name
                                array,           &  ! input: variable data
                                iStart,          &  ! input: start index
                                iCount,          &  ! input: length of vector
                                ierr, message)      ! output: error control
  implicit none
  ! input
  type(iosystem_desc_t), intent(inout) :: pioIoSystem
  character(*),          intent(in)    :: fname        ! filename
  character(*),          intent(in)    :: vname        ! variable name
  real(dp),              intent(in)    :: array(:)     ! variable data
  integer(i4b),          intent(in)    :: iStart(:)    ! start index
  integer(i4b),          intent(in)    :: iCount(:)    ! length of vector
  ! output
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message      ! error message
  ! local variables
  type(file_desc_t)                    :: pioFileDesc  ! pio file handle
  type(var_desc_t)                     :: pioVarId

  ierr=0; message='write_real_array1D/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, iStart, iCount, real(array,kind=sp))
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_closefile(pioFileDesc)

  end subroutine write_real_array1D

  ! ---------------------------------------------------------------
  ! write global integer 2D into 2D variable
  subroutine write_int_array2D(pioIoSystem,     &
                                fname,           &  ! input: filename
                                vname,           &  ! input: variable name
                                array,           &  ! input: variable data
                                iStart,          &  ! input: start index
                                iCount,          &  ! input: length of vector
                                ierr, message)      ! output: error control
  implicit none
  ! input
  type(iosystem_desc_t), intent(inout) :: pioIoSystem
  character(*),          intent(in)    :: fname        ! filename
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:,:)   ! variable data
  integer(i4b),          intent(in)    :: iStart(:)    ! start index
  integer(i4b),          intent(in)    :: iCount(:)    ! length of vector
  ! output
  integer(i4b),     intent(out)        :: ierr
  character(*),     intent(out)        :: message      ! error message
  ! local variables
  type(file_desc_t)                    :: pioFileDesc  ! pio file handle
  type(var_desc_t)                     :: pioVarId

  ierr=0; message='write_int_array2D/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, iStart, iCount, array)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_closefile(pioFileDesc)

  end subroutine write_int_array2D

  ! ---------------------------------------------------------------
  ! write global real 2D into 2D variable
  subroutine write_real_array2D(pioIoSystem,     &
                                fname,           &  ! input: filename
                                vname,           &  ! input: variable name
                                array,           &  ! input: variable data
                                iStart,          &  ! input: start index
                                iCount,          &  ! input: length of vector
                                ierr, message)      ! output: error control
  implicit none
  ! input
  type(iosystem_desc_t), intent(inout) :: pioIoSystem
  character(*),          intent(in)    :: fname        ! filename
  character(*),          intent(in)    :: vname        ! variable name
  real(dp),              intent(in)    :: array(:,:)   ! variable data
  integer(i4b),          intent(in)    :: iStart(:)    ! start index
  integer(i4b),          intent(in)    :: iCount(:)    ! length of vector
  ! output
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message      ! error message
  ! local variables
  type(file_desc_t)                    :: pioFileDesc  ! pio file handle
  type(var_desc_t)                     :: pioVarId

  ierr=0; message='write_real_array2D/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  ierr = pio_put_var(pioFileDesc, pioVarId, iStart, iCount, real(array,kind=sp))
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_closefile(pioFileDesc)

  end subroutine write_real_array2D


  ! ---------------------------------------------------------------
  ! write distributed integer vector into 1D variable
  subroutine write_int_darray1D(pioIoSystem,     &
                                fname,           &  ! input: filename
                                vname,           &  ! input: variable name
                                array,           &  ! input: variable data
                                iodesc,          &  ! input: ??? it is from initdecomp routine
                                ierr, message)      ! output: error control
  implicit none
  ! input variables
  type(iosystem_desc_t), intent(inout) :: pioIoSystem
  character(*),          intent(in)    :: fname        ! filename
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:)     ! variable data
  type(io_desc_t),       intent(inout) :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  ! output variables
  integer(i4b), intent(out)            :: ierr         ! error code
  character(*), intent(out)            :: message      ! error message
  ! local variables
  type(file_desc_t)                    :: pioFileDesc  ! pio file handle
  type(var_desc_t)                     :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_int_darray1D/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_syncfile(pioFileDesc)

  call pio_closefile(pioFileDesc)

  end subroutine write_int_darray1D

  ! ---------------------------------------------------------------
  ! write distributed 2D integer array into 2D variable
  subroutine write_int_darray2D(pioIoSystem,     &
                                fname,           &  ! input: filename
                                vname,           &  ! input: variable name
                                array,           &  ! input: variable data
                                iodesc,          &  ! input: ??? it is from initdecomp routine
                                ierr, message)      ! output: error control
  implicit none
  ! input variables
  type(iosystem_desc_t), intent(inout) :: pioIoSystem
  character(*),          intent(in)    :: fname        ! filename
  character(*),          intent(in)    :: vname        ! variable name
  integer(i4b),          intent(in)    :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout) :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  ! output variables
  integer(i4b),          intent(out)   :: ierr         ! error code
  character(*),          intent(out)   :: message      ! error message
  ! local variables
  type(file_desc_t)                    :: pioFileDesc  ! pio file handle
  type(var_desc_t)                     :: pioVarId     ! netCDF variable ID

  ! initialize error control
  ierr=0; message='write_int_darray2D/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_syncfile(pioFileDesc)

  call pio_closefile(pioFileDesc)

  end subroutine write_int_darray2D

  ! ---------------------------------------------------------------
  ! write distributed real vector into 1D data
  subroutine write_real_darray1D(pioIoSystem,     &  ! input:
                                 fname,           &  ! input: filename
                                 vname,           &  ! input: variable name
                                 array,           &  ! input: variable data
                                 iodesc,          &  ! input: ??? it is from initdecomp routine
                                 ierr, message)      ! output: error control
  implicit none
  ! input variables
  type(iosystem_desc_t), intent(inout) :: pioIoSystem
  character(*),          intent(in)    :: fname        ! filename
  character(*),          intent(in)    :: vname        ! variable name
  real(dp),              intent(in)    :: array(:)     ! variable data
  type(io_desc_t),       intent(inout) :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  ! output variables
  integer(i4b), intent(out)            :: ierr         ! error code
  character(*), intent(out)            :: message      ! error message
  ! local variables
  type(file_desc_t)                    :: pioFileDesc  ! pio file handle
  type(var_desc_t)                     :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_real_darray1D/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, real(array,kind=sp), ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_syncfile(pioFileDesc)

  call pio_closefile(pioFileDesc)

  end subroutine write_real_darray1D

  ! ---------------------------------------------------------------
  ! write distributed real 2D integer array into 2D variable
  subroutine write_real_darray2D(pioIoSystem,     &
                                 fname,           &  ! input: filename
                                 vname,           &  ! input: variable name
                                 array,           &  ! input: variable data
                                 iodesc,          &  ! input: ??? it is from initdecomp routine
                                 ierr, message)      ! output: error control
  implicit none
  ! input variables
  type(iosystem_desc_t), intent(inout)  :: pioIoSystem
  character(*),          intent(in)     :: fname        ! filename
  character(*),          intent(in)     :: vname        ! variable name
  real(dp),              intent(in)     :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  ! output variables
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(file_desc_t)                     :: pioFileDesc  ! pio file handle
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ierr=0; message='write_real_darray2D/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, real(array,kind=sp), ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_syncfile(pioFileDesc)

  call pio_closefile(pioFileDesc)

  end subroutine write_real_darray2D

  !
  ! Writing distributed data in nc variable with record dimension
  !
  ! ---------------------------------------------------------------
  ! write distributed integer vector into 2D variable with record dimension
  subroutine write_int_darray2D_recdim(pioIoSystem,     &
                                       fname,           &  ! input: filename
                                       vname,           &  ! input: variable name
                                       array,           &  ! input: variable data
                                       iodesc,          &  ! input: ??? it is from initdecomp routine
                                       nr,              &  ! input: index of record dimension
                                       ierr, message)      ! output: error control
  implicit none
  ! input variables
  type(iosystem_desc_t), intent(inout)  :: pioIoSystem
  character(*),          intent(in)     :: fname        ! filename
  character(*),          intent(in)     :: vname        ! variable name
  integer(i4b),          intent(in)     :: array(:)     ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  ! output variables
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(file_desc_t)                     :: pioFileDesc  ! pio file handle
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ! initialize error control
  ierr=0; message='write_int_darray2D_recdim/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr,kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_syncfile(pioFileDesc)

  call pio_closefile(pioFileDesc)

  end subroutine write_int_darray2D_recdim

  ! ---------------------------------------------------------------
  ! write distributed real vector into 2D variable with record dimension
  subroutine write_real_darray2D_recdim(pioIoSystem,     &
                                        fname,           &  ! input: filename
                                        vname,           &  ! input: variable name
                                        array,           &  ! input: variable data
                                        iodesc,          &  ! input: it is from initdecomp routine
                                        nr,              &  ! input: index of record dimension
                                        ierr, message)      ! output: error control
  implicit none
  ! input variables
  type(iosystem_desc_t), intent(inout)  :: pioIoSystem
  character(*),          intent(in)     :: fname        ! filename
  character(*),          intent(in)     :: vname        ! variable name
  real(dp),              intent(in)     :: array(:)     ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  ! output variables
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(file_desc_t)                     :: pioFileDesc  ! pio file handle
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ! initialize error control
  ierr=0; message='write_real_darray2D_recdim/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr, kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, real(array,kind=sp), ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_syncfile(pioFileDesc)

  call pio_closefile(pioFileDesc)

  end subroutine write_real_darray2D_recdim

  ! ---------------------------------------------------------------
  ! write distributed integer 2D array into 3D variable with record dimension
  subroutine write_int_darray3D_recdim(pioIoSystem,     &
                                       fname,           &  ! input: filename
                                       vname,           &  ! input: variable name
                                       array,           &  ! input: variable data
                                       iodesc,          &  ! input: ??? it is from initdecomp routine
                                       nr,              &  ! input: index of record dimension
                                       ierr, message)      ! output: error control
  implicit none
  ! input variables
  type(iosystem_desc_t), intent(inout)  :: pioIoSystem
  character(*),          intent(in)     :: fname        ! filename
  character(*),          intent(in)     :: vname        ! variable name
  integer(i4b),          intent(in)     :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  ! output variables
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(file_desc_t)                     :: pioFileDesc  ! pio file handle
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ! initialize error control
  ierr=0; message='write_int_darray3D_recdim/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr,kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, array, ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_syncfile(pioFileDesc)

  call pio_closefile(pioFileDesc)

  end subroutine write_int_darray3D_recdim

  ! ---------------------------------------------------------------
  ! write distributed real 2D array into 3D variable with record dimension
  subroutine write_real_darray3D_recdim(pioIoSystem,     &
                                        fname,           &  ! input: filename
                                        vname,           &  ! input: variable name
                                        array,           &  ! input: variable data
                                        iodesc,          &  ! input: ??? it is from initdecomp routine
                                        nr,              &  ! input: index of record dimension
                                        ierr, message)      ! output: error control
  implicit none
  ! input variables
  type(iosystem_desc_t), intent(inout)  :: pioIoSystem
  character(*),          intent(in)     :: fname        ! filename
  character(*),          intent(in)     :: vname        ! variable name
  real(dp),              intent(in)     :: array(:,:)   ! variable data
  type(io_desc_t),       intent(inout)  :: iodesc       ! io descriptor handle that is generated in PIO_initdecomp
  integer(i4b),          intent(in)     :: nr           ! index of record dimension
  ! output variables
  integer(i4b),          intent(out)    :: ierr         ! error code
  character(*),          intent(out)    :: message      ! error message
  ! local variables
  type(file_desc_t)                     :: pioFileDesc  ! pio file handle
  type(var_desc_t)                      :: pioVarId     ! netCDF variable ID

  ! initialize error control
  ierr=0; message='write_real_darray3D_recdim/'

  ierr=pio_openfile(pioIoSystem, pioFileDesc, iotype, trim(fname), mode=pio_write)
  if(ierr/=0)then; message=trim(message)//'ERROR: opening nc file'; return; endif

  ierr = pio_inq_varid(pioFileDesc, trim(vname), pioVarId)
  if(ierr/=0)then; message=trim(message)//'ERROR: getting variable id'; return; endif

  call pio_setframe(pioFileDesc, pioVarId, int(nr, kind=pio_offset_kind))

  call pio_write_darray(pioFileDesc, pioVarId, iodesc, real(array,kind=sp), ierr)
  if(ierr/=pio_noerr)then; message=trim(message)//'cannot write data'; return; endif

  call pio_syncfile(pioFileDesc)

  call pio_closefile(pioFileDesc)

  end subroutine write_real_darray3D_recdim


end module pio_utils
