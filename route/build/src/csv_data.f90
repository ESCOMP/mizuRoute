MODULE csv_data

! csv class
!
! usage example:
!
!  ! --initialize csv (must be done)
!  csv_data = csv(trim(csv_file), ierr, cmessage,  mode='r')
!
!  ! --Read csv file into the csv_data (here header is first row)
!  call csv_data%csv_read(ierr, cmessage, header_row=1)
!
!  ! --safely closed now
!  call csv_data%fclose()
!
!  ! -- Explore csv_data --
!
!  ! --Get header name
!  call csv_data%get_header(header, ierr, cmessage)
!
!  ! --Get data type for each column
!  call csv_data%get_dtype(dtype)
!
!  ! --Get the second column as integer
!  call csv_data%get_data(col=2, col_data_integer, ierr)

USE nrtype
USE public_var
USE ascii_utils

implicit none
integer(i4b), parameter :: type_string  = 1  ! a character string cell
integer(i4b), parameter :: type_double  = 2  ! a real64 cell
integer(i4b), parameter :: type_float   = 3  ! a real32 cell
integer(i4b), parameter :: type_integer = 4  ! an integer cell
integer(i4b), parameter :: type_logical = 5  ! a logical cell

private

type, public :: string_array
  character(len=:), allocatable :: str
end type string_array

type, public :: csv
  private
  character(len=strLen)           :: fname            ! Filename of the opened dataset
  character(1)                    :: mode ='r'        ! File open mode
  character(1)                    :: delimiter =','   !
  character(1)                    :: quote     = '"'  ! quotation character
  integer(i4b)                    :: fid       = 1    !
  integer(i4b)                    :: nrows = 0        ! number of rows in the file
  integer(i4b)                    :: ncols = 0        ! number of columns in the file
  type(string_array), allocatable :: header(:)        ! header array
  integer(i4b),       allocatable :: dtype(:)         ! datatype identification
  type(string_array), allocatable :: csv_data(:,:)    ! data array

CONTAINS

  procedure, public  :: csv_read => read_data
  procedure, public  :: get_dtype
  procedure, public  :: fclose
  procedure, private :: variable_types
  procedure, private :: split_csv_line
  procedure, private :: initCSV
  generic,   public  :: get_data => get_csv_data_as_str, &
                                    get_cell_value,      &
                                    get_real_sp_column,  &
                                    get_real_dp_column,  &
                                    get_integer_column,  &
                                    get_long_column,     &
                                    get_logical_column,  &
                                    get_character_column,  &
                                    get_string_column
  procedure, private :: get_csv_data_as_str
  procedure, private :: get_cell_value
  procedure, private :: get_column
  procedure, private :: get_real_sp_column
  procedure, private :: get_real_dp_column
  procedure, private :: get_integer_column
  procedure, private :: get_long_column
  procedure, private :: get_logical_column
  procedure, private :: get_character_column
  procedure, private :: get_string_column
  procedure, private :: get_icol_by_name
  generic,   public  :: get_header => get_header_str, &
                                     get_header_csv_str
  procedure, private :: get_header_str
  procedure, private :: get_header_csv_str
  final              :: clean

end type csv

private :: read_data

INTERFACE csv
  module procedure constructor
END INTERFACE csv

CONTAINS

  FUNCTION constructor(fname, ierr, message, mode, delimiter) RESULT(instCSV)
    implicit none
    type(csv)                            :: instCSV
    ! Argument variables:
    character(*),  intent(in)            :: fname
    integer(i4b),  intent(out)           :: ierr
    character(*),  intent(out)           :: message
    character(1),  intent(in), optional  :: mode
    character(1),  intent(in), optional  :: delimiter
    ! Local variables:
    character(len=strLen)                :: cmessage

    instCSV%fname     = fname
    if (present(mode)) then
      instCSV%mode = mode
    end if
    if (present(delimiter)) then
      instCSV%delimiter = delimiter
    end if

    call instCSV%initCSV(ierr, cmessage)
    if(ierr/=0)then; message="Initialization csv failed: "//trim(cmessage); return; endif

  END FUNCTION constructor

  SUBROUTINE initCSV(this, ierr, message)

    implicit none
    ! Argument variables:
    class(csv),           intent(inout) :: this
    integer(i4b),         intent(out)   :: ierr
    character(len=strLen),intent(out)   :: message

    ierr=0; message='initCSV/'

    select case(this%mode)
    case("w")
       open(this%fid, file=this%fname, status='replace', action='write', iostat=ierr)
    case("r")
       open(this%fid, file=this%fname, action='read', status='old', iostat=ierr)
    case("a")
    case default
       ierr=20; message=trim(message)//"Mode argument must be in 'w','r','a' !"
       return
    end select
    if(ierr/=0)then; message=trim(message)//"Error: open failed: "//trim(this%fname); return; endif

  END SUBROUTINE initCSV

  SUBROUTINE read_data(this, ierr, message, header_row, skip_rows)

    implicit none
    class(csv),       intent(inout)       :: this
    integer(i4b),     intent(in),optional :: header_row    ! the header row
    integer(i4b),     intent(in),optional :: skip_rows(:)  ! rows to skip
    integer(i4b),     intent(out)         :: ierr          !
    character(*),     intent(out)         :: message
    ! Local variables
    type(string_array), allocatable       :: row_data(:)      ! a split row
    character(len=strLen),allocatable     :: cLines(:)        ! vector of character strings
    integer(i4b)                          :: iLine            !
    integer(i4b)                          :: irow             ! counter
    integer(i4b)                          :: jcol             ! counter
    integer(i4b)                          :: nRows            ! number of rows in the output data matrix
    integer(i4b)                          :: nCols            ! number of columns in the file (and output data matrix)
    logical(lgt)                          :: arrays_allocated !
    integer(i4b)                          :: iheader          ! row number of header row
    character(len=strLen)                 :: cmessage         ! error message from subroutine

    ierr=0; message="read_data/"

    arrays_allocated = .false.

    !get character array (nrow)
    call get_vlines(this%fid, cLines, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage);return;endif

    !get number of lines in the data array
    if (present(skip_rows)) then
      nrows = size(cLines) - size(skip_rows)
    else
      nrows = size(cLines)
    end if
    if (present(header_row)) then
      iheader = max(0, header_row)
      nrows = nrows - 1
    else ! no header
      iheader = 0
    end if

    this%nrows = nrows

    ! we don't know the number of columns
    ! until we parse the first row (or the header)

    !read each line in the file, parse it, and populate data
    irow = 0
    do iLine=1,size(cLines)
      ! skip row if necessary
      if (present(skip_rows)) then
        if (any(iLine==skip_rows)) then
          cycle
        end if
      end if

      call this%split_csv_line(cLines(iLine), row_data)

      if (.not. arrays_allocated) then
        ! note: the number of columns is obtained
        ! from the first one read.
        ncols = size(row_data)
        this%ncols = ncols
        allocate(this%csv_data(nrows, ncols))
        if (iheader/=0) allocate(this%header(ncols))
        arrays_allocated = .true.
      end if

      if (iLine==iheader) then
        do jcol=1, this%ncols
          this%header(jcol)%str = row_data(jcol)%str
        end do
      else
        irow = irow + 1  !! row counter in data array
        do jcol=1, ncols
          this%csv_data(irow,jcol) = row_data(jcol)
        end do
      end if
    end do

    call this%variable_types()

  END SUBROUTINE read_data

  SUBROUTINE fclose(this)
    implicit none
    ! Argument variables
    class(csv),           intent(in)   :: this
    close(this%fid)
  END SUBROUTINE

  SUBROUTINE get_header_csv_str(this, header, ierr, message)
    ! Returns the header as a type(string_array) array.
    implicit none
    class(csv),                      intent(in)  :: this
    type(string_array), allocatable, intent(out) :: header(:)
    integer(i4b),                    intent(out) :: ierr          ! the actual rows to skip
    character(*),                    intent(out) :: message
    integer(i4b)                                 :: icol

    ierr=0; message="get_header_csv_str/"

    if (allocated(this%header)) then
      allocate(header(this%ncols))
      do icol=1, this%ncols
        header(icol) = this%header(icol)
      end do
    else
      ierr=10; message=trim(message)//'csv data is not populated. call "csv_read" first'
    end if

  END SUBROUTINE get_header_csv_str

  SUBROUTINE get_header_str(this, header, ierr, message)
    ! Returns the header as a character(len=*) array.
    implicit none
    class(csv),                    intent(in)    :: this
    character(len=*), allocatable, intent(out)   :: header(:)
    integer(i4b),                  intent(out)   :: ierr          ! the actual rows to skip
    character(*),                  intent(out)   :: message
    integer(i4b)                                 :: icol

    ierr=0; message="get_header_str/"

    if (allocated(this%header)) then
      allocate(header(this%ncols))
      do icol=1, this%ncols
        header(icol) = this%header(icol)%str
      end do
    else
      ierr=10; message=trim(message)//'csv data is not populated. call "csv_read" first'
    end if

  END SUBROUTINE get_header_str

  SUBROUTINE get_dtype(this, dtype, ierr, message)
    ! Returns datatype array as a integer array.
    implicit none
    class(csv),                intent(inout) :: this
    integer(i4b), allocatable, intent(out)   :: dtype(:)
    integer(i4b),              intent(out)   :: ierr          ! the actual rows to skip
    character(*),              intent(out)   :: message
    integer(i4b)                             :: icol

    if (allocated(this%dtype)) then
      allocate(dtype(this%ncols))
      do icol=1, this%ncols
        dtype(icol) = this%dtype(icol)
      end do
    else
      ierr=10; message=trim(message)//'csv data is not populated. call "csv_read" first'
    end if

  END SUBROUTINE get_dtype

  SUBROUTINE get_csv_data_as_str(this, csv_data, ierr, message)
    ! Returns a character(len=*) array containing the csv data
    implicit none
    class(csv),                    intent(inout) :: this
    character(len=*), allocatable, intent(out)   :: csv_data(:,:)  !! the data
    integer(i4b),                  intent(out)   :: ierr          ! the actual rows to skip
    character(*),                  intent(out)   :: message
    integer(i4b)                                 :: i              !! row counter
    integer(i4b)                                 :: j              !! column counter

    if (allocated(this%csv_data)) then
      ! size the output array:
      allocate(csv_data(this%nrows,this%ncols))
      ! convert each element to a string:
      do j=1, this%ncols
        do i=1, this%nrows
          csv_data(i,j) = this%csv_data(i,j)%str
        end do
      end do
    else
      ierr=10; message=trim(message)//'csv_data is not populated. call "csv_read" first'
    end if

  END SUBROUTINE get_csv_data_as_str

  SUBROUTINE variable_types(this)
    ! Returns an array indicating the variable type of each columns.
    implicit none
    class(csv), intent(inout) :: this
    integer(i4b)              :: icol

    if (allocated(this%csv_data)) then
      allocate(this%dtype(this%ncols))
      do icol=1,this%ncols
        call infer_variable_type(this%csv_data(1,icol)%str, this%dtype(icol))
      end do
    else
      write(*,'(A,1X,I5)') 'Error: csv_data is not populated. call "csv_read" first '
    end if

    CONTAINS

      SUBROUTINE infer_variable_type(str,itype)
        !  Infers the variable type, assuming the following precedence:
        !  * string  = 1
        !  * double  = 2
        !  * integer = 4
        !  * logical = 5
        implicit none
        character(len=*),intent(in) :: str
        integer,intent(out) :: itype
        real(dp)            :: rval      ! a real64 value
        integer(i4b)        :: ival      ! an iteger value
        logical(lgt)        :: lval      ! a logical value
        integer(i4b)        :: ierr      !

        call to_integer(str, ival, ierr)
        if (ierr==0) then
          itype = type_integer
          return
        end if

        call to_real_dp(str, rval, ierr)
        if (ierr==0) then
          itype = type_double
          return
        end if

        call to_logical(str, lval, ierr)
        if (ierr==0) then
          itype = type_logical
          return
        end if

        ! default is string:
        itype = type_string

      END SUBROUTINE infer_variable_type

  END SUBROUTINE variable_types

  pure elemental SUBROUTINE to_real_sp(str,val, istat)
    ! Convert a string to a `real(sp)`
    implicit none
    character(len=*),intent(in)  :: str
    real(sp),        intent(out) :: val
    integer(i4b),    intent(out) :: istat

    read(str,fmt=*,iostat=istat) val
    if (istat/=0) then
      val = 0._sp
    end if

  END SUBROUTINE to_real_sp

  pure elemental SUBROUTINE to_real_dp(str, val, istat)
    ! Convert a string to a `real(dp)`
    implicit none
    character(len=*),intent(in)  :: str
    real(dp),        intent(out) :: val
    integer(i4b),    intent(out) :: istat

    read(str,fmt=*,iostat=istat) val
    if (istat/=0) then
      val = 0._dp
    end if

  END SUBROUTINE to_real_dp

  pure elemental SUBROUTINE to_integer(str, val, istat)
    ! Convert a string to a `integer(ip)`
    implicit none
    character(len=*),intent(in)  :: str
    integer(i4b),    intent(out) :: val
    integer(i4b),    intent(out) :: istat

    read(str,fmt='(I256)',iostat=istat) val
    if (istat/=0) then
      val = 0
    end if

  END SUBROUTINE to_integer

  pure elemental SUBROUTINE to_logical(str, val, istat)
    ! The string match is not case sensitive.
    implicit none
    character(len=*),intent(in)    :: str
    logical(lgt),    intent(out)   :: val
    character(len=50), allocatable :: tmp
    integer(i4b),    intent(out) :: istat

    ! True and False options (all lowercase):
    character(len=*),dimension(3),parameter :: true_str  = ['t     ', &
                                                            'true  ', &
                                                            '.true.']
    character(len=*),dimension(3),parameter :: false_str = ['f      ', &
                                                            'false  ', &
                                                            '.false.']
    tmp = lower(str)
    if ( any(trim(tmp)==true_str) ) then
      val = .true.
      istat = 0
    else if ( any(trim(tmp)==false_str) ) then
      val = .false.
      istat = 0
    else
      val = .false.
      istat = 1
    end if

  END SUBROUTINE to_logical

  SUBROUTINE get_cell_value(this, row, col, val, istat)
    ! Get an individual value from the `csv_data` structure in the CSV class.
    ! The output `val` can be an `integer(ip)`, `real(wp)`,
    ! `logical`, or `character(len=*)` variable.

    implicit none
    class(csv),    intent(in)    :: this
    integer(i4b),  intent(in)    :: row   !! row number
    integer(i4b),  intent(in)    :: col   !! column number
    class(*),      intent(out)   :: val   !! the returned value
    integer(i4b),  intent(out)   :: istat !! status flag

    select type (val)
    type is (integer(i4b))
      call to_integer(this%csv_data(row,col)%str, val, istat)
    type is (real(sp))
      call to_real_sp(this%csv_data(row,col)%str, val, istat)
    type is (real(dp))
      call to_real_dp(this%csv_data(row,col)%str, val, istat)
    type is (logical(lgt))
      call to_logical(this%csv_data(row,col)%str, val, istat)
    type is (character(len=*))
      istat = 0
      if (allocated(this%csv_data(row,col)%str)) then
        val = this%csv_data(row,col)%str
      else
        val = ''
      end if
    type is (string_array)
      val = this%csv_data(row,col)
      istat = 0
    class default
      istat = 1
    end select

  END SUBROUTINE get_cell_value

  SUBROUTINE get_icol_by_name(this, col_name, icol, istat)
    implicit none
    class(csv),  intent(in)             :: this
    character(*),intent(in)             :: col_name     ! column name
    integer(i4b),intent(out)            :: icol
    integer(i4b), intent(out)           :: istat        ! status flag
    logical(lgt)                        :: found
    character(len=strLen), allocatable  :: header(:)
    character(len=strLen)               :: cmessage

    call this%get_header(header, istat, cmessage)

    found=.false.
    do icol=1,this%ncols
      if (trim(header(icol))==trim(col_name)) then
        found=.true.
        exit
      end if
    end do

    if (found) then
      istat = 0
    else
      write(*,'(A,1X,A)') 'Error: column name is not found: ',trim(col_name)
      istat = 1
    end if

  END SUBROUTINE get_icol_by_name

!@note This routine requires that the `r` array already be allocated.
!      This is because Fortran doesn't want to allow to you pass
!      a non-polymorphic variable into a routine with a dummy variable
!      with `class(*),dimension(:),allocatable,intent(out)` attributes.
  SUBROUTINE get_column(this, icol, col_data, istat)
    ! Return a column from a CSV file vector.
    implicit none
    class(csv),  intent(in)    :: this
    integer(i4b),intent(in)    :: icol        ! column number
    class(*),    intent(out)   :: col_data(:) ! assumed to have been allocated to
                                              ! the correct size by the caller.
                                              ! (`n_rows`)
    integer(i4b), intent(out) :: istat        ! status flag
    integer(i4b)              :: i !! counter
    character(len=:),allocatable :: tmp !! for gfortran workaround

    istat=0

    ! we know the data is allocated, since that
    ! was checked by the calling routines.

    if (this%ncols>=icol .and. icol>0) then
      do i= 1, this%nrows  ! row loop
        ! the following is a workaround for gfortran bugs:
        select type (col_data)
        type is (character(len=*))
          tmp = repeat(' ',len(col_data)) ! size the string
          call this%get_cell_value(i, icol, tmp, istat)
          col_data(i) = tmp
        class default
          call this%get_cell_value(i, icol, col_data(i), istat)
        end select

        if (istat/=0) then
          select type (col_data)
          ! note: character conversion can never fail, so not
          ! checking for that here. also we know it is real,
          ! integer, or logical at this point.
          type is (integer(i4b))
            write(*,'(A)') 'Error converting string to integer: '//trim(this%csv_data(i,icol)%str)
            col_data(i) = integerMissing
          type is (real(sp))
            write(*,'(A)') 'Error converting string to real(real32): '//trim(this%csv_data(i,icol)%str)
            col_data(i) = floatMissing
          type is (real(dp))
            write(*,'(A)') 'Error converting string to real(real64): '//trim(this%csv_data(i,icol)%str)
            col_data(i) = realMissing
          type is (logical(lgt))
            write(*,'(A)') 'Error converting string to logical: '//trim(this%csv_data(i,icol)%str)
            col_data(i) = .false.
          end select
        end if
      end do
    else
      write(*,'(A,1X,I5)') 'Error: invalid column number: ',icol
      istat = 1
    end if

  END SUBROUTINE get_column

  SUBROUTINE get_real_sp_column(this, col, col_data, istat)
    !  Return a column from a CSV file as a `real(sp)` vector.
    implicit none
    class(csv),            intent(in)    :: this
    class(*),              intent(in)    :: col          ! name:character or index:integer
    real(sp), allocatable, intent(out)   :: col_data(:)
    integer(i4b),          intent(out)   :: istat
    integer(i4b)                         :: icol         ! column number

    istat=0

    select type (col)
      type is (integer(i4b))
       icol=col
      type is (character(len=*))
        call this%get_icol_by_name(col, icol, istat)
    end select

    if (allocated(this%csv_data)) then
      allocate(col_data(this%nrows))  ! size the output vector
      call this%get_column(icol, col_data, istat)
    else
      write(*,'(A,1X,I5)') 'Error: class has not been initialized'
      istat = 1
    end if

  END SUBROUTINE get_real_sp_column

  SUBROUTINE get_real_dp_column(this, col, col_data, istat)
    !!!  Return a column from a CSV file as a `real(wp)` vector.
    implicit none
    class(csv),            intent(in)    :: this
    class(*),              intent(in)    :: col         ! name:character or index:integer
    real(dp), allocatable, intent(out)   :: col_data(:)
    integer(i4b),          intent(out)   :: istat
    integer(i4b)                         :: icol        ! column number

    istat=0

    select type (col)
      type is (integer(i4b))
       icol=col
      type is (character(len=*))
        call this%get_icol_by_name(col, icol, istat)
    end select

    if (allocated(this%csv_data)) then
      allocate(col_data(this%nrows))  ! size the output vector
      call this%get_column(icol, col_data, istat)
    else
      write(*,'(A,1X,I5)') 'Error: class has not been initialized'
      istat = 1
    end if

  END SUBROUTINE get_real_dp_column

  SUBROUTINE get_integer_column(this, col, col_data, istat)
    ! Return a column from a CSV file as a `integer(ip)` vector.
    implicit none
    class(csv),               intent(in)    :: this
    class(*),                 intent(in)    :: col         ! name:character or index:integer
    integer(i4b),allocatable, intent(out)   :: col_data(:)
    integer(i4b),             intent(out)   :: istat
    integer(i4b)                            :: icol        ! column number

    istat=0

    select type (col)
      type is (integer(i4b))
       icol=col
      type is (character(len=*))
        call this%get_icol_by_name(col, icol, istat)
    end select

    if (allocated(this%csv_data)) then
      allocate(col_data(this%nrows))  ! size the output vector
      call this%get_column(icol, col_data, istat)
    else
      write(*,'(A,1X,I5)') 'Error: class has not been initialized'
      istat = 1
    end if

  END SUBROUTINE get_integer_column

  SUBROUTINE get_long_column(this, col, col_data, istat)
    ! Return a column from a CSV file as a `integer(ip)` vector.
    implicit none
    class(csv),               intent(in)    :: this
    class(*),                 intent(in)    :: col         ! name:character or index:integer
    integer(i8b),allocatable, intent(out)   :: col_data(:)
    integer(i4b),             intent(out)   :: istat
    integer(i4b)                            :: icol        ! column number

    istat=0

    select type (col)
      type is (integer(i4b))
       icol=col
      type is (character(len=*))
        call this%get_icol_by_name(col, icol, istat)
    end select

    if (allocated(this%csv_data)) then
      allocate(col_data(this%nrows))  ! size the output vector
      call this%get_column(icol, col_data, istat)
    else
      write(*,'(A,1X,I5)') 'Error: class has not been initialized'
      istat = 1
    end if

  END SUBROUTINE get_long_column

  SUBROUTINE get_logical_column(this, col, col_data, istat)
    ! Convert a column from a `string_array` matrix to a `logical` vector.
    implicit none
    class(csv),                intent(in)    :: this
    class(*),                  intent(in)    :: col         ! name:character or index:integer
    logical(lgt), allocatable, intent(out)   :: col_data(:)
    integer(i4b),              intent(out)   :: istat
    integer(i4b)                             :: icol        ! column number

    istat=0

    select type (col)
      type is (integer(i4b))
       icol=col
      type is (character(len=*))
        call this%get_icol_by_name(col, icol, istat)
    end select

    if (allocated(this%csv_data)) then
      allocate(col_data(this%nrows))  ! size the output vector
      call this%get_column(icol, col_data, istat)
    else
      write(*,'(A,1X,I5)') 'Error: class has not been initialized'
      istat = 1
    end if

  END SUBROUTINE get_logical_column

  SUBROUTINE get_character_column(this, col, col_data, istat)
    ! Convert a column from a `string_array` matrix to a `character(len=*)` vector.
    implicit none
    class(csv),                    intent(in)    :: this
    class(*),                      intent(in)    :: col         ! name:character or index:integer
    character(len=*), allocatable, intent(out)   :: col_data(:)
    integer(i4b),                  intent(out)   :: istat
    integer(i4b)                                 :: icol        ! column number

    istat=0

    select type (col)
      type is (integer(i4b))
       icol=col
      type is (character(len=*))
        call this%get_icol_by_name(col, icol, istat)
    end select

    if (allocated(this%csv_data)) then
      allocate(col_data(this%nrows))  ! size the output vector
      call this%get_column(icol, col_data, istat)
    else
      write(*,'(A,1X,I5)') 'Error: class has not been initialized'
      istat = 1
    end if

  END SUBROUTINE get_character_column

  SUBROUTINE get_string_column(this, col, col_data, istat)
    ! Convert a column from a `string_array` matrix to a `type(string_array)` vector.
    implicit none
    class(csv),                      intent(in)    :: this
    class(*),                        intent(in)    :: col         ! name:character or index:integer
    type(string_array), allocatable, intent(out)   :: col_data(:)
    integer(i4b),                    intent(out)   :: istat
    integer(i4b)                                   :: icol        ! column number

    istat=0

    select type (col)
      type is (integer(i4b))
       icol=col
      type is (character(len=*))
        call this%get_icol_by_name(col, icol, istat)
    end select

    if (allocated(this%csv_data)) then
      allocate(col_data(this%nrows))  ! size the output vector
      call this%get_column(icol, col_data, istat)
    else
      write(*,'(A,1X,I5)') 'Error: class has not been initialized'
      istat = 1
    end if

  END SUBROUTINE get_string_column


  SUBROUTINE split_csv_line(this, line, cells)
    implicit none
    class(csv),                   intent(inout) :: this
    character(len=*),             intent(in)    :: line
    type(string_array), allocatable,intent(out)   :: cells(:)
    ! Local variables
    integer(i4b)                                :: i !! counter
    character(len=:),allocatable                :: tmp !! a temp string with whitespace removed
    integer(i4b)                                :: n !! length of compressed string

    call split(line, this%delimiter, cells)

    ! remove quotes if present:
    do i = 1, size(cells)
        ! remove whitespace from the string:
        tmp = trim(adjustl(cells(i)%str))
        n = len(tmp)

        if (n>1) then
          ! if the first and last non-blank character is
          ! a quote, then remove them and replace with what
          ! is inside the quotes. Otherwise, leave it as is.
          if (tmp(1:1)==this%quote .and. tmp(n:n)==this%quote) then
            if (n>2) then
              cells(i)%str = tmp(2:n-1)  ! remove the quotes
            else
              cells(i)%str = ''  ! empty string
            end if
          end if
        end if
    end do

  END SUBROUTINE split_csv_line

  SUBROUTINE clean(this)
    implicit none
    type(csv), intent(inout) :: this
    integer(i4b)              :: ierr
    character(len=strLen)     :: message

    if (allocated(this%header)) then
      deallocate(this%header, stat=ierr, errmsg=message)
    end if
    if (allocated(this%dtype)) then
      deallocate(this%dtype, stat=ierr, errmsg=message)
    end if
    if (allocated(this%csv_data)) then
      deallocate(this%csv_data, stat=ierr, errmsg=message)
    end if

  END SUBROUTINE clean

 ! *********************************************************************
 ! private subroutine: split a string into sub-strings at delimiter
 ! *********************************************************************
  pure SUBROUTINE split(str, delimiter, vals)
    !  This routine is inspired by the Python split function.
    !  e.g.
    !   character(len=:),allocatable :: s
    !   type(string_array),dimension(:),allocatable :: vals
    !   s = '1,2,3,4,5'
    !   call split(s,',',vals)

    implicit none
    ! Argument variables
    character(len=*),              intent(in)  :: str
    character(len=*),              intent(in)  :: delimiter
    type(string_array),allocatable,intent(out) :: vals(:)
    ! Local variables
    integer(i4b)                             :: ix             ! counter
    integer(i4b)                             :: len_str        ! significant length of `str`
    integer(i4b)                             :: len_delimiter  ! length of the delimiter
    integer(i4b)                             :: n_delimiters   ! number of delimiters
    integer(i4b)                             :: i1             ! index
    integer(i4b)                             :: i2             ! index
    integer(i4b)                             :: j              ! counters
    integer(i4b), allocatable                :: iDelimiters(:) ! start indices of the delimiter locations in `str`

    ! first, count the number of the delimiters
    ! in the string, and get the delimiter indices.
    ! Examples:
    !  ',         '    --> 1
    !  '1234,67,90'    --> 5,8
    !  '123,      '    --> 4

    len_delimiter = len(delimiter)  ! length of the delimiter, usually 1

    ! length of the string
    if (delimiter == ' ') then ! white space
      ! in this case, we can't ignore trailing space
      len_str = len(str)
    else
      ! safe to ignore trailing space when looking for delimiters
      len_str = len_trim(str)
    end if

    j = 1
    n_delimiters = 0
    do
      if (j>len_str) exit      ! end of string, finished
      ix = index(str(j:),delimiter) ! index of next delimiter in remaining string
      if (ix<=0) exit           ! no more delimiters found
      call pop_integer_array(iDelimiters, ix+j-1)  ! save the delimiter location
      n_delimiters = n_delimiters + 1
      j = j + ix + (len_delimiter - 1)
    end do

    allocate(vals(n_delimiters+1))

    if (n_delimiters>0) then
        len_str = len(str)

        i1 = 1
        i2 = idelimiters(1)-1
        if (i2>=i1) then
            vals(1)%str = str(i1:i2)
        else
            vals(1)%str = ''  !the first character is a delimiter
        end if

        do ix=2,n_delimiters
            i1 = idelimiters(ix-1)+len_delimiter
            i2 = idelimiters(ix)-1
            if (i2>=i1) then
                vals(ix)%str = str(i1:i2)
            else
                vals(ix)%str = ''  !empty element (e.g., 'abc,,def')
            end if
        end do

        i1 = idelimiters(n_delimiters) + len_delimiter
        i2 = len_str
        if (idelimiters(n_delimiters)+len_delimiter<=len_str) then
          vals(n_delimiters+1)%str = str(i1:i2)
        else
          vals(n_delimiters+1)%str = ''  !the last character was a delimiter
        end if
    else ! no delimiters present, so just return the original string:
        vals(1)%str = str
    end if

  END SUBROUTINE split

 ! *********************************************************************
 ! private subroutine: poping in array
 ! *********************************************************************
  pure SUBROUTINE pop_integer_array(vec, val)
    implicit none
    integer(i4b), allocatable, intent(inout) :: vec(:)
    integer(i4b),              intent(in)    :: val      ! adding this to `vec`
    integer(i4b)                             :: n
    integer(i4b),allocatable                 :: tmp(:)

    if (allocated(vec)) then
      allocate(tmp(size(vec)+1))
      tmp(1:size(vec)) = vec
      call move_alloc(tmp,vec)
      n = size(vec)
    else ! first element
      n = 1
      allocate(vec(n))
    end if
    vec(n) = val

  END SUBROUTINE pop_integer_array

END MODULE csv_data
