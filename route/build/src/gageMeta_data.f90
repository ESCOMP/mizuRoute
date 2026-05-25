module gageMeta_data

  use nrtype,        only: i4b, lgt, strLen, gageStrLen
  use public_var,    only: integerMissing
  use csv_data,      only: csv
  use nr_utils,      only: match_index

  implicit none

  public:: gageMeta

  type :: gageMeta
    private
    integer(i4b)                   :: nGage
    character(len=gageStrLen), allocatable :: gageID(:)
    integer(i4b),      allocatable :: reachID(:)

    contains

      procedure, public  :: get_reach_index
      procedure, public  :: gage_num => fn_get_gage_num
      generic,   public  :: gage_id  => fn_get_gageID, fn_get_gageID_vec, fn_get_gageID_scalar
      generic,   public  :: reach_id => fn_get_reachID, fn_get_reachID_vec, fn_get_reachID_scalar
      procedure, private :: fn_get_gageID
      procedure, private :: fn_get_gageID_vec
      procedure, private :: fn_get_gageID_scalar
      procedure, private :: fn_get_reachID
      procedure, private :: fn_get_reachID_vec
      procedure, private :: fn_get_reachID_scalar

  end type gageMeta

  private :: fn_get_gage_num

  interface gageMeta
    module procedure constructor
  end interface gageMeta

  contains

    ! -----------------------------------------------------
    ! constructor
    ! -----------------------------------------------------
    function constructor(gageCsvName, ierr, message) RESULT(instGageMeta)

      implicit none
      type(gageMeta)                    :: instGageMeta
      character(*),       intent(in)    :: gageCsvName
      integer(i4b),       intent(out)   :: ierr              ! error code
      character(*),       intent(out)   :: message           ! error message
      ! local variables
      type(csv)                  :: gage_meta_csv
      character(len=strLen)      :: cmessage          ! error message from subroutine

      ierr=0; message='gageMeta'

      ! import gauge meta data
      gage_meta_csv = csv(trim(gageCsvName), ierr, cmessage,  mode='r')
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call gage_meta_csv%csv_read(ierr, cmessage, header_row=1)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call gage_meta_csv%fclose()

      ! populate gauge data
      call gage_meta_csv%get_data('gage_id', instGageMeta%gageID, ierr)
      call gage_meta_csv%get_data('reach_id', instGageMeta%reachID, ierr)
      instGageMeta%nGage = size(instGageMeta%gageID)

    end function constructor

    ! -----------------------------------------------------
    ! function: get gage site number
    ! -----------------------------------------------------
    function fn_get_gage_num(this) result(nGage)
      implicit none
      class(gageMeta), intent(in) :: this
      integer(i4b), allocatable   :: nGage
      nGage = this%nGage
    end function fn_get_gage_num

    ! -----------------------------------------------------
    ! functions: get gage ID
    ! -----------------------------------------------------
    function fn_get_gageID(this) result(gageID)
      implicit none
      class(gageMeta), intent(in) :: this
      character(len=gageStrLen), allocatable :: gageID(:)
      allocate(gageID(this%nGage))
      gageID = this%gageID
    end function fn_get_gageID

    function fn_get_gageID_vec(this, ix) result(gageID)
      implicit none
      class(gageMeta), intent(in) :: this
      integer(i4b)                :: ix(:)
      character(len=gageStrLen), allocatable :: gageID(:)
      allocate(gageID(size(ix)))
      gageID = this%gageID(ix)
    end function fn_get_gageID_vec

    function fn_get_gageID_scalar(this, ix) result(gageID)
      implicit none
      class(gageMeta), intent(in) :: this
      integer(i4b)                :: ix
      character(len=gageStrLen)   :: gageID
      gageID = this%gageID(ix)
    end function fn_get_gageID_scalar

    ! -----------------------------------------------------
    ! functions: get reach ID
    ! -----------------------------------------------------
    function fn_get_reachID(this) result(reachID)
      implicit none
      class(gageMeta), intent(in) :: this
      integer(i4b), allocatable   :: reachID(:)
      allocate(reachID(this%nGage))
      reachID = this%reachID
    end function fn_get_reachID

    function fn_get_reachID_vec(this, ix) result(reachID)
      implicit none
      class(gageMeta), intent(in) :: this
      integer(i4b)                :: ix(:)
      integer(i4b), allocatable   :: reachID(:)
      allocate(reachID(size(ix)))
      reachID = this%reachID(ix)
    end function fn_get_reachID_vec

    function fn_get_reachID_scalar(this, ix) result(reachID)
      implicit none
      class(gageMeta), intent(in) :: this
      integer(i4b)                :: ix
      integer(i4b)                :: reachID
      reachID = this%reachID(ix)
    end function fn_get_reachID_scalar

    ! -----------------------------------------------------
    ! subroutine: index of reach ID array
    ! -----------------------------------------------------
    subroutine get_reach_index(this, reachID_in, reach_index)

      implicit none
      ! Argument variables
      class(gageMeta),                 intent(in)  :: this
      integer(i4b), allocatable,       intent(in)  :: reachID_in(:)
      integer(i4b), allocatable,       intent(out) :: reach_index(:)     ! index in the same size as reachID_in
      ! local variables
      integer(i4b)                                 :: ierr              ! error code
      integer(i4b), allocatable                    :: index1(:)
      integer(i4b)                                 :: nGage
      logical(lgt), allocatable                    :: mask(:)
      character(len=strLen)                        :: cmessage          ! error message from subroutine

      ! array size = number of gauges in each mpi core
      allocate(index1(this%nGage), mask(this%nGage))

      ! Find index of matching reachID in reachID_in array
      index1 = match_index(reachID_in, this%reachID, ierr, cmessage)
      mask=(index1/=integerMissing)
      nGage = count(mask)

      if (nGage>0) then
        allocate(reach_index(nGage))
        reach_index = pack(index1, mask)
      else
        allocate(reach_index(1))
        reach_index = integerMissing
      end if

    end subroutine get_reach_index

end module gageMeta_data
