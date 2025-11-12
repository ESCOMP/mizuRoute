MODULE process_gage_meta

  USE nrtype
  USE public_var, only: integerMissing

  implicit none

  private
  public :: read_gage_meta
  public :: reach_subset

  CONTAINS

    SUBROUTINE read_gage_meta(gageCsvName, ierr, message)

      USE globalData, ONLY: gage_meta_data
      USE gageMeta_data, ONLY: gageMeta            ! gage meta class

      implicit none
      ! Argument variables
      character(*),intent(in)    :: gageCsvName       ! gauge meta csv name
      integer(i4b),intent(out)   :: ierr              ! error code
      character(*),intent(out)   :: message           ! error message
      ! local variables
      character(len=strLen)      :: cmessage          ! error message from subroutine

      ierr=0; message='read_gage_meta'

      gage_meta_data = gageMeta(gageCsvName, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    END SUBROUTINE read_gage_meta

    SUBROUTINE reach_subset(reachID_in, gage_data_in, ierr, message, compdof, index2)

      USE nr_utils,  ONLY: match_index
      USE nr_utils,  ONLY: arth
      USE gageMeta_data, ONLY: gageMeta

      implicit none
      ! Argument variables
      integer(i4b), allocatable,           intent(in)  :: reachID_in(:)
      type(gageMeta),                      intent(in)  :: gage_data_in
      integer(i4b),                        intent(out) :: ierr
      character(len=strLen),               intent(out) :: message
      integer(i4b), allocatable, optional, intent(out) :: compdof(:)    ! gindex in gauge_data space
      integer(i4b), allocatable, optional, intent(out) :: index2(:)     ! index in the same size as reachID_in
      ! local variables
      integer(i4b), allocatable                        :: index1(:)
      integer(i4b)                                     :: nGage
      logical(lgt), allocatable                        :: mask(:)
      character(len=strLen)                            :: cmessage

      ierr=0; message='reach_subset/'

      ! array size = number of gauges in each mpi core
      allocate(index1(gage_data_in%gage_num()), mask(gage_data_in%gage_num()))

      ! Find index of matching reachID in reachID_in array
      index1 = match_index(reachID_in, gage_data_in%reach_id(), ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); endif

      mask=(index1/=integerMissing)
      nGage = count(mask)

      if (nGage>0) then
        if (present(compdof)) then
          allocate(compdof(nGage))
          compdof = pack(arth(1,1,gage_data_in%gage_num()), mask)
        end if
        if (present(index2)) then
          allocate(index2(nGage))
          index2  = pack(index1, mask)
        end if
      else
        if (present(compdof)) then
          allocate(compdof(1))
          compdof = 0
        end if
        if (present(index2)) then
          allocate(index2(1))
          index2 = integerMissing
        end if
      end if

    END SUBROUTINE reach_subset

END MODULE process_gage_meta
