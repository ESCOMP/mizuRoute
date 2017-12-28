module remap
  use nrtype
  use public_var
  use data_remap, only:remap_data            ! data structures holding the data for remapping runoff hru to river network hru

  implicit none
  private
  public ::remap_runoff

  contains

  subroutine remap_runoff(rn_hru_ids, runoff_in, runoff_hru_in, qsim_remapped, err, message)
    implicit none
    ! input
    integer(i4b),         intent(in)  :: rn_hru_ids(:)
    real(dp),             intent(in)  :: runoff_in(:)
    integer(i4b),         intent(in)  :: runoff_hru_in(:)
    ! output
    real(dp),             intent(out) :: qsim_remapped(:)
    integer(i4b),         intent(out) :: err
    character(len=strLen),intent(out) :: message
    ! local
    real(dp),allocatable              :: map_weights(:)
    real(dp),allocatable              :: runoff_sub(:)
    real(dp)                          :: runoff_out
    integer(i4b),allocatable          :: map_runoff_hrus(:)
    integer(i4b)                      :: rn_hru_id
    integer(i4b)                      :: iHru,jHru,kHru ! hru loop indices
    integer(i4b)                      :: idx
    integer(i4b)                      :: idx_start
    integer(i4b)                      :: idx_end
    integer(i4b)                      :: nHRU
    character(len=strLen)             :: cmessage

    err=0; message="remap_runoff/"

    do iHru = 1, size(rn_hru_ids)
      rn_hru_id = rn_hru_ids(iHru) ! river network hru id where average runoff is computed
      call find_index(rn_hru_id, remap_data%hru_id, idx, err, cmessage)
      if(err/=0)then; message=trim(cmessage)//cmessage; return; endif

      nHRU = remap_data%num_qhru(idx)  ! number of runoff hrus contributing to a river network hru

      allocate(map_weights(nHRU),stat=err)
      if(err/=0)then;message=message//'error allocating map_weights';return;endif
      allocate(map_runoff_hrus(nHRU),stat=err)
      if(err/=0)then;message=message//'error allocating runoff_hru_id';return;endif
      allocate(runoff_sub(nHRU), stat=err)
      if(err/=0)then;message=message//'error allocating runoff_sub';return;endif

      idx_end     = sum(remap_data%num_qhru(1:idx));
      idx_start   = idx_end-nHRU+1

      map_runoff_hrus = remap_data%qhru_id(idx_start:idx_end) ! ids of runoff hrus overlapping a river network hru
      map_weights = remap_data%weight(idx_start:idx_end)      ! weights of runoff hrus

      do jHru = 1,nHRU
        do kHru = 1, size(runoff_hru_in)
          if (runoff_hru_in(kHru) == map_runoff_hrus(jHru)) then
            runoff_sub(jHru) = runoff_in(kHru)
          endif
        end do
      end do

      call aggregate_runoff(runoff_sub, map_weights, runoff_out, err, cmessage)
      if(err/=0)then;message=trim(message)//trim(cmessage);return;endif

      qsim_remapped(iHru) = runoff_out

      deallocate(map_weights,stat=err)
      if(err/=0)then;message=message//'error deallocating map_weights';return;endif
      deallocate(map_runoff_hrus,stat=err)
      if(err/=0)then;message=message//'error deallocating map_runoff_hrus';return;endif
      deallocate(runoff_sub,stat=err)
      if(err/=0)then;message=message//'error deallocating runoff_sub';return;endif

    end do
    return
  end subroutine


  subroutine aggregate_runoff(qsim_in, weight, qsim_remapped, err, message)
    implicit none
    ! input
    real(dp),               intent(in)  :: qsim_in(:)         ! runoff data for all runoff hrus
    real(dp),               intent(in)  :: weight(:)          ! weight of runoff hrus
    ! output
    real(dp),               intent(out) :: qsim_remapped      ! weighted (scaled) value
    integer(i4b),           intent(out) :: err                ! error code
    character(len=strLen),  intent(out) :: message            ! error message for current routine
    ! local
    real(dp),parameter                  :: wgtMin=1.e-50_dp   ! minimum value for weight
    logical(lgt),allocatable            :: mask(:)            ! maks
    real(dp),allocatable                :: weight_packed(:)   ! packed weight vector
    real(dp),allocatable                :: qsim_packed(:)     ! packed data value vector
    real(dp)                            :: weight_sum         ! packed weight vector
    integer(i4b)                        :: nElm_org           ! number of vector elements -original vector
    integer(i4b)                        :: nElm               ! number of vector elements -packed vector
    integer(i4b)                        :: iElm               ! index of vector

    err=0; message="aggregate_runoff/"

    ! Compute size of dimension
    nElm_org=size(weight)
    if (nElm_org /= size(qsim_in))then; err=20;message=trim(message)//'data and weight are different size';return;endif
    ! Create mask
    allocate(mask(nElm_org),stat=err); if(err/=0)then;message=trim(message)//'problem allocating mask';return;endif
    mask=(weight > wgtMin)
    ! Pack vector
    allocate(weight_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'problem allocating weight_packed';return;endif
    allocate(qsim_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'problem allocating qsim_packed';return;endif
    weight_packed=pack(weight, mask)
    qsim_packed=pack(qsim_in,mask)
    ! Re-check size of packed vector
    nElm=size(weight_packed)
    if (nElm /= size(qsim_packed)) then; err=20;message=trim(message)//'packed weight and packed qsim are different size';return;endif
    if (nElm > 0) then
      ! Recompute weight
      weight_sum = sum(weight_packed)
      if (weight_sum /= 1.0) weight_packed = weight_packed / weight_sum

      ! compute weighted runoff
      qsim_remapped = 0._dp
      do iElm = 1, nElm
        qsim_remapped = qsim_remapped + qsim_packed(iElm)*weight_packed(iElm)
      end do
    else
      qsim_remapped=dmiss
    endif

    return
  end subroutine


  subroutine find_index(scl, vec, iSelect, ierr, message)
    ! Find index where the value match up with scl
    implicit none
    !input
    integer(i4b),intent(in)              :: scl
    integer(i4b),intent(in)              :: vec(:)
    integer(i4b),intent(out)             :: iSelect
    integer(i4b), intent(out)            :: ierr      ! error code
    character(*), intent(out)            :: message   ! error message
    integer(i4b)                         :: i(1)

    ! initialize error control
    ierr=0; message='findix/'

    i = minloc(abs(vec - scl))
    iSelect = i(1)  ! de-vectorize the found index
    if(vec(iSelect) /= scl)&
      ierr=60; message=trim(message)//'unable to find matched value'; return
    return
  end subroutine

end module
