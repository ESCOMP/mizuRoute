MODULE remap_hru
 USE nrtype

 implicit none

 save

 ! data to remap runoff hru to river network hrus
 type remap_data
   integer(i4b),dimension(:),allocatable  :: hru_id    ! Id of hrus associated with river network (="hru")
   integer(i4b),dimension(:),allocatable  :: qhru_id   ! id of hrus associated with runoff simulation (="qhru")
   real(dp),    dimension(:),allocatable  :: weight    ! area weight of "qhru" intersecting "hru"
   integer(i4b),dimension(:),allocatable  :: num_qhru  ! number of "qhru" within "hru"
 end type remap_data

end module
