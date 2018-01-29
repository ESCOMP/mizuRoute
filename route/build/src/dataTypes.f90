module dataTypes
USE nrtype
! used to create specific data types
! --
! data type containing a name and a variable (double precision)
type namepvar
 character(len=256)    :: varName
 real(dp), pointer     :: varData(:) => null()
endtype namepvar

! data type containing a name and a variable (integer)
type nameivar
 character(len=256)    :: varName
 integer(i4b), pointer :: varData(:) => null()
endtype nameivar

 ! data to remap runoff hru to river network hrus
 type remap
   integer(i4b),dimension(:),allocatable  :: hru_id    ! Id of hrus associated with river network (="hru")
   integer(i4b),dimension(:),allocatable  :: qhru_id   ! id of hrus associated with runoff simulation (="qhru")
   real(dp),    dimension(:),allocatable  :: weight    ! area weight of "qhru" intersecting "hru"
   integer(i4b),dimension(:),allocatable  :: num_qhru  ! number of "qhru" within "hru"
 end type remap

 ! simulated runoff data
 type runoff
   integer(i4b),dimension(:),  allocatable  :: hru_id    ! id of hrus at which runoff is simulated
   real(dp),    dimension(:),  allocatable  :: time      ! time
   real(dp),    dimension(:),  allocatable  :: qsim      ! runoff(hru) at one time step
 end type runoff


end module dataTypes
