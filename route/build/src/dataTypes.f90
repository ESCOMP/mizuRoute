module dataTypes
USE nrtype
USE nrtype, only: strLen   ! string length
! used to create specific data types
! --
implicit none

 ! everything private unless specified otherwise
 private

 ! ---------- information on the data structures -----------------------------------------------------------

 ! data structure information
 type,public :: struct_info
  character(len=strLen)  :: structName = 'empty'          ! name of the data structure
  integer(i4b)           :: varType    = integerMissing   ! type of variable in each data structure
  integer(i4b)           :: nSpace     = integerMissing   ! length of the spatial dimension in each data structure
  integer(i4b)           :: nVars      = integerMissing   ! number of variables in each data structure
 end type struct_info

 ! ---------- metadata structures --------------------------------------------------------------------------

 ! define derived type for model variables, including name, description, and units
 type,public :: var_info
  character(len=strLen)  :: varName  = 'empty'         ! variable name
  character(len=strLen)  :: varDesc  = 'empty'         ! variable description
  character(len=strLen)  :: varUnit  = 'empty'         ! variable units
  integer(i4b)           :: varType  = integerMissing  ! variable type (vectors of different size)
  logical(lgt)           :: varFile  = .true.          ! .true. if the variable should be read from a file
 endtype var_info

 ! ---------- general data structures ----------------------------------------------------------------------

 ! ** double precision type
 ! dat vector
 type, public :: dlength
  real(dp),allocatable                :: dat(:)    ! dat(:)
 endtype dlength
 ! var vector
 type, public :: var_dlength
  type(dlength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_dlength

 ! ** integer type
 ! dat vector
 type, public :: ilength
  integer(i4b),allocatable            :: dat(:)    ! dat(:)
 endtype ilength
 ! var vector
 type, public :: var_ilength
  type(ilength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_ilength

 ! ---------- specific data structures ----------------------------------------------------------------------

 ! data to remap runoff hru to river network hrus
 type, public :: remap
   integer(i4b),dimension(:),  allocatable  :: hru_id    ! Id of hrus associated with river network (="hru")
   integer(i4b),dimension(:),  allocatable  :: qhru_id   ! id of hrus associated with runoff simulation (="qhru")
   real(dp),    dimension(:),  allocatable  :: weight    ! area weight of "qhru" intersecting "hru"
   integer(i4b),dimension(:),  allocatable  :: num_qhru  ! number of "qhru" within "hru"
 end type remap

 ! simulated runoff data
 type, public :: runoff
   integer(i4b),dimension(:),  allocatable  :: hru_id    ! id of hrus at which runoff is simulated
   real(dp),    dimension(:),  allocatable  :: time      ! time
   real(dp),    dimension(:),  allocatable  :: qsim      ! runoff(hru) at one time step
 end type runoff

 ! ---------- unused data structures ------------------------------------------------------------------------

 ! data type containing a name and a variable (double precision)
 type, public :: namepvar
  character(len=256)    :: varName
  real(dp), pointer     :: varData(:) => null()
 endtype namepvar

 ! data type containing a name and a variable (integer)
 type, public :: nameivar
  character(len=256)    :: varName
  integer(i4b), pointer :: varData(:) => null()
 endtype nameivar

end module dataTypes
