MODULE data_runoff
 USE nrtype

 implicit none

 ! simulated runoff data
 type runoff
   integer(i4b),dimension(:),  allocatable  :: hru_id    ! id of hrus at which runoff is simulated
   real(dp),    dimension(:),  allocatable  :: time      ! time
   real(dp),    dimension(:),  allocatable  :: qsim      ! runoff(hru) at one time step
 end type runoff

 type(runoff), save :: runoff_data

end module
