MODULE globalData

  ! General rules
  ! This module include data/variables that can be accessed from any other modules
  ! User do not directly interface the data/variables
  ! Data/values can be altered throughout the runtime
  ! See public_var.f90 for difference

  USE nrtype

  implicit none

  save

  ! ---------- MPI/OMP/PIO variables ----------------------------------------------------------------

  integer(i4b),                    public :: mpicom_route                       ! communicator for this program
  integer(i4b),                    public :: pid                                ! process id
  integer(i4b),                    public :: nNodes                             ! number of MPI processors
  integer(i4b),                    public :: nThreads                           ! number of threads
  logical(lgt),                    public :: masterproc                         ! root logical. root processor => true, other => false
  logical(lgt),                    public :: multiProcs                         ! MPI multi-processors logical. => number of processors>1 true, other => false

END MODULE globalData
