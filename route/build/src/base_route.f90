MODULE base_route

  ! Description: Definition of base (or template) reach routing method class.
  ! this abstract class needs to be extended to specific routing method types for
  ! implementation and instantiation.

  implicit none

  private
  public:: base_route_rch   ! base (abstract) reach routing method class (to be extended to specific)
  public:: routeContainer   ! a holder of instantiated reach routing method object


  ! --- base (abstract or template) reach routing method
  type, abstract :: base_route_rch
    CONTAINS
    procedure(sub_route_rch), deferred :: route
  end type

  ! --- routing method container
  ! This container (holder) include instantiated reach routing method
  type :: routeContainer
    class(base_route_rch), allocatable :: rch_route
  end type


  ABSTRACT INTERFACE

    SUBROUTINE sub_route_rch(this,          & ! object to bound the procedure
                             segIndex,      & ! input: reach indice
                             ixDesire,      & ! input: index of verbose reach
                             T0, T1,        & ! input: start and end of simulation time step since start [sec]
                             NETOPO_in,     & ! input: reach topology data structure
                             RPARAM_in,     & ! input: reach parameter data structure
                             RCHSTA_out,    & ! inout: reach state data structure
                             RCHFLX_out,    & ! inout: reach flux data structure
                             ierr, message)   ! output: error control

      ! Description: template interfade for reach routing subroutine
      !   to perform a routing (after instantiated) at a given reasch (segIndex) and time step
      !   reach parameters (RPARAM), river network topology (NETOPO) to get upstream location,
      !   state (RCHSTA) and flux (RCHFLX) are required for a set of input
      !   ixDesire is index of reach where more information is writting in log along the computation

      USE nrtype
      USE dataTypes, ONLY: STRFLX       ! fluxes in each reach
      USE dataTypes, ONLY: STRSTA       ! states in each reach
      USE dataTypes, ONLY: RCHTOPO      ! Network topology
      USE dataTypes, ONLY: RCHPRP       ! Reach parameter

      import base_route_rch
      ! Arguments
      class(base_route_rch)                     :: this             ! object to bound the procedure
      integer(i4b),  intent(in)                 :: segIndex         ! segment index
      integer(i4b),  intent(in)                 :: ixDesire         ! index of the reach for verbose output
      real(dp),      intent(in)                 :: T0, T1           ! start and end of the time step (seconds)
      type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)     ! River Network topology
      type(RCHPRP),  intent(inout), allocatable :: RPARAM_in(:)     ! River reach parameter
      type(STRSTA),  intent(inout)              :: RCHSTA_out(:)    ! reach state data
      type(STRFLX),  intent(inout)              :: RCHFLX_out(:)    ! Reach fluxes (space [reaches]) for decomposed domains
      integer(i4b),  intent(out)                :: ierr             ! error code
      character(*),  intent(out)                :: message          ! error message

    END SUBROUTINE sub_route_rch

  END INTERFACE

END MODULE base_route
