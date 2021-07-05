module local_arrays

  use mfe_constants, only: wp
  use mfe_types, only: NodeVar, NodeMtx
  private

  type, public :: TwoVec
    real(kind=wp) :: x, u
  end type TwoVec

  integer, save, public :: nelt

  ! Local solution arrays.
  type(NodeVar), dimension(:,:),   allocatable, save, public :: u, udot

  ! Local result arrays.
  type(NodeVar), dimension(:,:),   allocatable, save, public :: r
  type(NodeMtx), dimension(:,:,:), allocatable, save, public :: mtx

  ! Intermediate data arrays.
  real(kind=wp), dimension(:,:),   allocatable, save, public :: l
  type(TwoVec),  dimension(:,:),   allocatable, save, public :: n
  real(kind=wp), dimension(:,:),   allocatable, save, public :: dudx
  real(kind=wp), dimension(:),     allocatable, save, public :: dx

end module local_arrays
