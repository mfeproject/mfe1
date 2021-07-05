
module bc_data

  use mfe_constants, only: wp
  private

  type, public :: NodeBC
    integer       :: x_type
    integer       :: u_type
    real(kind=wp) :: x_value
    real(kind=wp) :: u_value
  end type NodeBC

  integer, parameter, public :: FREE  = 0
  integer, parameter, public :: FIXED = 1

  type(NodeBC), save, public :: bc_left, bc_right

end module bc_data

module bc_procs

  use mfe_constants, only: wp
  use bc_data
  use mfe_types, only: NodeVar, NodeMtx
  private

  public :: residual_bc, matrix_bc

  real(kind=wp), parameter, private :: BIG = 1.0e20_wp   !!! FIX THIS !!!

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! RESIDUAL_BC
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine residual_bc (u, r)

      type(NodeVar), dimension(:), intent(in)    :: u
      type(NodeVar), dimension(:), intent(inout) :: r

      integer :: n

      n = size(r)

      !!! LEFT ENDPOINT !!!

      if (bc_left % x_type == FIXED) then
        r(1) % x = r(1) % x + BIG * (u(1) % x - bc_left % x_value)
      end if

      if (bc_left % u_type == FIXED) then
        r(1) % u = r(1) % u + BIG * (u(1) % u - bc_left % u_value)
      end if

      !!! RIGHT ENDPOINT !!!

      if (bc_right % x_type == FIXED) then
        r(n) % x = r(n) % x + BIG * (u(n) % x - bc_right % x_value)
      end if

      if (bc_right % u_type == FIXED) then
        r(n) % u = r(n) % u + BIG * (u(n) % u - bc_right % u_value)
      end if

    end subroutine residual_bc

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! MATRIX_BC
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine matrix_bc (diag)

      type(NodeMtx), dimension(:), intent(inout) :: diag

      integer :: n

      n = size(diag)

      !!! LEFT ENDPOINT !!!

      if (bc_left % x_type == FIXED) then
        diag(1) % xx = diag(1) % xx + BIG
      end if

      if (bc_left % u_type == FIXED) then
        diag(1) % uu = diag(1) % uu + BIG
      end if

      !!! RIGHT ENDPOINT !!!

      if (bc_right % x_type == FIXED) then
        diag(n) % xx = diag(n) % xx + BIG
      end if

      if (bc_right % u_type == FIXED) then
        diag(n) % uu = diag(n) % uu + BIG
      end if

    end subroutine matrix_bc

end module bc_procs
