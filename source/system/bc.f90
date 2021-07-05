
module bc_data

  use mfe_constants, only: wp, neq
  private

  type, public :: NodeBC
    integer                       :: x_type
    integer,       dimension(neq) :: u_type
    real(kind=wp)                 :: x_value
    real(kind=wp), dimension(neq) :: u_value
  end type NodeBC

  integer, parameter, public :: FREE  = 0
  integer, parameter, public :: FIXED = 1

  type(NodeBC), save, public :: bc_left, bc_right

end module bc_data

module bc_procs

  use mfe_constants, only: wp, neq
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

      where (bc_left % u_type == FIXED)
        r(1) % u = r(1) % u + BIG * (u(1) % u - bc_left % u_value)
      end where

      !!! RIGHT ENDPOINT !!!

      if (bc_right % x_type == FIXED) then
        r(n) % x = r(n) % x + BIG * (u(n) % x - bc_right % x_value)
      end if

      where (bc_right % u_type == FIXED)
        r(n) % u = r(n) % u + BIG * (u(n) % u - bc_right % u_value)
      end where

    end subroutine residual_bc

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! MATRIX_BC
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine matrix_bc (diag)

      type(NodeMtx), dimension(:), intent(inout) :: diag

      integer :: n, k

      n = size(diag)

      !!! LEFT ENDPOINT !!!

      if (bc_left % x_type == FIXED) then
        diag(1) % xx = diag(1) % xx + BIG
      end if

      do k = 1, neq
        if (bc_left % u_type(k) == FIXED) then
          diag(1) % uu(k,k) = diag(1) % uu(k,k) + BIG
        end if
      end do

      !!! RIGHT ENDPOINT !!!

      if (bc_right % x_type == FIXED) then
        diag(n) % xx = diag(n) % xx + BIG
      end if

      do k = 1, neq
        if (bc_right % u_type(k) == FIXED) then
          diag(n) % uu(k,k) = diag(n) % uu(k,k) + BIG
        end if
      end do

    end subroutine matrix_bc

end module bc_procs
