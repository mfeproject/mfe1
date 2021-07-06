!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bc_data

  use mfe_constants, only: wp
  implicit none
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

  public  :: force_bc
  private :: bc_res, bc_diag, bc_jac, bc_udot

  interface force_bc
    module procedure bc_res, bc_diag, bc_jac, bc_udot
  end interface

  contains

    subroutine bc_res (u, r)

      type(NodeVar), dimension(:), intent(in)    :: u
      type(NodeVar), dimension(:), intent(inout) :: r

      integer :: n

      n = size(r)

      !!! LEFT ENDPOINT !!!

      if (bc_left % x_type == FIXED) then
        r(1) % x = u(1) % x - bc_left % x_value
      end if

      if (bc_left % u_type == FIXED) then
        r(1) % u = u(1) % u - bc_left % u_value
      end if

      !!! RIGHT ENDPOINT !!!

      if (bc_right % x_type == FIXED) then
        r(n) % x = u(n) % x - bc_right % x_value
      end if

      if (bc_right % u_type == FIXED) then
        r(n) % u = u(n) % u - bc_right % u_value
      end if

    end subroutine bc_res

    subroutine bc_diag (diag)

      type(NodeMtx), dimension(:), intent(inout) :: diag

      integer :: n

      n = size(diag)

      !!! LEFT ENDPOINT !!!

      if (bc_left % x_type == FIXED) then
        diag(1) % xu = 0.0_wp
        diag(1) % ux = 0.0_wp
        diag(1) % xx = 1.0_wp
      end if

      if (bc_left % u_type == FIXED) then
        diag(1) % uu = 1.0_wp
        diag(1) % ux = 0.0_wp
        diag(1) % xu = 0.0_wp
      end if

      !!! RIGHT ENDPOINT !!!

      if (bc_right % x_type == FIXED) then
        diag(n) % xu = 0.0_wp
        diag(n) % ux = 0.0_wp
        diag(n) % xx = 1.0_wp
      end if

      if (bc_right % u_type == FIXED) then
        diag(n) % uu = 1.0_wp
        diag(n) % ux = 0.0_wp
        diag(n) % xu = 0.0_wp
      end if

    end subroutine bc_diag

    subroutine bc_jac (jac_l, jac_d, jac_u)

      type(NodeMtx), dimension(:), intent(inout) :: jac_l, jac_d, jac_u

      integer :: n

      n = size(jac_d)

      !!! LEFT ENDPOINT !!!

      if (bc_left % x_type == FIXED) then
        jac_d(1) % xu = 0.0_wp
        jac_d(1) % xx = 1.0_wp
        jac_u(1) % xu = 0.0_wp
        jac_u(1) % xx = 0.0_wp
      end if

      if (bc_left % u_type == FIXED) then
        jac_d(1) % uu = 1.0_wp
        jac_d(1) % ux = 0.0_wp
        jac_u(1) % uu = 0.0_wp
        jac_u(1) % ux = 0.0_wp
      end if

      !!! RIGHT ENDPOINT !!!

      if (bc_right % x_type == FIXED) then
        jac_l(n) % xu = 0.0_wp
        jac_l(n) % xx = 0.0_wp
        jac_d(n) % xu = 0.0_wp
        jac_d(n) % xx = 1.0_wp
      end if

      if (bc_right % u_type == FIXED) then
        jac_l(n) % uu = 0.0_wp
        jac_l(n) % ux = 0.0_wp
        jac_d(n) % uu = 1.0_wp
        jac_d(n) % ux = 0.0_wp
      end if

    end subroutine bc_jac

    subroutine bc_udot (a_l, a_d, a_u, g)

      type(NodeMtx), dimension(:), intent(inout) :: a_l, a_d, a_u
      type(NodeVar), dimension(:), intent(inout) :: g

      integer :: n

      n = size(g)

      !!! LEFT ENDPOINT !!!

      if (bc_left % x_type == FIXED) then
        a_d(1) % xu = 0.0_wp
        a_d(1) % xx = 1.0_wp
        a_u(1) % xu = 0.0_wp
        a_u(1) % xx = 0.0_wp
        g(1) % x = 0.0_wp
      end if

      if (bc_left % u_type == FIXED) then
        a_d(1) % uu = 1.0_wp
        a_d(1) % ux = 0.0_wp
        a_u(1) % uu = 0.0_wp
        a_u(1) % ux = 0.0_wp
        g(1) % u = 0.0_wp
      end if

      !!! RIGHT ENDPOINT !!!

      if (bc_right % x_type == FIXED) then
        a_l(n) % xu = 0.0_wp
        a_l(n) % xx = 0.0_wp
        a_d(n) % xu = 0.0_wp
        a_d(n) % xx = 1.0_wp
        g(n) % x = 0.0_wp
      end if

      if (bc_right % u_type == FIXED) then
        a_l(n) % uu = 0.0_wp
        a_l(n) % ux = 0.0_wp
        a_d(n) % uu = 1.0_wp
        a_d(n) % ux = 0.0_wp
        g(n) % u = 0.0_wp
      end if

    end subroutine bc_udot

end module bc_procs
