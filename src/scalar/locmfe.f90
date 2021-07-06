!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module local_mfe

  use mfe_constants
  use mfe_types, only: NodeVar, NodeMtx
  use mfe_data
  use local_arrays
  use problem_pde
  implicit none
  private

  public :: preprocessor, res_mass_matrix, eval_mass_matrix, reg_rhs, eval_dfdy

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! PREPROCESSOR --
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine preprocessor (errc)

      integer, intent(out), optional :: errc

      integer :: j
      logical :: check_dx

      if (present(errc)) then
        check_dx = .true.
        errc = 0
      else
        check_dx = .false.
      end if

      do j = 1, ncell

        ! CELL LENGTH
        dx(j) = u(2,j) % x - u(1,j) % x

        if (check_dx .and. dx(j) < dxmin) then
          errc = j
          return
        end if

        ! U DIFFERENCE ACROSS CELL
        du(j) = u(2,j) % u - u(1,j) % u

        ! CELL LENGTH ON SOLUTION MANIFOLD
        l(j) = sqrt( dx(j)**2 + du(j)**2 )

        ! UNIT NORMAL TO THE SOLUTION MANIFOLD
        n(j) % x = - du(j) / l(j)
        n(j) % u =   dx(j) / l(j)

        ! SOLUTION GRADIENT
        dudx(j) = du(j) / dx(j)

      end do

    end subroutine preprocessor

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! RES_MASS_MATRIX
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine res_mass_matrix ()

      integer :: j
      real(kind=wp) :: c, term2, dxdot, dudot
      real(kind=wp), dimension(NVERT) :: ndot, term

     !!!
     !!!  Pure MFE mass matrix terms.

      c = 1.0_wp / 6.0_wp

      do j = 1, ncell

        ! Normal velocity at each vertex.
        ndot(:) = n(j) % x * udot(:,j) % x + n(j) % u * udot(:,j) % u
        term(:) = (c * l(j)) * (sum(ndot) + ndot(:))

        r(:,j) % x = r(:,j) % x - term * n(j) % x
        r(:,j) % u = r(:,j) % u - term * n(j) % u

      end do

     !!!
     !!!  Regularization contribution to the mass matrix.

      select case (kreg)

        case (1)

         !!!
         !!! RATE OF DEFORMATION PENALIZATION

          do j = 1, ncell

            dxdot = udot(2,j) % x - udot(1,j) % x
            dudot = udot(2,j) % u - udot(1,j) % u
            term2 = (eltvsc / l(j)) * (n(j) % u * dxdot - n(j) % x * dudot)

            r(1,j) % x = r(1,j) % x + (term2 * n(j) % u)
            r(2,j) % x = r(2,j) % x - (term2 * n(j) % u)

            r(1,j) % u = r(1,j) % u - (term2 * n(j) % x)
            r(2,j) % u = r(2,j) % u + (term2 * n(j) % x)

          end do

        case (2)

         !!!
         !!! TOTAL GRADIENT PENALIZATION

          do j = 1, ncell

            dxdot = udot(2,j) % x - udot(1,j) % x
            dudot = udot(2,j) % u - udot(1,j) % u
            term2 = eltvsc / l(j)

            r(1,j) % x = r(1,j) % x + (term2 * dxdot)
            r(2,j) % x = r(2,j) % x - (term2 * dxdot)

            r(1,j) % u = r(1,j) % u + (term2 * dudot)
            r(2,j) % u = r(2,j) % u - (term2 * dudot)

          end do

      end select

    end subroutine res_mass_matrix

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! EVAL_MASS_MATRIX
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_mass_matrix (factor, diag_only)

      real(kind=wp), intent(in), optional :: factor
      logical, intent(in), optional :: diag_only

      integer :: j
      logical :: diagonal
      type(NodeMtx) :: block
      real(kind=wp) :: c, fac, term, term_xx, term_xu, term_uu

      if (present(factor)) then
        fac = factor
      else
        fac = 1.0_wp
      end if

      if (present(diag_only)) then
        diagonal = diag_only
      else
        diagonal = .false.
      end if

     !!!
     !!! PURE MFE MASS MATRIX

      c = fac / 3.0_wp

      do j = 1, ncell

        block % xx =  (c * l(j)) * n(j) % x * n(j) % x
        block % xu =  (c * l(j)) * n(j) % x * n(j) % u
        block % ux =  (c * l(j)) * n(j) % x * n(j) % u
        block % uu =  (c * l(j)) * n(j) % u * n(j) % u

        ! COPY THE BASIC BLOCK

        mtx(1,1,j) = block
        mtx(2,2,j) = block

        if (diagonal) then
          cycle
        end if

        block % xx = 0.5_wp * block % xx
        block % xu = 0.5_wp * block % xu
        block % ux = 0.5_wp * block % ux
        block % uu = 0.5_wp * block % uu

        mtx(2,1,j) = block
        mtx(1,2,j) = block

      end do

     !!!
     !!!  Regularization contribution to the mass matrix.

      select case (kreg)

        case (1)

         !!!
         !!! RATE OF DEFORMATION PENALIZATION

          c = fac * eltvsc

          do j = 1, ncell

            term = c / l(j)
            term_xx =   term * n(j) % u * n(j) % u
            term_xu = - term * n(j) % x * n(j) % u
            term_uu =   term * n(j) % x * n(j) % x

            mtx(1,1,j) % xx = mtx(1,1,j) % xx + term_xx
            mtx(1,1,j) % xu = mtx(1,1,j) % xu + term_xu
            mtx(1,1,j) % ux = mtx(1,1,j) % ux + term_xu
            mtx(1,1,j) % uu = mtx(1,1,j) % uu + term_uu

            mtx(2,2,j) % xx = mtx(2,2,j) % xx + term_xx
            mtx(2,2,j) % xu = mtx(2,2,j) % xu + term_xu
            mtx(2,2,j) % ux = mtx(2,2,j) % ux + term_xu
            mtx(2,2,j) % uu = mtx(2,2,j) % uu + term_uu

            if (diagonal) then
              cycle
            end if

            mtx(2,1,j) % xx = mtx(2,1,j) % xx - term_xx
            mtx(2,1,j) % xu = mtx(2,1,j) % xu - term_xu
            mtx(2,1,j) % ux = mtx(2,1,j) % ux - term_xu
            mtx(2,1,j) % uu = mtx(2,1,j) % uu - term_uu

            mtx(1,2,j) % xx = mtx(1,2,j) % xx - term_xx
            mtx(1,2,j) % xu = mtx(1,2,j) % xu - term_xu
            mtx(1,2,j) % ux = mtx(1,2,j) % ux - term_xu
            mtx(1,2,j) % uu = mtx(1,2,j) % uu - term_uu

          end do

        case (2)

         !!!
         !!! TOTAL GRADIENT PENALIZATION

          c = fac * eltvsc

          do j = 1, ncell

            term = c / l(j)

            mtx(1,1,j) % xx = mtx(1,1,j) % xx + term
            mtx(1,1,j) % uu = mtx(1,1,j) % uu + term

            mtx(2,2,j) % xx = mtx(2,2,j) % xx + term
            mtx(2,2,j) % uu = mtx(2,2,j) % uu + term

            if (diagonal) then
              cycle
            end if

            mtx(2,1,j) % xx = mtx(2,1,j) % xx - term
            mtx(2,1,j) % uu = mtx(2,1,j) % uu - term

            mtx(1,2,j) % xx = mtx(1,2,j) % xx - term
            mtx(1,2,j) % uu = mtx(1,2,j) % uu - term

          end do

      end select

    end subroutine eval_mass_matrix

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! REG_RHS --
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reg_rhs ()

      integer :: j
      real(kind=wp) :: term, term_x, term_u

      do j = 1, ncell

        term = segspr / l(j)**2
        term_x =   term * n(j) % u
        term_u = - term * n(j) % x

        r(1,j) % x = r(1,j) % x - term_x
        r(1,j) % u = r(1,j) % u - term_u

        r(2,j) % x = r(2,j) % x + term_x
        r(2,j) % u = r(2,j) % u + term_u

      end do

    end subroutine reg_rhs

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! EVAL_DFDY --
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_dfdy (t, errc)

      real(kind=wp), intent(in) :: t
      integer, intent(out) :: errc

      integer :: j, k, l
      real(kind=wp) :: rh
      type(NodeVar), dimension(NVERT,ncell) :: r_save

      rh = 1.0_wp / fdinc

     !!!
     !!!  Compute and save the unperturbed residual.

      call pde_rhs (t)
      call reg_rhs ()
      call res_mass_matrix ()

      r_save = r

      do k = 1, NVERT

       !!!
       !!! PARTIALS WITH RESPECT TO X AT VERTEX K

        u(k,:) % x = u(k,:) % x + fdinc

        call preprocessor (errc)
        if (errc /= 0) then
          return
        end if
        call pde_rhs (t)
        call reg_rhs ()
        call res_mass_matrix ()

        u(k,:) % x = u(k,:) % x - fdinc

        do j = 1, ncell
          do l = 1, NVERT
            mtx(l,k,j) % xx = mtx(l,k,j) % xx + rh * (r(l,j) % x - r_save(l,j) % x)
            mtx(l,k,j) % ux = mtx(l,k,j) % ux + rh * (r(l,j) % u - r_save(l,j) % u)
          end do
        end do

       !!!
       !!! PARTIALS WITH RESPECT TO THE DEPENDENT VARIABLES AT VERTEX K

        u(k,:) % u = u(k,:) % u + fdinc

        call preprocessor (errc)
        if (errc /= 0) then
          return
        end if
        call pde_rhs (t)
        call reg_rhs ()
        call res_mass_matrix ()

        u(k,:) % u = u(k,:) % u - fdinc

        do j = 1, ncell
          do l = 1, NVERT
            mtx(l,k,j) % xu = mtx(l,k,j) % xu + rh * (r(l,j) % x - r_save(l,j) % x)
            mtx(l,k,j) % uu = mtx(l,k,j) % uu + rh * (r(l,j) % u - r_save(l,j) % u)
          end do
        end do

      end do

      errc = 0

    end subroutine eval_dfdy

end module local_mfe

