
module local_mfe

  use mfe_constants, only: wp, neq
  use mfe_types, only: NodeVar, NodeMtx
  use mfe_data
  use local_arrays
  use problem_pde
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

      integer :: i, k
      logical :: check_dx
      real(kind=wp) :: du

      if (present(errc)) then
        check_dx = .true.
        errc = 0
      else
        check_dx = .false.
      end if

      do i = 1, nelt

        ! Element width.
        dx(i) = u(2,i) % x - u(1,i) % x

        if (check_dx) then
          if (dx(i) < dxmin) then
            errc = i
            return
          end if
        end if

        do k = 1, neq

          ! Solution difference across the element.
          du = u(2,i) % u(k) - u(1,i) % u(k)

          ! Solution manifold element length.
          l(k,i) = sqrt (du**2 + dx(i)**2)

          ! Unit normal to the solution manifold
          n(k,i) % x = - du * (1.0_wp / l(k,i))
          n(k,i) % u = dx(i) * (1.0_wp / l(k,i))

          ! Solution gradient
          dudx(k,i) = du * (1.0_wp / dx(i))

        end do

      end do

    end subroutine preprocessor

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! RES_MASS_MATRIX
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine res_mass_matrix ()

      integer :: i, k
      real(kind=wp), dimension(neq) :: c
      real(kind=wp) :: r1x, r2x, ndot1, ndot2, sum_ndot, term, term1, term2, &
                       dxdot, dudot

     !!!
     !!!  Pure MFE mass matrix terms.

      c = eqw / 6.0_wp

      do i = 1, nelt

        r1x = r(1,i) % x
        r2x = r(2,i) % x

        do k = 1, neq

          ! Normal velocity of each node.
          ndot1 = n(k,i) % x * udot(1,i) % x + n(k,i) % u * udot(1,i) % u(k)
          ndot2 = n(k,i) % x * udot(2,i) % x + n(k,i) % u * udot(2,i) % u(k)

          sum_ndot = ndot1 + ndot2

          term1 = c(k) * l(k,i)
          term2 = term1 * (sum_ndot + ndot1)
          r1x = r1x - term2 * n(k,i) % x
          r(1,i) % u(k) = r(1,i) % u(k) - term2 * n(k,i) % u

          term2 = term1 * (sum_ndot + ndot2)
          r2x = r2x - term2 * n(k,i) % x
          r(2,i) % u(k) = r(2,i) % u(k) - term2 * n(k,i) % u

        end do

        r(1,i) % x = r1x
        r(2,i) % x = r2x

      end do

     !!!
     !!!  Regularization contribution to the mass matrix.

      select case (kreg)

        case (1)

         !!!
         !!! REGULARIZATION ON THE GRADIENT OF THE TANGENTIAL MOTION

          c = eqw * eltvsc

          do i = 1, nelt

            r1x = r(1,i) % x
            r2x = r(2,i) % x

            dxdot = udot(2,i) % x - udot(1,i) % x

            do k = 1, neq

              dudot = udot(2,i) % u(k) - udot(1,i) % u(k)
              term = (c(k) / l(k,i)) * (n(k,i) % u * dxdot - n(k,i) % x * dudot)

              r1x = r1x + (term * n(k,i) % u)
              r2x = r2x - (term * n(k,i) % u)

              r(1,i) % u(k) = r(1,i) % u(k) - (term * n(k,i) % x)
              r(2,i) % u(k) = r(2,i) % u(k) + (term * n(k,i) % x)

            end do

            r(1,i) % x = r1x
            r(2,i) % x = r2x

          end do

        case (2)

         !!!
         !!! REGULARIZATION ON THE GRADIENT OF THE X AND U MOTION

          c = eqw * eltvsc

          do i = 1, nelt

            r1x = r(1,i) % x
            r2x = r(2,i) % x

            dxdot = udot(2,i) % x - udot(1,i) % x

            do k = 1, neq

              dudot = udot(2,i) % u(k) - udot(1,i) % u(k)
              term = c(k) / l(k,i)

              r1x = r1x + (term * dxdot)
              r2x = r2x - (term * dxdot)

              r(1,i) % u(k) = r(1,i) % u(k) + (term * dudot)
              r(2,i) % u(k) = r(2,i) % u(k) - (term * dudot)

            end do

            r(1,i) % x = r1x
            r(2,i) % x = r2x

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

      integer :: i, k
      logical :: diagonal_only
      type(NodeMtx) :: block
      real(kind=wp) :: fact, term, term_xx, term_xu, term_uu
      real(kind=wp), dimension(neq) :: c

      if (present(factor)) then
        fact = factor
      else
        fact = 1.0_wp
      end if

      if (present(diag_only)) then
        diagonal_only = diag_only
      else
        diagonal_only = .false.
      end if

     !!!
     !!! Load the pure MFE mass matrix.

      c = (fact / 3.0_wp) * eqw

      block % uu = 0.0_wp

      do i = 1, nelt

        term_xx = 0.0_wp
        do k = 1, neq
          term = c(k) * l(k,i)
          term_xx = term * n(k,i) % x * n(k,i) % x + term_xx
          term_xu = term * n(k,i) % x * n(k,i) % u
          term_uu = term * n(k,i) % u * n(k,i) % u
          block % xu(k) = term_xu
          block % ux(k) = term_xu
          block % uu(k,k) = term_uu
        end do
        block % xx = term_xx

        ! Copy the basic block
        mtx(1,1,i) = block
        mtx(2,2,i) = block

        if (diagonal_only) then
          cycle
        end if

        block % xx = 0.5_wp * block % xx
        block % xu = 0.5_wp * block % xu
        block % ux = 0.5_wp * block % ux
        block % uu = 0.5_wp * block % uu

        mtx(2,1,i) = block
        mtx(1,2,i) = block

      end do

     !!!
     !!!  Regularization contribution to the mass matrix.

      select case (kreg)

        case (1)

         !!!
         !!! REGULARIZATION ON THE GRADIENT OF THE TANGENTIAL MOTION

          c = fact * eqw * eltvsc

          do i = 1, nelt

            do k = 1, neq

              term = c(k) / l(k,i)
              term_xx =   term * n(k,i) % u * n(k,i) % u
              term_xu = - term * n(k,i) % x * n(k,i) % u
              term_uu =   term * n(k,i) % x * n(k,i) % x

              mtx(1,1,i) % xx      = mtx(1,1,i) % xx      + term_xx
              mtx(1,1,i) % xu(k)   = mtx(1,1,i) % xu(k)   + term_xu
              mtx(1,1,i) % ux(k)   = mtx(1,1,i) % ux(k)   + term_xu
              mtx(1,1,i) % uu(k,k) = mtx(1,1,i) % uu(k,k) + term_uu

              mtx(2,2,i) % xx      = mtx(2,2,i) % xx      + term_xx
              mtx(2,2,i) % xu(k)   = mtx(2,2,i) % xu(k)   + term_xu
              mtx(2,2,i) % ux(k)   = mtx(2,2,i) % ux(k)   + term_xu
              mtx(2,2,i) % uu(k,k) = mtx(2,2,i) % uu(k,k) + term_uu

              if (diagonal_only) then
                cycle
              end if

              mtx(2,1,i) % xx      = mtx(2,1,i) % xx      - term_xx
              mtx(2,1,i) % xu(k)   = mtx(2,1,i) % xu(k)   - term_xu
              mtx(2,1,i) % ux(k)   = mtx(2,1,i) % ux(k)   - term_xu
              mtx(2,1,i) % uu(k,k) = mtx(2,1,i) % uu(k,k) - term_uu

              mtx(1,2,i) % xx      = mtx(1,2,i) % xx      - term_xx
              mtx(1,2,i) % xu(k)   = mtx(1,2,i) % xu(k)   - term_xu
              mtx(1,2,i) % ux(k)   = mtx(1,2,i) % ux(k)   - term_xu
              mtx(1,2,i) % uu(k,k) = mtx(1,2,i) % uu(k,k) - term_uu

            end do

          end do

        case (2)

         !!!
         !!! REGULARIZATION ON THE GRADIENT OF THE X AND U MOTION

          c = fact * eqw * eltvsc

          do i = 1, nelt

            do k = 1, neq

              term = c(k) / l(k,i)

              mtx(1,1,i) % xx      = mtx(1,1,i) % xx      + term
              mtx(1,1,i) % uu(k,k) = mtx(1,1,i) % uu(k,k) + term

              mtx(2,2,i) % xx      = mtx(2,2,i) % xx      + term
              mtx(2,2,i) % uu(k,k) = mtx(2,2,i) % uu(k,k) + term

              if (diagonal_only) then
                cycle
              end if

              mtx(2,1,i) % xx      = mtx(2,1,i) % xx      - term
              mtx(2,1,i) % uu(k,k) = mtx(2,1,i) % uu(k,k) - term

              mtx(1,2,i) % xx      = mtx(1,2,i) % xx      - term
              mtx(1,2,i) % uu(k,k) = mtx(1,2,i) % uu(k,k) - term

            end do

          end do

      end select

    end subroutine eval_mass_matrix

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! REG_RHS --
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reg_rhs ()

      integer :: i, k
      real(kind=wp) :: term, term_x, term_u
      real(kind=wp), dimension(neq) :: c

      c = eqw * segspr

      do i = 1, nelt

        do k = 1, neq

          term = c(k) / l(k,i)**2
          term_x =   term * n(k,i) % u
          term_u = - term * n(k,i) % x

          r(1,i) % x    = r(1,i) % x    - term_x
          r(1,i) % u(k) = r(1,i) % u(k) - term_u

          r(2,i) % x    = r(2,i) % x    + term_x
          r(2,i) % u(k) = r(2,i) % u(k) + term_u

        end do

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

      integer :: i, j, k
      real(kind=wp) :: rh
      type(NodeVar), dimension(2,nelt) :: r_save

      rh = 1.0_wp / fdinc

     !!!
     !!!  Compute and save the unperturbed residual.

      call pde_rhs (t)
      call reg_rhs ()
      call res_mass_matrix ()

      r_save = r

      do j = 1, 2

       !!!
       !!! PARTIALS WITH RESPECT TO X AT VERTEX J

        u(j,:) % x = u(j,:) % x + fdinc

        call preprocessor (errc)
        if (errc /= 0) then
          return
        end if
        call pde_rhs (t)
        call reg_rhs ()
        call res_mass_matrix ()

        u(j,:) % x = u(j,:) % x - fdinc

        do i = 1, nelt
          mtx(1,j,i) % xx = mtx(1,j,i) % xx + rh * (r(1,i) % x - r_save(1,i) % x)
          mtx(1,j,i) % ux = mtx(1,j,i) % ux + rh * (r(1,i) % u - r_save(1,i) % u)
          mtx(2,j,i) % xx = mtx(2,j,i) % xx + rh * (r(2,i) % x - r_save(2,i) % x)
          mtx(2,j,i) % ux = mtx(2,j,i) % ux + rh * (r(2,i) % u - r_save(2,i) % u)
        end do

       !!!
       !!! PARTIALS WITH RESPECT TO THE DEPENDENT VARIABLES AT VERTEX VK

        do k = 1, neq

          u(j,:) % u(k) = u(j,:) % u(k) + fdinc

          call preprocessor (errc)
          if (errc /= 0) then
            return
          end if
          call pde_rhs (t)
          call reg_rhs ()
          call res_mass_matrix ()

          u(j,:) % u(k) = u(j,:) % u(k) - fdinc

          do i = 1, nelt
            mtx(1,j,i) % xu(k)   = mtx(1,j,i) % xu(k)   + rh * (r(1,i) % x - r_save(1,i) % x)
            mtx(1,j,i) % uu(:,k) = mtx(1,j,i) % uu(:,k) + rh * (r(1,i) % u - r_save(1,i) % u)
            mtx(2,j,i) % xu(k)   = mtx(2,j,i) % xu(k)   + rh * (r(2,i) % x - r_save(2,i) % x)
            mtx(2,j,i) % uu(:,k) = mtx(2,j,i) % uu(:,k) + rh * (r(2,i) % u - r_save(2,i) % u)
          end do

        end do

      end do

      errc = 0

    end subroutine eval_dfdy

end module local_mfe

