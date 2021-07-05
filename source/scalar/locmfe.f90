
module local_mfe

  use mfe_constants, only: wp
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

      integer :: i
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

        ! Solution difference across the element.
        du = u(2,i) % u - u(1,i) % u

        ! Solution manifold element length.
        l(i) = sqrt (du**2 + dx(i)**2)

        ! Unit normal to the solution manifold
        n(i) % x = - du * (1.0_wp / l(i))
        n(i) % u = dx(i) * (1.0_wp / l(i))

        ! Solution gradient
        dudx(i) = du * (1.0_wp / dx(i))

      end do

    end subroutine preprocessor

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! RES_MASS_MATRIX
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine res_mass_matrix ()

      integer :: i
      real(kind=wp) :: c, ndot1, ndot2, sum_ndot, term, term1, term2, &
                       dxdot, dudot

     !!!
     !!!  Pure MFE mass matrix terms.

      c = 1.0_wp / 6.0_wp

      do i = 1, nelt

        ! Normal velocity of each node.
        ndot1 = n(i) % x * udot(1,i) % x + n(i) % u * udot(1,i) % u
        ndot2 = n(i) % x * udot(2,i) % x + n(i) % u * udot(2,i) % u

        sum_ndot = ndot1 + ndot2

        term1 = c * l(i)
        term2 = term1 * (sum_ndot + ndot1)
        r(1,i) % x = r(1,i) % x - term2 * n(i) % x
        r(1,i) % u = r(1,i) % u - term2 * n(i) % u

        term2 = term1 * (sum_ndot + ndot2)
        r(2,i) % x = r(2,i) % x - term2 * n(i) % x
        r(2,i) % u = r(2,i) % u - term2 * n(i) % u

      end do

     !!!
     !!!  Regularization contribution to the mass matrix.

      select case (kreg)

        case (1)

         !!!
         !!! REGULARIZATION ON THE GRADIENT OF THE TANGENTIAL MOTION

          do i = 1, nelt

            dxdot = udot(2,i) % x - udot(1,i) % x
            dudot = udot(2,i) % u - udot(1,i) % u
            term = (eltvsc / l(i)) * (n(i) % u * dxdot - n(i) % x * dudot)

            r(1,i) % x = r(1,i) % x + (term * n(i) % u)
            r(2,i) % x = r(2,i) % x - (term * n(i) % u)

            r(1,i) % u = r(1,i) % u - (term * n(i) % x)
            r(2,i) % u = r(2,i) % u + (term * n(i) % x)

          end do

        case (2)

         !!!
         !!! REGULARIZATION ON THE GRADIENT OF THE X AND U MOTION

          do i = 1, nelt

            dxdot = udot(2,i) % x - udot(1,i) % x
            dudot = udot(2,i) % u - udot(1,i) % u
            term = eltvsc / l(i)

            r(1,i) % x = r(1,i) % x + (term * dxdot)
            r(2,i) % x = r(2,i) % x - (term * dxdot)

            r(1,i) % u = r(1,i) % u + (term * dudot)
            r(2,i) % u = r(2,i) % u - (term * dudot)

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

      integer :: i
      logical :: diagonal_only
      type(NodeMtx) :: block
      real(kind=wp) :: fact, c, term, term_xx, term_xu, term_uu

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

      c = fact / 3.0_wp

      do i = 1, nelt

        term = c * l(i)
        block % xx = term * n(i) % x * n(i) % x
        block % xu = term * n(i) % x * n(i) % u
        block % ux = block % xu
        block % uu = term * n(i) % u * n(i) % u

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

          c = fact * eltvsc

          do i = 1, nelt

            term = c / l(i)
            term_xx =   term * n(i) % u * n(i) % u
            term_xu = - term * n(i) % x * n(i) % u
            term_uu =   term * n(i) % x * n(i) % x

            mtx(1,1,i) % xx = mtx(1,1,i) % xx + term_xx
            mtx(1,1,i) % xu = mtx(1,1,i) % xu + term_xu
            mtx(1,1,i) % ux = mtx(1,1,i) % ux + term_xu
            mtx(1,1,i) % uu = mtx(1,1,i) % uu + term_uu

            mtx(2,2,i) % xx = mtx(2,2,i) % xx + term_xx
            mtx(2,2,i) % xu = mtx(2,2,i) % xu + term_xu
            mtx(2,2,i) % ux = mtx(2,2,i) % ux + term_xu
            mtx(2,2,i) % uu = mtx(2,2,i) % uu + term_uu

            if (diagonal_only) then
              cycle
            end if

            mtx(2,1,i) % xx = mtx(2,1,i) % xx - term_xx
            mtx(2,1,i) % xu = mtx(2,1,i) % xu - term_xu
            mtx(2,1,i) % ux = mtx(2,1,i) % ux - term_xu
            mtx(2,1,i) % uu = mtx(2,1,i) % uu - term_uu

            mtx(1,2,i) % xx = mtx(1,2,i) % xx - term_xx
            mtx(1,2,i) % xu = mtx(1,2,i) % xu - term_xu
            mtx(1,2,i) % ux = mtx(1,2,i) % ux - term_xu
            mtx(1,2,i) % uu = mtx(1,2,i) % uu - term_uu

          end do

        case (2)

         !!!
         !!! REGULARIZATION ON THE GRADIENT OF THE X AND U MOTION

          c = fact * eltvsc

          do i = 1, nelt

            term = c / l(i)

            mtx(1,1,i) % xx = mtx(1,1,i) % xx + term
            mtx(1,1,i) % uu = mtx(1,1,i) % uu + term

            mtx(2,2,i) % xx = mtx(2,2,i) % xx + term
            mtx(2,2,i) % uu = mtx(2,2,i) % uu + term

            if (diagonal_only) then
              cycle
            end if

            mtx(2,1,i) % xx = mtx(2,1,i) % xx - term
            mtx(2,1,i) % uu = mtx(2,1,i) % uu - term

            mtx(1,2,i) % xx = mtx(1,2,i) % xx - term
            mtx(1,2,i) % uu = mtx(1,2,i) % uu - term

          end do

      end select

    end subroutine eval_mass_matrix

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! REG_RHS --
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reg_rhs ()

      integer :: i
      real(kind=wp) :: term, term_x, term_u

      do i = 1, nelt

        term = segspr / l(i) ** 2
        term_x =   term * n(i) % u
        term_u = - term * n(i) % x

        r(1,i) % x = r(1,i) % x - term_x
        r(1,i) % u = r(1,i) % u - term_u

        r(2,i) % x = r(2,i) % x + term_x
        r(2,i) % u = r(2,i) % u + term_u

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

      integer :: i, vk
      real(kind=wp) :: rh
      type(NodeVar), dimension(2,nelt) :: r0

      rh = 1.0_wp / fdinc

     !!!
     !!!  Compute and save the unperturbed residual.

      call pde_rhs (t)
      call reg_rhs ()
      call res_mass_matrix ()

      r0 = r

      do vk = 1, 2

       !!!
       !!! PARTIALS WITH RESPECT TO X AT VERTEX J

        u(vk,:) % x = u(vk,:) % x + fdinc

        call preprocessor (errc)
        if (errc /= 0) then
          return
        end if
        call pde_rhs (t)
        call reg_rhs ()
        call res_mass_matrix ()

        u(vk,:) % x = u(vk,:) % x - fdinc

        do i = 1, nelt
          mtx(1,vk,i) % xx = mtx(1,vk,i) % xx + rh * (r(1,i) % x - r0(1,i) % x)
          mtx(1,vk,i) % ux = mtx(1,vk,i) % ux + rh * (r(1,i) % u - r0(1,i) % u)
          mtx(2,vk,i) % xx = mtx(2,vk,i) % xx + rh * (r(2,i) % x - r0(2,i) % x)
          mtx(2,vk,i) % ux = mtx(2,vk,i) % ux + rh * (r(2,i) % u - r0(2,i) % u)
        end do

       !!!
       !!! PARTIALS WITH RESPECT TO U AT VERTEX J

        u(vk,:) % u = u(vk,:) % u + fdinc

        call preprocessor (errc)
        if (errc /= 0) then
          return
        end if
        call pde_rhs (t)
        call reg_rhs ()
        call res_mass_matrix ()

        u(vk,:) % u = u(vk,:) % u - fdinc

        do i = 1, nelt
          mtx(1,vk,i) % xu = mtx(1,vk,i) % xu + rh * (r(1,i) % x - r0(1,i) % x)
          mtx(1,vk,i) % uu = mtx(1,vk,i) % uu + rh * (r(1,i) % u - r0(1,i) % u)
          mtx(2,vk,i) % xu = mtx(2,vk,i) % xu + rh * (r(2,i) % x - r0(2,i) % x)
          mtx(2,vk,i) % uu = mtx(2,vk,i) % uu + rh * (r(2,i) % u - r0(2,i) % u)
        end do

      end do

      errc = 0

    end subroutine eval_dfdy

end module local_mfe

