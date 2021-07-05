module local_laplacian

  use mfe_constants, only: wp
  use local_arrays
  private

  public  :: laplacian
  private :: lapl_const_coef, lapl_var_coef

  interface laplacian
    module procedure lapl_const_coef, lapl_var_coef
  end interface

  real(kind=wp), parameter, private :: ETA = 0.01_wp, &
    C3 = 1.0_wp / 3.0_wp, C5 = 1.0_wp / 5.0_wp, C7 = 1.0_wp / 7.0_wp

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! LAPL_CONST_COEF (LAPLACIAN)
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lapl_const_coef (coef)

      real(kind=wp), intent(in) :: coef

      integer :: i
      real(kind=wp) :: r1, m, e, s1, s2

      do i = 1, nelt

        r1 = 1.0_wp / n(i) % u
        m = dudx(i)

       !!!
       !!! FIRST- AND SECOND-KIND INTEGRALS

        e = m / r1
        if (abs(e) > ETA) then
          s1 = - sign ( log(abs(m) + r1), m )
        else
          s1 = - e * (1.0_wp + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
        end if

        s2 = m ** 2 / (1.0_wp + r1)

       !!!
       !!! LOAD THE INNER PRODUCTS

        r(1,i) % x = r(1,i) % x - (coef * s2)
        r(1,i) % u = r(1,i) % u - (coef * s1)

        r(2,i) % x = r(2,i) % x + (coef * s2)
        r(2,i) % u = r(2,i) % u + (coef * s1)

      end do

    end subroutine lapl_const_coef

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! LAPL_VAR_COEF (LAPLACIAN)
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lapl_var_coef (coef)

      real(kind=wp), dimension(:,:), intent(in) :: coef

      integer :: i
      real(kind=wp) :: r1, m, e, s1, s2

      do i = 1, nelt

        r1 = 1.0_wp / n(i) % u
        m = dudx(i)

       !!!
       !!! FIRST- AND SECOND-KIND INTEGRALS

        e = m / r1
        if (abs(e) > ETA) then
          s1 = - sign ( log(abs(m) + r1), m )
        else
          s1 = - e * (1.0_wp + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
        end if

        s2 = m ** 2 / (1.0_wp + r1)

       !!!
       !!! LOAD THE INNER PRODUCTS

        r(1,i) % x = r(1,i) % x - coef(1,i) * s2
        r(1,i) % u = r(1,i) % u - coef(1,i) * s1

        r(2,i) % x = r(2,i) % x + coef(2,i) * s2
        r(2,i) % u = r(2,i) % u + coef(2,i) * s1

      end do

    end subroutine lapl_var_coef

end module local_laplacian
