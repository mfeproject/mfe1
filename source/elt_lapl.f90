!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module element_laplacian

  use precision
  use mfe_extents
  use mfe_data, only: wpde
  use bc_data, only: bcl, bcr
  use element_arrays, only: du, r0, r1

  implicit none
  private

  interface lapl
    module procedure lapl_const_coef, lapl_var_coef
  end interface

  public :: lapl

  real(wp), parameter :: eta = 0.01_wp, c3 = 1.0_wp / 3.0_wp, &
                         c5 = 1.0_wp / 5.0_wp, c7 = 1.0_wp / 7.0_wp

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! LAPL (generic)
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lapl_const_coef (k, visc)

      integer, intent(in)  :: k
      real(wp), intent(in) :: visc

      integer  :: j
      real(wp) :: c, dudx, dS, e, i1(nelt), i2(nelt)

      c = wpde(k) * visc

      do j = 1, nelt

        dudx = du(k,j) / du(x,j)
        dS = sqrt (1.0_wp + dudx**2)
        e = dudx / dS

        i2(j) = c * dudx**2 / (1.0_wp + dS)

        if (abs(e) > eta) then
          i1(j) = c * sign ( log(abs(dudx) + dS), dudx )
        else
          i1(j) = c * e * (1.0_wp + (e**2) * (c3 + (e**2) * (c5 + c7 * (e**2))))
        end if

      end do

      r1(x,1) = r1(x,1) + i2(1)
      r1(k,1) = r1(k,1) - i1(1)

      do j = 2, nelt - 1

        r0(x,j) = r0(x,j) - i2(j)
        r0(k,j) = r0(k,j) + i1(j)

        r1(x,j) = r1(x,j) + i2(j)
        r1(k,j) = r1(k,j) - i1(j)

      end do

      r0(x,nelt) = r0(x,nelt) - i2(nelt)
      r0(k,nelt) = r0(k,nelt) + i1(nelt)

      if (bcl(k) == 2) r0(k,1)    = r0(k,1)    + i1(1)
      if (bcr(k) == 2) r1(k,nelt) = r1(k,nelt) - i1(nelt)

    end subroutine lapl_const_coef


    subroutine lapl_var_coef (k, a0, a1)

      integer, intent(in)  :: k
      real(wp), intent(in) :: a0(:), a1(:)

      integer  :: j
      real(wp) :: dudx, dS, e, i1(nelt), i2(nelt), term

      do j = 1, nelt

        dudx = du(k,j) / du(x,j)
        dS = sqrt (1.0_wp + dudx**2)
        e = dudx / dS

        i2(j) = dudx**2 / (1.0_wp + dS)

        if (abs(e) > eta) then
          i1(j) = sign ( log(abs(dudx) + dS), dudx )
        else
          i1(j) = e * (1.0_wp + (e**2) * (c3 + (e**2) * (c5 + c7 * (e**2))))
        end if

      end do

      term = wpde(k) * a1(1)
      r1(k,1) = r1(k,1) - term * i1(1)
      r1(x,1) = r1(x,1) + term * i2(1)

      do j = 2, nelt - 1

        term = wpde(k) * a0(j)
        r0(x,j) = r0(x,j) - term * i2(j)
        r0(k,j) = r0(k,j) + term * i1(j)

        term = wpde(k) * a1(j)
        r1(x,j) = r1(x,j) + term * i2(j)
        r1(k,j) = r1(k,j) - term * i1(j)

      end do

      term = wpde(k) * a0(nelt)
      r0(k,nelt) = r0(k,nelt) + term * i1(nelt)
      r0(x,nelt) = r0(x,nelt) - term * i2(nelt)

      if (bcl(k) == 2) r0(k,1)    = r0(k,1)    + (wpde(k) * a0(1)) * i1(1)
      if (bcr(k) == 2) r1(k,nelt) = r1(k,nelt) - (wpde(k) * a1(nelt)) * i1(nelt)

    end subroutine lapl_var_coef

end module element_laplacian
