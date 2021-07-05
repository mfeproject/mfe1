!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module element_mfe

  use precision
  use mfe_extents
  use element_arrays
  use problem_procs, only: eltrhs

  implicit none

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! ELTFY
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eltfy (t, errc)

      use mfe_data, only: fdinc

      real(wp), intent(in) :: t
      integer, intent(out) :: errc

      integer k
      real(wp), dimension(nepn,nelt) :: r0save, r1save

     !!!
     !!!  Compute and save the unperturbed residual.

      call eltrhs (t)      ! locpp already called by caller.
      call regrhs
      call cres

      r0save = r0
      r1save = r1

      do k = 1, nepn

       !!!
       !!! Partials with respect to the kth unknown at node 0.

        u0(k,:) = u0(k,:) + fdinc

        call eltpp (errc)
        if (errc /= 0) return
        call eltrhs (t)
        call regrhs
        call cres

        u0(k,:) = u0(k,:) - fdinc

        a00(:,k,:) = a00(:,k,:) + (r0 - r0save) / fdinc
        a10(:,k,:) = a10(:,k,:) + (r1 - r1save) / fdinc

       !!!
       !!! Partials with respect to the kth unknown at node 1.

        u1(k,:) = u1(k,:) + fdinc

        call eltpp (errc)
        if (errc /= 0) return
        call eltrhs (t)
        call regrhs
        call cres

        u1(k,:) = u1(k,:) - fdinc

        a01(:,k,:) = a01(:,k,:) + (r0 - r0save) / fdinc
        a11(:,k,:) = a11(:,k,:) + (r1 - r1save) / fdinc

      end do

    end subroutine eltfy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! ELTPP
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eltpp (errc)

      use mfe_data, only: dxmin

      integer, intent(out) :: errc

      integer j

      do j = 1, nelt

        du(:,j) = u1(:,j) - u0(:,j)

        if (du(nepn,j) < dxmin) then
           errc = j
           return
        end if

        l(:,j) = sqrt(du(1:npde,j)**2 + (du(nepn,j)**2))
        n1(:,j) = -du(1:npde,j) / l(:,j)
        n2(:,j) = du(nepn,j) / l(:,j)

      end do

      errc = 0

    end subroutine eltpp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! CRES
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cres

      use mfe_data, only: wpde, kreg, segvsc

      integer j, k
      real(wp) c(npde), nU0, nU1, term, term1, term2

     !!!
     !!!  Pure MFE mass matrix terms.

      c = wpde / 6.0_wp

      do j = 1, nelt
      do k = 1, npde

        nU0 = n1(k,j)*u0dot(x,j) + n2(k,j)*u0dot(k,j)
        nU1 = n1(k,j)*u1dot(x,j) + n2(k,j)*u1dot(k,j)

        term1 = c(k) * l(k,j)
        term2 = term1 * (2_wp * nU0 + nU1)
        r0(x,j) = r0(x,j) - term2 * n1(k,j)
        r0(k,j) = r0(k,j) - term2 * n2(k,j)

        term2 = term1 * (nU0 + 2_wp * nU1)
        r1(x,j) = r1(x,j) - term2 * n1(k,j)
        r1(k,j) = r1(k,j) - term2 * n2(k,j)

      end do
      end do

     !!!
     !!!  Regularization contribution to the mass matrix.

      select case (kreg)

        case (1)

          c = wpde * segvsc

          do j = 1, nelt       ! Old regularization on grad of tangential motion.
          do k = 1, npde

             term = (c(k) / l(k,j)) * (  n2(k,j) * (u1dot(x,j) - u0dot(x,j)) &
                                       - n1(k,j) * (u1dot(k,j) - u0dot(k,j)) )

             r0(x,j) = r0(x,j) + (term * n2(k,j))
             r0(k,j) = r0(k,j) - (term * n1(k,j))

             r1(x,j) = r1(x,j) - (term * n2(k,j))
             r1(k,j) = r1(k,j) + (term * n1(k,j))

          end do
          end do

        case (2)

          c = wpde * segvsc

          do j = 1, nelt       ! New regularization on grad of x and u motion.
          do k = 1, npde

             r0(x,j) = r0(x,j) - ((c(k) / l(k,j)) * (u0dot(x,j) - u1dot(x,j)))
             r0(k,j) = r0(k,j) - ((c(k) / l(k,j)) * (u0dot(k,j) - u1dot(k,j)))

             r1(x,j) = r1(x,j) + ((c(k) / l(k,j)) * (u0dot(x,j) - u1dot(x,j)))
             r1(k,j) = r1(k,j) + ((c(k) / l(k,j)) * (u0dot(k,j) - u1dot(k,j)))

          end do
          end do

      end select

    end subroutine cres

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! CMTX
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cmtx (factor)

      use mfe_data, only: wpde, kreg, segvsc

      real(wp), intent(in), optional :: factor

      integer j, k
      real(wp) fac, c(npde), aa, ab, bb, term

      fac = 1.0_wp
      if (present(factor)) fac = factor

     !!!
     !!! Load the pure MFE mass matrix.

      c = (fac / 6.0_wp) * wpde

      do j = 1, nelt

        a00(:,:,j) = 0.0_wp
        a10(:,:,j) = 0.0_wp
        a01(:,:,j) = 0.0_wp
        a11(:,:,j) = 0.0_wp

        bb = 0.0_wp

        do k = 1, npde

          aa = (c(k) * l(k,j)) * n2(k,j)**2
          ab = (c(k) * l(k,j)) * n1(k,j) * n2(k,j)
          bb = bb + (c(k) * l(k,j)) * n1(k,j)**2

          a00(k,k,j) = 2.0_wp * aa      ! Load the (alpha,alpha) inner products.
          a10(k,k,j) = aa
          a01(k,k,j) = aa
          a11(k,k,j) = 2.0_wp * aa

          a00(k,x,j) = 2.0_wp * ab      ! Load the (alpha,beta) inner products.
          a00(x,k,j) = 2.0_wp * ab
          a10(k,x,j) = ab
          a10(x,k,j) = ab
          a01(k,x,j) = ab
          a01(x,k,j) = ab
          a11(k,x,j) = 2.0_wp * ab
          a11(x,k,j) = 2.0_wp * ab

        end do

        a00(x,x,j) = 2.0_wp * bb         ! Load the (beta,beta) inner products.
        a10(x,x,j) = bb
        a01(x,x,j) = bb
        a11(x,x,j) = 2.0_wp * bb

      end do

     !!!
     !!!  Regularization contribution to the mass matrix.

      select case (kreg)

        case (1)    ! Old regularization on the grad of the tangential motion.

          c = fac * wpde * segvsc

          do j = 1, nelt
          do k = 1, npde

            term = (c(k)/l(k,j)) * n1(k,j)**2
            a00(k,k,j) = a00(k,k,j) + term
            a10(k,k,j) = a10(k,k,j) - term
            a01(k,k,j) = a01(k,k,j) - term
            a11(k,k,j) = a11(k,k,j) + term

            term = (c(k)/l(k,j)) * n1(k,j) * n2(k,j)
            a00(k,x,j) = a00(k,x,j) - term
            a00(x,k,j) = a00(x,k,j) - term
            a10(k,x,j) = a10(k,x,j) + term
            a10(x,k,j) = a10(x,k,j) + term
            a01(k,x,j) = a01(k,x,j) + term
            a01(x,k,j) = a01(x,k,j) + term
            a11(k,x,j) = a11(k,x,j) - term
            a11(x,k,j) = a11(x,k,j) - term

            term = (c(k)/l(k,j)) * n2(k,j)**2
            a00(x,x,j) = a00(x,x,j) + term
            a10(x,x,j) = a10(x,x,j) - term
            a01(x,x,j) = a01(x,x,j) - term
            a11(x,x,j) = a11(x,x,j) + term

          end do
          end do

        case (2)    ! New regularization on the grad of x and u motions.

          c = fac * wpde * segvsc

          do j = 1, nelt
          do k = 1, npde

            term = c(k) / l(k,j)

            a00(k,k,j) = a00(k,k,j) + term
            a10(k,k,j) = a10(k,k,j) - term
            a01(k,k,j) = a01(k,k,j) - term
            a11(k,k,j) = a11(k,k,j) + term

            a00(x,x,j) = a00(x,x,j) + term
            a10(x,x,j) = a10(x,x,j) - term
            a01(x,x,j) = a01(x,x,j) - term
            a11(x,x,j) = a11(x,x,j) + term

          end do
          end do

      end select

    end subroutine cmtx

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! CDGL
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cdgl

      use mfe_data, only: wpde, kreg, segvsc

      integer j, k
      real(wp) c(npde), aa, ab, bb, term

     !!!
     !!! Load the pure MFE mass matrix diagonal.

      c = wpde / 3.0_wp

      do j = 1, nelt

         a00(:,:,j) = 0.0_wp
         a11(:,:,j) = 0.0_wp

         bb = 0.0_wp

         do k = 1, npde

            aa = (c(k) * l(k,j)) * n2(k,j)**2
            ab = (c(k) * l(k,j)) * n1(k,j) * n2(k,j)
            bb = bb + (c(k) * l(k,j)) * n1(k,j)**2

            a00(k,k,j) = aa      ! Load the (alpha,alpha) inner products.
            a11(k,k,j) = aa

            a00(k,x,j) = ab      ! Load the (alpha,beta) inner products.
            a00(x,k,j) = ab
            a11(k,x,j) = ab
            a11(x,k,j) = ab

         end do

         a00(x,x,j) = bb         ! Load the (beta,beta) inner products.
         a11(x,x,j) = bb

      end do

     !!!
     !!!  Regularization contribution to the mass matrix diagonal.

      select case (kreg)

        case (1)    ! Old regularization on the grad of the tangential motion.

          c = wpde * segvsc

          do j = 1, nelt
          do k = 1, npde

            term = (c(k)/l(k,j)) * n1(k,j)**2
            a00(k,k,j) = a00(k,k,j) + term
            a11(k,k,j) = a11(k,k,j) + term

            term = (c(k)/l(k,j)) * n1(k,j) * n2(k,j)
            a00(k,x,j) = a00(k,x,j) - term
            a00(x,k,j) = a00(x,k,j) - term
            a11(k,x,j) = a11(k,x,j) - term
            a11(x,k,j) = a11(x,k,j) - term

            term = (c(k)/l(k,j)) * n2(k,j)**2
            a00(x,x,j) = a00(x,x,j) + term
            a11(x,x,j) = a11(x,x,j) + term

          end do
          end do

        case (2)    ! New regularization on the grad of x and u motions.

          c = wpde * segvsc

          do j = 1, nelt
          do k = 1, npde

            term = c(k) / l(k,j)

            a00(k,k,j) = a00(k,k,j) + term
            a11(k,k,j) = a11(k,k,j) + term

            a00(x,x,j) = a00(x,x,j) + term
            a11(x,x,j) = a11(x,x,j) + term

          end do
          end do

      end select

    end subroutine cdgl

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! REGRHS
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine regrhs

      use mfe_data, only: wpde, segspr

      integer j, k
      real(wp) c(npde), term

      c = wpde * segspr

      do j = 1, nelt
      do k = 1, npde

        term = c(k) / l(k,j)**2

        r0(x,j) = r0(x,j) - term * n2(k,j)
        r0(k,j) = r0(k,j) + term * n1(k,j)

        r1(x,j) = r1(x,j) + term * n2(k,j)
        r1(k,j) = r1(k,j) - term * n1(k,j)

      end do
      end do

    end subroutine regrhs

end module element_mfe
