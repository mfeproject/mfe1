!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bls_solver

  use precision

  implicit none
  private

  public :: btfct, btslv, vfct, vslv, vmslv

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  BTFCT -- BLOCK TRIDIAGONAL MATRIX FACTORIZATION.
   !!  BTSLV -- BLOCK TRIDIAGONAL MATRIX SOLVE.
   !!
   !!  THESE ROUTINES SOLVE UNIFORMLY PARTITIONED LINEAR SYSTEMS A*X=B
   !!  FOR THE VECTOR X WHEN A IS A BLOCK TRIDIAGONAL MATRIX.  THE ROUTINE
   !!  BTFCT COMPUTES THE BLOCK LU-FACTORIZATION OF THE COEFFICIENT
   !!  MATRIX A.  IT IS ASSUMED THAT GAUSSIAN ELIMINATION APPLIED TO A
   !!  IS STABLE WITHOUT PIVOTING.  THE MATRIX A IS STORED BY BLOCK
   !!  DIAGONALS.  THE ELEMENTS OF THE FACTORS OVERWRITE THOSE OF A
   !!  IN THE STANDARD MANNER.  SUBSEQUENT CALLS TO THE ROUTINE BTSLV
   !!  SOLVES THE SYSTEM, CARRYING OUT THE FORWARD AND BACKWARD
   !!  SUBSTITUTIONS.
   !!
   !!  ARGUMENT  DESCRIPTION
   !!  --------  -----------
   !!   D,U,L    ON ENTRY, D, U, AND L CONTAIN THE ELEMENTS OF THE
   !!            DIAGONAL, UPPER DIAGONAL, AND LOWER DIAGONAL BLOCKS,
   !!            RESPECTIVELY, OF THE COEFFICIENT MATRIX A.  THEY ARE
   !!            OVERWRITTEN WITH THE ELEMENTS OF THE LOWER AND UPPER
   !!            FACTORS.
   !!
   !!   B        ON ENTRY, THE ARRAY CONTAINS THE ELEMENTS OF THE RHS
   !!            VECTOR AND IS OVERWRITTEN WITH THE SOLUTION VECTOR.
   !!
   !!  SUBROUTINES:  BTFCT -- CMAB, FCT, MSLV.
   !!                BTSLV -- SLV, YMAX.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine btfct (d, u, l)

      real(wp), dimension(:,:,:), intent(inout) :: d, u, l

      ! local variables
      integer  :: m, n, j

      n = size(d,1)  ! Must = size(d,2)
      m = size(d,3)  ! Must have shape(d) = shape(u) = shape(l)

      select case (n)

       case (:1)

         ! Do nothing.

       case default

         do j = 1, m-1

           call fct (d(:,:,j))
           call mslv (d(:,:,j), u(:,:,j))
           call cmab (l(:,:,j+1),  u(:,:,j),  d(:,:,j+1))

         end do

         call fct (d(:,:,m))

      end select

    end subroutine btfct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine btslv (d, u, l, b)

      real(wp), dimension(:,:,:), intent(in)  :: d, u, l
      real(wp), dimension(:,:), intent(inout) :: b

      integer m, n, j

      n = size(d,1)   ! Must = size(d,2) = size(b,1).
      m = size(d,3)   ! Must = size(b,2) and also have
                      !   shape(d) = shape(u) = shape(l).

      select case (n)

        case (:1)

          ! Do nothing.

        case default

          call slv (d(:,:,1),  b(:,1))     ! Forward substitution.

          do j = 2, m

            call ymax (l(:,:,j),  b(:,j-1),  b(:,j))
            call slv (d(:,:,j),  b(:,j))

          end do

          do j = m-1, 1, -1                ! Backward substitution

            call ymax (u(:,:,j),  b(:,j+1),  b(:,j))

          end do

      end select

    end subroutine btslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  VFCT  -- VECTOR OF MATRIX FACTORIZATIONS.
   !!  VSLV  -- VECTOR OF FORWARD AND BACKWARD SOLVES.
   !!  VMSLV -- VECTOR OF FORWARD AND BACKWARD MATRIX SOLVES.
   !!
   !!  THESE ROUTINES SOLVE "VECTORS" OF LINEAR SYSTEMS A*X=B FOR X WHERE
   !!  A IS A SQUARE MATRIX AND EITHER B AND X ARE BOTH VECTORS OR BOTH
   !!  MATRICES.  THE ROUTINE VFCT COMPUTES THE LU-FACTORIZATION OF EACH
   !!  COEFFICIENT MATRIX A.  IT IS ASSUMED THAT GAUSSIAN ELIMINATION
   !!  APPLIED TO EACH A IS STABLE WITHOUT PIVOTING.  THE ELEMENTS OF THE
   !!  FACTORS OVERWRITE THOSE OF A IN THE STANDARD MANNER.  A SUBSEQUENT
   !!  CALL TO THE ROUTINE VSLV THEN SOLVES THE SYSTEMS WHEN THE B'S ARE
   !!  VECTORS, CARRYING OUT THE FORWARD AND BACKWARD SUBSTITUTIONS.  THE
   !!  ROUTINE VMSLV SOLVES THE SYSTEMS WHEN THE B'S ARE MATRICIES.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vfct (a)

      real(wp), intent(inout) :: a(:,:,:)

      ! local variables
      integer  :: m, n, i, j, k, l
      real(wp) :: lkk, lkj, ujk

      n = size(a,1)  ! Must = size(a,2)
      m = size(a,3)

      do l = 1, m

        do k = 1, n

          lkk = a(k,k,l)

          do j = 1, k-1

            lkj = a(k,j,l)
            ujk = a(j,k,l)

            do i = 1, j-1

              lkj = lkj - a(k,i,l) * a(i,j,l)
              ujk = ujk - a(j,i,l) * a(i,k,l)

            end do

            ujk = a(j,j,l) * ujk
            lkk = lkk - lkj * ujk

            a(k,j,l) = lkj
            a(j,k,l) = ujk

          end do

          a(k,k,l) = 1.0_wp / lkk

        end do

      end do

    end subroutine vfct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vslv (a, b)

      real(wp), intent(in)    :: a(:,:,:)
      real(wp), intent(inout) :: b(:,:)

      ! local variables
      integer  :: m, n, j, k, l
      real(wp) :: bk

      n = size(a,1)  ! Must = size(a,2) = size(b,1)
      m = size(a,3)  ! Must = size(b,2)

      do l = 1, m

        do k = 1, n     ! Forward Substition.
          bk = b(k,l)
          do j = 1, k-1
            bk = bk - a(k,j,l) * b(j,l)
          end do
           b(k,l) = bk * a(k,k,l)
        end do

        do k = n-1, 1, -1       ! Backward Substition.
          bk = b(k,l)
          do j = k+1, n
            bk = bk - a(k,j,l) * b(j,l)
          end do
          b(k,l) = bk
        end do

      end do

    end subroutine vslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vmslv (a, b)

      real(wp), intent(in)    :: a(:,:,:)
      real(wp), intent(inout) :: b(:,:,:)

      ! local variables
      integer  :: m, n, p, i, j, k, l
      real(wp) :: bk

      n = size(a,1)  ! Must = size(a,2) = size(b,1)
      m = size(a,3)  ! Must = size(b,3)
      p = size(b,2)

      do l = 1, m

        do i = 1, p

          do k = 1, n           ! Forward Substition.
            bk = b(k,i,l)
            do j = 1, k-1
              bk = bk - a(k,j,l) * b(j,i,l)
            end do
            b(k,i,l) = bk * a(k,k,l)
          end do

          do k = n-1, 1, -1     ! Backward Substition.
            bk = b(k,i,l)
            do j = k+1, n
              bk = bk - a(k,j,l) * b(j,i,l)
            end do
            b(k,i,l) = bk
          end do

        end do

      end do

    end subroutine vmslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  FCT  -- MATRIX FACTORIZATION.
   !!  SLV  -- FORWARD AND BACKWARD SOLVE.
   !!  MSLV -- FORWARD AND BACKWARD MATRIX SOLVE.
   !!
   !!  THESE ROUTINES SOLVE LINEAR SYSTEMS OF THE FORM A*X=B FOR X WHERE
   !!  A IS A SQUARE MATRIX AND B AND X ARE EITHER BOTH VECTORS OR BOTH
   !!  MATRICES.  THE ROUTINE FCT COMPUTES THE LU-FACTORIZATION OF THE
   !!  COEFFICIENT MATRIX A.  IT IS ASSUMED THAT GAUSSIAN ELIMINATION
   !!  APPLIED TO A IS STABLE WITHOUT PIVOTING.  THE ELEMENTS OF THE
   !!  FACTORS OVERWRITE THOSE OF A IN THE STANDARD MANNER.  A SUBSEQUENT
   !!  CALL TO THE ROUTINE SLV THEN SOLVES THE SYSTEM WHEN B IS A VECTOR,
   !!  CARRYING OUT THE FORWARD AND BACKWARD SUBSTITUTIONS.  THE ROUTINE
   !!  MSLV SOLVES THE SYSTEM WHEN B IS A MATRIX.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine fct (a)

      real(wp), intent(inout) :: a(:,:)

      ! local variables
      integer  :: n, i, j, k
      real(wp) :: lkk, lkj, ujk

      n = size(a,1)  ! Must = size(a,2)

      select case (n)

        case (:1)

          ! Do nothing.

        case (2)

          a(1,1) = 1.0_wp / a(1,1)
          a(1,2) = a(1,1) * a(1,2)
          a(2,2) = 1.0_wp / (a(2,2) - a(2,1) * a(1,2))

        case default

          do k = 1, n

            lkk = a(k,k)

            do j = 1, k-1

              lkj = a(k,j)
              ujk = a(j,k)

              do i = 1, j-1

                lkj = lkj - a(k,i) * a(i,j)
                ujk = ujk - a(j,i) * a(i,k)

              end do

              ujk = a(j,j) * ujk
              lkk = lkk - lkj * ujk

              a(k,j) = lkj
              a(j,k) = ujk

            end do

            a(k,k) = 1.0_wp / lkk

          end do

      end select

    end subroutine fct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine slv (a, b)

      real(wp), intent(in)    :: a(:,:)
      real(wp), intent(inout) :: b(:)

      ! local variables
      integer  :: n, j, k
      real(wp) :: bk

      n = size(a,1)  ! Must = size(a,2) = size(b)

      select case (n)

        case (:1)

          ! do nothing.

        case (2)

          b(1) = a(1,1) * b(1)
          b(2) = a(2,2) * (b(2) - a(2,1) * b(1))
          b(1) = b(1) - a(1,2) * b(2)

        case default

          do k = 1, n
            bk = b(k)
            do j = 1, k-1
              bk = bk - a(k,j) * b(j)
            end do
            b(k) = bk * a(k,k)
          end do

          do k = n-1, 1, -1
            bk = b(k)
            do j = k+1, n
              bk = bk - a(k,j) * b(j)
            end do
            b(k) = bk
          end do

      end select

    end subroutine slv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mslv (a, b)

      real(wp), intent(in)    :: a(:,:)
      real(wp), intent(inout) :: b(:,:)

      integer  :: m, n, i, j, k
      real(wp) :: bk

      n = size(a,1)  ! Must = size(a,2) = size(b,1)
      m = size(b,2)

      select case (n)

        case (:1)

          ! do nothing.

        case (2)

           b(1,:) = a(1,1) * b(1,:)
           b(2,:) = a(2,2) * (b(2,:) - a(2,1) * b(1,:))
           b(1,:) = b(1,:) - a(1,2) * b(2,:)

        case default

          do i = 1, m

            do k = 1, n
              bk = b(k,i)
              do j = 1, k-1
                bk = bk - a(k,j) * b(j,i)
              end do
              b(k,i) = bk * a(k,k)
            end do

            do k = n-1, 1, -1
              bk = b(k,i)
              do j = k+1, n
                bk = bk - a(k,j) * b(j,i)
              end do
              b(k,i) = bk
            end do

          end do

      end select

    end subroutine mslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  YMAX -- VECTOR MINUS MATRIX-VECTOR PRODUCT.
   !!  CMAB -- MATRIX MINUS MATRIX-MATRIX PRODUCT.
   !!
   !!  THESE ROUTINES COMPUTE THE TERMS Y-A*X AND C-A*B RESPECTIVELY
   !!  WHERE A, B, AND C ARE MATRICES AND X AND Y ARE VECTORS.  THESE
   !!  ARE NOT GENERAL ROUTINES AND ARE INTENDED TO BE USED ONLY BY
   !!  BBFCT AND BBSLV.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ymax (a, x, y)

      real(wp), intent(in)    :: a(:,:), x(:)
      real(wp), intent(inout) :: y(:)

      ! local variables
      integer  :: m, n, j, k
      real(wp) :: yj

      m = size(a,1)  ! Must = size(y)
      n = size(a,2)  ! Must = size(x)

      do j = 1, m
        yj = y(j)
        do k = 1, n
          yj = yj - a(j,k) * x(k)
        end do
        y(j) = yj
      end do

    end subroutine ymax

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cmab (a, b, c)

      real(wp), intent(in)    :: a(:,:), b(:,:)
      real(wp), intent(inout) :: c(:,:)

      ! local variables
      integer  :: m, n, p, i, j, k
      real(wp) :: cij

      m = size(c,1)  ! Must = size(a,1)
      n = size(c,2)  ! Must = size(b,2)
      p = size(a,2)  ! Must = size(b,1)

      do j = 1, n
        do i = 1, m
          cij = c(i,j)
          do k = 1, p
            cij = cij - a(i,k) * b(k,j)
          end do
          c(i,j) = cij
        end do
      end do

    end subroutine cmab

end module bls_solver
