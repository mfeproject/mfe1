!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! FP_ACCELERATOR -- Fixed Point Iteration Accelerator
!!
!! This module implements a Krylov subspace type algorithm that accelerates the
!! convergence of the standard fixed point iteration x_{n+1} = x_{n} - f(x_{n}),
!! for the nonlinear system f(x) = 0.  A modified Newton iteration, for example,
!! gives rise to such an iteration.
!!
!! The accelerator must be initialized by calling
!!
!!    fpa_init (f, maxv, vtol)
!!
!!    f     -- Function value vector; it is used ONLY to glean its shape
!!             in order to allocate working storage.
!!    maxv  -- Maximum number of vectors used in the algorithm (>0).
!!    vtol  -- Tolerance for dropping vectors.  Optional; defaults to 0.1.
!!
!! The accelerated correction is computed by
!!
!!    fpa_correction (itr, f)
!!
!!    itr   -- Iteration count.  Necessary initialization is done for ITR = 1,
!!             but otherwise ITR is ignored.
!!    f     -- Vector containing f at the current iterate.  It is overwritten
!!             with the accelerated correction.  Note that the unaccelerated
!!             correction would be simply f itself.
!!
!! At the end of the iteration, fpa_finish may be called (no arguments)
!! to free the working storage used by the module (this is not required).
!! If it is called, then fpa_init will need to be called before calling
!! fpa_correction again.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fp_accelerator

  use prec

!  implicit none
  private

  public :: fpa_init, fpa_correction, fpa_finish

  real(kind=wp), save, private  :: tol
  integer, save, private   :: mvec

  real(kind=wp), dimension(:,:,:), allocatable, private :: v, w
  real(kind=wp), dimension(:,:), allocatable, private :: h
  integer, dimension(:), allocatable, private  :: next, prev

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! FPA_INIT
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine fpa_init (f, maxv, vtol)

      real(kind=wp), dimension(:,:), intent(in) :: f
      integer, intent(in)            :: maxv
      real(kind=wp), intent(in), optional :: vtol

      integer :: n

      mvec = maxv
      n = mvec + 1

      allocate (v(size(f,1),size(f,2),n), w(size(f,1),size(f,2),n))
      allocate (h(n,n), next(n), prev(n))

      if (present(vtol)) then
        tol = vtol
      else
        tol = 0.1_wp
      end if

    end subroutine fpa_init

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! FPA_FINISH
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine fpa_finish ()

      deallocate (v, w, h, next, prev)

    end subroutine fpa_finish

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! FPA_CORRECTION
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine fpa_correction (itr, f)

      integer, intent(in)     :: itr
      real(kind=wp), dimension(:,:), intent(inout) :: f

      ! local variables.
      real(kind=wp)                :: s, hkk, hkj, cj
      real(kind=wp), dimension(mvec+1) :: c
      integer                 :: i, j, k, new, nvec
      integer, save           :: free, first, last

      if (itr == 1) then

       !!!
       !!! FIRST ITERATION

        v(:,:,1) = f      ! Save the (unaccelerated) correction.
        w(:,:,1) = f      ! Save f to compute the difference on the next call.

        first   = 1       ! Initialize the linked list.
        last    = 1
        next(1) = 0
        prev(1) = 0

        free = 2          ! Initialize the free storage linked list.
        do k = 2, mvec
           next(k) = k + 1
        end do
        next(mvec+1) = 0

      else

       !!!
       !!! NEXT DIFFERENCE W

        w(:,:,first) = w(:,:,first) - f
        !s = 1.0 / sqrt (dot_product (w(:,:,first), w(:,:,first))) ! Ackk...
        s = 1_wp / sqrt (sum (w(:,:,first)**2))

        ! Normalize w_1 and apply same factor to v_1.
        v(:,:,first) = s * v(:,:,first)
        w(:,:,first) = s * w(:,:,first)

        k = next(first)   ! Update H.
        do
          if (k == 0) then
            exit
          end if
          !h(first,k) = dot_product (w(:,:,first), w(:,:,k)) ! Ackk...
          h(first,k) = sum (w(:,:,first) * w(:,:,k))
          k = next(k)
        end do

       !!!
       !!! CHOLESKI FACTORIZATION OF H

        h(first,first) = 1_wp
        k = next(first)
        nvec = 1

        do

          if (k == 0) then
            exit            ! No more vectors.
          end if

          if (nvec == mvec) then      ! Retain at most MVEC vectors:

            next(last) = free           ! truncate the list and
            free = k                    ! update the free storage list.
            last = prev(k)
            next(last) = 0
            exit

          end if

          hkk = 1_wp                  ! Single stage of Choleski factorization.
          j = first
          do
            if (j == k) then
              exit
            end if
            hkj = h(j,k)
            i = first
            do
              if (i == j) then
                exit
              end if
              hkj = hkj - h(k,i)*h(j,i)
              i = next(i)
            end do
            hkj = hkj / h(j,j)
            hkk = hkk - hkj**2
            h(k,j) = hkj
            j = next(j)
          end do

          if (hkk > tol**2) then

            h(k,k) = sqrt(hkk)
            nvec = nvec + 1

          else    ! The current w nearly lies in the span of the previous set.

            next(prev(k)) = next(k)   ! Drop the current vector,
            if (next(k) == 0) then
               last = prev(k)
            else
               prev(next(k)) = prev(k)
            end if

            next(k) = free            ! update the free storage list,
            free = k

            k = prev(k)               ! and back-up.

          endif

          k = next(k)

        end do

       !!!
       !!! PROJECT F ONTO THE SPAN OF THE W VECTORS.

        j = first
        do                            ! Forward substitution.
          if (j == 0) then
            exit
          end if
          !cj = dot_product (f, w(:,:,j))  ! Ackk...
          cj = sum (f * w(:,:,j))
          i = first
          do
            if (i == j) then
              exit
            end if
            cj = cj - h(j,i)*c(i)
            i = next(i)
          end do
          c(j) = cj / h(j,j)
          j = next(j)
        end do

        j = last
        do                            ! Backward substitution.
          if (j == 0) then
            exit
          end if
          cj = c(j)
          i = last
          do
            if (i == j) then
              exit
            end if
            cj = cj - h(i,j)*c(i)
            i = prev(i)
          end do
          c(j) = cj / h(j,j)
          j = prev(j)
        end do

       !!!
       !!! ACCELERATED CORRECTION

        new = free                    ! Find storage for the new vectors.
        free = next(free)

        w(:,:,new) = f                ! Save f for next call.

        k = first                     ! Compute the next correction,
        do
          f = f - c(k)*w(:,:,k) + c(k)*v(:,:,k)
          k = next(k)
          if (k == 0) then
            exit
          end if
        end do

        v(:,:,new) = f                ! and save it for the next call.

        prev(new) = 0                 ! Prepend the vectors to the list.
        next(new) = first
        prev(first) = new
        first = new

      end if

    end subroutine fpa_correction

end module fp_accelerator
