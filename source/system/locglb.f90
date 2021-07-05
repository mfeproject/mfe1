module local_global

  use mfe_constants, only: wp, neq
  use mfe_types
  use local_arrays
  private

  public :: gather_local_solution, free_local_arrays, assemble_vector, &
            assemble_matrix, assemble_diagonal

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  GATHER_LOCAL_SOLUTION
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine gather_local_solution (global_u, global_udot)

      type(NodeVar), dimension(:), intent(in)           :: global_u
      type(NodeVar), dimension(:), intent(in), optional :: global_udot

      integer :: i

      nelt = size(global_u) - 1

     !!!
     !!! ALLOCATE WORKING STORAGE

      allocate (u(2,nelt), r(2,nelt), mtx(2,2,nelt), &
                l(neq,nelt), n(neq,nelt), dudx(neq,nelt), dx(nelt))

      if (present(global_udot)) then
        allocate (udot(2,nelt))
      end if

     !!!
     !!! GATHER-UP LOCAL COPIES OF THE SOLUTION VECTORS

      do i = 1, nelt

        u(1,i) % x = global_u(i) % x
        u(1,i) % u = global_u(i) % u

        u(2,i) % x = global_u(i+1) % x
        u(2,i) % u = global_u(i+1) % u

      end do

      if (present(global_udot)) then

        do i = 1, nelt

          udot(1,i) % x = global_udot(i) % x
          udot(1,i) % u = global_udot(i) % u

          udot(2,i) % x = global_udot(i+1) % x
          udot(2,i) % u = global_udot(i+1) % u

        end do

      end if

    end subroutine gather_local_solution

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  FREE_LOCAL_ARRAYS
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine free_local_arrays ()

      if (allocated(u)) then
        deallocate (u, r, mtx, l, n, dudx, dx)
      end if

      if (allocated(udot)) then
        deallocate (udot)
      end if

    end subroutine free_local_arrays

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ASSEMBLE_VECTOR
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine assemble_vector (global_r)

      type(NodeVar), dimension(:), intent(out) :: global_r

      integer :: i

      global_r(1) % x = r(1,1) % x
      global_r(1) % u = r(1,1) % u

      do i = 2, nelt
        global_r(i) % x = r(1,i) % x + r(2,i-1) % x
        global_r(i) % u = r(1,i) % u + r(2,i-1) % u
      end do

      global_r(nelt+1) % x = r(2,nelt) % x
      global_r(nelt+1) % u = r(2,nelt) % u

    end subroutine assemble_vector

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ASSEMBLE_MATRIX
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine assemble_matrix (lowr, diag, uppr)

      type(NodeMtx), dimension(:), intent(out) :: lowr, diag, uppr

      integer :: i

      diag(1) = mtx(1,1,1)
      do i = 2, nelt
        diag(i) % xx = mtx(1,1,i) % xx + mtx(2,2,i-1) % xx
        diag(i) % xu = mtx(1,1,i) % xu + mtx(2,2,i-1) % xu
        diag(i) % ux = mtx(1,1,i) % ux + mtx(2,2,i-1) % ux
        diag(i) % uu = mtx(1,1,i) % uu + mtx(2,2,i-1) % uu
      end do
      diag(nelt+1) = mtx(2,2,nelt)

      do i = 1, nelt
        uppr(i) = mtx(1,2,i)
      end do
      uppr(nelt+1) = 0.0_wp

      lowr(1) = 0.0_wp
      do i = 2, nelt + 1
        lowr(i) = mtx(2,1,i-1)
      end do

    end subroutine assemble_matrix

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ASSEMBLE_DIAGONAL
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine assemble_diagonal (diag)

      type(NodeMtx), dimension(:), intent(out) :: diag

      integer :: i

      diag(1) = mtx(1,1,1)
      do i = 2, nelt
        diag(i) % xx = mtx(1,1,i) % xx + mtx(2,2,i-1) % xx
        diag(i) % xu = mtx(1,1,i) % xu + mtx(2,2,i-1) % xu
        diag(i) % ux = mtx(1,1,i) % ux + mtx(2,2,i-1) % ux
        diag(i) % uu = mtx(1,1,i) % uu + mtx(2,2,i-1) % uu
      end do
      diag(nelt+1) = mtx(2,2,nelt)

    end subroutine assemble_diagonal

end module local_global
