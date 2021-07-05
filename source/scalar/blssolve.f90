module bls_solver

  use mfe_constants
  use mfe_types
  private

  public  :: factor, solve
  private :: fct, vfct, slv, mslv, vslv, vmslv, cmab, ymax, btfct, btslv, update

  interface factor
    module procedure fct, vfct, btfct
  end interface

  interface solve
    module procedure slv, mslv, vslv, vmslv, btslv
  end interface

  interface update
    module procedure cmab, ymax
  end interface

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! BTFCT (FACTOR) -- Factor a (2x2 block) tridiagonal matrix
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine btfct (l, d, u)

      type(NodeMtx), dimension(:), intent(inout) :: l, d, u

      integer :: i
      real(kind=wp) :: tmp_uu, tmp_xu, tmp_ux, tmp_xx

      ! CALL FACTOR (D(1))
      d(1) % uu = 1.0_wp / d(1) % uu
      d(1) % ux = d(1) % uu * d(1) % ux
      d(1) % xx = 1.0_wp / (d(1) % xx - d(1) % xu * d(1) % ux)

      do i = 2, size(d)

        ! CALL SOLVE (D(I-1), U(I-1))
        tmp_uu = d(i-1) % uu * u(i-1) % uu
        tmp_xu = d(i-1) % xx * (u(i-1) % xu - d(i-1) % xu * tmp_uu)
        u(i-1) % xu = tmp_xu
        u(i-1) % uu = tmp_uu - d(i-1) % ux * tmp_xu

        tmp_ux = d(i-1) % uu * u(i-1) % ux
        tmp_xx = d(i-1) % xx * (u(i-1) % xx - d(i-1) % xu * tmp_ux)
        u(i-1) % xx = tmp_xx
        u(i-1) % ux = tmp_ux - d(i-1) % ux * tmp_xx

        ! CALL UPDATE (D(I), L(I), U(I-1))
        d(i) % xx = d(i) % xx - l(i) % xx * u(i-1) % xx - l(i) % xu * u(i-1) % ux
        d(i) % xu = d(i) % xu - l(i) % xx * u(i-1) % xu - l(i) % xu * u(i-1) % uu
        d(i) % ux = d(i) % ux - l(i) % ux * u(i-1) % xx - l(i) % uu * u(i-1) % ux
        d(i) % uu = d(i) % uu - l(i) % ux * u(i-1) % xu - l(i) % uu * u(i-1) % uu

        ! CALL FACTOR (D(I))
        d(i) % uu = 1.0_wp / d(i) % uu
        d(i) % ux = d(i) % uu * d(i) % ux
        d(i) % xx = 1.0_wp / (d(i) % xx - d(i) % xu * d(i) % ux)

      end do

    end subroutine btfct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! BTSLV (SOLVE) -- Solve a (2x2 block) tridiagonal linear system
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine btslv (l, d, u, b)

      type(NodeMtx), dimension(:), intent(in)    :: l, d, u
      type(NodeVar), dimension(:), intent(inout) :: b

      integer :: i
      real(kind=wp) :: bu, bx

      ! CALL SOLVE (D(1),  B(1))      !!! FORWARD SUBSTITUTION !!!
      bu = d(1) % uu * b(1) % u
      bx = d(1) % xx * (b(1) % x - d(1) % xu * bu)
      b(1) % x = bx
      b(1) % u = bu - d(1) % ux * bx


      do i = 2, size(d)

        ! CALL UPDATE (B(I), L(I), B(I-1))
        b(i) % x = b(i) % x - l(i) % xx * b(i-1) % x - l(i) % xu * b(i-1) % u
        b(i) % u = b(i) % u - l(i) % ux * b(i-1) % x - l(i) % uu * b(i-1) % u

        ! CALL SOLVE (D(I), B(I))
        bu = d(i) % uu * b(i) % u
        bx = d(i) % xx * (b(i) % x - d(i) % xu * bu)
        b(i) % x = bx
        b(i) % u = bu - d(i) % ux * bx

      end do

      do i = size(d) - 1, 1, -1     !!! BACKWARD SUBSTITUTION !!!

        ! CALL UPDATE (B(I), U(I), B(I+1))
        b(i) % x = b(i) % x - u(i) % xx * b(i+1) % x - u(i) % xu * b(i+1) % u
        b(i) % u = b(i) % u - u(i) % ux * b(i+1) % x - u(i) % uu * b(i+1) % u

      end do

    end subroutine btslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! FCT (FACTOR) -- LU-factor a matrix
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine fct (a)

      type(NodeMtx), intent(inout) :: a

      a % uu = 1.0_wp / a % uu
      a % ux = a % uu * a % ux
      a % xx = 1.0_wp / (a % xx - a % xu * a % ux)

    end subroutine fct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! VFCT (FACTOR) -- LU-factor a rank-1 array of matrices
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vfct (a)

      type(NodeMtx), dimension(:), intent(inout) :: a

      integer :: l

      do l = 1, size(a)

        a(l) % uu = 1.0_wp / a(l) % uu
        a(l) % ux = a(l) % uu * a(l) % ux
        a(l) % xx = 1.0_wp / (a(l) % xx - a(l) % xu * a(l) % ux)

      end do

    end subroutine vfct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! SLV (SOLVE) -- Solve a matrix-vector linear system
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine slv (a, b)

      type(NodeMtx), intent(in)    :: a
      type(NodeVar), intent(inout) :: b

      real(kind=wp) :: bu, bx

      bu = a % uu * b % u
      bx = a % xx * (b % x - a % xu * bu)
      b % x = bx
      b % u = bu - a % ux * bx

    end subroutine slv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! MSLV (SOLVE) -- Solve a nodal matrix-matrix linear system
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mslv (a, b)

      type(NodeMtx), intent(in)    :: a
      type(NodeMtx), intent(inout) :: b

      real(kind=wp) :: buu, bxu, bux, bxx

      buu = a % uu * b % uu
      bxu = a % xx * (b % xu - a % xu * buu)
      b % xu = bxu
      b % uu = buu - a % ux * bxu

      bux = a % uu * b % ux
      bxx = a % xx * (b % xx - a % xu * bux)
      b % xx = bxx
      b % ux = bux - a % ux * bxx

    end subroutine mslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! VSLV (SOLVE) -- Solve a rank-1 array of matrix-vector linear systems
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vslv (a, b)

      type(NodeMtx), dimension(:), intent(in)    :: a
      type(NodeVar), dimension(:), intent(inout) :: b

      integer :: l
      real(kind=wp) :: bu, bx

      do l = 1, size(a)

        bu = a(l) % uu * b(l) % u
        bx = a(l) % xx * (b(l) % x - a(l) % xu * bu)
        b(l) % x = bx
        b(l) % u = bu - a(l) % ux * bx

      end do

    end subroutine vslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! VMSLV (SOLVE) -- Solve a rank-1 array of nodal matrix-matrix linear systems
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vmslv (a, b)

      type(NodeMtx), dimension(:), intent(in)    :: a
      type(NodeMtx), dimension(:), intent(inout) :: b

      integer :: l
      real(kind=wp) :: buu, bxu, bux, bxx

      do l = 1, size(a)

        buu = a(l) % uu * b(l) % uu
        bxu = a(l) % xx * (b(l) % xu - a(l) % xu * buu)
        b(l) % xu = bxu
        b(l) % uu = buu - a(l) % ux * bxu

        bux = a(l) % uu * b(l) % ux
        bxx = a(l) % xx * (b(l) % xx - a(l) % xu * bux)
        b(l) % xx = bxx
        b(l) % ux = bux - a(l) % ux * bxx

      end do

    end subroutine vmslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! YMAX (UPDATE) -- Perform a nodal "vector-minus-matrix-vector" update step
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ymax (c, a, b)

      type(NodeMtx), intent(in)    :: a
      type(NodeVar), intent(in)    :: b
      type(NodeVar), intent(inout) :: c

      c % x = c % x - a % xx * b % x - a % xu * b % u
      c % u = c % u - a % ux * b % x - a % uu * b % u

    end subroutine ymax

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! CMAB (UPDATE) -- Perform a nodal "matrix-minus-matrix-matrix" update step
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cmab (c, a, b)

      type(NodeMtx), intent(in)    :: a, b
      type(NodeMtx), intent(inout) :: c

      c % xx = c % xx - a % xx * b % xx - a % xu * b % ux
      c % xu = c % xu - a % xx * b % xu - a % xu * b % uu
      c % ux = c % ux - a % ux * b % xx - a % uu * b % ux
      c % uu = c % uu - a % ux * b % xu - a % uu * b % uu

    end subroutine cmab

end module bls_solver
