!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module problem_procs

  use precision

  implicit none
  private

  ! Publically accessible procedures.
  public :: read_prob_data, eltrhs

  ! Problem parameters.
  real(kind=wp), dimension(2), save :: scale
  real(kind=wp), save :: lambda0, eps

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  READ_PROB_DATA
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_prob_data

      use common_io

      call read_tagged_data (scale, "Hole and voltage scale factors")
      call read_tagged_data (lambda0, "Initial diffusion coefficient")
      call read_tagged_data (eps,"Voltage relaxation coefficient")

    end subroutine read_prob_data

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ELTRHS -- Inner products for gasdynamics.
   !!
   !!        u1(1,:) == mass density,
   !!        u2(2,:) == momentum density,
   !!        u3(3,:) == total energy density (internal plus kinetic).
   !!
   !!     Each equation is in conservation law form
   !!
   !!        du/dt = -df/dx + visc* d2u/dx2.
   !!
   !!  FLUX -- Gas dynamics flux f.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine eltrhs (t)

       use mfe_extents
       use mfe_data, only: wpde
       use element_laplacian
       use element_arrays, only: u0, u1, r0, r1, du, l, n1, n2

       real(kind=wp), intent(in) :: t

       ! Local variables
       integer :: j
       real(kind=wp) :: c1, c2, c3, c4, p, pa0, pa1, term1, term2, lambda

       ! 3-point Gaussian quadrature weights and points on [0,1].
       real(kind=wp), dimension(3), parameter :: gw = &
       (/ .277777777777778_wp, .444444444444444_wp, .277777777777778_wp /)
       real(kind=wp), dimension(3), parameter :: gp = &
       (/ .112701665379259_wp, .500000000000000_wp, .887298334620741_wp /)

       if (t < 40.0_wp) then
         lambda = lambda0
       else
         lambda = lambda0 / 10.0_wp**((t-40.0_wp)/100.0_wp)
       end if

       c1 = 0.5_wp * wpde(1) * scale(2)
       c2 = 0.5_wp * wpde(1) * scale(1) * lambda
       c3 = wpde(1) / scale(1)
       c4 = -wpde(2) / (eps * scale(2))

       do j = 1, nelt

         ! Average of p*alphas by 3 pt Gauss quadrature.
           p = exp (scale(1) * (gp(3) * u0(1,j) + gp(1) * u1(1,j)))
         pa0 = gw(1) * gp(3) * p
         pa1 = gw(1) * gp(1) * p

           p = exp (scale(1) * (gp(2) * u0(1,j) + gp(2) * u1(1,j)))
         pa0 = pa0 + gw(2) * gp(2) * p
         pa1 = pa1 + gw(2) * gp(2) * p

           p = exp (scale(1) * (gp(1) * u0(1,j) + gp(3) * u1(1,j)))
         pa0 = pa0 + gw(3) * gp(1) * p
         pa1 = pa1 + gw(3) * gp(3) * p

         term1 = du(1,j) * (c1*du(2,j) + c2*du(1,j)) / du(x,j)
         term2 = du(x,j) * (0.5_wp - pa0)
         r0(x,j) = (term1 + c3 * term2) * n1(1,j) + (c4 * term2) * n1(2,j)
         r0(1,j) = (term1 + c3 * term2) * n2(1,j)
         r0(2,j) = (c4 * term2) * n2(2,j)

         term2 = du(x,j) * (0.5_wp - pa1)
         r1(x,j) = (term1 + c3 * term2) * n1(1,j) + (c4 * term2) * n1(2,j)
         r1(1,j) = (term1 + c3 * term2) * n2(1,j)
         r1(2,j) = (c4 * term2) * n2(2,j)

       end do

       call lapl (1, lambda)
       call lapl (2, 1.0_wp / eps)

     end subroutine eltrhs

end module problem_procs
