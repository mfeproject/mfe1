module problem_data

  use mfe_constants, only: wp
  private
  
  real(kind=wp), dimension(2), save, public :: scf
  real(kind=wp), save, public :: lambda0, eps
  
end module problem_data
  
module problem_init

  use problem_data
  use common_io
  private
  
  public :: read_problem_data
  
  contains
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  READ_PROBLEM_DATA
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    subroutine read_problem_data ()

      call read_tagged_data (scf, "Hole and voltage scale factors")
      call read_tagged_data (lambda0, "Initial diffusion coefficient")
      call read_tagged_data (eps,"Voltage relaxation coefficient")

    end subroutine read_problem_data

end module problem_init

module problem_pde

  use problem_data
  use mfe_constants, only: wp
  use mfe_data, only: eqw
  use local_arrays
  use local_laplacian
  private
  
  public  :: pde_rhs
  
  contains
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  LOAD_PDE_RHS -- Inner products for the drift-diffusion equations
   !!                  of single-carrier semiconductor device simulation.
   !!  
   !!  The hole density is given by p = exp(SCF*U) and the voltage is
   !!  v = VSCF*V.  The equations for U and V are 
   !!
   !!    dU/dt = Grad(U)*(VSCF*Grad(V) + LAMBDA*USCF*Grad(U))
   !!
   !!             + (1 - exp(USCF*U))/USCF + LAMBDA*Lapl(U)
   !!
   !!    dV/dt = (1/EPS)*(Lapl(V) - (1 - exp(USCF*U))/VSCF)
   !!
   !!  LAMBDA and EPS are parameters and USCF and VSCF are scaling factors.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t

      ! Local variables
      integer :: i
      real(kind=wp) :: c1, c2, c3, c4, p, pa1, pa2, term1, term2, lambda

      ! 3-point Gaussian quadrature weights and points on [0,1].
      real(kind=wp), dimension(3), parameter :: gw = &
      (/ 0.277777777777778_wp, 0.444444444444444_wp, 0.277777777777778_wp /)
      real(kind=wp), dimension(3), parameter :: gp = &
      (/ 0.112701665379259_wp, 0.500000000000000_wp, 0.887298334620741_wp /)

      if (t < 40.0_wp) then
        lambda = lambda0
      else
        lambda = lambda0 / 10.0_wp ** ((t - 40.0_wp) / 100.0_wp)
      end if

      c1 = 0.5_wp * eqw(1) * scf(2)
      c2 = 0.5_wp * eqw(1) * scf(1) * lambda
      c3 = eqw(1) / scf(1)
      c4 = - eqw(2) / (eps * scf(2))

      do i = 1, nelt

       !!!
       !!! AVERAGE OF P * ALPHAS BY 3-PT GAUSS QUADRATURE
        
        p = exp (scf(1) * (gp(3) * u(1,i) % u(1) + gp(1) * u(2,i) % u(1)))
        pa1 = gw(1) * gp(3) * p
        pa2 = gw(1) * gp(1) * p

        p = exp (scf(1) * (gp(2) * u(1,i) % u(1) + gp(2) * u(2,i) % u(1)))
        pa1 = pa1 + gw(2) * gp(2) * p
        pa2 = pa2 + gw(2) * gp(2) * p

        p = exp (scf(1) * (gp(1) * u(1,i) % u(1) + gp(3) * u(2,i) % u(1)))
        pa1 = pa1 + gw(3) * gp(1) * p
        pa2 = pa2 + gw(3) * gp(3) * p

       !!!
       !!! LOAD THE INNER PRODUCTS

        term1 = dx(i) * (c1 * dudx(2,i) + c2 * dudx(1,i)) * dudx(1,i)
        term2 = dx(i) * (0.5_wp - pa1)
        r(1,i) % x    = (term1 + c3 * term2) * n(1,i) % x + &
                        (c4 * term2) * n(2,i) % x
        r(1,i) % u(1) = (term1 + c3 * term2) * n(1,i) % u
        r(1,i) % u(2) = (c4 * term2) * n(2,i) % u

        term2 = dx(i) * (0.5_wp - pa2)
        r(2,i) % x    = (term1 + c3 * term2) * n(1,i) % x + &
                        (c4 * term2) * n(2,i) % x
        r(2,i) % u(1) = (term1 + c3 * term2) * n(1,i) % u
        r(2,i) % u(2) = (c4 * term2) * n(2,i) % u

      end do

      call laplacian (eqno=1, coef=lambda)
      call laplacian (eqno=2, coef=1.0_wp / eps)

    end subroutine pde_rhs

end module problem_pde
