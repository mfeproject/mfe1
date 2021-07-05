module problem_data

  use mfe_constants, only: wp
  private
  
  ! Change of variable functions
  public  :: c_of_u, u_of_c
  private :: c_of_u_scalar, c_of_u_vector, u_of_c_scalar, u_of_c_vector
  
  real(kind=wp), save, public :: mu, u_scale = 1.0_wp, u_shift = 0.0_wp, t_scale
  
  interface c_of_u
    module procedure c_of_u_scalar, c_of_u_vector
  end interface
  
  interface u_of_c
    module procedure u_of_c_scalar, u_of_c_vector
  end interface
  
  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  C_OF_U -- Compute c(u).
   !!
   !!  Using a Newton iteration, this routine inverts the function
   !!
   !!	U(C) = (C/MU + log(C/MU) - U_SHIFT) / U_SCALE.
   !!
   !!  This is done by writing
   !!
   !!	U_SCALE * U + U_SHIFT = (C/MU) + log(C/MU)
   !!
   !!  and then inverting the function y = x + log(x).
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    function c_of_u_scalar (u) result (c)
    
      real(kind=wp), intent(in) :: u
      real(kind=wp) :: c
      
      integer :: k
      real(kind=wp) :: y
      integer, parameter :: NITR = 4

      y = u_scale * u + u_shift
      
      if (y < 1.0_wp) then
        c = exp(y - 1.0_wp)
      else
        c = y
      end if
      
      y = y + 1.0_wp
      
      do k = 1, NITR
        c = c * (y - log(c)) / (1.0_wp + c)
      end do
      
      c = mu * c
      
    end function c_of_u_scalar
    
    function c_of_u_vector (u) result (c)
    
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(size(u)) :: c
      
      integer :: k
      real(kind=wp), dimension(size(u)) :: y
      integer, parameter :: NITR = 4

      y = u_scale * u + u_shift
      
      where (y < 1.0_wp)
        c = exp(y - 1.0_wp)
      elsewhere
        c = y
      end where
      
      y = y + 1.0_wp
      
      do k = 1, NITR
        c = c * (y - log(c)) / (1.0_wp + c)
      end do
      
      c = mu * c
      
    end function c_of_u_vector
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  U_OF_C -- Compute u(c).
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    function u_of_c_scalar (c) result (u)
    
      real(kind=wp), intent(in) :: c
      real(kind=wp) :: u
      
      u = ((c/mu) + log(c/mu) - u_shift) / u_scale
      
    end function u_of_c_scalar

    function u_of_c_vector (c) result (u)
    
      real(kind=wp), dimension(:), intent(in) :: c
      real(kind=wp), dimension(size(c)) :: u
      
      u = ((c/mu) + log(c/mu) - u_shift) / u_scale
      
    end function u_of_c_vector

end module problem_data

module problem_init

  use mfe_constants, only: wp
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

      call read_tagged_data (mu, "Transformation parameter MU")
      call read_tagged_data (t_scale, "Time scale factor")

    end subroutine read_problem_data

end module problem_init

module problem_pde

  use mfe_constants, only: wp
  use problem_data
  use local_arrays
  use local_laplacian
  private
  
  public  :: pde_rhs
  private :: b_of_u
  
  contains
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  PDE_RHS --
   !!  
   !!  RHS inner products for the (transformed) nonlinear diffusion problem
   !!
   !!     dc/dt = Div((1 + c) Grad(c)).
   !!
   !!  The equation is transformed via the change of variable
   !!
   !!     U_SCALE * u + U_SHIFT = (c / MU) + log(c / MU)
   !!
   !!  With a scaled time T = t / T_SCALE the final equation is
   !!
   !!     du/dT = T_SCALE * a(u) * Lapl(u) + (T_SCALE * U_SCALE) * b(u) * |Grad(u)|**2
   !!
   !!  where
   !!
   !!     a(u) = 1 + c(u),   b(u) = MU + (1 - MU) / (1 + c(u)/MU)**2.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t

      real(kind=wp) :: fac
      real(kind=wp), dimension(2,nelt) :: visc
      real(kind=wp), dimension(nelt) :: b_phi_1, b_phi_2
      
      real(kind=wp), parameter ::      & ! 3PT GAUSSIAN QUADRATURE PARAMETERS
          w1 = 0.2464717596168727_wp,  & ! (5 / 18) * (Sqrt[5] + Sqrt[3]) / (2 Sqrt[5])
          w2 = 0.2222222222222222_wp,  & ! (4 / 9)  * (1 / 2)
          w3 = 0.03130601816090511_wp, & ! (5 / 18) * (Sqrt[5] - Sqrt[3]) / (2 Sqrt[5])
          c1 = 0.8872983346207417_wp,  & ! (Sqrt[5] + Sqrt[3]) / (2 Sqrt[5])
          c2 = 0.1127016653792583_wp     ! (Sqrt[5] - Sqrt[3]) / (2 Sqrt[5])
      
      fac = u_scale * t_scale

      b_phi_1 = w1 * b_of_u (c1 * u(1,:) % u + c2 * u(2,:) % u)  + &
                w2 * b_of_u (0.5_wp * (u(1,:) % u + u(2,:) % u)) + &
                w3 * b_of_u (c2 * u(1,:) % u + c1 * u(2,:) % u)

      b_phi_2 = w3 * b_of_u (c1 * u(1,:) % u + c2 * u(2,:) % u)  + &
                w2 * b_of_u (0.5_wp * (u(1,:) % u + u(2,:) % u)) + &
                w1 * b_of_u (c2 * u(1,:) % u + c1 * u(2,:) % u)

      r(1,:) % x = ((fac * dx * dudx**2) * b_phi_1) * n % x
      r(1,:) % u = ((fac * dx * dudx**2) * b_phi_1) * n % u

      r(2,:) % x = ((fac * dx * dudx**2) * b_phi_2) * n % x
      r(2,:) % u = ((fac * dx * dudx**2) * b_phi_2) * n % u

      ! Diffusion coefficients at the cell endpoints.
      visc(1,:) = t_scale * (1.0_wp + c_of_u (u(1,:) % u))
      visc(2,:) = t_scale * (1.0_wp + c_of_u (u(2,:) % u))

      call laplacian (coef=visc)

    end subroutine pde_rhs
    
    function b_of_u (u) result (b)
    
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(size(u)) :: b
      
      b = mu + (1.0_wp - mu) / (1.0_wp + c_of_u (u) / mu)**2
      
    end function b_of_u

end module problem_pde
