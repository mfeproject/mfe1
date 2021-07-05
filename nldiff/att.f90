module problem_data

  use mfe_constants, only: wp
  private
  
  public :: c_of_u, u_of_c
  
  real(kind=wp), parameter, public :: D0 = 7.3e-8, D1 = 1.17e-7, ni = 7.12e18

  real(kind=wp), save, public :: mu, t_scale, u_scale, u_shift, c_min, c_max
  
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
    
    function c_of_u (u) result (c)
    
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(size(u)) :: c, y
      
      integer :: k
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
      
    end function c_of_u
    
    function u_of_c (c) result (u)
    
      real(kind=wp), dimension(:), intent(in) :: c
      real(kind=wp), dimension(size(c)) :: u
      
      u = ((c/mu) + log(c/mu) - u_shift) / u_scale
      
    end function u_of_c

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
      call read_tagged_data (c_min, "Minimum concentration (for transf)")
      call read_tagged_data (c_max, "Maximum concentration (for transf)")

      c_min = c_min / (2.0_wp * ni)
      c_max = c_max / (2.0_wp * ni)

    end subroutine read_problem_data

end module problem_init

module problem_pde

  use mfe_constants, only: wp
  use problem_data
  use local_arrays
  use local_laplacian
  private
  
  public :: pde_rhs
  
  contains
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  PDE_RHS --
   !!  
   !!  RHS inner products for Kent Smith's arsenic diffusion problem,
   !!  using the transformation
   !!
   !!     USCF*U(x,T) + GAMMA = C(x,t)/MU + log(C(x,t)/MU)
   !!
   !!  where C(x,t) is the scaled concentration of arsenic, USCF is
   !!  the scaling factor for U and the scaled time T = t/TSCF.  The
   !!  equation for U is
   !!
   !!     dU/dT = TSCF*A(C)*Lapl(U) + (TSCF*USCF)*B(C)*|Grad(U)|**2
   !!
   !!  where
   !!
   !!     A(C) = N(C)*(D0 + D1*N(C))/sqrt(1 + C**2)
   !!
   !!              C*A'(C)        A(C)
   !!     B(C) = ---------- + -------------
   !!            (1 + C/MU)   (1 + C/MU)**2
   !!
   !!     N(C) = C + sqrt(1 + C**2)
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t

      real(kind=wp) :: fac
      real(kind=wp), dimension(2,nelt) :: visc
      real(kind=wp), dimension(nelt) :: b_phi_1, b_phi_2
      
      ! 3pt gaussian quadrature parameters
      real(kind=wp), parameter :: w1 = 0.2464717596168727_wp,  &
                                  w2 = 0.2222222222222222_wp,  &
                                  w3 = 0.03130601816090511_wp, &
                                  c1 = 0.8872983346207417_wp,  &
                                  c2 = 0.1127016653792583_wp
      
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

      visc(1,:) = t_scale * a_of_u (u(1,:) % u)
      visc(2,:) = t_scale * a_of_u (u(2,:) % u)

      call laplacian (coef=visc)

    end subroutine pde_rhs
    
    function a_of_u (u) result (a)
    
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(size(u)) :: a
      
      real(kind=wp), dimension(size(u)) :: c, r
      
      c = c_of_u (u)
      r = sqrt(1.0_wp + c**2)
      a = (c + r) * (D0 + D1 * (c + r)) / r
      
    end function a_of_u
    
    function b_of_u (u) result (b)
    
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(size(u)) :: b
      
      real(kind=wp), dimension(size(u)) :: c, r
      
      c = c_of_u (u)
      r = sqrt(1.0_wp + c**2)
      b = ( (c + r) * (D0 + D1 * (c + r)) / (1.0_wp + c/mu) + &
            (D0 + D1 * (c + r)**2 * (2.0_wp * r - c)) / (1.0_wp + c**2) ) &
          / (r * (1.0_wp + c/mu))
      
    end function b_of_u

end module problem_pde
