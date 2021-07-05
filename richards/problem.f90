module problem_data

  use mfe_constants, only: wp
  private
  
  real(kind=wp), public :: t_scale, x_scale, u_scale
  
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

      call read_tagged_data (t_scale, "Time scale factor")
      call read_tagged_data (x_scale, "Space scale factor")
      call read_tagged_data (u_scale, "U scale factor")

    end subroutine read_problem_data

end module problem_init

module problem_pde

  use mfe_constants, only: wp
  use problem_data
  use local_arrays
  use local_laplacian
  private
  
  public  :: pde_rhs
  private :: diffusivity, conductivity
  
  real(kind=wp), parameter, private :: THETA_R = 0.075_wp
  real(kind=wp), parameter, private :: THETA_S = 0.287_wp
  real(kind=wp), parameter, private :: ALPHA   = 1.611e6_wp
  real(kind=wp), parameter, private :: BETA    = 3.96_wp
  
  real(kind=wp), parameter, private :: K_S     = 0.00944_wp
  real(kind=wp), parameter, private :: A       = 1.175e6_wp
  real(kind=wp), parameter, private :: GAMMA   = 4.74_wp
  
  contains
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  PDE_RHS --
   !!  
   !!  RHS inner products for Richards' equation
   !!
   !!     ds/dt = d/dx (D(s) ds/dx + K(s)).
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t

      real(kind=wp), dimension(2,nelt) :: visc
      
      real(kind=wp), parameter ::      & ! 3PT GAUSSIAN QUADRATURE PARAMETERS
          w1 = 5.0_wp / 18.0_wp,       & ! (5 / 18)
          w2 = 4.0_wp / 9.0_wp,        & ! (4 / 9)
          c1 = 0.8872983346207417_wp,  & ! (Sqrt[5] + Sqrt[3]) / (2 Sqrt[5])
          c2 = 0.1127016653792583_wp     ! (Sqrt[5] - Sqrt[3]) / (2 Sqrt[5])
          
      integer :: j
      real(kind=wp) :: u1, u2, u3, f_avg, f_1, f_2
      
      do j = 1, nelt
      
        u1 = c1 * u(1,j) % u + c2 * u(2,j) % u
        u2 = 0.5_wp * (u(1,j) % u + u(2,j) % u)
        u3 = c2 * u(1,j) % u + c1 * u(2,j) % u
        
        f_avg = w1 * (diffusivity(u1) * dudx(j) - conductivity(u1)) + &
                w2 * (diffusivity(u2) * dudx(j) - conductivity(u2)) + &
                w1 * (diffusivity(u3) * dudx(j) - conductivity(u3))
                
        f_1 = diffusivity(u(1,j) % u) * dudx(j) - conductivity(u(1,j) % u)
        f_2 = diffusivity(u(2,j) % u) * dudx(j) - conductivity(u(2,j) % u)
      
        ! Diffusion coefficients at the cell endpoints.
        visc(1,j) = diffusivity(u(1,j) % u)
        visc(2,j) = diffusivity(u(2,j) % u)

        r(1,j) % x = (f_avg - f_1) * n(j) % x
        r(1,j) % u = (f_avg - f_1) * n(j) % u
        
        r(2,j) % x = (f_2 - f_avg) * n(j) % x
        r(2,j) % u = (f_2 - f_avg) * n(j) % u
        
      end do

      call laplacian (coef=visc)

    end subroutine pde_rhs
    
    function conductivity (theta) result (value)
    
      real(kind=wp), intent(in) :: theta
      real(kind=wp) :: value
      real(kind=wp) :: p, uu
      
      uu = u_scale * theta
      p = - (ALPHA * (THETA_S - uu) / (uu - THETA_R)) ** (1.0_wp / BETA)
      
      value = (t_scale/ (x_scale * u_scale )) * (K_S * A) / (A + abs(p) ** GAMMA)
      
    end function conductivity
    
    function diffusivity (theta) result (value)
    
      real(kind=wp), intent(in) :: theta
      real(kind=wp) :: value
      real(kind=wp) :: p, uu
      
      uu = u_scale * theta
      
      p = - (ALPHA * (THETA_S - uu) / (uu - THETA_R)) ** (1.0_wp / BETA)
      
      value = (t_scale / x_scale**2) * (K_S * A * (THETA_S - THETA_R) * ALPHA / BETA) / &
              (((uu-THETA_R) ** 2) * (abs(p) ** (BETA-1.0_wp)) * (A + abs(p) ** GAMMA))
                    
    end function diffusivity
    
end module problem_pde
