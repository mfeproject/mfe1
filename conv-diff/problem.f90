module problem_data

  use mfe_constants, only: wp
  private
  
  real(kind=wp), save, public :: speed, visc
  
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

      call read_tagged_data (speed, "Convection speed")
      call read_tagged_data (visc, "PDE viscosity coefficient")

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
   !!  PDE_RHS -- Inner products for a simple convection-diffusion equation.
   !!  
   !!        du/dt = - speed * du/dx + visc * d2u/dx2.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t

      integer :: i
      real(kind=wp) :: c, term, term_x, term_u

      c = - 0.5_wp * speed

      do i = 1, nelt

        term = c * dx(i) * dudx(i)
        term_x = term * n(i) % x
        term_u = term * n(i) % u
        
        r(1,i) % x = term_x
        r(1,i) % u = term_u

        r(2,i) % x = term_x
        r(2,i) % u = term_u

      end do

      call laplacian (coef=visc)

    end subroutine pde_rhs

end module problem_pde
