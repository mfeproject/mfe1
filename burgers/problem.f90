module problem_data

  use mfe_constants, only: wp
  private
  
  real(kind=wp), save, public :: visc, u_scale
  
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

      call read_tagged_data (visc, "PDE viscosity coefficient")
      call read_tagged_data (u_scale, "Dependent variable scale factor")

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
   !!  PDE_RHS -- Inner products for Burgers' equation.
   !!  
   !!        du/dt = - u * du/dx + visc * d2u/dx2.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t

      integer :: i
      real(kind=wp) :: c, term1, term2

      c = - u_scale / 6.0_wp

      do i = 1, nelt

        term1 = c * dx(i) * dudx(i)
        term2 = term1 * (2.0_wp * u(1,i) % u + u(2,i) % u)
        r(1,i) % x = term2 * n(i) % x
        r(1,i) % u = term2 * n(i) % u

        term2 = term1 * (u(1,i) % u + 2.0_wp * u(2,i) % u)
        r(2,i) % x = term2 * n(i) % x
        r(2,i) % u = term2 * n(i) % u

      end do

      call laplacian (coef=visc)

    end subroutine pde_rhs

end module problem_pde
