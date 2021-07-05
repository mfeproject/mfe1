module problem_data

  use mfe_constants, only: wp
  private

  !real(kind=wp), save, public :: visc

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

      !call read_tagged_data (visc, "PDE viscosity coefficient")

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
   !!  PDE_RHS -- Stub procedure for the RHS inner products
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t

      integer :: i

      do i = 1, nelt

        r(1,i) % x = 0.0_wp
        r(1,i) % u = 0.0_wp

        r(2,i) % x = 0.0_wp
        r(2,i) % u = 0.0_wp

      end do
      
      !call laplacian (coef=visc)

    end subroutine pde_rhs

end module problem_pde
