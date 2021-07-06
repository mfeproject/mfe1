!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module problem_data

  use mfe_constants
  implicit none
  private

  !real(kind=wp), save, public :: visc

end module problem_data

module problem_init

  use problem_data
  use common_io
  implicit none
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

  use problem_data
  use mfe_constants
  use mfe_data, only: eqw
  use local_arrays
  use local_laplacian
  implicit none
  private

  public  :: pde_rhs

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  PDE_RHS -- Stub routine for the RHS inner products of a PDE system.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t

      integer :: i, k

      do i = 1, nelt

        do k = 1, NEQNS
          r(1,i) % u(k) = 0.0_wp
          r(2,i) % u(k) = 0.0_wp
        end do

        r(1,i) % x = 0.0_wp
        r(2,i) % x = 0.0_wp

      end do

      do k = 1, NEQNS
       !call laplacian (eqno=k, coef=visc)
      end do

    end subroutine pde_rhs

end module problem_pde
