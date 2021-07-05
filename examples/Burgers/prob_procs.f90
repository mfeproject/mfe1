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
  real(wp), save :: visc, scale

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  READ_PROB_DATA
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_prob_data

      use common_io

      call read_tagged_data (visc, "PDE viscosity coefficient")
      call read_tagged_data (scale, "Dependent variable scale factor")

    end subroutine read_prob_data

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ELTRHS -- Inner products for Burgers' equation.
   !!
   !!        du/dt = - u * du/dx + visc* d2u/dx2.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine eltrhs (t)

       use mfe_extents
       use element_laplacian
       use element_arrays, only: u0, u1, r0, r1, du, l, n1, n2

       real(wp), intent(in) :: t

       ! Local variables
       integer  :: j
       real(wp) :: c, term

       c = -scale / 6.0_wp

       do j = 1, nelt

         term = ((c * du(1,j)) * (2.0_wp * u0(1,j) + u1(1,j)))
         r0(x,j) = term * n1(1,j)
         r0(1,j) = term * n2(1,j)

         term = ((c * du(1,j)) * (u0(1,j) + 2.0_wp * u1(1,j)))
         r1(x,j) = term * n1(1,j)
         r1(1,j) = term * n2(1,j)

       end do

       call lapl (1, visc)

     end subroutine eltrhs

end module problem_procs
