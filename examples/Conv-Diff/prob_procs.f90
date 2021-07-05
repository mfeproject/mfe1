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
  
  ! Publically accesible procedures.
  public :: read_prob_data, eltrhs
  
  ! Problem parameters.
  real(wp), save :: speed, visc

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  READ_PROB_DATA
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    subroutine read_prob_data 

      use common_io

      call read_tagged_data (speed, "Convection speed")
      call read_tagged_data (visc, "Diffusion coefficient")

    end subroutine read_prob_data

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ELTRHS -- Inner products for convection-diffusion equation.
   !! 
   !!     du/dt = - speed * du/dx + visc * d2u/dx2
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    subroutine eltrhs (t)

      use mfe_extents
      use element_laplacian
      use element_arrays, only: u0, u1, r0, r1, du, l, n1, n2

      real(wp), intent(in) :: t

      ! local variables
      real(wp) c, term
      integer j

      c = -0.5_wp * speed

      do j = 1, nelt

        term = c * du(1,j)

        r0(2,j) = term * n1(1,j)
        r0(1,j) = term * n2(1,j)

        r1(2,j) = term * n1(1,j)
        r1(1,j) = term * n2(1,j)

      end do

      if (visc > 0.0_wp) call lapl (1, visc)

    end subroutine eltrhs

end module problem_procs
