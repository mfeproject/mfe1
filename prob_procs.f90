!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module problem_procs

  use prec
  
  !implicit none
  private

  ! Publically accessible procedures.
  public :: read_prob_data, eltrhs
  private :: flux
  
  ! Problem parameters.  
  real(kind=wp), save, public :: visc
  
  contains
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  READ_PROB_DATA
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    subroutine read_prob_data ()

      use common_io

      call read_tagged_data (visc, "PDE viscosity coefficient")

    end subroutine read_prob_data

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ELTRHS -- Inner products for gasdynamics.
   !!  
   !!        u1(1,:) == mass density,
   !!        u2(2,:) == momentum density,
   !!        u3(3,:) == total energy density (internal plus kinetic).
   !!
   !!     Each equation is in conservation law form
   !!
   !!        du/dt = -df/dx + visc* d2u/dx2.
   !!
   !!  FLUX -- Gas dynamics flux f.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine eltrhs (t)

       use mfe_extents
       use mfe_data, only: wpde
       use element_laplacian
       use element_arrays, only: u0, u1, r0, r1, du, l, n1, n2

       real(kind=wp), intent(in) :: t

       ! Local variables
       integer  :: j, k
       real(kind=wp) :: rx0, rx1, term
       real(kind=wp), dimension(3) :: umid, f0, f1, favg
       real(kind=wp), parameter :: c1 = 1.0_wp / 6.0_wp, c2 = 4.0_wp / 6.0_wp

       do j = 1, nelt

         f0 = flux (u0(1:3,j))                       ! Flux at left endpoint.
         f1 = flux (u1(1:3,j))                       ! Flux at right endpoint.
         umid = 0.5_wp * (u0(1:3,j) + u1(1:3,j))     ! Variables at midpoint.
         favg = c1 * (f0 + f1) + c2 * flux (umid)    ! Average flux (Simpson).

         rx0 = 0.0_wp
         rx1 = 0.0_wp

         do k = 1, 3

           term = (wpde(k) * (favg(k) - f0(k)))         
           rx0 = rx0 - term * n1(k,j)
           r0(k,j) = - term * n2(k,j)

           term = (wpde(k) * (f1(k) - favg(k)))
           rx1 = rx1 - term * n1(k,j)
           r1(k,j) = - term * n2(k,j)

         end do

         r0(x,j) = rx0
         r1(x,j) = rx1

       end do

       call lapl (1, visc)
       call lapl (2, visc)
       call lapl (3, visc)

     end subroutine eltrhs


     function flux (u) result (f)

       real(kind=wp), dimension(:), intent(in) :: u
       real(kind=wp), dimension(3) :: f

       real(kind=wp), parameter :: gamma = 1.4_wp,                        &
                                 c1 = 0.5_wp * (3.0_wp - gamma),     &
                                 c2 = gamma - 1.0_wp,                &
                                 c3 = 0.5_wp * (1.0_wp - gamma),     &
                                 c4 = gamma

       ! Mass flux (equation 1).
       f(1) = u(2)

       ! Momentum flux (equation 2).
       f(2) = c1 * (u(2) * (u(2) / u(1))) + c2 * u(3)

       ! Energy flux (equation 3).
       f(3) = (c3 * (u(2) * (u(2) / u(1))) + c4 * u(3)) * (u(2) / u(1))

     end function flux

end module problem_procs
