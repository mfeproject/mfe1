!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1999 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module problem_data

  use mfe_constants, only: wp
  implicit none
  private
  
  !! Universal constants
  real(kind=wp), parameter, private :: F = 96484.6_wp     ! Faraday constant, C/mol
  real(kind=wp), parameter, private :: R = 8.3143_wp      ! Ideal gas constant, J/(mol-K)
  
  !! Material constants
  real(kind=wp), parameter, private :: e = 7.0832e-10_wp  ! Dielectric constant of H20, C/(volt-m)
  real(kind=wp), parameter, private :: T = 300.0_wp       ! Temperature, K
  
  !! Hardwired scaling factors
  real(kind=wp), parameter, private :: scf_x = 0.01_wp    ! length units are cm
  real(kind=wp), parameter, private :: scf_t = 3600.0_wp  ! time units are hours
  real(kind=wp), parameter, private :: scf_c = 1.0_wp     ! Don't scale concentration
  
  !! Material parameters for the species SO_4(2-), Mg(2+), Na(+), K(+):
  integer, parameter, public :: NSPECIES = 4
  real(kind=wp), dimension(NSPECIES), parameter, public :: &
      D = (/ 3.0_wp, 3.0_wp, 3.0_wp, 5.0_wp /) * 1.0e-10 * scf_t / scf_x**2
  real(kind=wp), dimension(NSPECIES), parameter, public :: &
      valence = (/ -2.0_wp, 2.0_wp, 1.0_wp, 1.0_wp /)
  
  !! Physical value for lambda^2 -- it's TINY!  We'll read in a value instead.
  !! real(kind=wp), parameter, public :: lambda_sq = e * R * T / (F**2 * scf_x**2 * scf_c)
  
  real(kind=wp), save, public :: lambda_sq, rho0, eps

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

      call read_tagged_data (rho0, "Fixed charge density")
      call read_tagged_data (lambda_sq, "lambda^2 coefficient")
      call read_tagged_data (eps,"Voltage relaxation coefficient")

    end subroutine read_problem_data

end module problem_init

module problem_pde

  use problem_data
  use mfe_constants, only: wp
  use mfe_data, only: eqw
  use local_arrays
  use local_laplacian
  implicit none
  private
  
  public  :: pde_rhs
  
  contains
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  LOAD_PDE_RHS -- Inner products for the poisson-nernst-planck equations
   !!  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pde_rhs (t)

      real(kind=wp), intent(in) :: t
      
      ! Local variables
      integer :: i, k
      real(kind=wp) :: c2, c3
      real(kind=wp), dimension(NSPECIES) :: u1, u2, gradu, Lphi1, Lphi2, c1
      real(kind=wp) :: v1, v2, gradv, Mphi1, Mphi2, rho1, rho2

      c1 = eqw(1:NSPECIES) * valence * D
      c2 = 1.0_wp / (12.0_wp * lambda_sq)
      c3 = eqw(NSPECIES+1) / (6.0_wp * eps * lambda_sq)
      
      do i = 1, ncell
      
        u1 = u(1,i) % u(1:NSPECIES)  ! Concentrations at left
        u2 = u(2,i) % u(1:NSPECIES)  ! Concentrations at right
        v1 = u(1,i) % u(NSPECIES+1)  ! Voltage at left
        v2 = u(2,i) % u(NSPECIES+1)  ! Voltage at right
        gradu = dudx(1:NSPECIES,i)   ! Concentration gradients
        gradv = dudx(NSPECIES+1,i)   ! Voltage gradient
        
        rho1 = rho0 + dot_product(valence, u1)  ! Charge density at left
        rho2 = rho0 + dot_product(valence, u2)  ! Charge density at right
        
        Lphi1 = c1 * ( ((0.5_wp*dx(i)*gradv)*gradu) - &
                       (c2*dx(i)) * ((3.0_wp * rho1 + rho2) * u1 + (rho1 + rho2) * u2) )
        Lphi2 = c1 * ( ((0.5_wp*dx(i)*gradv)*gradu) - &
                       (c2*dx(i)) * ((rho1 + rho2) * u1 + (rho1 + 3.0_wp * rho2) * u2) )
                 
        Mphi1 = c3 * dx(i) * (2.0_wp * rho1 + rho2)
        Mphi2 = c3 * dx(i) * (rho1 + 2.0_wp * rho2)

       !!!
       !!! LOAD THE INNER PRODUCTS
       
        r(1,i) % x = dot_product(Lphi1, n(1:NSPECIES,i) % x) + Mphi1 * n(NSPECIES+1,i) % x
        r(1,i) % u(1:NSPECIES) = Lphi1 * n(1:NSPECIES,i) % u
        r(1,i) % u(NSPECIES+1) = Mphi1 * n(NSPECIES+1,i) % u

        r(2,i) % x = dot_product(Lphi2, n(1:NSPECIES,i) % x) + Mphi2 * n(NSPECIES+1,i) % x
        r(2,i) % u(1:NSPECIES) = Lphi2 * n(1:NSPECIES,i) % u
        r(2,i) % u(NSPECIES+1) = Mphi2 * n(NSPECIES+1,i) % u

      end do

      do k = 1, NSPECIES
        call laplacian (eqno=k, coef=D(k))
      end do
      call laplacian (eqno=NSPECIES+1, coef=1.0_wp / eps)

    end subroutine pde_rhs

end module problem_pde
