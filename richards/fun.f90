module richards_functions

  use kind_parameters
  private
  
  public :: conductivity, pressure, pressure_deriv, saturation, cond_theta, diffusivity
  
  real(kind=r8), parameter, private :: THETA_R = 0.075_r8
  real(kind=r8), parameter, private :: THETA_S = 0.287_r8
  real(kind=r8), parameter, private :: ALPHA   = 1.611e6_r8
  real(kind=r8), parameter, private :: BETA    = 3.96_r8
  
  real(kind=r8), parameter, private :: K_S     = 0.00944_r8
  real(kind=r8), parameter, private :: A       = 1.175e6_r8
  real(kind=r8), parameter, private :: GAMMA   = 4.74_r8
  
  contains
  
    function cond_theta (theta) result (value)
    
      real(kind=r8), intent(in) :: theta
      real(kind=r8) :: value
      real(kind=r8) :: p
      
      p = - (ALPHA * (THETA_S - theta) / (theta - THETA_R)) ** (1.0_r8 / BETA)
      
      value = (K_S * A) / (A + abs(p) ** GAMMA)
      
    end function cond_theta
    
    function diffusivity (theta) result (value)
    
      real(kind=r8), intent(in) :: theta
      real(kind=r8) :: value
      real(kind=r8) :: p
      
      p = - (ALPHA * (THETA_S - theta) / (theta - THETA_R)) ** (1.0_r8 / BETA)
      
      value = (K_S * A * (THETA_S - THETA_R) * ALPHA / BETA) / &
              (((theta-THETA_R) ** 2) * (abs(p) ** (BETA-1.0_r8)) * (A + abs(p) ** GAMMA))
                    
    end function diffusivity
    
    function conductivity (p) result (value)
    
      real(kind=r8), intent(in) :: p
      real(kind=r8) :: value
      
      value = (K_S * A) / (A + abs(p) ** GAMMA)
      
    end function conductivity
    
    function pressure (s) result (value)
    
      real(kind=r8), intent(in) :: s
      real(kind=r8) :: value
      
      value = - (ALPHA * (THETA_S - s) / (s - THETA_R)) ** (1.0_r8 / BETA)
      
    end function pressure
    
    function pressure_deriv (s) result (value)
    
      real(kind=r8), intent(in) :: s
      real(kind=r8) :: value
      
      value = (ALPHA ** (1.0_r8 / BETA) * (THETA_S - THETA_R) / BETA) * &
              ((THETA_S - s) / (s - THETA_R)) ** ((1.0_r8 - BETA)/ BETA) / &
              (s - THETA_R) ** 2
      
    end function pressure_deriv
      
    function saturation (p) result (value)
    
      real(kind=r8), intent(in) :: p
      real(kind=r8) :: value
      
      value = ALPHA * (THETA_S - THETA_R) / (ALPHA + abs(p)**BETA) + THETA_R
      
    end function saturation
      
end module richards_functions

module write_functions

  use kind_parameters
  use richards_functions
  private
  
  public :: write_conductivity, write_diffusivity, write_saturation
  
  contains
  
    subroutine write_conductivity (s0, s1, ds, filename)
    
      real(kind=r8), intent(in) :: s0, s1, ds
      character(len=*), intent(in) :: filename
      
      integer :: j
      real(kind=r8) :: s
      
      open(unit=10, file=filename, action="write", position="rewind", status="replace")
        
      j = 0
      
      do
        s = s0 + j * ds
        if (s > s1) then
          exit
        end if
        write(unit=10, fmt=*) s, cond_theta(s)
        j = j + 1
      end do
      
      close(unit=10)
      
    end subroutine write_conductivity
    
    subroutine write_diffusivity (s0, s1, ds, filename)
    
      real(kind=r8), intent(in) :: s0, s1, ds
      character(len=*), intent(in) :: filename
      
      integer :: j
      real(kind=r8) :: s
      
      open(unit=10, file=filename, action="write", position="rewind", status="replace")
        
      j = 0
      
      do
        s = s0 + j * ds
        if (s > s1) then
          exit
        end if
        write(unit=10, fmt=*) s, diffusivity(s)
        j = j + 1
      end do
      
      close(unit=10)
      
    end subroutine write_diffusivity
    
    subroutine write_saturation (p0, p1, dp, filename)
    
      real(kind=r8), intent(in) :: p0, p1, dp
      character(len=*), intent(in) :: filename
      
      integer :: j
      real(kind=r8) :: p
      
      open(unit=10, file=filename, action="write", position="rewind", status="replace")
        
      j = 0
      
      do
        p = p0 + j * dp
        if (p > p1) then
          exit
        end if
        write(unit=10, fmt=*) p, saturation(p)
        j = j + 1
      end do
      
      close(unit=10)
      
    end subroutine write_saturation
    
end module write_functions

program plot_em

  use kind_parameters
  use richards_functions
  use write_functions
  
  real(kind=r8), parameter :: s0 = 0.076_r8, s1 = 0.286_r8, ds = 0.001_r8
  
  call write_diffusivity (s0=0.1_r8, s1=0.27_r8, ds=0.001_r8, filename="diff")
  call write_conductivity (s0=0.1_r8, s1=0.27_r8, ds=0.001_r8, filename="cond")
  !call write_saturation (p0=-61.5_r8, p1=-20.7_r8, dp=0.1_r8, filename="sat")
  
end program plot_em
