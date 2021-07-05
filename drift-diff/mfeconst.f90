module mfe_constants

  use kind_parameters, only: r8
  private
  
  integer, parameter, public :: wp = r8   ! Working precision for reals
  integer, parameter, public :: neq = 2   ! Number of PDEs
  
  integer, parameter, public :: nvar = neq + 1

end module mfe_constants
