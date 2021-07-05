module mfe_constants

  use kind_parameters, only: r8
  private

  integer, parameter, public :: wp = r8   ! Working precision for reals

  !!! DO NOT CHANGE THE FOLLOWING PARAMETERS !!!
  integer, parameter, private :: neq = 1
  integer, parameter, public  :: nvar = neq + 1

end module mfe_constants
