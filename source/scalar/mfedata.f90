module mfe_data

  use mfe_constants, only: wp
  private

  integer, save, public :: kreg
  real(kind=wp), save, public :: eltvsc, segspr
  real(kind=wp), save, public :: fdinc, dxmin

end module mfe_data
