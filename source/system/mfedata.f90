module mfe_data

  use mfe_constants, only: wp, neq
  private

  integer, save, public :: kreg
  real(kind=wp), dimension(neq), save, public :: eqw, eltvsc, segspr
  real(kind=wp), save, public :: fdinc, dxmin

end module mfe_data
