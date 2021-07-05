!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_data

  use prec
  use mfe_extents
  !use mfe_extents, only: npde

  private

  real(kind=wp), dimension(npde), save, public :: wpde    ! Relative PDE weights.
  integer, public, save  :: kreg          ! Switch for segment viscosity type.
  real(kind=wp), dimension(npde), save, public :: segvsc  ! Internodal viscosity coefficients.
  real(kind=wp), dimension(npde), save, public :: segspr  ! Internodal spring coefficients.
  real(kind=wp), save, public :: fdinc         ! Increment for FD calc of Jacobian.
  real(kind=wp), save, public :: dxmin         ! Minimum allowed element length.

end module mfe_data
