!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_data

  use precision
  use mfe_extents, only: npde

  private

  real(wp), public, save :: wpde(npde)    ! Relative PDE weights.
  integer, public, save  :: kreg          ! Switch for segment viscosity type.
  real(wp), public, save :: segvsc(npde)  ! Internodal viscosity coefficients.
  real(wp), public, save :: segspr(npde)  ! Internodal spring coefficients.
  real(wp), public, save :: fdinc         ! Increment for FD calc of Jacobian.
  real(wp), public, save :: dxmin         ! Minimum allowed element length.

end module mfe_data
