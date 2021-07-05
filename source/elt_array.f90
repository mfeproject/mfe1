!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module element_arrays

  use precision

  ! Local solution arrays.
  real(wp), dimension(:,:), allocatable   :: u0, u1
  real(wp), dimension(:,:), allocatable   :: u0dot, u1dot

  ! Local result arrays.
  real(wp), dimension(:,:), allocatable   :: r0, r1
  real(wp), dimension(:,:,:), allocatable :: a00, a10, a01, a11

  ! Intermediate data arrays.
  real(wp), dimension(:,:), allocatable   :: du
  real(wp), dimension(:,:), allocatable   :: l, n1, n2

end module element_arrays
