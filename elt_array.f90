!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module element_arrays

  use prec

  private

  ! Local solution arrays.
  real(kind=wp), dimension(:,:), allocatable, public   :: u0, u1
  real(kind=wp), dimension(:,:), allocatable, public   :: u0dot, u1dot

  ! Local result arrays.
  real(kind=wp), dimension(:,:), allocatable, public   :: r0, r1
  real(kind=wp), dimension(:,:,:), allocatable, public :: a00, a10, a01, a11

  ! Intermediate data arrays.
  real(kind=wp), dimension(:,:), allocatable, public   :: du
  real(kind=wp), dimension(:,:), allocatable, public   :: l, n1, n2

end module element_arrays
