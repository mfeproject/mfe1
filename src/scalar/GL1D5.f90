!!
!! Fifth order Gauss-Labatto quadrature on an interval.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module GL1D4

  implicit none
  private

  integer, parameter :: r8 = selected_real_kind(10,50)

  integer, parameter :: NUM_QUAD_PTS = 4

  real(kind=r8), parameter :: a0 = 1.0_r8
  real(kind=r8), parameter :: b0 = 0.0_r8
  real(kind=r8), parameter :: a1 = 0.27639320225002103_r8  ! (5 - Sqrt(5))/10
  real(kind=r8), parameter :: b1 = 0.72360679774997897_r8  ! (5 + Sqrt(5))/10
  real(kind=r8), parameter :: w0 = 1.0_r8 / 12.0_r8
  real(kind=r8), parameter :: w1 = 5.0_r8 / 12.0_r8

  real(kind=r8), dimension(2,NUM_QUAD_PTS), parameter, public :: phi = &
    reshape( shape=(/2,NUM_QUAD_PTS/), source=(/ a0, b0, b0, a0, a1, b1, b1, a1 /) )

  real(kind=r8), dimension(NUM_QUAD_PTS), parameter, public :: wgt = (/ w0, w0, w1, w1 /)

end module GL1D4

!!
!!  USAGE TEMPLATES
!!
!!  Assume X = (/ X1, X2 /) are the cell endpoints and U = (/ U1, U2 /) are the values
!!  of a primary variable u at the cell endpoints.  Suppose F and G are defined as
!!  elemental functions of position x and u, respectively.
!!
!!    !! Quadrature points.
!!    q_pts = matmul(X, phi)
!!
!!    !! Average of f.
!!    f_avg = dot_product( wgt, f(matmul(X, phi)) )
!!
!!    !! Basis-weighted averages of f.
!!    f_phi = matmul( phi, wgt * f(matmul(X, phi)) )
!!
!!    !! Average of g.
!!    g_avg = dot_product( wgt, g(matmul(U,phi)) )
!!
!!    !! Basis-weighted averages of g.
!!    f_phi = matmul( phi, wgt * f(matmul(U,phi)) )
