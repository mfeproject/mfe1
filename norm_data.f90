!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module norm_data

  use prec
  !use mfe_extents, only: nepn
  use mfe_extents

  private

  real(kind=wp), dimension(nepn), save, public  :: ptol
  real(kind=wp), save, public  :: rtol
  real(kind=wp), dimension(:), allocatable, public :: dx

end module norm_data
