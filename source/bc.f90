!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bc_data

  use precision
  use mfe_extents, only: nepn

  implicit none
  private

  public :: bcl, bcr, bvl, bvr

  integer, dimension(nepn), save  :: bcl, bcr  ! Left/right boundary conditions.
  real(wp), dimension(nepn), save :: bvl, bvr  ! Left/right boundary values.

end module bc_data


module bc_procs

  use precision
  use mfe_extents, only: nepn, nnod
  use bc_data

  implicit none
  private

  public :: resbc, jacbc

  real(wp), parameter :: BIG = 1.e20_wp	!!! FIX THIS !!!

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! RESBC
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine resbc (u, t, r)

      real(wp), intent(in)    :: t
      real(wp), intent(in)    :: u(:,:)
      real(wp), intent(inout) :: r(:,:)

      integer k

      do k = 1, nepn

        if (bcl(k) == 1) r(k,1)    = r(k,1)    + BIG * (u(k,1)    - bvl(k))
        if (bcr(k) == 1) r(k,nnod) = r(k,nnod) + BIG * (u(k,nnod) - bvr(k))

      end do

    end subroutine resbc

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! JACBC
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine jacbc (u, t, a)

      real(wp), intent(in)    :: t
      real(wp), intent(in)    :: u(:,:)
      real(wp), intent(inout) :: a(:,:,:)

      integer k

      do k = 1, nepn

        if (bcl(k) == 1) a(k,k,1)    = a(k,k,1)    + BIG
        if (bcr(k) == 1) a(k,k,nnod) = a(k,k,nnod) + BIG

      end do

    end subroutine jacbc

end module bc_procs
