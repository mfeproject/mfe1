!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bc_data

  use prec
  use mfe_extents
  !use mfe_extents, only: nepn

  !implicit none
  private

  integer, dimension(nepn), save, public  :: bcl, bcr  ! Left/right boundary conditions.
  real(kind=wp), dimension(nepn), save, public :: bvl, bvr  ! Left/right boundary values.

end module bc_data


module bc_procs

  use prec
  use mfe_extents
  !use mfe_extents, only: nepn, nnod
  use bc_data

  !implicit none
  private

  public :: resbc, jacbc

  real(kind=wp), parameter, private :: BIG = 1.0e20_wp        !!! FIX THIS !!!

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! RESBC
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine resbc (u, t, r)

      real(kind=wp), intent(in)    :: t
      real(kind=wp), dimension(:,:), intent(in)    :: u
      real(kind=wp), dimension(:,:), intent(inout) :: r

      integer :: k

      do k = 1, nepn

        if (bcl(k) == 1) then
          r(k,1)    = r(k,1)    + BIG * (u(k,1)    - bvl(k))
        end if
        if (bcr(k) == 1) then
          r(k,nnod) = r(k,nnod) + BIG * (u(k,nnod) - bvr(k))
        end if

      end do

    end subroutine resbc

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! JACBC
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine jacbc (u, t, a)

      real(kind=wp), intent(in)    :: t
      real(kind=wp), dimension(:,:), intent(in)    :: u
      real(kind=wp), dimension(:,:,:), intent(inout) :: a

      integer :: k

      do k = 1, nepn

        if (bcl(k) == 1) then
          a(k,k,1)    = a(k,k,1)    + BIG
        end if
        if (bcr(k) == 1) then
          a(k,k,nnod) = a(k,k,nnod) + BIG
        end if

      end do

    end subroutine jacbc

end module bc_procs
