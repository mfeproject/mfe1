
module output

  use mfe_constants
  use mfe_types
  use common_io
  private

  public :: write_soln

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  WRITE_SOLN
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine write_soln (u, t)

      type(NodeVar), dimension(:), intent(in) :: u
      real(kind=wp), intent(in) :: t

      integer :: j
      character(len=16) :: fmt

      write (unit=out_unit, fmt="(a,es13.5)") "TIME = ", t
      write (unit=fmt,fmt="(a,i2,a)") "(", nvar, "es17.8)"
      write (unit=out_unit,fmt=fmt) (u(j) % x, u(j) % u, j = 1, size(u))

    end subroutine write_soln

end module output
