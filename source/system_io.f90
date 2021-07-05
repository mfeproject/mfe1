!  Rudimentary (evolving) system-dependent I/O module, for IBM XLF Version 3.2
!
!  Neil N. Carlson, June 1996

module system_io

implicit none

! IOSTAT values... (there are a great many more!)
integer, parameter :: END_OF_FILE = -1
integer, parameter :: END_OF_INTERNAL_FILE = -2
integer, parameter :: END_OF_RECORD = -4

! Default preconnected units...
integer, parameter :: STDERR = 0
integer, parameter :: STDIN  = 5
integer, parameter :: STDOUT = 6
integer, parameter :: PRECONNECTED_UNITS(3) = (/ STDERR, STDIN, STDOUT /)

! Largest possible unit number.
integer, parameter :: MAX_UNIT_NUMBER = 2147483647

contains

 !!!
 !!!  Find an unconnected unit number.

   function new_unit ()  result (unit)

   integer unit
   logical exists, opened
   integer ios

   do unit = 1, MAX_UNIT_NUMBER

      if (any (unit == PRECONNECTED_UNITS)) cycle

      inquire (unit, exist=exists, opened=opened, iostat=ios)
      if (exists .and. .not. opened .and. ios == 0) return

   end do

   unit = -1

   end function new_unit

end module system_io
