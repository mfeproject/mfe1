
module kind_parameters     ! Temporary hack.

  ! This should give 8-byte reals (double precision)
  integer, parameter, public :: r8 = selected_real_kind (10, 50)

end module kind_parameters

