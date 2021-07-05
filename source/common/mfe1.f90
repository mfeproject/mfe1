
program mfe1

  use mfe_constants
  use mfe_types
  use common_io
  use system_io
  use output
  use initialize
  use mfe_ode_solver
  use mfe_procs, only: eval_udot
  !use xlfutility, only: mclock

  real(kind=wp) :: t
  real(kind=wp), dimension(6) :: rvar
  integer :: mode, rtype, j, errc, debug_unit
  integer, dimension(3) :: ivar
  type(NodeVar), dimension(:), pointer  :: u, udot
  type(soln_profile) :: profile

  ! Open the input and output files.

  call new_unit (input_unit)
  open (unit=input_unit, file="mfein", position="rewind", action="read", status="old")

  call new_unit (log_unit)
  open (unit=log_unit, file="mfelog", position="rewind", action="write", status="replace")

  call new_unit (out_unit)
  open (unit=out_unit, file="mfegrf", position="rewind", action="write", status="replace")

  call read_soln (u, udot)
  call read_data ()

  t = tout(1)
  call write_soln (u, t)

  if (mstep <= 0) then
    stop
  end if

  call eval_udot (u, t, udot, errc)

  if (errc /= 0) then
    call abort ( (/log_unit,STDERR/), "Bad initial solution.")
  end if

  if (debug /= 0) then

    call new_unit (debug_unit)
    open (unit=debug_unit, file="bdfout", action="write", status="replace")
    call set_solver_messages (debug, debug_unit)

  end if

  rvar(1) = h
  rvar(2) = hlb
  rvar(3) = hub
  rvar(4) = ntol
  rvar(5) = margin
  rvar(6) = vtol

  ivar(1) = mtry
  ivar(2) = mitr
  ivar(3) = mvec

  write (unit=log_unit, fmt="(//a,es11.3)") "** BEGINNING SOLUTION AT T =", t

  mode = START_SOLN
  if (ofreq <= 0) then
    ofreq = mstep
  end if

  j = 2

  do

    call bdf2_solver (mode, rvar, ivar, tout(j), ofreq, rtype, u, udot, t)

    mode = RESUME_SOLN
    profile = get_soln_profile ()

    write (unit=log_unit, fmt="(/a,es11.3)") "** SOLUTION OUTPUT AT T =", t
    call write_soln_profile ( (/log_unit/) )
    call write_soln (u, t)

    select case (rtype)

      case (SOLN_AT_TOUT)       ! Integrated to TOUT.

        j = j + 1

      case (SOLN_AT_STEP)       ! Integrated OFREQ more steps.

      case (FAIL_ON_STEP)

        call abort ( (/ log_unit, STDERR /), "Repeated failure at a step.")

      case (SMALL_H_FAIL)

        call abort ( (/ log_unit, STDERR /), "Next time step is too small.")

      case (FAIL_ON_START)

        call abort ( (/ log_unit, STDERR /), "Starting procedure failed.")

      case default

        call abort ( (/ log_unit, STDERR /), "Unknown return type!")

    end select

    if (j > size(tout)) then

      write (unit=log_unit,fmt="(/a)") "** Integrated to final TOUT.  Done."
      exit

    end if

    if (profile % nstep >= mstep) then

      write (unit=log_unit,fmt="(/a)") "** Maximum number of steps taken.  Done."
      exit

    end if

  end do

end program mfe1
