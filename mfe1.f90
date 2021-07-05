!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module output

  use prec
  use mfe_extents
  use common_io
  use mfe_ode_solver, only: soln_profile
  !private

  public :: write_log, write_soln
  
  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  WRITE_LOG
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine write_log (profile, t, cpusec)

      type (soln_profile), intent(in) :: profile
      real(kind=wp), intent(in) :: t
      real(kind=wp), intent(in) :: cpusec
      
      write (unit=log_unit,fmt="(/, i5, a, es11.3, a, es10.2, a, f8.2, a,/, 8x, a, i4.3, ':', i3.3, a, i4.3, 4(':', i3.3))")     &
                                         profile % nstep,                  &
                              ":  T =",  t,                                &
                              ",  H =",  profile % hlast,                  &
                            ",  CPU =",  cpusec, "(SEC)",                  &
                            "NRE:NJE=",  profile % nre,                    &
                                         profile % nje,                    &
            ",  NPF:NJF:NNR:NNF:NPEF=",  profile % npf,                    &
                                         profile % njf,                    &
                                         profile % nnr,                    &
                                         profile % nnf,                    &
                                         profile % npef

    end subroutine write_log

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  WRITE_SOLN
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine write_soln (u, profile, t, cpusec)

      real(kind=wp), dimension(:,:), intent(in) :: u
      type (soln_profile), intent(in) :: profile
      real(kind=wp), intent(in) :: t
      real(kind=wp), intent(in) :: cpusec
      
      integer :: j
      character(len=16) :: fmt

      write (unit=out_unit, fmt="(a, es13.5, a, i4, a, f8.1)") &
            "t=", t, ", nstep=", profile % nstep, ", cpu=", cpusec
      write (unit=fmt,fmt="(a,i2,a)") "(", nepn, "es17.8)"
      write (unit=out_unit,fmt=fmt) (u(x,j), u(1:npde,j), j = 1, nnod)

    end subroutine write_soln

end module output


program mfe1s

  use prec
  use mfe_extents
  use common_io
  use system_io
  use output
  use init_procs
  use mfe_ode_solver
  use mfe_procs, only: eval_udot
  !use xlfutility, only: mclock

  !implicit none

  real(kind=wp) :: t, cpusec
  real(kind=wp), dimension(6) :: rvar
  integer :: mode, rtype, j, errc, debug_unit
  integer, dimension(3) :: ivar
  real(kind=wp), dimension(:,:), pointer  :: u, udot
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
  call write_soln (u, profile, t, cpusec)

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
    !cpusec = mclock() / 100.0_wp

    call write_log (profile, t, cpusec)
    call write_soln (u, profile, t, cpusec)

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

end program mfe1s
