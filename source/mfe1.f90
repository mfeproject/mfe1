!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program mfe1s

  use precision
  use mfe_extents
  use common_io
  use init_procs
  use mfe_ode_solver
  use mfe_procs, only: eval_udot
  use xlfutility, only: mclock

  implicit none

  real(wp) :: t, rvar(6), cpusec
  integer  :: mode, rtype, j, errc, ivar(3), debug_unit
  real(wp), pointer  :: u(:,:), udot(:,:)
  type(soln_profile) :: profile

  ! Open the input and output files.

  input_unit = new_unit ()
  open (input_unit, file="mfein", action="read", status="old")

  log_unit = new_unit ()
  open (log_unit, file="mfelog", action="write", status="replace")

  out_unit = new_unit ()
  open (out_unit, file="mfegrf", action="write", status="replace")

  call read_soln (u, udot)
  call read_data

  t = tout(1)
  call write_soln

  if (mstep <= 0) stop

  call eval_udot (u, t, udot, errc)

  if (errc /= 0) call abort ( (/log_unit,stderr/), "Bad initial solution.")

  if (debug /= 0) then

    debug_unit = new_unit ()
    open (debug_unit, file="bdfout", action="write", status="replace")
    call set_solver_messages (debug, debug_unit)

  endif

  rvar(1) = h
  rvar(2) = hlb
  rvar(3) = hub
  rvar(4) = ntol
  rvar(5) = margin
  rvar(6) = vtol

  ivar(1) = mtry
  ivar(2) = mitr
  ivar(3) = mvec

  write (log_unit, "(//a,es11.3)") "** BEGINNING SOLUTION AT T =", t

  mode = START_SOLN
  if (ofreq <= 0) ofreq = mstep
  j = 2

  do

    call bdf2_solver (mode, rvar, ivar, tout(j), ofreq, rtype, u, udot, t)

    mode = RESUME_SOLN
    profile = get_soln_profile ()
    cpusec = mclock() / 100.0_wp

    call write_log
    call write_soln

    select case (rtype)

      case (SOLN_AT_TOUT)       ! Integrated to TOUT.

        j = j + 1

      case (SOLN_AT_STEP)       ! Integrated OFREQ more steps.

      case (FAIL_ON_STEP)

        call abort ( (/ log_unit, stderr /), "Repeated failure at a step.")

      case (SMALL_H_FAIL)

        call abort ( (/ log_unit, stderr /), "Next time step is too small.")

      case (FAIL_ON_START)

        call abort ( (/ log_unit, stderr /), "Starting procedure failed.")

      case default

        call abort ( (/ log_unit, stderr /), "Unknown return type!")

    end select

    if (j > size(tout)) then

      write (log_unit,"(/a)") "** Integrated to final TOUT.  Done."
      exit

    end if

    if (profile % nstep >= mstep) then

      write (log_unit,"(/a)") "** Maximum number of steps taken.  Done."
      exit

    end if

  end do

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  WRITE_LOG
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine write_log

      write (log_unit,'(/, i5, a, es11.3, a, es10.2, a, f8.2, a,           &
                  &/, 8x, a, i4.3, ":", i3.3, a, i4.3, 4(":", i3.3))')     &
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

    subroutine write_soln

      integer j
      character fmt*16

      write (out_unit, '(a, es13.5, a, i4, a, f8.1)') &
            't=', t, ', nstep=', profile % nstep, ', cpu=', cpusec
      write (fmt,'(a,i2,a)') '(', nepn, 'es17.8)'
      write (out_unit,fmt) (u(x,j), u(1:npde,j), j = 1, nnod)

    end subroutine write_soln

end program
