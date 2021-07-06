!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  MFE1 --
!!
!!  This is a Fortran 90 implementation of the piecewise-linear
!!  Gradient-Weighted Moving Finite Element Method (GWMFE) for
!!  time-dependent systems of partial differential equations in
!!  one space dimension.
!!
!!  Version 0.3, 12 August 1997
!!
!!  Neil N. Carlson, Dept of Math, Purdue University
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997, Neil N. Carlson
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program mfe1

  use mfe_constants
  use mfe_types
  use common_io
  use system_io
  use output
  use initialize
  use mfe_ode_solver
  use mfe_procs, only: eval_udot
  implicit none

  real(kind=wp) :: t
  real(kind=wp), dimension(6) :: rvar
  integer :: mode, rtype, j, errc, debug_unit, nstep
  integer, dimension(3) :: ivar
  type(NodeVar), dimension(:), pointer :: u, udot

  character(len=16) :: string

  real :: cpusec, cpusec0

  call cpu_time (cpusec0)

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
  mode = START_SOLN
  call eval_udot (u, t, udot, errc)
  call write_soln (u, t)
  j = 2

  if (errc /= 0) then
    call abort ( (/log_unit,STDERR/), "Bad initial solution.")
  end if

  if (mstep <= 0) then
    stop
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

  write (unit=string, fmt="(es11.3)" ) t
  call info ( (/ log_unit, STDOUT /), "BEGINNING SOLUTION AT T = " // string )

  if (ofreq <= 0) then
    ofreq = mstep
  end if

  do

    call bdf2_solver (mode, rvar, ivar, tout(j), ofreq, rtype, u, udot, t)

    mode = RESUME_SOLN

    call cpu_time (cpusec)
    cpusec = cpusec - cpusec0
    call write_soln (u, t)

    select case (rtype)

      case (SOLN_AT_TOUT)       ! Integrated to TOUT.
        call write_status( (/ log_unit, STDOUT /), cpusec, t )
        j = j + 1

      case (SOLN_AT_STEP)       ! Integrated OFREQ more steps.
        call write_status( (/ log_unit, STDOUT /), cpusec )

      case (FAIL_ON_STEP)
        call write_status( (/ log_unit, STDOUT /), cpusec )
        call abort ( (/ log_unit, STDERR /), "Repeated failure at a step.")

      case (SMALL_H_FAIL)
        call write_status( (/ log_unit, STDOUT /), cpusec )
        call abort ( (/ log_unit, STDERR /), "Next time step is too small.")

      case (FAIL_ON_START)
        call write_status( (/ log_unit, STDOUT /), cpusec )
        call abort ( (/ log_unit, STDERR /), "Starting procedure failed.")

      case default
        call write_status( (/ log_unit, STDOUT /), cpusec )
        call abort ( (/ log_unit, STDERR /), "Unknown return type!")

    end select

    if (j > size(tout)) then
      call info( (/ log_unit, STDOUT /), "Integrated to final TOUT.  Done." )
      exit
    end if

    call bdf2_inquire( nstep=nstep )
    if (nstep >= mstep) then
      call info( (/ log_unit, STDOUT /), "Maximum number of steps taken.  Done." )
      exit
    end if

  end do

end program mfe1
