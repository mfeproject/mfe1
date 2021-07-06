!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_ode_solver

  use mfe_constants, only: wp
  use mfe_types
  use mfe_procs, only: eval_residual, eval_jacobian
  use norm_procs, only: eval_norm, check_soln
  implicit none
  private

  public  :: bdf2_solver, set_solver_messages, &
             bdf2_write_checkpoint, bdf2_read_checkpoint, bdf2_inquire

  ! Mode codes.
  integer, parameter, public :: START_SOLN  = 0
  integer, parameter, public :: RESUME_SOLN = 1

  ! Return codes.
  integer, parameter, public :: FAIL_ON_STEP = -1
  integer, parameter, public :: SMALL_H_FAIL = -2
  integer, parameter, public :: FAIL_ON_START = -3
  integer, parameter, public :: SOLN_AT_TOUT = 2
  integer, parameter, public :: SOLN_AT_STEP = 3

  ! Aliases.
  logical, parameter :: YES = .true.
  logical, parameter :: NO = .false.

  ! Diagnostic message formats.
  integer, save :: logfile=6, msgs=0
  character(len=*), parameter ::                                  &
        FMT_N = "('N:',i4.4,':',i1,':',e12.6,':',e12.6,':')", &
        FMT_J = "('J:')",                                     &
        FMT_A = "('A:',e12.6,':')",                           &
        FMT_H = "('H:',e12.6,':',e12.6,':')",                 &
       FMT_PF = "('PF:',e12.6,':')",                          &
       FMT_JF = "('JF:',e12.6,':')",                          &
       FMT_NR = "('NR:')",                                    &
       FMT_NF = "('NF:',e12.6,':')",                          &
      FMT_PEF = "('PEF:',e12.6,':',e12.6,':')",               &
       FMT_NI = "('NI:',i1,':',e12.6,':',e12.6,':')",         &
      FMT_NIF = "('NIF:',i1,':',e12.6,':')"

 !!!
 !!! SOLVER PARAMETERS

  real(kind=wp), save :: hlb, hub      ! Upper and lower bounds on the step size
  real(kind=wp), save :: ntol          ! BCE step error tolerance (relative to 1)
  real(kind=wp), save :: margin        ! Jacobian step size slack.
  integer, save       :: mtry          ! Maximum number of attempts at a step
  integer, save       :: mitr          ! Maximum number of BCE step iterations

  real(kind=wp), parameter :: RMIN = 0.25_wp
  real(kind=wp), parameter :: RMAX = 4.0_wp

  ! FP accelerator parameters
  real(kind=wp), save :: vtol    ! Dependent vector tolerance
  integer,       save :: mvec    ! Maximum number of vectors used

 !!!
 !!! SOLVER STATE VARIABLES

  ! Last solution and its backward differences
  type(NodeVar), pointer, dimension(:), save :: u0, u1, u2, u3

  real(kind=wp), save :: hlast, h2, h3 ! Last three time steps

  real(kind=wp), save :: h             ! Next time step
  real(kind=wp), save :: hjac          ! BCE time step built into current Jacobian
  logical,       save :: RenewJacobian ! Update jacobian at the next opportunity
  logical,       save :: StaleJacobian ! Jacobian evaluated at a previous step
  integer,       save :: FreezeCount   ! Don't increase step size for this number of steps

  real(kind=wp), save :: tlast         ! Time at last converged solution step.
  integer,       save :: step          ! Step number

  integer,      save :: ReturnState    ! What triggered the (last) return
  integer, parameter :: RET_ON_TIME = 1
  integer, parameter :: RET_ON_STEP = 2
  integer, parameter :: RET_ON_FAIL = 4

 !!!
 !!! SOLVER PROFILE VARIABLES

  ! These counters are informative only and are not used by the solver.
  ! Their values can be read with the subroutine BDF2_INQUIRE.

  integer, save :: residual_calls      ! Number of calls to EVAL_RESIDUAL
  integer, save :: jacobian_calls      ! Number of calls to EVAL_JACOBIAN
  integer, save :: bad_predictors      ! Number of bad predicted solutions (CHECK_SOLN)
  integer, save :: failed_jacobians    ! Number of failed Jacobian evaluations
  integer, save :: retried_bce_steps   ! Number of retried BCE steps
  integer, save :: failed_bce_steps    ! Number of completely failed BCE steps
  integer, save :: rejected_steps      ! Number of steps rejected on error tolerance

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  BDF2_SOLVER -- Second order backward difference solver.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bdf2_solver (mode, rvar, ivar, tout, rfreq, ReturnType, u, udot, t)

      integer, dimension(:), intent(in)     :: ivar
      integer, intent(in) :: rfreq, mode
      integer, intent(out)    :: ReturnType
      real(kind=wp), dimension(:), intent(in)    :: rvar
      real(kind=wp), intent(in) :: tout
      type(NodeVar), dimension(:), intent(inout) :: u, udot
      real(kind=wp), intent(inout) :: t

      ! local variables.
      logical  :: ReturnCheck
      integer  :: errc, try
      real(kind=wp) :: perr, eta, etah, dt
      type(NodeVar), dimension(:), pointer :: ptr

      select case (mode)

        case (START_SOLN)

          allocate (u0(size(u)), u1(size(u)), u2(size(u)), u3(size(u)))

          mtry = ivar(1)
          mitr = ivar(2)
          mvec = ivar(3)
          h    = rvar(1)
          hlb  = rvar(2)
          hub  = rvar(3)
          ntol = rvar(4)
          margin = rvar(5)
          vtol = rvar(6)

          step = 0
          jacobian_calls = 0
          residual_calls = 0
          bad_predictors = 0
          failed_jacobians = 0
          retried_bce_steps = 0
          failed_bce_steps = 0
          rejected_steps = 0

          call start_bdf2 (u, udot, t, errc)

          if (errc /= 0) then           ! Starting procedure failed.

            t = tlast
            u = u0
            ReturnType = FAIL_ON_START
            ReturnState = RET_ON_FAIL
            return

          end if

          h3 = 3.0_wp * h
          h2 = 2.0_wp * h
          hlast = h
          hjac = (2.0_wp / 3.0_wp) * h
          RenewJacobian = NO
          FreezeCount = 0
          ReturnCheck = YES

        case (RESUME_SOLN)

          if (ReturnState == RET_ON_TIME) then
            ReturnCheck = YES
          else
            ReturnCheck = NO
          end if

        case default

          stop ! ERROR!

      end select



      bdf2_step: do

       !!!
       !!! Return Solution

        if (ReturnCheck) then     ! Check for a normal return before proceding.

          if (tout <= tlast) then                ! Cubic interpolation to TOUT.

            dt = tout - tlast
            t = tout
            u = u0 + dt*(u1 + (dt+hlast)*(u2 + (dt+h2)*u3))
            ReturnType  = SOLN_AT_TOUT
            ReturnState = RET_ON_TIME
            exit

          else if (modulo(step,rfreq) == 0) then   ! Return current solution.

            t = tlast
            u = u0
            ReturnType  = SOLN_AT_STEP
            ReturnState = RET_ON_STEP
            exit

          end if

        end if

        step = step + 1
        try = 0

       !!!
       !!! Atempt a BDF2 Step

        attempt: do

          try = try + 1

          if (try > mtry) then          ! Repeated failure at a single step.

            t = tlast
            u = u0
            ReturnType  = FAIL_ON_STEP
            ReturnState = RET_ON_FAIL
            exit bdf2_step

          end if

          if (h < hlb) then             ! Time step is too small.

            t = tlast
            u = u0
            ReturnType  = SMALL_H_FAIL
            ReturnState = RET_ON_FAIL
            exit bdf2_step

          end if

          if (msgs > 0) then
            write (unit=logfile,fmt=FMT_N) step, try, tlast, h
          end if

          t = tlast + h
          eta = (hlast + h) / (hlast + 2.0_wp * h)
          etah = eta * h
          StaleJacobian = YES
          if (hjac/etah > 1.0_wp + margin) then
            RenewJacobian = YES
          end if
          if (etah/hjac > 1.0_wp + margin) then
            RenewJacobian = YES
          end if

         !!!
         !!! Predicted Solution and Backward Difference

          u = u0 + h * (u1 + (h + hlast) * u2)
          udot = u1 + ((h + hlast) / eta) * u2

          call check_soln (u, 0, errc)  ! Check the predicted solution.

          if (errc /= 0) then             ! It's bad; cut h and retry.

            h = 0.5_wp * h
            FreezeCount = 1
            bad_predictors = bad_predictors + 1
            if (msgs > 0) then
              write (unit=logfile,fmt=FMT_PF) h
            end if
            cycle attempt

          end if

          bce: do

           !!!
           !!! Jacobian Evaluation

            if (RenewJacobian) then     ! Reevaluate the Jacobian.

              jacobian_calls = jacobian_calls + 1
              if (msgs > 0) then
                write (unit=logfile,fmt=FMT_J)
              end if

              call eval_jacobian (u, udot, t, etah, errc)

              if (errc /= 0) then         ! Evaluation failed; cut h and retry.

                h = 0.5_wp * h
                RenewJacobian = YES
                failed_jacobians = failed_jacobians + 1
                if (msgs > 0) then
                  write (unit=logfile,fmt=FMT_JF) h
                end if
                cycle attempt

              end if

              hjac = etah
              RenewJacobian = NO
              StaleJacobian = NO

            end if

           !!!
           !!! BCE Step

            call bce_step (u, udot, t, etah, errc)

            if (errc == 0) then     ! The BCE step was successful.
              exit bce
            end if

            if (StaleJacobian) then     ! Update jacobian and retry BCE step.

              RenewJacobian = YES
              u = u0 + h * (u1 + (h + hlast) * u2)
              udot = u1 + ((h + hlast) / eta) * u2
              retried_bce_steps = retried_bce_steps + 1
              if (msgs > 0) then
                write (unit=logfile,fmt=FMT_NR)
              end if
              cycle bce

            else                        ! Jacobian was fresh; cut h and retry.

              h = 0.5_wp * h
              FreezeCount = 1
              ! Should we insist upon reevaluating the Jacobian?
              failed_bce_steps = failed_bce_steps + 1
              if (msgs > 0) then
                write (unit=logfile,fmt=FMT_NF) h
              end if
              cycle attempt

            end if

          end do bce

         !!!
         !!! Local Truncation Error Estimate

          ! Predictor error.
          u3 = u - (u0+ h * (u1 + (h + hlast) * u2))
          perr = eval_norm (u3, 2)

          if (perr < 4.0_wp) then       ! Accept the step.

            if (msgs > 0) then
              write (unit=logfile,fmt=FMT_A) perr
            end if
            exit attempt

          else                          ! Reject the step; cut h and retry.

            h = 0.5_wp * h
            FreezeCount = 1
            ! Should we insist upon reevaluating the Jacobian?
            rejected_steps = rejected_steps + 1
            if (msgs > 0) then
              write (unit=logfile,fmt=FMT_PEF) perr, h
            end if
            cycle attempt

          end if

        end do attempt

       !!!
       !!! Clean-Up

        ptr => u3                       ! Shift the solution history.
        u3 => u2
        u2 => u1
        u1 => u0
        u0 => ptr
        u0 = u

        h3 = h2 + h                     ! Shift the time-step history.
        h2 = hlast + h
        hlast = h
        tlast = t

        u1 = (1.0_wp / hlast) * (u0 - u1)     ! Update the divided differences.
        u2 = (1.0_wp / h2)    * (u1 - u2)
        u3 = (1.0_wp / h3)    * (u2 - u3)

        ReturnCheck = YES

       !!!
       !!! Next Time Step

        call ChooseH (h, perr)
        h = max ( RMIN * hlast, min ( hub, RMAX * hlast, h ) )
        if (FreezeCount /= 0) then
          h = min ( hlast, h )
        end if
        FreezeCount = max ( 0, FreezeCount - 1 )
        if (msgs > 0) then
          write (unit=logfile,fmt=FMT_H) hlast, h
        end if

      end do bdf2_step

    end subroutine bdf2_solver

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ChooseH -- Choose a new time step.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ChooseH (h, perr)

      real(kind=wp), intent(inout) :: h
      real(kind=wp), intent(inout) :: perr
      real(kind=wp), parameter     :: tol = 0.001_wp
      real(kind=wp)                :: a, dh, phi, dphi

      a = 0.5_wp * hlast * h2 * h3 / max(perr, 0.001_wp)

      do ! until converged -- DANGEROUS!

        phi  = h * (h + hlast) * (h + h2) - a
        dphi = (2.0_wp * h + hlast) * (h + h2) + h * (h + hlast)

        dh = phi / dphi
        h = h - dh
        if (abs(dh) / h < tol) then
          exit
        end if

      end do

    end subroutine ChooseH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  BCE_STEP -- Backward Cauchy-Euler Step.
   !!
   !!  The Backward Cauchy-Euler method applied to the ODE system F(u,u',t) = 0
   !!  yields a nonlinear system of equations of the form
   !!
   !!     R(u) = F(u,udot,t) = 0,   udot = (u - ubar)/h,
   !!
   !!  for the unknown u,  where udot is the backward difference with time
   !!  step h.  This routine solves this nonlinear system by a modified Newton
   !!  iteration (aka fixed point iteration):
   !!
   !!     du := Solution of R' du = F(u_k,udot_k,t),
   !!     u_(k+1) := u_k - du,
   !!     udot_(k+1) := udot_k - du/h,
   !!
   !!  where the Jacobian matrix R' has been evaluated at some (u,udot,t) and
   !!  is used for every iteration.  The routine EVAL_RESIDUAL computes the
   !!  correction du (that is, the residual of the preconditioned system
   !!  R'^{-1} F(u,udot,t)).  The routine EVAL_NORM computes the norm of the
   !!  corrections du.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bce_step (u, udot, t, h, errc)

      use fp_accelerator

      real(kind=wp), intent(in)    :: t, h
      type(NodeVar), dimension(:), intent(inout) :: u, udot
      integer, intent(out)    :: errc

      ! local variables.
      integer  :: itr
      real(kind=wp) :: error
      type(NodeVar), dimension(size(u)) :: du

      if (mvec > 0) then
        call fpa_init (u, mvec, vtol)
      end if

      itr = 0

      do

        itr = itr + 1

        if (itr > mitr) then              ! Failed to converge.

          if (msgs > 1) then
            write (unit=logfile,fmt=FMT_NI) itr, error
          end if
          errc = 1
          exit

        end if

        residual_calls = residual_calls + 1
        call eval_residual (u, udot, t, du)  ! Solve for the correction.

        if (mvec > 0) then                ! Compute the accelerated correction.
          call fpa_correction (itr, du)
        end if

        u = u - du                        ! Next solution iterate
        udot = udot - ((1.0_wp / h) * du) ! and the backward difference.

        call check_soln (u, 1, errc)      ! Check the solution.

        if (errc /= 0) then               ! Solution is bad; exit.

          errc = 2
          if (msgs > 1) then
            write (unit=logfile,fmt=FMT_NIF) itr
          end if
          exit

        end if

        error = eval_norm (du, 1)         ! Estimated error.

        ! Check for convergence.
        ! We require 100 times greater accuracy on the first iteration.

        if ((error < 0.01_wp * ntol) .or. ((error < ntol) .and. (itr > 1))) then

          if (msgs > 1) then
            write (unit=logfile,fmt=FMT_NI) itr, error
          end if
          errc = 0
          exit

        end if

      end do

      if (mvec > 0) then
        call fpa_finish ()
      end if

    end subroutine bce_step

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  START_BDF2 -- Starting procedure for BDF2.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine start_bdf2 (u, udot, t, errc)

      type(NodeVar), dimension(:), intent(inout) :: u, udot
      real(kind=wp), intent(inout) :: t
      integer, intent(out)    :: errc

      ! local variable.
      real(kind=wp) :: etah
      type(NodeVar), dimension(:), pointer :: ptr

      u0 = u
      u1 = udot
      tlast = t

     !!!
     !!!  Step 1:  Trapezoid method with FCE as the predictor.

      if (msgs > 0) then
        write (unit=logfile,fmt=FMT_N) step, 1, tlast, h
      end if

      step = step + 1
      t = t + h
      etah = 0.5_wp * h

      ! Predicted solution and backward difference for the BCE step.
      u = u0 + h * u1
      udot = u1

      ! Check the predicted solution.
      call check_soln (u, 0, errc)
      if (errc /= 0) then
        return
      end if

      ! Evaluate the Jacobian at the predicted solution.
      jacobian_calls = jacobian_calls + 1
      call eval_jacobian (u, udot, t, etah, errc)
      if (errc /= 0) then
        return
      end if

      call bce_step (u, udot, t, etah, errc)
      if (errc /= 0) then
        return
      end if

      ! Shift the solution history.
      ptr => u2
      u2 => u1
      u1 => u0
      u0 => ptr
      u0 = u
      tlast = t

      ! Update the divided differences.
      u1 = (1.0_wp / h)            * (u0 - u1)
      u2 = (1.0_wp / (2.0_wp * h)) * (u1 - u2)

     !!!
     !!! Step 2:  BDF2 with FCE as the predictor.

      if (msgs > 0) then
        write (unit=logfile,fmt=FMT_N) step, 1, tlast, h
      end if

      step = step + 1
      t = tlast + h
      etah = (2.0_wp / 3.0_wp) * h

      ! Predicted solution and backward difference for the BCE step.
      u = u0 + h * (u1 + (2.0_wp * h) * u2)
      udot = u1 + (3.0_wp * h) * u2

      ! Check the predicted solution.
      call check_soln (u, 0, errc)
      if (errc /= 0) then
        return
      end if

      ! Evaluate the Jacobian at the predicted solution.
      jacobian_calls = jacobian_calls + 1
      call eval_jacobian (u, udot, t, etah, errc)
      if (errc /= 0) then
        return
      end if

      call bce_step (u, udot, t, etah, errc)
      if (errc /= 0) then
        return
      end if

      ! Shift the solution history.
      ptr => u3
      u3 => u2
      u2 => u1
      u1 => u0
      u0 => ptr
      u0 = u
      tlast = t

      ! Update the divided differences.
      u1 = (1.0_wp / h)            * (u0 - u1)
      u2 = (1.0_wp / (2.0_wp * h)) * (u1 - u2)
      u3 = (1.0_wp / (3.0_wp * h)) * (u2 - u3)

      errc = 0

    end subroutine start_bdf2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  SET_SOLVER_MESSAGES -- Set message level and output unit for the solver.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine set_solver_messages (level, unit)

      integer, intent(in) :: level, unit

      if (level >= 0) then
        msgs = level
      end if

      if (unit  >= 0) then
        logfile = unit
      end if

    end subroutine set_solver_messages

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  BDF2_WRITE_CHECKPOINT -- Read internal solver variables.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bdf2_write_checkpoint (cp)

      integer, intent(in) :: cp

      integer :: j

      write(unit=cp) hlb, hub, ntol, margin, vtol
      write(unit=cp) mtry, mitr, mvec
      write(unit=cp) h, hlast, h2, h3, hjac, tlast
      write(unit=cp) RenewJacobian, StaleJacobian, FreezeCount, ReturnState
      write(unit=cp) step, jacobian_calls, residual_calls, bad_predictors, failed_jacobians, &
                     retried_bce_steps, failed_bce_steps, rejected_steps

      write(unit=cp) shape(u0)

      do j = 1, size(u0)
        write(unit=cp) u0(j), u1(j), u2(j), u3(j)
      end do

    end subroutine bdf2_write_checkpoint

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  BDF2_READ_CHECKPOINT -- Read internal solver variables.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bdf2_read_checkpoint (cp)

      integer, intent(in) :: cp

      integer :: d1, j

      read(unit=cp) hlb, hub, ntol, margin, vtol
      read(unit=cp) mtry, mitr, mvec
      read(unit=cp) h, hlast, h2, h3, hjac, tlast
      read(unit=cp) RenewJacobian, StaleJacobian, FreezeCount, ReturnState
      read(unit=cp) step, jacobian_calls, residual_calls, bad_predictors, failed_jacobians, &
                    retried_bce_steps, failed_bce_steps, rejected_steps

      read(unit=cp) d1
      allocate (u0(d1), u1(d1), u2(d1), u3(d1))

      do j = 1, d1
        read(unit=cp) u0(j), u1(j), u2(j), u3(j)
      end do

      ! Don't assume the Jacobian was restored -- force an update;
      RenewJacobian = YES

    end subroutine bdf2_read_checkpoint

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  BDF2_INQUIRE -- Read internal solver variables.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bdf2_inquire (t, h_next, h_last, nstep, nres, njac, nbp, njf, nnr, nnf, nsr)

      real(kind=wp), intent(out), optional :: t, h_next, h_last
      integer,       intent(out), optional :: nstep, nres, njac, nbp, njf, nnr, nnf, nsr

      if (present(t)) then
        t = tlast
      end if

      if (present(h_next)) then
        h_next = h
      end if

      if (present(h_last)) then
        h_last = hlast
      end if

      if (present(nstep)) then
        nstep = step
      end if

      if (present(nres)) then
        nres = residual_calls
      end if

      if (present(njac)) then
        njac = jacobian_calls
      end if

      if (present(nbp)) then
        nbp = bad_predictors
      end if

      if (present(njf)) then
        njf = failed_jacobians
      end if

      if (present(nnr)) then
        nnr = retried_bce_steps
      end if

      if (present(nnf)) then
        nnf = failed_bce_steps
      end if

      if (present(nsr)) then
        nsr = rejected_steps
      end if

    end subroutine bdf2_inquire

end module mfe_ode_solver
