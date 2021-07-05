!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_ode_solver

  use prec
  use mfe_procs, only: check_soln, eval_jacobian, eval_residual, eval_norm

  !implicit none
  private

  public  :: bdf2_solver, get_soln_profile, set_solver_messages
  private :: start_bdf2, bce_step, ChooseH
  
  ! Mode codes.
  integer, parameter, public :: START_SOLN  = 0
  integer, parameter, public :: RESUME_SOLN = 1

  ! Return codes.
  integer, parameter, public :: FAIL_ON_STEP = -1
  integer, parameter, public :: SMALL_H_FAIL = -2
  integer, parameter, public :: FAIL_ON_START = -3
  integer, parameter, public :: SOLN_AT_TOUT = 2
  integer, parameter, public :: SOLN_AT_STEP = 3

  ! State codes.
  integer, parameter, private :: RET_ON_TIME = 2
  integer, parameter, private :: RET_ON_STEP = 3
  integer, parameter, private :: RET_ON_FAIL = 4

  ! Aliases.
  logical, parameter, private :: YES = .true.
  logical, parameter, private :: NO = .false.

  ! Diagnostic message formats.
  integer, save, private :: logfile=6, msgs=0
  character(len=*), parameter, private ::                                  &
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

  ! Type for bundling the solution profile.
  type, public :: soln_profile
    integer  :: nstep, nje, nre, npf, njf, nnr, nnf, npef
    real(kind=wp) :: tlast, hlast
  end type soln_profile

  ! Solver parameters.
  real(kind=wp), save, private :: hlb, hub, ntol, margin, vtol
  integer, save, private  :: mtry, mitr, mvec

  real(kind=wp), parameter, private :: RMIN = 0.25_wp, RMAX = 4.0_wp

  ! Solution state and history variables.
  real(kind=wp), save, private :: h, hlast, h2, h3, hjac, tlast
  logical, save, private  :: RenewJacobian, StaleJacobian
  integer, save, private  :: FreezeCount, ReturnState
  integer, save, private  :: nstep, nje, nre, npf, njf, nnr, nnf, npef

  real(kind=wp), pointer, dimension(:,:), private :: u0, u1, u2, u3, ptr

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  BDF2_SOLVER -- Second order backward difference solver.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bdf2_solver (mode, rvar, ivar, tout, rfreq, ReturnType, u, udot, t)

      integer, dimension(3), intent(in)     :: ivar
      integer, intent(in) :: rfreq, mode
      integer, intent(out)    :: ReturnType
      real(kind=wp), dimension(6), intent(in)    :: rvar
      real(kind=wp), intent(in) :: tout
      real(kind=wp), dimension(:,:), intent(inout) :: u, udot
      real(kind=wp), intent(inout) :: t

      ! local variables.
      logical  :: ReturnCheck
      integer  :: errc, try, d1, d2
      real(kind=wp) :: perr, eta, etah, dt

      select case (mode)

        case (START_SOLN)

          d1 = size(u,1)
          d2 = size(u,2)
          allocate (u0(d1,d2), u1(d1,d2), u2(d1,d2), u3(d1,d2))

          mtry = ivar(1)
          mitr = ivar(2)
          mvec = ivar(3)
          h    = rvar(1)
          hlb  = rvar(2)
          hub  = rvar(3)
          ntol = rvar(4)
          margin = rvar(5)
          vtol = rvar(6)

          nstep = 0
          nje = 0
          nre = 0
          npf = 0
          njf = 0
          nnr = 0
          nnf = 0
          npef = 0

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



      step: do

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

          else if (mod(nstep,rfreq) == 0) then   ! Return current solution.

            t = tlast
            u = u0
            ReturnType  = SOLN_AT_STEP
            ReturnState = RET_ON_STEP
            exit

          end if

        end if

        nstep = nstep + 1
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
            exit step

          end if

          if (h < hlb) then             ! Time step is too small.

            t = tlast
            u = u0
            ReturnType  = SMALL_H_FAIL
            ReturnState = RET_ON_FAIL
            exit step

          end if

          if (msgs > 0) then
            write (unit=logfile,fmt=FMT_N) nstep, try, tlast, h
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
            npf = npf + 1
            if (msgs > 0) then
              write (unit=logfile,fmt=FMT_PF) h
            end if
            cycle attempt

          end if

          bce: do

           !!!
           !!! Jacobian Evaluation

            if (RenewJacobian) then     ! Reevaluate the Jacobian.

              nje = nje + 1
              if (msgs > 0) then
                write (unit=logfile,fmt=FMT_J)
              end if

              call eval_jacobian (u, udot, t, etah, errc)

              if (errc /= 0) then         ! Evaluation failed; cut h and retry.

                h = 0.5_wp * h
                RenewJacobian = YES
                njf = njf + 1
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
              nnr = nnr + 1
              if (msgs > 0) then
                write (unit=logfile,fmt=FMT_NR)
              end if
              cycle bce

            else                        ! Jacobian was fresh; cut h and retry.

              h = 0.5_wp * h
              FreezeCount = 1
              ! Should we insist upon reevaluating the Jacobian?
              nnf = nnf + 1
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
            npef = npef + 1
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

        u1 = (u0 - u1) / hlast          ! Update the divided differences.
        u2 = (u1 - u2) / h2
        u3 = (u2 - u3) / h3

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

      end do step

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

      a = 0.5_wp * hlast * h2 * h3 / perr

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
      real(kind=wp), dimension(:,:), intent(inout) :: u, udot
      integer, intent(out)    :: errc

      ! local variables.
      integer  :: itr
      real(kind=wp) :: error, last_error, conv_rate
      real(kind=wp), dimension(size(u,1),size(u,2)) :: du

      if (mvec > 0) then
        call fpa_init (u, mvec, vtol)
      end if
      
      itr = 0
      last_error = 1.0_wp

      do

        itr = itr + 1

        if (itr > mitr) then              ! Failed to converge.

          if (msgs > 1) then
            write (unit=logfile,fmt=FMT_NI) itr, last_error, conv_rate
          end if
          errc = 1
          exit

        end if

        nre = nre + 1
        call eval_residual (u, udot, t, du)  ! Solve for the correction.

        if (mvec > 0) then                ! Compute the accelerated correction.
          call fpa_correction (itr, du)
        end if

        u = u - du                        ! Next solution iterate
        udot = udot - du / h              ! and the backward difference.

        call check_soln (u, 1, errc)      ! Check the solution.

        if (errc /= 0) then               ! Solution is bad; exit.

          errc = 2
          if (msgs > 1) then
            write (unit=logfile,fmt=FMT_NIF) itr
          end if
          exit

        end if

        error = eval_norm (du, 1)         ! Estimated error.
        conv_rate = error / last_error
        last_error = error

        ! Check for convergence.
        ! We require 100 times greater accuracy on the first iteration.

        if ((error < 0.01_wp * ntol) .or. ((error < ntol) .and. (itr > 1))) then

          if (msgs > 1) then
            write (unit=logfile,fmt=FMT_NI) itr, error, conv_rate
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

      real(kind=wp), dimension(:,:), intent(inout) :: u, udot
      real(kind=wp), intent(inout) :: t
      integer, intent(out)    :: errc

      ! local variable.
      real(kind=wp) :: etah

      u0 = u
      u1 = udot
      tlast = t

     !!!
     !!!  Step 1:  Trapezoid method with FCE as the predictor.

      if (msgs > 0) then
        write (unit=logfile,fmt=FMT_N) nstep, 1, tlast, h
      end if

      nstep = nstep + 1
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
      nje = nje + 1
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
      u1 = (u0 - u1) / h
      u2 = (u1 - u2) / (2.0_wp * h)

     !!!
     !!! Step 2:  BDF2 with FCE as the predictor.

      if (msgs > 0) then
        write (unit=logfile,fmt=FMT_N) nstep, 1, tlast, h
      end if

      nstep = nstep + 1
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
      nje = nje + 1
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
      u1 = (u0 - u1) / h
      u2 = (u1 - u2) / (2.0_wp * h)
      u3 = (u2 - u3) / (3.0_wp * h)

      errc = 0

    end subroutine start_bdf2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  GET_SOLN_PROFILE -- Get the current solution profile.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function get_soln_profile () result (p)

      type(soln_profile) :: p

      p % nstep = nstep    ! Number of steps completed.
      p % nje = nje        ! Number of Jacobian evaluations.
      p % nre = nre        ! Number of residual evaluations.
      p % npf = npf        ! Number of predicted solution failures.
      p % njf = njf        ! Number of Jacobian evaluation failures.
      p % nnr = nnr        ! Number of retries of a BCE step.
      p % nnf = nnf        ! Number of complete BCE step failures.
      p % npef = npef      ! Number of predictor error failures.
      p % hlast = hlast    ! Last successful time step
      p % tlast = tlast    ! Time of last solution.

    end function get_soln_profile

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

end module mfe_ode_solver
