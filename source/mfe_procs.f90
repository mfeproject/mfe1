!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_procs

  use precision
  use mfe_extents
  use common_io, only: log_unit

  implicit none
  private

  ! Publically accessible procedures.
  public :: eval_norm, check_soln, eval_residual, eval_jacobian, eval_udot

  ! Storage for the Jacobian matrix
  real(wp), dimension(:,:,:), allocatable, save :: jaca, jacb, jacc

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  EVAL_NORM
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function eval_norm (du, type) result (norm)

      use norm_data

      real(wp)             :: norm
      real(wp), intent(in) :: du(:,:)
      integer, intent(in)  :: type

      integer j

      ! Weighted max-norm.
      norm = maxval(maxval(abs(du),2)/ptol)

      ! Relative norm on element lengths.
      do j = 1, nelt
        norm = max ( norm, abs(du(x,j+1) - du(x,j)) / (dx(j) * rtol) )
      end do

    end function eval_norm

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  CHECK_SOLN
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine check_soln (u, type, errc)

      use mfe_data,  only: dxmin
      use norm_data, only: dx

      real(wp), intent(in) :: u(:,:)
      integer, intent(in)  :: type
      integer, intent(out) :: errc

      integer j

     !!!
     !!! Check for bad elements, saving their lengths for ENORM.

      do j = 1, nelt

        dx(j) = u(x,j+1) - u(x,j)
        if (dx(j) > dxmin) cycle

        write (log_unit,"(/a,i4)") "** CHECK_SOLN: Bad element -- no. ", j
        errc = 1
        return

      end do

      errc = 0

    end subroutine check_soln

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  EVAL_RESIDUAL
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_residual (u, udot, t, r)

      use element_mfe
      use bc_procs
      use bls_solver

      real(wp), intent(in)  :: u(:,:), udot(:,:)
      real(wp), intent(out) :: r(:,:)
      real(wp), intent(in)  :: t

      ! Local variables.
      integer  :: errc
      real(wp) :: diag(nepn,nepn,nnod)

      ! Allocate storage for the local arrays.
      allocate (u0(nepn,nelt), u1(nepn,nelt), u0dot(nepn,nelt), u1dot(nepn,nelt))
      allocate (r0(nepn,nelt), r1(nepn,nelt), a00(nepn,nepn,nelt), a11(nepn,nepn,nelt))
      allocate (du(nepn,nelt), l(npde,nelt), n1(npde,nelt), n2(npde,nelt))

     !!!
     !!! EVALUATE THE RESIDUAL

      u0 = u(:,1:nelt)       ! Localize the solution arrays.
      u1 = u(:,2:nnod)
      u0dot = udot(:,1:nelt)
      u1dot = udot(:,2:nnod)

      call eltpp (errc)                ! Compute intermediate local terms.

      if (errc /= 0) then

        write (log_unit, "(/,a,i4)") "** EVAL_RESIDUAL: Bad element -- no. ", errc
        errc = 1
        return

      end if

      call eltrhs (t)                  ! Evaluate the residual locally.
      call regrhs
      call cres

      r(:,1)      = r0(:,1)            ! Assemble the residual.
      r(:,2:nelt) = r0(:,2:nelt) + r1(:,1:nelt-1)
      r(:,nnod)   = r1(:,nelt)

      call resbc (u, t, r)             ! BC contribution to the residual.

     !!!
     !!! Diagonal Preconditioning.

      call cdgl                        ! Evaluate the matrix diagonal locally.

      diag(:,:,1)      = a00(:,:,1)    ! Assemble the mass matrix block diagonal.
      diag(:,:,2:nelt) = a00(:,:,2:nelt) + a11(:,:,1:nelt-1)
      diag(:,:,nnod)   = a11(:,:,nelt)

      call jacbc (u, t, diag)          ! BC contribution to the diagonal.

      call vfct (diag)
      call vslv (diag, r)

     !!!
     !!! Solve.

      call btslv (jaca, jacb, jacc, r)

      errc = 0

      ! Free the storage for the local arrays.
      deallocate (u0, u1, u0dot, u1dot, r0, r1, a00, a11, du, l, n1, n2)

    end subroutine eval_residual

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  EVAL_JACOBIAN
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_jacobian (u, udot, t, h, errc)

      use element_mfe
      use bc_procs
      use bls_solver

      real(wp), intent(in) :: u(:,:), udot(:,:)
      real(wp), intent(in) :: t, h
      integer, intent(out) :: errc

      ! local variables.
      real(wp) diag(nepn,nepn,nnod)

      ! Allocate storage for the local arrays.
      allocate (u0(nepn,nelt), u1(nepn,nelt), u0dot(nepn,nelt), u1dot(nepn,nelt))
      allocate (r0(nepn,nelt), r1(nepn,nelt))
      allocate (a00(nepn,nepn,nelt), a10(nepn,nepn,nelt), a01(nepn,nepn,nelt), a11(nepn,nepn,nelt))
      allocate (du(nepn,nelt), l(npde,nelt), n1(npde,nelt), n2(npde,nelt))

      ! Localize the solution arrays.
      u0 = u(:,1:nelt)
      u1 = u(:,2:nnod)
      u0dot = udot(:,1:nelt)
      u1dot = udot(:,2:nnod)

      call eltpp (errc)

      if (errc /= 0) then

        write (log_unit, "(/,a,i4)") "** EVAL_JACOBIAN: Bad element -- no. ", errc
        errc = 1
        return

      end if

      call cmtx (-1.0_wp/h)

      ! Assemble the (unscaled) mass matrix block diagonal.
      diag(:,:,1)      = -h * a00(:,:,1)
      diag(:,:,2:nelt) = -h * (a00(:,:,2:nelt) + a11(:,:,1:nelt-1))
      diag(:,:,nnod)   = -h * a11(:,:,nelt)

      call jacbc (u, t, diag)

      call eltfy (t, errc)

      if (errc /= 0) then

        write (log_unit, "(/,a,i4)") "** ELTFY: Bad element -- no. ", errc
        errc = 2
        return

      end if

      if (.not. allocated(jaca)) &
      allocate (jaca(nepn,nepn,nnod), jacb(nepn,nepn,nnod), jacc(nepn,nepn,nnod))

      ! Assemble the Jacobian.
      jaca(:,:,1)      = a00(:,:,1)
      jaca(:,:,2:nelt) = a00(:,:,2:nelt) + a11(:,:,1:nelt-1)
      jaca(:,:,nnod)   = a11(:,:,nelt)

      jacb(:,:,1:nelt) = a01

      jacc(:,:,2:nnod) = a10

      call jacbc (u, t, jaca)

      ! Diagonal preconditioning.
      call vfct (diag)
      call vmslv (diag, jaca)
      call vmslv (diag, jacb)
      call vmslv (diag, jacc)

      ! Factorize the Jacobian.
      call btfct (jaca, jacb, jacc)

      errc = 0

      ! Free the storage for the local arrays.
      deallocate (u0, u1, u0dot, u1dot, r0, r1, a00, a10, a01, a11, du, l, n1, n2)

    end subroutine eval_jacobian

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  EVAL_UDOT
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_udot (u, t, udot, errc)

      use element_mfe
      use bc_procs
      use bls_solver

      real(wp), intent(in)  :: u(:,:)
      real(wp), intent(out) :: udot(:,:)
      real(wp), intent(in)  :: t
      integer, intent(out)  :: errc

      ! Allocate storage for the local arrays.
      allocate (u0(nepn,nelt), u1(nepn,nelt))
      allocate (r0(nepn,nelt), r1(nepn,nelt))
      allocate (a00(nepn,nepn,nelt), a10(nepn,nepn,nelt), a01(nepn,nepn,nelt), a11(nepn,nepn,nelt))
      allocate (du(nepn,nelt), l(npde,nelt), n1(npde,nelt), n2(npde,nelt))

      ! Localize the solution array.
      u0 = u(:,1:nelt)
      u1 = u(:,2:nnod)

      call eltpp (errc)

      if (errc /= 0) then

        write (log_unit, "(/,a,i4)") "** EVAL_UDOT: Bad element -- no. ", errc
        errc = 1
        return

      end if

      call cmtx

      if (.not. allocated(jaca)) &
      allocate (jaca(nepn,nepn,nnod), jacb(nepn,nepn,nnod), jacc(nepn,nepn,nnod))

      ! Assemble the mass matrix.
      jaca(:,:,1)      = a00(:,:,1)
      jaca(:,:,2:nelt) = a00(:,:,2:nelt) + a11(:,:,1:nelt-1)
      jaca(:,:,nnod)   = a11(:,:,nelt)

      jacb(:,:,1:nelt) = a01

      jacc(:,:,2:nnod) = a10

      call jacbc (u, t, jaca)

      call eltrhs (t)                  ! Evaluate the residual locally.
      call regrhs

      udot(:,1)      = r0(:,1)         ! Assemble the right hand side vector.
      udot(:,2:nelt) = r0(:,2:nelt) + r1(:,1:nelt-1)
      udot(:,nnod)   = r1(:,nelt)

      ! Solve for du/dt.
      call btfct (jaca, jacb, jacc)
      call btslv (jaca, jacb, jacc, udot)

      deallocate (jaca, jacb, jacc)
      errc = 0

      ! Free the storage for the local arrays.
      deallocate (u0, u1, r0, r1, a00, a10, a01, a11, du, l, n1, n2)

    end subroutine eval_udot

end module mfe_procs
