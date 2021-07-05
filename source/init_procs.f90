!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module init_procs

  use precision

  implicit none
  private

  public :: read_soln, read_data

  integer, public :: ofreq, mstep, debug
  real(wp), allocatable, public :: tout(:)

  ! MFE ODE solver parameters.
  integer, public  :: mtry, mitr, mvec
  real(wp), public :: h, hlb, hub, ntol, margin, vtol

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !!  READ_SOLN
    !!
    !!    This procedure should:
    !!    1) complete the initialization of the data in mfe_extents;
    !!    2) generate the initial solution u and allocate storage for udot; and
    !!    3) initialize all data in the module bc_data.
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_soln (u, udot)

      use common_io
      use mfe_extents
      use bc_data

      real(wp), pointer :: u(:,:), udot(:,:)

      integer nseg
      integer, allocatable  :: niseg(:)
      real(wp), allocatable :: useg(:,:)

      write (log_unit, "(//a/)") "I N I T I A L   C O N D I T I O N S"

      call read_tagged_data (nseg, "Mesh segments (NSEG)")

      allocate (niseg(nseg), useg(nepn,nseg+1))

      call read_tagged_data (niseg, "Elements per segment (NISEG)")
      call read_tagged_data (useg, "Initial solution (USEG)")

      nelt = sum (niseg)
      nnod = nelt + 1

      allocate (u(nepn,nnod), udot(nepn, nnod))

      call refine (useg, niseg, u)

      deallocate (niseg, useg)

      write (log_unit, *)
      call read_tagged_data (bcl, "Left-end BC (BCL)")
      call read_tagged_data (bcr, "Right-end BC (BCR)")

      ! Save the initial boundary values.
      bvl = u(:,1)
      bvr = u(:,nnod)

    end subroutine read_soln

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  READ_DATA
   !!
   !!    This procedure should:
   !!    1) initialize the parameters in module mfe_data;
   !!    2) initialize the parameters in module norm_data;
   !!    3) ...
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_data

      use common_io
      use mfe_extents
      use mfe_data
      use norm_data
      use problem_procs, only: read_prob_data

      integer nout

      write (log_unit, "(//a/)") "M F E   P A R A M E T E R S"

      if (npde > 1) then
        call read_tagged_data (wpde, "PDE weights (WPDE)")
      else
        wpde(1) = 1.0_wp
      end if

      call read_tagged_data (kreg, "Segment visc type (KREG)")
      call read_tagged_data (segvsc, "Segment visc coef (SEGVSC)")
      call read_tagged_data (segspr, "Segment spring coef (SEGSPR)")
      call read_tagged_data (dxmin, "Minimum element size (DXMIN)")
      call read_tagged_data (fdinc, "Jacob FD increment (FDINC)")

      write (log_unit, *)
      call read_tagged_data (debug, "Solver debug level (DEBUG)")
      call read_tagged_data (ofreq, "Output frequency (OFREQ)")
      call read_tagged_data (mstep, "Maximum steps (MSTEP)")

      call read_tagged_data (nout)
      allocate (tout(nout))
      call read_tagged_data (tout, "Output times (TOUT)")

      write (log_unit, *)
      call read_tagged_data (rtol, "Relative dx tolerance (RTOL)")
      call read_tagged_data (ptol, "Predictor tolerance (PTOL)")
      allocate (dx(nelt))

      write (log_unit, *)
      call read_tagged_data (h, "Initial time step (H)")
      call read_tagged_data (hlb, "Time step lower bound (HLB)")
      call read_tagged_data (hub, "Time step upper bound (HUB)")
      call read_tagged_data (margin, "Jacob update margin (MARGIN)")
      call read_tagged_data (ntol, "FPI tolerance (NTOL)")
      call read_tagged_data (mitr, "Max FP iterations (MITR)")

      write (log_unit, *)
      call read_tagged_data (vtol, "FPA vector tolerance (VTOL)")
      call read_tagged_data (mvec, "Max FPA vectors (MVEC)")

      mtry = 9

      write (log_unit, "(//a/)") "P R O B L E M   S P E C I F I C   P A R A M E T E R S"

      call read_prob_data

    end subroutine read_data

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  REFINE
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine refine (x0, n, x)

      integer, intent(in)   :: n(:)
      real(wp), intent(in)  :: x0(:,:)
      real(wp), intent(out) :: x(:,:)

      integer node, j, k
      real(wp) dx(size(x0,1))

      node = 1

      do k = 1, size(n)

        x(:,node) = x0(:,k)
        node = node + 1

        dx = (x0(:,k+1) - x0(:,k)) / n(k)

        do j = 1, n(k) - 1

          x(:,node) = x(:,node-1) + dx
          node = node + 1

        end do

      end do

      x(:,node) = x0(:,size(n)+1)

    end subroutine refine

end module init_procs
