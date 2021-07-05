      program gridgen
************************************************************************
*
************************************************************************
      integer n,itr,j,ios
      real*8 x(0:100),xtraj(0:100,0:500),w1(100),w2(100),tol,e
      double precision xl,xr,h
      external f,simp,gauss3
      character*8 label
      real*8 dose, sig, pos, mu, c0
      common /funcdata/ dose, sig, pos, mu, c0

  800 format(a$)
  810 format('*** Crossed knots at iteration ',i3)
  820 format('*** Failed to converge.')
  830 format('*** Converged.')

      read(5,*) tol
      read(5,*) dose
      read(5,*) sig
      read(5,*) pos
      read(5,*) mu
      read(5,*) c0
      
      open(unit=1,file="soln",position="rewind")
      do j = 0, 100
        read(unit=1,fmt=*,iostat=ios) x(j)
        if(ios /= 0) then
          n = j - 2
          exit
        end if
      end do
      close(unit=1)

      itr = 0
      
      do 10 j=0,n+1
       xtraj(j,itr) = x(j)
   10 continue

  100 itr = itr + 1
      call sglstp(n,x,f,gauss3,e)

*     Record knot positions.
      do 110 j=0,n+1
  110 xtraj(j,itr) = x(j)

*     Check for crossed knots.
      do 120 j=1,n+1
       if(x(j) .le. x(j-1)) then
         write(6,810) itr
         call dplfit(n,x,f,gauss3,  w1,w2)
         call output(n,itr,xtraj,w1,w2)
         stop
         endif
  120 continue

*     Check for convergence.
      if(e .le. tol) then
        write(6,830)
        call dplfit(n,x,f,gauss3,  w1,w2)
        call output(n,itr,xtraj,w1,w2)
        stop
        endif

      if(itr .lt. 500) go to 100

      write(6,820)
      call dplfit(n,x,f,gauss3,  w1,w2)
      call output(n,itr,xtraj,w1,w2)

      stop
      end
      subroutine output(n,nitr,xtraj,w1,w2)

      integer n,nitr,j,itr
      real*8 xtraj(0:100,0:500),w1(100),w2(100)

      open(7,file="traj")

      do 10 j=0,n+1
       write(7,810) "Node ", j
       write(7,820) (xtraj(j,itr), itr, itr=0,nitr)
   10 continue

      open(8,file="soln",status="replace",position="rewind")

      write(6,810) "Best fit after iteration ", nitr
      write(8,fmt="(2es13.5)") xtraj(0,nitr), w1(1)
      do 20 j=1,n
       write(8,fmt="(2es13.5)") xtraj(j,nitr), .5*(w2(j) + w1(j+1))
   20 continue
      write(8,fmt="(2es13.5)") xtraj(n+1,nitr), w2(n+1)

  810 format(a,i3)
  820 format((e12.5,i5))
  830 format(2e13.5)

      return
      end



      subroutine dplfit(n,x,f,quad,  w1,w2)

      integer n,j
      real*8 x(0:n+1),w1(n+1),w2(n+1),b1,b2
      external f,quad

      do 10 j=1,n+1
       call quad(x(j-1),x(j),f,  b1,b2)
       w1(j) = 4.d0*b1 - 2.d0*b2
       w2(j) = 4.d0*b2 - 2.d0*b1
   10 continue

      return
      end


      real*8 function f(x)

      real*8 x
      real*16 xx
      real*8 dose, sig, pos, mu, c0
      common /funcdata/ dose, sig, pos, mu, c0

      real(kind=8), parameter :: twopi = 6.283185307179587_8

      xx = (x - pos) / (sqrt(2.0) * sig)
      f = c0 + (dose / (2.0 * sqrt(twopi) * sig)) * erfc (xx)

      f = f/mu + log(f/mu)

      return
      end
      
c      double precision function f(x)
c      double precision x
c
c      real*8 dose, sig, pos, mu, c0
c      common /funcdata/ dose, sig, pos, mu, c0
c
c      real(kind=8), parameter :: twopi = 6.283185307179587_8
c
c      f = c0 + (dose / (sqrt(twopi) * sig)) *
c     &        exp (- (x - pos) ** 2 / (2.0 * sig ** 2 ))
c
c      f = f/mu + log(f/mu)
c
c      return
c      end
      subroutine sglstp(n,x,f,quad,e)
************************************************************************
*
*   SGLSTP -- Performs a single simultaneous-step of the algorithm.
*
*   Argument  I/O/M/T  Description
*   --------  -------  -----------
*     N          I     Number of free knots.
*
*     X          M     The coordinates of the knots.  X(0) and X(N+1)
*                      are the coordinates of the fixed endpoints and
*                      are not modified.
*
*     E          O     The max norm of the correction.
*
*     F                External function name to be fitted.
*
*     QUAD             External subroutine name of the quadrature
*                      method to be used.
*
************************************************************************
      integer n,j
      double precision x(0:n+1),e,b1,b2,w1l,w2l,w1r,w2r,ml,mr,dx
      external f,quad
      
*     Best linear fit on [x_0, x_1].
      call quad(x(0),x(1),f,  b1,b2)
      w1r = 4.d0*b1 - 2.d0*b2
      w2r = 4.d0*b2 - 2.d0*b1
      mr = (w2r - w1r)/(x(1) - x(0))

      e = 0.d0
      do 10 j=1,n
*      Save the previously computed best fit on [x_(j-1),x_j].
       w1l = w1r
       w2l = w2r
       ml = mr

*      Best linear fit on [x_j, x_(j+1)].
       call quad(x(j),x(j+1),f,  b1,b2)
       w1r = 4.d0*b1 - 2.d0*b2
       w2r = 4.d0*b2 - 2.d0*b1
       mr = (w2r - w1r)/(x(j+1) - x(j))

*      Next iterate for x_j.
       dx = (w1r - w2l)/(mr - ml)
       x(j) = x(j) - dx
       e = max( e, abs(dx) )
   10 continue

      return
      end
      subroutine sweep(n,x,f,quad,e)
************************************************************************
*
*   SWEEP -- Performs a single sweep-step of the algorithm.
*
*   Argument  I/O/M/T  Description
*   --------  -------  -----------
*     N          I     Number of free knots.
*
*     X          M     The coordinates of the knots.  X(0) and X(N+1)
*                      are the coordinates of the fixed endpoints and
*                      are not modified.
*
*     E          O     The max norm of the correction.
*
*     F                External function name to be fitted.
*
*     QUAD             External subroutine name of the quadrature
*                      method to be used.
*
************************************************************************
      integer n,j
      double precision x(0:n+1),e,b1,b2,w1l,w2l,w1r,w2r,ml,mr,dx
      external f,quad
      
      e = 0.d0
      do 10 j=1,n
*      Best linear fit on [x_(j-1), x_j].
       call quad(x(j-1),x(j),f,  b1,b2)
       w1l = 4.d0*b1 - 2.d0*b2
       w2l = 4.d0*b2 - 2.d0*b1
       ml  = (w2l - w1l)/(x(j) - x(j-1))

*      Best linear fit on [x_j, x_(j+1)].
       call quad(x(j),x(j+1),f,  b1,b2)
       w1r = 4.d0*b1 - 2.d0*b2
       w2r = 4.d0*b2 - 2.d0*b1
       mr  = (w2r - w1r)/(x(j+1) - x(j))

*      Next iterate for x_j.
       dx = (w1r - w2l)/(mr - ml)
       x(j) = x(j) - dx
       e = max( e, abs(dx) )
   10 continue

      return
      end
      subroutine simp(x1,x2,f,  fphi1,fphi2)
************************************************************************
*
*   SIMP -- Simpson's rule quadrature.
*
************************************************************************
      double precision x1,x2,f,fphi1,fphi2
      external f

      fphi1 = (f(x1) + 2.d0*f(.5d0*(x1+x2)))/6.d0
      fphi2 = (f(x2) + 2.d0*f(.5d0*(x1+x2)))/6.d0

      return
      end
      subroutine gauss3(x1,x2,f,  fphi1,fphi2)
************************************************************************
*
*   Gauss3 -- 3 point Gauss quadrature.
*
************************************************************************
      double precision x1,x2,f,fphi1,fphi2,a,b1,b2,b3
      external f

      parameter(a=.112701665379259d0, b1=.31306018160905d-1,
     +         b2=.222222222222222d0, b3=.24647175961687d0)

      fphi1 = b3*f(x1+a*(x2-x1))+b2*f(.5d0*(x1+x2))+b1*f(x2-a*(x2-x1))
      fphi2 = b1*f(x1+a*(x2-x1))+b2*f(.5d0*(x1+x2))+b3*f(x2-a*(x2-x1))

      return
      end
      block data
************************************************************************
*
*   Quadrature formula:
*
*   G1D3 -- 3 point 1D Gaussian quadrature formula.  Fifth order
*   accurate.  The weights are g3w(.) and the corresponding quadrature
*   points are (g3a0(.), g3a1(.)) given as the isoparametric coordinates
*   (alpha_0, alpha_1) of the interval.
*
************************************************************************
      real*8 g3w,g3a0,g3a1
      common /g1d3/ g3w(3),g3a0(3),g3a1(3)

      integer k

      data g3w/
     + 2.77777777777778d-1,4.44444444444444d-1,2.77777777777778d-1/
      data (g3a0(k),g3a1(k),k=1,3)/
     + 1.12701665379259d-1, 8.87298334620741d-1,
     + 5.00000000000000d-1, 5.00000000000000d-1,
     + 8.87298334620741d-1, 1.12701665379259d-1/
 
      end
