July 2021

This code dates from the Fall of 1996. It solves the nonlinear diffusion of
a dopant (e.g. arsenic) in semiconductor crystal. This particular formulation
uses a transformed variable to better capture the exponential tails.

The particular problem setup here starts with the dopant implanted at depth
below the surface of the semiconductor with a gaussian profile, which then
diffuses.

The directory originally contained some working files which I decided weren't
important to save. In particular there was a subdirectory 2d which contained
an older F77 version of GWMFE2DS, which appears to be there only for reference.

The program `fit.f` is a standalone program that generates the initial
gaussian profile solution, using the best L2-fit algorithm I developed in
the late 1980's while visiting Mike Baines at the University of Reading in
the UK. **This may be the only copy of this program that I still have.**

The source file `att.f90` is not used. It's a version of `problem.f90` that
uses a different, or differently parameterized, transformation. I chose to
hang onto this in case there is something important here.
