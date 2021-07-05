VPATH = .:../source/common:../source/system

VERSION = 0.2
PROBLEM = NAVIER-STOKES

###############################################################################
#
# Select your Fortran 90 compiler.
#
#   Imagine1's F compiler (F90 subset language):
#
#FC = F
#FSUFFIX = F
#FFLAGS = -O
#LIB =
#
#   IBM"s XLF compiler:
#
 FC = f90
 FSUFFIX = f90
 FFLAGS = -u -O3 -qarch=604 #-Q -qalias=noaryovrlp
 LIB = -lmass
#
###############################################################################

OBJ = mfeconst.o comio.o mfedata.o mfetypes.o locarray.o loclapl.o \
      problem.o locmfe.o locglb.o locmfeprocs.o blssolve.o bc.o \
      mfeprocs.o norm.o fpa.o mfesolve.o init.o output.o mfe1.o

LIBOBJ = kindpar.o sysio.o

SRC = $(LIBOBJ:.o=.$(FSUFFIX)) $(OBJ:.o=.$(FSUFFIX))

go: $(LIBOBJ) $(OBJ)
	$(FC) $(FFLAGS) -o go $(LIBOBJ) $(OBJ) $(LIB)

runem: ex1
ex1: go
	@cp -f Data/mfein-ex1 mfein
	@chmod u+w mfein
	time go
	@mv -f mfelog mfelog-ex1
	@mv -f mfegrf mfegrf-ex1
	@echo "Results moved to mfelog-ex1 and mfegrf-ex1;"
	@echo "compare with those archived the directory Data."
	
source: title.txt $(SRC)
	cat $? > MFE1-$(VERSION)-$(PROBLEM).$(FSUFFIX)

clean:;	-rm -f *.o *.mod f_files/* F.err mod/* go mfelog* mfegrf* bdfout

# Module dependencies

bc.o: mfeconst.o mfetypes.o
blssolve.o: mfeconst.o mfetypes.o
comio.o: mfeconst.o sysio.o
fpa.o: mfeconst.o mfetypes.o
init.o: bc.o comio.o mfeconst.o mfedata.o mfetypes.o norm.o problem.o
locarray.o: mfeconst.o mfetypes.o
locglb.o: locarray.o mfeconst.o mfetypes.o
loclapl.o: locarray.o mfeconst.o
locmfe.o: locarray.o mfeconst.o mfedata.o mfetypes.o problem.o
locmfeprocs.o: locglb.o locmfe.o problem.o
mfe1.o: comio.o init.o mfeconst.o mfeprocs.o mfesolve.o mfetypes.o output.o \
        sysio.o
mfeconst.o: kindpar.o
mfedata.o: mfeconst.o
mfeprocs.o: bc.o blssolve.o comio.o locmfeprocs.o mfeconst.o mfetypes.o
mfesolve.o: fpa.o mfeconst.o mfeprocs.o mfetypes.o norm.o
mfetypes.o: mfeconst.o
norm.o: comio.o mfeconst.o mfedata.o mfetypes.o
output.o: comio.o mfeconst.o mfetypes.o
problem.o: comio.o locarray.o loclapl.o mfeconst.o

.SUFFIXES: .F .f90

.F:
	$(FC) $(FFLAGS) $(LDFLAGS) $< -o $@
.F.o:
	$(FC) $(FFLAGS) -c $<
.f90:
	$(FC) $(FFLAGS) $(LDFLAGS) $< -o $@
.f90.o:
	$(FC) $(FFLAGS) -c $<
