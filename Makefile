F	= f90
F90	= f90

SRC =	precision.f90 mfe_extents.f90 mfe_data.f90 norm_data.f90 bc.f90 \
	system_io.f90 common_io.f90 elt_array.f90 elt_lapl.f90 prob_procs.f90 \
	elt_mfe.f90 bls_solver.f90 mfe_procs.f90 fpa.f90 mfe_solver.f90 \
	init_procs.f90 mfe1.f90

OBJ =	precision.o mfe_extents.o mfe_data.o norm_data.o bc.o \
	system_io.o common_io.o elt_array.o elt_lapl.o prob_procs.o \
	elt_mfe.o bls_solver.o mfe_procs.o fpa.o mfe_solver.o \
	init_procs.o mfe1.o

LIB =

FFLAGS	= -O

go: $(OBJ)
	$(F) $(FFLAGS) -o go $(OBJ) $(LIB)

runem: ex1
ex1: go
	cp -f Data/mfein-ex1 mfein
	go
	mv -f mfelog mfelog-ex1
	mv -f mfegrf mfegrf-ex1
	echo "Results moved to mfelog-ex1 and mfegrf-ex2;"
	echo "compare with those archived the directory Data."

clean:;	rm -f *.o mod/* go mfelog mfegrf bdfout

source:;	cat $(SRC) > MFE1-all.f90

# Dependencies...

mfe1.o:		precision.o common_io.o mfe_extents.o norm_data.o init_procs.o \
		mfe_solver.o mfe_procs.o
init_procs.o:	precision.o common_io.o mfe_extents.o bc.o mfe_data.o \
		norm_data.o prob_procs.o
mfe_solver.o:	precision.o mfe_procs.o fpa.o
mfe_procs.o:	precision.o common_io.o norm_data.o mfe_extents.o mfe_data.o \
		bc.o elt_mfe.o bls_solver.o
elt_mfe.o:	precision.o mfe_extents.o mfe_data.o elt_array.o prob_procs.o
prob_procs.o:	precision.o common_io.o mfe_extents.o mfe_data.o elt_lapl.o \
		elt_array.o
elt_lapl.o:	precision.o mfe_extents.o mfe_data.o bc.o elt_array.o
bc.o:		precision.o mfe_extents.o
mfe_data.o:	precision.o mfe_extents.o
elt_array.o:	precision.o
norm_data.o:	precision.o mfe_extents.o
bls_solver.o:	precision.o
fpa.o:		precision.o
common_io.o:	precision.o system_io.o
mfe_extents.o:	
precision.o:
system_io:

.SUFFIXES: .f90 .F

.F:
	$(F) $(FFLAGS) $(LDFLAGS) $< -o $@

.F.o:
	$(F) $(FFLAGS) -c $<

.f90:
	$(F90) $(FFLAGS) $(LDFLAGS) $< -o $@

.f90.o:
	$(F90) $(FFLAGS) -c $<

