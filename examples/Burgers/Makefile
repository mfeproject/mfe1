VPATH = .:../../source

OBJ =	precision.o \
	mfe_extents.o \
	mfe_data.o \
	norm_data.o \
	bc.o \
	system_io.o \
	common_io.o \
	elt_array.o \
	elt_lapl.o \
	prob_procs.o \
	elt_mfe.o \
	bls_solver.o \
	mfe_procs.o \
	fpa.o \
	mfe_solver.o \
	init_procs.o \
	mfe1.o

SRC =	$(OBJ:.o=.f90)

LIB	= -lmass
F90	= f90
FFLAGS	= -O3 -qhot -qalias=noaryovrlp -qmoddir=./mod -I./mod $(EXTRA_FLAGS)

go: $(OBJ)
	$(F90) $(FFLAGS) -o go $(OBJ) $(LIB)

runem: ex1 ex2
ex1: go
	cp -f Data/mfein-ex1 mfein
	go
	mv -f mfelog mfelog-ex1
	mv -f mfegrf mfegrf-ex1
	@echo "Results moved to mfelog-ex1 and mfegrf-ex2;"
	@echo "compare with those archived the directory Data."

ex2: go
	cp -f Data/mfein-ex2 mfein
	go
	mv -f mfelog mfelog-ex2
	mv -f mfegrf mfegrf-ex2
	@echo "Results moved to mfelog-ex2 and mfegrf-ex2;"
	@echo "compare with those archived the directory Data."

clean:;	rm -f *.o mod/* go mfelog* mfegrf* bdfout

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

.SUFFIXES: .f90

.f90:
	$(F90) $(FFLAGS) $(LDFLAGS) $< -o $@

.f90.o:
	$(F90) $(FFLAGS) -c $<

