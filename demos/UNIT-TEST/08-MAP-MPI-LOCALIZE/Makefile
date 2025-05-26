
topdir := $(shell i=. ; for c in 0 1 2 3 4 5 6 7 8; do [ ! -f $$i/.chimas_top ] && i=$$i/..; done; echo $$i)#
TOPDIR=$(shell realpath $(topdir))

include $(topdir)/share/makefile/Makefile.demos

export CHDEBUG=0

run-%: set-% .force
	mpirun -n 4 $(topdir)/build/$*/map-mpi-localize

gdb-dbg: set-dbg .force
	gdb --ex 'set breakpoint pending on' --ex 'b MPI_Abort' --ex 'b breakpoint' --ex 'run' --args $(topdir)/build/dbg/map-mpi-localize $(PETSC_OPT)

r:
	 g++ -I../../../include -o map-mpi-localize map-mpi-localize.cpp
	./map-mpi-localize
