
topdir := $(shell i=. ; for c in 0 1 2 3 4 5 6 7 8; do [ ! -f $$i/.chimas_top ] && i=$$i/..; done; echo $$i)#
TOPDIR=$(shell realpath $(topdir))

include $(topdir)/share/makefile/Makefile.demos

export CHDEBUG=1

run-%: set-% .force
	 mpirun -n 20 $(topdir)/build/$*/ut-netcdf
#     $(topdir)/build/$*/ut-netcdf

gdb-devel: set-devel .force
	gdb --ex 'set breakpoint pending on' --ex 'b MPI_Abort' --ex 'b breakpoint' --ex 'run' --args $(topdir)/build/devel/ut-netcdf

