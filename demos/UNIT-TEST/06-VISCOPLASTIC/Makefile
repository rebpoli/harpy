
topdir := $(shell i=. ; for c in 0 1 2 3 4 5 6 7 8; do [ ! -f $$i/.chimas_top ] && i=$$i/..; done; echo $$i)#
TOPDIR=$(shell realpath $(topdir))

include $(topdir)/share/makefile/Makefile.demos

export CHDEBUG=1
export N=1

PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -ksp_monitor -pc_type lu -ksp_type cg -snes_converged_reason -ksp_converged_reason


run-%: set-% .force
	rm -rf run/
	mpirun -n $(N) $(topdir)/build/$*/ut-viscoplastic $(PETSC_OPT)

gdb-dbg: set-dbg .force
	gdb --ex 'set breakpoint pending on' --ex 'b MPI_Abort' --ex 'b breakpoint' --ex 'run' --args $(topdir)/build/dbg/ut-viscoplastic $(PETSC_OPT)

gdb: gdb-dbg



#     rm -f run/exo/*   # Avoid deleting dir to keep dir reference (Paraview likes)
#     rm -f run/log/*
#     rm -f run/csv/*
#     rm -f run/msh/*
