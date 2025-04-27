
topdir := $(shell i=. ; for c in 0 1 2 3 4 5 6 7 8; do [ ! -f $$i/.chimas_top ] && i=$$i/..; done; echo $$i)#
TOPDIR=$(shell realpath $(topdir))

include $(topdir)/share/makefile/Makefile.demos

export CHDEBUG=1
export N=1

PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -pc_type lu -ksp_type gmres -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-10 -snes_max_it 10 -snes_atol 1e-7 -snes_linesearch_type basic -snes_divergence_tolerance 1e10 -ksp_rtol 1e-20
# PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -ksp_monitor -pc_type lu -ksp_type gmres -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-10 -snes_max_it 10 -snes_atol 1e-7 -snes_linesearch_type basic -snes_divergence_tolerance 1e10 -ksp_rtol 1e-20
# PETSC_OPT=-snes_monitor -log_view -ksp_monitor -pc_type hypre -ksp_type bcgsl -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-10 -snes_max_it 10 -snes_atol 1e-7

run-%: set-% .force
	rm -rf run/
	mpirun -n $(N) $(topdir)/build/$*/ut-viscoplastic $(PETSC_OPT) 

gdb-mpi: set-devel
	mpiexec -np $(N) xterm -e gdb --ex 'set breakpoint pending on' --ex 'b MPI_Abort' --ex 'b breakpoint' --ex 'run' --args $(topdir)/build/devel/ut-viscoplastic $(PETSC_OPT)

gdb-%: set-% .force
	gdb --ex 'set breakpoint pending on' --ex 'b MPI_Abort' --ex 'b breakpoint' --ex 'run' --args $(topdir)/build/$*/ut-viscoplastic $(PETSC_OPT)

## MESH
mesh: model/cube-hole.e
model/cube-hole.e: model/cube-hole.msh set-opt .force
	$(topdir)/build/opt/msh2exo --in model/cube-hole.msh --out model/cube-hole.e
model/cube-hole.msh: model/mesh_single_well.py
	model/mesh_single_well.py 

model/cube-sphere.msh: model/cube-sphere.py
	model/cube-sphere.py 

model/cube-sphere.e: model/cube-sphere.msh
	$(topdir)/build/opt/msh2exo --in model/cube-sphere.msh --out model/cube-sphere.e

gdb: gdb-dbg

.force:
