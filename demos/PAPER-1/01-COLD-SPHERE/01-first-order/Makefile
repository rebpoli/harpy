
topdir := $(shell i=. ; for c in 0 1 2 3 4 5 6 7 8; do [ ! -f $$i/.chimas_top ] && i=$$i/..; done; echo $$i)#
TOPDIR=$(shell realpath $(topdir))

include $(topdir)/share/makefile/Makefile.demos

export CHDEBUG=0
export N=10

# PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -pc_type hypre -ksp_type bcgsl -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-10 -snes_max_it 10 -snes_atol 1e-7 -snes_linesearch_type basic -snes_divergence_tolerance 1e10 -ksp_rtol 1e-20 -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7 -pc_hypre_boomeramg_agg_nl 4 -pc_hypre_boomeramg_agg_num_paths 5 -pc_hypre_boomeramg_max_levels 25 -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_P_max 2 -pc_hypre_boomeramg_truncfactor 0.3 
PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -pc_type hypre -ksp_type bcgsl -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-4 -snes_max_it 20 -ksp_rtol 1e-15 -snes_linesearch_type basic -snes_atol 1e-7 -snes_divergence_tolerance 1e10 -ksp_rtol 1e-20 -ksp_atol 1e-15 -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7 -pc_hypre_boomeramg_agg_nl 4 -pc_hypre_boomeramg_agg_num_paths 5 -pc_hypre_boomeramg_max_levels 25 -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_P_max 2 -pc_hypre_boomeramg_truncfactor 0.3 -ksp_monitor

# PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -pc_type lu -ksp_type cg -ksp_gmres_restart 50 -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-4 -snes_max_it 8 -ksp_rtol 1e-15 -ksp_dtol 1e50 -ksp_atol 1e-15 -snes_linesearch_type basic -snes_atol 1e-2 -snes_divergence_tolerance 1e10 -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7 -pc_hypre_boomeramg_agg_nl 4 -pc_hypre_boomeramg_agg_num_paths 5 -pc_hypre_boomeramg_max_levels 25 -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_P_max 2 -pc_hypre_boomeramg_truncfactor 0.3 -ksp_max_it 5000 -ksp_monitor_true_residual 
# PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -pc_type hypre -ksp_type gmres -ksp_gmres_restart 50 -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-4 -snes_max_it 8 -ksp_rtol 1e-15 -ksp_dtol 1e50 -ksp_atol 1e-12 -snes_linesearch_type basic -snes_atol 1e-2 -snes_divergence_tolerance 1e10 -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7 -pc_hypre_boomeramg_agg_nl 4 -pc_hypre_boomeramg_agg_num_paths 5 -pc_hypre_boomeramg_max_levels 25 -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_P_max 2 -pc_hypre_boomeramg_truncfactor 0.3 -ksp_max_it 5000 -ksp_monitor_true_residual 

# PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -ksp_monitor -pc_type lu -ksp_type gmres -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-10 -snes_max_it 10 -snes_atol 1e-7 -snes_linesearch_type basic -snes_divergence_tolerance 1e10 -ksp_rtol 1e-20
# PETSC_OPT=-snes_monitor -log_view -ksp_monitor -pc_type hypre -ksp_type bcgsl -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-10 -snes_max_it 10 -snes_atol 1e-7

run-%: set-% .force
	rm -rf run/
	mpirun -n $(N) $(topdir)/build/$*/harpy $(PETSC_OPT) 

gdb-mpi: set-devel
	mpiexec -np $(N) xterm -e gdb --ex 'set breakpoint pending on' --ex 'b MPI_Abort' --ex 'b breakpoint' --ex 'run' --args $(topdir)/build/devel/harpy $(PETSC_OPT)

gdb-%: set-% .force
	gdb --ex 'set breakpoint pending on' --ex 'b MPI_Abort' --ex 'b breakpoint' --ex 'run' --args $(topdir)/build/$*/harpy $(PETSC_OPT)

## MESH
mesh: model/cube.e

# Mesh targets
model/cube-hole.e: model/cube-hole.msh set-opt .force
	$(topdir)/build/opt/msh2exo --in model/cube-hole.msh --out model/cube-hole.e
model/cube-hole.msh: model/mesh_single_well.py
	model/mesh_single_well.py 

model/cube-sphere.msh: model/cube-sphere.py
	model/cube-sphere.py 
model/cube-sphere.e: model/cube-sphere.msh
	$(topdir)/build/opt/msh2exo --in model/cube-sphere.msh --out model/cube-sphere.e

model/cube.e: model/cube.msh set-opt .force
	$(topdir)/build/opt/msh2exo --in model/cube.msh --out model/cube.e
model/cube.msh: model/cube.py
	model/cube.py 

gdb: gdb-dbg

plot: .force
	py/plot_lin0.py

.force:
