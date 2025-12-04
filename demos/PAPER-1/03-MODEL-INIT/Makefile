
topdir := $(shell i=. ; for c in 0 1 2 3 4 5 6 7 8; do [ ! -f $$i/.chimas_top ] && i=$$i/..; done; echo $$i)#
TOPDIR=$(shell realpath $(topdir))

include $(topdir)/share/makefile/Makefile.demos

export CHDEBUG=0
export N=8

# PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -pc_type hypre -ksp_type bcgsl -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-10 -snes_max_it 10 -snes_atol 1e-7 -snes_linesearch_type basic -snes_divergence_tolerance 1e10 -ksp_rtol 1e-20 -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7 -pc_hypre_boomeramg_agg_nl 4 -pc_hypre_boomeramg_agg_num_paths 5 -pc_hypre_boomeramg_max_levels 25 -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_P_max 2 -pc_hypre_boomeramg_truncfactor 0.3 
PETSC_OPT=-snes_monitor -snes_force_iteration -log_view -pc_type hypre -ksp_type bcgsl -snes_converged_reason -ksp_converged_reason -snes_rtol 1e-5 -snes_max_it 12 -ksp_rtol 1e-15 -snes_linesearch_type basic -snes_atol 1e-7 -snes_divergence_tolerance 1e10 -ksp_rtol 1e-20 -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7 -pc_hypre_boomeramg_agg_nl 4 -pc_hypre_boomeramg_agg_num_paths 5 -pc_hypre_boomeramg_max_levels 25 -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_P_max 2 -pc_hypre_boomeramg_truncfactor 0.3 -ksp_max_it 1000 # -ksp_monitor_true_residual 

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

video:
	ffmpeg -framerate 10 -pattern_type glob -i 'png-linears/*.png' -c:v libx264 -pix_fmt yuv420p -crf 18 output.mp4

launch-athena: export HARPY_CMD=$(topdir)/build/opt/harpy $(PETSC_OPT) 
launch-athena: .force
	sbatch ../py/run_atena.sh


.force:
