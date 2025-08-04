#!/usr/bin/env -S python3 -i

import gmsh
import sys, os

# Initialize and set up
gmsh.initialize()
gmsh.model.add("cube-hole")

# Create geometry
cube = gmsh.model.occ.addBox(-0.1, -0.5, -0.5, .2, 1, 1)
cylinder = gmsh.model.occ.addCylinder(-0.6, 0, 0, 1.2, 0, 0, 0.2)
out, _ = gmsh.model.occ.cut([(3, cube)], [(3, cylinder)])
gmsh.model.occ.synchronize()

# Get all surfaces from the resulting volume
vol_tags = [tag for dim, tag in out if dim == 3]
boundary_dimtags = gmsh.model.getBoundary([(3, vol_tags[0])], oriented=False)

# Identify and name the surfaces based on their center coordinates
for dimtag in boundary_dimtags:
    # Get center of mass of the surface
    com = gmsh.model.occ.getCenterOfMass(dimtag[0], dimtag[1])

    # Assign physical name based on position
    if abs(com[0] + 0.1) < 1e-6:
        gmsh.model.addPhysicalGroup(2, [dimtag[1]], name="XN")  # X-negative (left)
    elif abs(com[0] - 0.1) < 1e-6:
        gmsh.model.addPhysicalGroup(2, [dimtag[1]], name="XP")  # X-positive (right)
    elif abs(com[1] + 0.5) < 1e-6:
        gmsh.model.addPhysicalGroup(2, [dimtag[1]], name="YN")  # Y-negative (front)
    elif abs(com[1] - 0.5) < 1e-6:
        gmsh.model.addPhysicalGroup(2, [dimtag[1]], name="YP")  # Y-positive (back)
    elif abs(com[2] + 0.5) < 1e-6:
        gmsh.model.addPhysicalGroup(2, [dimtag[1]], name="ZN")  # Z-negative (bottom)
    elif abs(com[2] - 0.5) < 1e-6:
        gmsh.model.addPhysicalGroup(2, [dimtag[1]], name="ZP")  # Z-positive (top)
    else:
        # This is the cylindrical surface
        gmsh.model.addPhysicalGroup(2, [dimtag[1]], name="WELL")
        cylinder_surface = dimtag[1]

# Add physical group for the volume
gmsh.model.addPhysicalGroup(3, vol_tags, name="TEST_FRAME")

min_mesh_size = 0.04  # Fine mesh size near cylinder
max_mesh_size = 0.2   # Coarse mesh size far from cylinder
max_dist = .2  # Controls how quickly mesh size increases with distance

print (f"Cylinder surface: {cylinder_surface}")
distance_field = gmsh.model.mesh.field.add("Distance")
# gmsh.model.mesh.field.setNumbers(distance_field, "SurfacesList", [cylinder_surface])
threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
gmsh.model.mesh.field.setNumber(threshold_field, "LcMin", min_mesh_size)
gmsh.model.mesh.field.setNumber(threshold_field, "LcMax", max_mesh_size)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 0)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", max_dist)

math_field = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(math_field, "F", "sqrt(y*y+z*z)-.25") # .2 is the cylinder radius

# Use this for the threshold instead
gmsh.model.mesh.field.setNumber(threshold_field, "IField", math_field)


## USE THIS FOR REFINEMENT
# gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

# Generate mesh and save
gmsh.model.mesh.generate(3)
mshfile = os.path.dirname(__file__) + "/cube-hole.msh"
print(f"Writing file: {mshfile}");

gmsh.write(mshfile)

# Optional visualization
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
