#!/usr/bin/env -S python3 -i

import gmsh
import sys

# Initialize Gmsh
gmsh.initialize(sys.argv)
gmsh.model.add("cube")

# Parameters
cube_size = 10.0

# Create geometry
cube_tag = gmsh.model.occ.addBox(-cube_size/2, -cube_size/2, -cube_size/2, cube_size, cube_size, cube_size)
gmsh.model.occ.synchronize()

gmsh.model.addPhysicalGroup(3, [cube_tag], name="TEST_FRAME")

# Get and label cube faces
cube_boundaries = gmsh.model.getBoundary([(3, cube_tag)], oriented=False)
face_labels = {"XN": lambda c: abs(c[0] + cube_size/2) < 1e-6, "XP": lambda c: abs(c[0] - cube_size/2) < 1e-6,
               "YN": lambda c: abs(c[1] + cube_size/2) < 1e-6, "YP": lambda c: abs(c[1] - cube_size/2) < 1e-6,
               "ZN": lambda c: abs(c[2] + cube_size/2) < 1e-6, "ZP": lambda c: abs(c[2] - cube_size/2) < 1e-6}

# Apply labels to faces
for dim, tag in cube_boundaries:
    center = gmsh.model.occ.getCenterOfMass(2, tag)
    for name, check in face_labels.items():
        if check(center):
            gmsh.model.addPhysicalGroup(2, [tag], name=name)
            break

# Generate mesh
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1)
gmsh.model.mesh.generate(3)

# Save mesh and show GUI if requested
gmsh.write("model/cube.msh")
gmsh.fltk.run()
gmsh.finalize()

