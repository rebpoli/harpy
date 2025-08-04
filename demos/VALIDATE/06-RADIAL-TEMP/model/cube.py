#!/usr/bin/env -S python3 -i

import gmsh
import sys

# Initialize Gmsh
gmsh.initialize(sys.argv)
gmsh.model.add("cube")

# Parameters
cube_width  = 200
cube_height = 200

# Create geometry
cube_tag = gmsh.model.occ.addBox(0, 0, -cube_height/2, cube_width/2, cube_width/2, cube_height)
gmsh.model.occ.synchronize()

gmsh.model.addPhysicalGroup(3, [cube_tag], name="TEST_FRAME")

# Get and label cube faces
cube_boundaries = gmsh.model.getBoundary([(3, cube_tag)], oriented=False)


face_labels = {
    "XN": lambda c: abs(c[0] - 0) < 1e-6,                       # X = 0 plane (Symmetry)
    "XP": lambda c: abs(c[0] - cube_width / 2) < 1e-6,          # X = +cube_width/2
    "YN": lambda c: abs(c[1] - 0) < 1e-6,                       # Y = 0 plane (Symmetry)
    "YP": lambda c: abs(c[1] - cube_width / 2) < 1e-6,          # Y = +cube_width/2
    "ZN": lambda c: abs(c[2] + cube_height / 2) < 1e-6,         # Z = -cube_height/2
    "ZP": lambda c: abs(c[2] - cube_height / 2) < 1e-6          # Z = +cube_height/2
}

# Apply labels to faces
for dim, tag in cube_boundaries:
    center = gmsh.model.occ.getCenterOfMass(2, tag)
    for name, check in face_labels.items():
        if check(center):
            gmsh.model.addPhysicalGroup(2, [tag], name=name)
            break

# Generate mesh
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 10)
gmsh.model.mesh.generate(3)

# Save mesh and show GUI if requested
gmsh.write("model/cube.msh")
gmsh.fltk.run()
gmsh.finalize()

