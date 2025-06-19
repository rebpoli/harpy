#!/usr/bin/env -S python3 -i

import gmsh
import sys

# Initialize Gmsh
gmsh.initialize(sys.argv)
gmsh.model.add("cube_with_sphere")

# Parameters
cube_size, sphere_radius = 10.0, 1
coarse_size, fine_size = 0.5, 0.5

# Create geometry
cube_tag = gmsh.model.occ.addBox(-cube_size/2, -cube_size/2, -cube_size/2, cube_size, cube_size, cube_size)
sphere_tag = gmsh.model.occ.addSphere(0, 0, 0, sphere_radius)

# Fragment to embed sphere in cube
gmsh.model.occ.fragment([(3, cube_tag)], [(3, sphere_tag)])
gmsh.model.occ.synchronize()

# Get volumes after boolean operation
volumes = gmsh.model.getEntities(3)

# Calculate volumes and sort to identify sphere and cube
volumes_data = []
for dim, tag in volumes:
    mass = gmsh.model.occ.getMass(dim, tag)
    volumes_data.append((tag, mass))

volumes_data.sort(key=lambda x: x[1])
sphere_tag, cube_tag = volumes_data[0][0], volumes_data[1][0]

# Assign physical groups to volumes
gmsh.model.addPhysicalGroup(3, [sphere_tag], name="TEST_FRAME")
gmsh.model.addPhysicalGroup(3, [cube_tag], name="SIDEBURDEN")

# Get the sphere-cube interface
sphere_boundaries = gmsh.model.getBoundary([(3, sphere_tag)], oriented=False)
interface_tags = [tag for dim, tag in sphere_boundaries]

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

# Simple sizing approach: use a ball field centered at origin
ball_field = gmsh.model.mesh.field.add("Ball")
gmsh.model.mesh.field.setNumber(ball_field, "Radius", sphere_radius * 1.3)  # Slightly larger than sphere
gmsh.model.mesh.field.setNumber(ball_field, "VIn", fine_size)  # Fine mesh size inside
gmsh.model.mesh.field.setNumber(ball_field, "VOut", coarse_size)  # Coarse mesh size outside
gmsh.model.mesh.field.setNumber(ball_field, "XCenter", 0)
gmsh.model.mesh.field.setNumber(ball_field, "YCenter", 0)
gmsh.model.mesh.field.setNumber(ball_field, "ZCenter", 0)

# Use this field to control mesh size
gmsh.model.mesh.field.setAsBackgroundMesh(ball_field)

# Set mesh size on the sphere-cube interface to ensure refinement
interface_vertices = []
for tag in interface_tags:
    boundary_vertices = gmsh.model.getBoundary([(2, tag)], oriented=False, recursive=True)
    interface_vertices.extend([t for d, t in boundary_vertices])

gmsh.model.mesh.setSize([(0, tag) for tag in set(interface_vertices)], fine_size)

# Generate mesh
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.model.mesh.generate(3)

# Save mesh and show GUI if requested
gmsh.write("model/cube-sphere.msh")
gmsh.fltk.run()
gmsh.finalize()

