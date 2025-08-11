#!/usr/bin/env -S python3 

import gmsh
import sys

gmsh.initialize()
gmsh.model.add("three_layers")

lx, ly = 500.0, 500.0
z_zero = -200.0
tol = 1e-6

layers = [
    {"name": "CAPROCK",         "thickness": -z_zero},
    {"name": "RESERVOIR",       "thickness": 300.0},
#     {"name": "UNDERBURDEN",     "thickness": 20.0},
]

z_base = z_zero
volume_tags = []

for l in layers:
    tag = gmsh.model.occ.addBox(0, 0, z_base, lx, ly, l['thickness'])
    volume_tags.append((3, tag))
    l["z_top"] = z_base
    l["z_bottom"] = z_base + l['thickness']
    z_base += l['thickness']

gmsh.model.occ.fragment(volume_tags, [])
gmsh.model.occ.synchronize()

for l in layers:
    vols = gmsh.model.getEntitiesInBoundingBox(
        - tol, - tol, l["z_top"] - tol,
         lx + tol,  ly + tol, l["z_bottom"] + tol, 3)
    if not vols:
        raise RuntimeError(f"Failed to find volume for {l['name']}")
    gmsh.model.addPhysicalGroup(3, [v[1] for v in vols], name=l["name"])

faces = gmsh.model.getBoundary(gmsh.model.getEntities(3), oriented=False, recursive=False)

face_labels = {
    "XN": lambda c: abs(c[0]) < tol,
    "XP": lambda c: abs(c[0] - lx) < tol,
    "YN": lambda c: abs(c[1]) < tol,
    "YP": lambda c: abs(c[1] - ly) < tol,
    "ZN": lambda c: abs(c[2] - z_zero) < tol,
    "ZP": lambda c: abs(c[2] - z_base) < tol
}

face_groups = {name: [] for name in face_labels}

for dim, tag in faces:
    center = gmsh.model.occ.getCenterOfMass(dim, tag)
    for name, check in face_labels.items():
        if check(center):
            face_groups[name].append(tag)
            break

for name, tags in face_groups.items():
    if not tags:
        raise RuntimeError(f"Failed to find face group {name}")
    gmsh.model.addPhysicalGroup(2, tags, name=name)

# ADD LINE FOR THE SIZING FIELD
z1 = z_zero
z2 = z_base
pt1 = gmsh.model.occ.addPoint(0, 0, z1)
pt2 = gmsh.model.occ.addPoint(0, 0, z2)
line = gmsh.model.occ.addLine(pt1, pt2)

## SYNC ALL UPDATES
gmsh.model.occ.synchronize()

# gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 10)

target_size = 3
min_size = target_size
max_size = 50
dmax = 20.0

# Field 1: Refine near (0, 0, z) in a column
f_dist = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(f_dist, "EdgesList", [line])
f1 = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(
    f1,
    "F",
    f"{min_size} + ({max_size} - {min_size}) * (F{f_dist}/{dmax})^2"
)

# Field 2: Thin slab at caprock-reservoir interface
interface_z = 0
f2 = gmsh.model.mesh.field.add("Box")
gmsh.model.mesh.field.setNumber(f2, "VIn", 3)
gmsh.model.mesh.field.setNumber(f2, "VOut", 1000)
gmsh.model.mesh.field.setNumber(f2, "XMin", 0)
gmsh.model.mesh.field.setNumber(f2, "XMax",  100)
gmsh.model.mesh.field.setNumber(f2, "YMin", 0)
gmsh.model.mesh.field.setNumber(f2, "YMax",  100)
gmsh.model.mesh.field.setNumber(f2, "ZMin", interface_z - 5)
gmsh.model.mesh.field.setNumber(f2, "ZMax", interface_z + 5)

f3 = gmsh.model.mesh.field.add("Constant")
gmsh.model.mesh.field.setNumber(f3, "VIn", max_size)
gmsh.model.mesh.field.setNumber(f3, "VOut", max_size)

fmin = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", [f1, f2, f3])
gmsh.model.mesh.field.setAsBackgroundMesh(fmin)

gmsh.model.mesh.generate(3)

gmsh.write("model/cube.msh")
# gmsh.fltk.run()
# gmsh.finalize()
