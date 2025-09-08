#!/usr/bin/env -S python3 

import gmsh
import sys

gmsh.initialize()
gmsh.model.add("cube")

# Parameters
lc = 1.0   # mesh size characteristic length
L = 1.0    # cube side length

# Define corner points of the cube
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(L, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(L, L, 0, lc)
p4 = gmsh.model.geo.addPoint(0, L, 0, lc)

p5 = gmsh.model.geo.addPoint(0, 0, L, lc)
p6 = gmsh.model.geo.addPoint(L, 0, L, lc)
p7 = gmsh.model.geo.addPoint(L, L, L, lc)
p8 = gmsh.model.geo.addPoint(0, L, L, lc)

# Define lines
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

l5 = gmsh.model.geo.addLine(p5, p6)
l6 = gmsh.model.geo.addLine(p6, p7)
l7 = gmsh.model.geo.addLine(p7, p8)
l8 = gmsh.model.geo.addLine(p8, p5)

l9 = gmsh.model.geo.addLine(p1, p5)
l10 = gmsh.model.geo.addLine(p2, p6)
l11 = gmsh.model.geo.addLine(p3, p7)
l12 = gmsh.model.geo.addLine(p4, p8)

# Define surfaces (as curve loops)
bottom = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
top    = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
front  = gmsh.model.geo.addCurveLoop([l1, l10, -l5, -l9])
back   = gmsh.model.geo.addCurveLoop([l3, l12, -l7, -l11])
left   = gmsh.model.geo.addCurveLoop([l4, l9, -l8, -l12])
right  = gmsh.model.geo.addCurveLoop([l2, l11, -l6, -l10])

s_bottom = gmsh.model.geo.addPlaneSurface([bottom])
s_top    = gmsh.model.geo.addPlaneSurface([top])
s_front  = gmsh.model.geo.addPlaneSurface([front])
s_back   = gmsh.model.geo.addPlaneSurface([back])
s_left   = gmsh.model.geo.addPlaneSurface([left])
s_right  = gmsh.model.geo.addPlaneSurface([right])

# Define surface loop and volume
sl = gmsh.model.geo.addSurfaceLoop([s_bottom, s_top, s_front, s_back, s_left, s_right])
vol = gmsh.model.geo.addVolume([sl])

# Synchronize before adding physical groups
gmsh.model.geo.synchronize()

# Assign physical groups (labels)
gmsh.model.addPhysicalGroup(2, [s_left], name="XN")
gmsh.model.addPhysicalGroup(2, [s_right], name="XP")
gmsh.model.addPhysicalGroup(2, [s_front], name="YN")
gmsh.model.addPhysicalGroup(2, [s_back], name="YP")
gmsh.model.addPhysicalGroup(2, [s_bottom], name="ZN")
gmsh.model.addPhysicalGroup(2, [s_top], name="ZP")

gmsh.model.addPhysicalGroup(3, [vol], name="RESERVOIR")

# gmsh.option.setNumber( "Mesh.MeshSizeMax", 2 )
gmsh.option.setNumber( "Mesh.CharacteristicLengthMax", 0.3 )
# Generate 3D mesh
gmsh.model.mesh.generate(3)

# Count number of 3D elements
eleTypes, eleTags, eleNodes = gmsh.model.mesh.getElements(dim=3)
numElems = sum(len(tags) for tags in eleTags)
print(f"Number of 3D elements: {numElems}")

gmsh.write("model/cube.msh")

gmsh.fltk.run()
gmsh.finalize()
