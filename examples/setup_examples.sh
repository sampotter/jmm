#!/usr/bin/env sh

TETGEN=~/Build/tetgen1.6.0/tetgen

# Make surface meshes
./make_box_mesh.py
./make_L_shaped_room_mesh.py

# Since tetgen can't handle comments in OFF files, filter them out
sed -i '' '/^#/d' box.off
sed -i '' '/^#/d' L.off

# Create a tetrahedron mesh for each
$TETGEN -q -k -C -B -N -E -F -a0.01 box.off > box_tetgen_output.txt
$TETGEN -q -k -C -B -N -E -F -a0.01 L.off > L_tetgen_output.txt
