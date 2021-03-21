#!/usr/bin/env sh

TETGEN=~/Build/tetgen1.6.0/tetgen

# Make surface meshes
echo "Making surface meshes"
./make_box_mesh.py
./make_L_shaped_room_mesh.py

# Since tetgen can't handle comments in OFF files, filter them out
echo "Stripping comments from OFF files"
./strip_off_comments.sh box.off
./strip_off_comments.sh L.off

# Create a tetrahedron mesh for each
echo "Generating tetrahedron meshes"
$TETGEN -q -k -C -B -N -E -F -a0.01 box.off > box_tetgen_output.txt
$TETGEN -q -k -C -B -N -E -F -a0.01 L.off > L_tetgen_output.txt

# Extract vertices and cells as binary files
echo "Writing vertex and cell arrays to disk as binary files"
./extract_verts_and_cells_from_tet_mesh.py --root=box > box_extract_arrays_output.txt
./extract_verts_and_cells_from_tet_mesh.py --root=L > L_extract_arrays_output.txt
