#!/usr/bin/env sh

OFF_PATH=L.off

MAX_VOLUME=0.01

cd ../common

./make_L_shaped_room_mesh.py
./strip_off_comments.sh $OFF_PATH
./tetrahedralize.sh $OFF_PATH $MAX_VOLUME > tetgen_output.txt
./extract_verts_and_cells_from_tet_mesh.py --path=L.1.vtk

rm L.1.smesh L.1.vtk L.obj L.off

cd ~-

mv ../common/tetgen_output.txt .
mv ../common/L_{verts,cells}.bin .
