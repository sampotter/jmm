#!/usr/bin/env bash

TETGEN=~/Build/tetgen1.6.0/tetgen

TETGEN_FLAGS="-q -k -C -B -N -E -F"

TETRAHEDRON_VOLUMES="a.txt"

OUTDIR=$1

MASK=$2

while IFS= read -r line
do
    echo "Running box test with a = ${line}:"

    OUTPATH=$OUTDIR/$line
    mkdir -p $OUTPATH

    $TETGEN $TETGEN_FLAGS -a$line box.off > $OUTPATH/tetgen_output.txt
    rm box.1.smesh
    mv box.1.vtk $OUTPATH/box.1.vtk

    ROOT=$OUTPATH/box

    echo "  ./extract_verts_and_cells_from_tet_mesh.py --root=$ROOT"
    ./extract_verts_and_cells_from_tet_mesh.py --root=$ROOT | sed 's/^/    /'

    echo "  ./solve_box.py --root=$ROOT --mask=$MASK"
    ./solve_box.py --root=$ROOT --mask=$MASK | sed 's/^/    /'

    echo ""
done < $TETRAHEDRON_VOLUMES
