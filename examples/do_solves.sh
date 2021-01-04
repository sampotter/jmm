#!/usr/bin/env bash

TETGEN=~/Build/tetgen1.6.0/tetgen

TETGEN_FLAGS="-q -k -C -B -N -E -F"

TETRAHEDRON_VOLUMES="a.txt"

SCENE=$1
OUTDIR=$2
MASK=$3

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <scene> <outdir> <mask>"
    exit 0
fi

while IFS= read -r line
do
    echo "Running ${SCENE} test with a = ${line}:"

    OUTPATH=$OUTDIR/$line
    mkdir -p $OUTPATH

    $TETGEN $TETGEN_FLAGS -a$line $SCENE.off > $OUTPATH/tetgen_output.txt
    rm $SCENE.1.smesh
    mv $SCENE.1.vtk $OUTPATH/$SCENE.1.vtk

    ROOT=$OUTPATH/$SCENE

    echo "  ./extract_verts_and_cells_from_tet_mesh.py --root=$ROOT"
    ./extract_verts_and_cells_from_tet_mesh.py --root=$ROOT | sed 's/^/    /'

    echo "  ./basic_solve.py --path=$OUTPATH --scene=$SCENE --mask=$MASK"
    ./basic_solve.py --path=$OUTPATH --scene=$SCENE --mask=$MASK | sed 's/^/    /'

    echo ""
done < $TETRAHEDRON_VOLUMES
