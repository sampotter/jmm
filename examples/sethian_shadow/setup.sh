#!/bin/bash

SCENE=room

MAX_VOLUME=0.1

cd ../common

./make_bounding_mesh.py $SCENE.obj
./strip_off_comments.sh ${SCENE}_bd.off
./tetrahedralize.sh ${SCENE}_bd.off $MAX_VOLUME > tetgen_output.txt
./sort_bd_vtu.py ${SCENE}_bd.1.vtk ${SCENE}.obj ${SCENE}.vtu ${SCENE}.json

rm ${SCENE}_bd.off ${SCENE}_bd.1.smesh ${SCENE}_bd.1.vtk

cd ~-

cp ../common/$SCENE.obj .
mv ../common/$SCENE.{vtu,json} .
mv ../common/tetgen_output.txt .
