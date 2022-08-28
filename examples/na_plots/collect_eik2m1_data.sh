#!/bin/bash

RUN_EIK2M1_PATH=../../build/run_eik2m1
RFAC=0.1
OUTDIR=./out

mkdir -p $OUTDIR

exe() { echo "$@"; "$@"; }

while read MAX_AREA
do
    exe ./make_mesh.py \
        $MAX_AREA \
        $OUTDIR/verts_$MAX_AREA.bin \
        $OUTDIR/faces_$MAX_AREA.bin \
        $OUTDIR/triplot_$MAX_AREA.pdf

    exe $RUN_EIK2M1_PATH \
        $OUTDIR/verts_$MAX_AREA.bin \
        $OUTDIR/faces_$MAX_AREA.bin \
        $RFAC \
        $OUTDIR/eik2m1_jets_gt_$MAX_AREA.bin \
        $OUTDIR/eik2m1_jets_$MAX_AREA.bin \
        $OUTDIR/eik2m1_l_$MAX_AREA.bin \
        $OUTDIR/eik2m1_lam_$MAX_AREA.bin
done < max_areas.txt
