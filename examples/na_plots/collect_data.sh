#!/bin/sh

RUN_EIK2G1_PATH=../../build/run_eik2g1
RFAC=0.3
OUTDIR=./out

mkdir -p $OUTDIR

exe() { echo "$@"; "$@"; }

for P in 7 8 9 10 11
do
    N=$((2**P + 1))
    exe $RUN_EIK2G1_PATH $N $RFAC \
        $OUTDIR/jet_gt_$N.bin \
        $OUTDIR/jet_$N.bin \
        $OUTDIR/l_$N.bin \
        $OUTDIR/lam_$N.bin \
	$OUTDIR/sol_infos_for_par_$N.txt
done
