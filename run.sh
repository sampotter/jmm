#!/usr/bin/env sh

OUTPUT_DIR=./data

mkdir -p $OUTPUT_DIR

for p in `seq 7 12`; do
    N=$((2**$p + 1))
    echo "N = 2^$p + 1 = $N"
    mkdir -p $OUTPUT_DIR/N$N
    ./scratch $N
    mv {T,Tx,Ty,Txy}.npy $OUTPUT_DIR/N$N
done
