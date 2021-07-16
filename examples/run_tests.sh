#!/usr/bin/env bash

L_FIRST=0
L_LAST=19

run_mult_arr_job() {
	echo $1
	(python testing_multiple_arrivals.py $1) > ./out/$1.txt 2>&1
}

export -f run_mult_arr_job

mkdir -p ./out

parallel --linebuffer run_mult_arr_job ::: $(seq $L_FIRST $L_LAST)

echo "done"
