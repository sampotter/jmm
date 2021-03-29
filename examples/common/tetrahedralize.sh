#!/bin/bash

TETGEN=~/Build/tetgen1.6.0/tetgen

$TETGEN -q -k -C -B -N -E -F -a$2 $1
