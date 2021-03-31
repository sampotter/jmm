#!/bin/bash

SCENE=L

NUM_VERTS=`python -c "import numpy as np; V = np.fromfile('${SCENE}_verts.bin'); print(V.size//3)"`

for (( l = 0 ; l < $NUM_VERTS; l++ )) ; do
    echo "l = $l"
    ./visualize_cutset.py -l $l -p srcs/${SCENE}/$l -s ${SCENE} # || exit 1
done
