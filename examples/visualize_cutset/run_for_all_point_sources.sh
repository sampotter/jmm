#!/usr/bin/env sh

NUM_VERTS=`python -c "import numpy as np; V = np.fromfile('L_verts.bin'); print(V.size//3)"`

for (( l = 0 ; l < $NUM_VERTS; l++ )) ; do
    echo "l = $l"
    ./visualize_cutset.py -l$l --outpath=srcs/$l || exit 1
done
