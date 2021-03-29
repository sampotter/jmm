#!/bin/bash

NUM_VERTS=`python -c "import numpy as np; V = np.fromfile('L_verts.bin'); print(V.size//3)"`

for (( l = 0 ; l < $NUM_VERTS; l++ )) ; do
    echo "l = $l"
    ./make_L_video.py -l$l -m0.5 -M2.0 -n7 --outpath=srcs/$l || exit 1
done
