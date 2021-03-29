#!/usr/bin/env sh

if [ -z "$REBUILD" ] ; then
   while true; do
       read -p "Do you want to rebuild jmm first? (y/N) " yn
       case $yn in
           [Yy]* )
               rm -rf build
               SANITIZE_ADDRESS=ON ./quick_build.sh
               break;;
           * )
               break;;
       esac
   done
else
    rm -rf build
    SANITIZE_ADDRESS=ON ./quick_build.sh
fi

cd examples/make_L_video

./setup.sh

NUM_VERTS=`python -c "import numpy as np; V = np.fromfile('L_verts.bin'); print(V.size//3)"`

cd ~-

mkdir -p stress_test_output

for (( l = 0 ; l < $NUM_VERTS; l++ )) ; do
    echo "l = $l"
    ./build/basic_solve \
        examples/make_L_video/L_verts.bin \
        examples/make_L_video/L_cells.bin \
        $l \
        stress_test_output/L_jets.$l.bin \
        || exit 1
done
