#!/usr/bin/env sh

set -v
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON .. || exit 1
make -j || exit 1
mv compile_commands.json ..
cd ..
python setup.py build_ext || exit 1
pip uninstall -y jmm
python setup.py install || exit 1
