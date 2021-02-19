#!/usr/bin/env sh

C_COMPILER=clang
CXX_COMPILER=clang++
BUILD_TYPE=Debug
SANITIZE_ADDRESS=OFF
SANITIZE_MEMORY=OFF
SANITIZE_UNDEFINED=ON

echo "C_COMPILER:         ${C_COMPILER}"
echo "CXX_COMPILER:       ${CXX_COMPILER}"
echo "BUILD_TYPE:         ${BUILD_TYPE}"
echo "SANITIZE_ADDRESS:   ${SANITIZE_ADDRESS}"
echo "SANITIZE_MEMORY:    ${SANITIZE_MEMORY}"
echo "SANITIZE_UNDEFINED: ${SANITIZE_UNDEFINED}"

set -v
mkdir -p build
cd build
cmake -DCMAKE_C_COMPILER=$C_COMPILER \
      -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      -DSANITIZE_ADDRESS=$SANITIZE_ADDRESS \
      -DSANITIZE_MEMORY=$SANITIZE_MEMORY \
      -DSANITIZE_UNDEFINED=$SANITIZE_UNDEFINED \
      -DBUILD_TESTS=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      .. || exit 1
make -j || exit 1
mv compile_commands.json ..
cd ..
python setup.py build_ext || exit 1
pip uninstall -y jmm
python setup.py install || exit 1
