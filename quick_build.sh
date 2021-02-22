#!/usr/bin/env bash

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    C_COMPILER=clang
else
    C_COMPILER=" "
fi

BUILD_TYPE=Debug
SANITIZE_ADDRESS=OFF
SANITIZE_MEMORY=OFF
SANITIZE_UNDEFINED=OFF

echo "C_COMPILER:         ${C_COMPILER}"
echo "BUILD_TYPE:         ${BUILD_TYPE}"
echo "SANITIZE_ADDRESS:   ${SANITIZE_ADDRESS}"
echo "SANITIZE_MEMORY:    ${SANITIZE_MEMORY}"
echo "SANITIZE_UNDEFINED: ${SANITIZE_UNDEFINED}"

CMAKE_ARGS="-DCMAKE_BUILD_TYPE=$BUILD_TYPE"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    CMAKE_ARGS="${CMAKE_ARGS} -DCMAKE_C_COMPILER=$C_COMPILER"
fi
CMAKE_ARGS="${CMAKE_ARGS} -DSANITIZE_ADDRESS=$SANITIZE_ADDRESS"
CMAKE_ARGS="${CMAKE_ARGS} -DSANITIZE_MEMORY=$SANITIZE_MEMORY"
CMAKE_ARGS="${CMAKE_ARGS} -DSANITIZE_UNDEFINED=$SANITIZE_UNDEFINED"
CMAKE_ARGS="${CMAKE_ARGS} -DBUILD_TESTS=ON"
CMAKE_ARGS="${CMAKE_ARGS} -DCMAKE_EXPORT_COMPILE_COMMANDS=ON"

echo $CMAKE_ARGS

set -v
mkdir -p build
cd build
cmake $CMAKE_ARGS .. || exit 1
make -j || exit 1
mv compile_commands.json ..
cd ..
python setup.py build_ext || exit 1
if pip list | grep jmm > /dev/null; then
    pip uninstall -y jmm
fi
python setup.py install || exit 1
