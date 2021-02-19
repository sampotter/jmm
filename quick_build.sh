#!/usr/bin/env sh

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    C_COMPILER=clang
    CXX_COMPILER=clang++
else
    C_COMPILER=" "
    CXX_COMPILER=" "
fi

BUILD_TYPE=Debug
SANITIZE_ADDRESS=ON
SANITIZE_MEMORY=OFF
SANITIZE_UNDEFINED=OFF

echo "C_COMPILER:         ${C_COMPILER}"
echo "CXX_COMPILER:       ${CXX_COMPILER}"
echo "BUILD_TYPE:         ${BUILD_TYPE}"
echo "SANITIZE_ADDRESS:   ${SANITIZE_ADDRESS}"
echo "SANITIZE_MEMORY:    ${SANITIZE_MEMORY}"
echo "SANITIZE_UNDEFINED: ${SANITIZE_UNDEFINED}"

CMAKE_ARGS="-DCMAKE_BUILD_TYPE=$BUILD_TYPE"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    CMAKE_ARGS="${CMAKE_ARGS} -DCMAKE_C_COMPILER=$C_COMPILER"
    CMAKE_ARGS="${CMAKE_ARGS} -DCMAKE_CXX_COMPILER=$CXX_COMPILER"
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
pip uninstall -y jmm
python setup.py install || exit 1
