#!/bin/bash

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    C_COMPILER=clang
else
    C_COMPILER=" "
fi

: ${BUILD_TYPE:=Debug}
: ${SANITIZE_ADDRESS:=OFF}
: ${SANITIZE_MEMORY:=OFF}
: ${SANITIZE_UNDEFINED:=OFF}

is_on_or_off() {
    case $2 in
        (ON|OFF) ;;
        (*) printf >&2 "$1 should be ON or OFF, got $2 instead\n";
            exit 1 ;;
    esac
}

is_on_or_off SANITIZE_ADDRESS "${SANITIZE_ADDRESS}"
is_on_or_off SANITIZE_MEMORY "${SANITIZE_MEMORY}"
is_on_or_off SANITIZE_UNDEFINED "${SANITIZE_UNDEFINED}"

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
python setup.py build -j`nproc` || exit 1
if pip list | grep jmm > /dev/null; then
    pip uninstall -y jmm
fi
python setup.py install || exit 1
