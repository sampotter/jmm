#!/usr/bin/env bash

set -v
cd build
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cgreen-runner libjmm_tests.so
elif [[ "$OSTYPE" == "darwin"* ]]; then
    cgreen-runner libjmm_tests.dylib
else
    echo "Unknown platform: $OSTYPE" 1>&2
    exit 1
fi
