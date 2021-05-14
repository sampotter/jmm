#!/bin/bash

# This is a simple one-off script to quickly rebuild everything after
# making some changes (anywhere!) in the code. The idea here is to
# waste a little time at the expense of not having to think too hard
# about what I need to do to ensure my changes propagate after
# rebuilding. Right now, this is platform-specific and really only
# appropriate for debugging.

set -v
cd build
make -j || exit 1
rm -r lib.macosx-10.15-x86_64-3.7
rm -r temp.macosx-10.15-x86_64-3.7
cd -
python setup.py build || exit 1
pip uninstall -y jmm
python setup.py install || exit 1
