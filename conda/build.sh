#!/bin/bash

set -euxo pipefail

mkdir build
pushd build

CLP_ROOT=${PREFIX} \
ARMADILLO_ROOT=${BUILD_PREFIX} \
CXXFLAGS="-I${BUILD_PREFIX}/include $CXXFLAGS" \
CMAKE_PREFIX_PATH=${BUILD_PREFIX} \
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DCMAKE_PLATFORM_FLAGS=$CMAKE_PLATFORM_FLAGS \
      ..

make -j${CPU_COUNT}
make install
make clean

popd
