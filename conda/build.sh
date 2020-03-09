#!/bin/bash

set -euxo pipefail

declare -a CMAKE_PLATFORM_FLAGS
if [[ ${HOST} =~ .*darwin.* ]]; then
  CMAKE_PLATFORM_FLAGS+=(-DCMAKE_OSX_SYSROOT="${CONDA_BUILD_SYSROOT}")
else
  CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${RECIPE_DIR}/cross-linux.cmake")
fi

mkdir build
pushd build

GLPK_ROOT=${BUILD_PREFIX} \
ARMADILLO_ROOT=${BUILD_PREFIX} \
CXXFLAGS="-I${BUILD_PREFIX}/include $CXXFLAGS" \
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DCMAKE_OSX_SYSROOT=${CONDA_BUILD_SYSROOT} ..

make -j${CPU_COUNT}
make install
make clean

popd