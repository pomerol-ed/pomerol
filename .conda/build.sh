#!/usr/bin/env bash

mkdir build
cd build

# Build pomerol
cmake \
    -DCMAKE_CXX_COMPILER=${BUILD_PREFIX}/bin/$(basename ${CXX}) \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DUSE_OPENMP=ON \
    -DDocumentation=OFF \
    -DProgs=ON \
    -Dprogs_list='anderson;hubbard2d' \
    ..

make -j2 VERBOSE=1

# Run unit tests
PATH="$(pwd)/test:$PATH" ctest --output-on-failure

# Install pomerol
make install
