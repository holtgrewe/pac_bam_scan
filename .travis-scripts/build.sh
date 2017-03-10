#!/bin/bash

set -exuo pipefail

export CXX="g++-4.9" CC="gcc-4.9"

pushd htslib
sed -i.bak -e 's~^prefix.*~prefix=./inst~' Makefile
make install
popd

mkdir -p build && pushd build
export HTSLIB_ROOT=${PWD}/../htslib/inst
cmake -DCMAKE_CXX_COMPILER=$(which $CXX) -DCMAKE_C_COMPILER=$(which $CC) -DCMAKE_BUILD_TYPE=Release ..
make VERBOSE=1

{ ./src/pac_bam_scan 2>&1 || true; } | grep Usage
