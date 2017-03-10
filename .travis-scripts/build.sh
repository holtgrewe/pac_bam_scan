#!/bin/bash

set -exuo pipefail

pushd htslib
sed -i.bak -e 's~^prefix.*~prefix=./inst~' Makefile
make install
popd

mkdir -p build && pushd build
export HTSLIB_ROOT=${PWD}/../htslib/inst
cmake -DCMAKE_BUILD_TYPE=Release ..
make

{ ./src/pac_bam_scan || true; } | grep Usage
