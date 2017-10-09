#!/bin/bash

mkdir -p bin
cd bin
ICC_PATH="$(which icc)"
ICPC_PATH="$(which icpc)"
export CC=${ICC_PATH}
export CXX=${ICPC_PATH}
cmake ..
make
