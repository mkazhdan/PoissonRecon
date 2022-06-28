#!/bin/bash

source common.sh

set -x

cmake \
    -S $REPO_ROOT \
    -B $BUILD_DIR \
    -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE
