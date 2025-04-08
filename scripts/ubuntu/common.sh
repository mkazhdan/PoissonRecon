#!/bin/bash

set -e
set -u

SCRIPTS_DIR=$PWD
PROJECT_ROOT=$(realpath $SCRIPTS_DIR/../../)
PLATFORM=ubuntu
BUILD_DIR=$PROJECT_ROOT/build/$PLATFORM
INSTALL_DIR=$PROJECT_ROOT/install/$PLATFORM
BUILD_TYPE=Release
NUM_CPU=$(nproc)

mkdir -p $BUILD_DIR $INSTALL_DIR
