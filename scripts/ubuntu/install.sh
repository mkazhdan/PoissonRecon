#!/bin/bash

source common.sh

cmake --build $BUILD_DIR --config $BUILD_TYPE --target install -j $NUM_CPU
