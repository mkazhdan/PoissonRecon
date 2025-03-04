#!/bin/bash

source common.sh

cmake --build $BUILD_DIR --config $BUILD_TYPE -j $NUM_CPU
