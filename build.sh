#!/bin/bash

echo "Start build"
if [[ -d "./build/" ]]; then
    rm -rf ./build/
fi
mkdir ./build/
echo "Build programme in build"
eval "cmake -B build"
echo "Compile programme in build"
cd build
eval "make" 2> "../build.log"
cd ..
echo "Build finish"
echo "See build.log"
