#!/bin/bash

if [ -d "build" ]; then
  rm -rf "build"
  echo "Deleting the build directory"
fi

if [ -f "build.log" ]; then
  rm "build.log"
fi
