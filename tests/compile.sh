#!/bin/env sh

# TODO: Remove me!!!

CFLAGS="-I../src -I/usr/include/eigen3 -I../jsoncpp/include -L../src -L../jsoncpp/src -lpomerol -ljsoncpp -lhdf5 -lhdf5_cpp -fopenmp"
g++ green.cpp -o green ${CFLAGS}
g++ gamma4.cpp -o gamma4 ${CFLAGS}
