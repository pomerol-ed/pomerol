#!/usr/bin/env sh

# TODO: tutorial/

SOURCES="include/*.hpp          \
         include/*/*.hpp        \
         include/*/*.hpp.in     \
         src/*/*.cpp            \
         prog/*.hpp             \
         prog/*.cpp             \
         test/*.hpp             \
         test/*.cpp             \
         test/catch2/catch2-*"

clang-format --verbose -i $SOURCES
