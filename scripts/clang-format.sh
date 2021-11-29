#!/usr/bin/env sh
#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# TODO: tutorial/

SOURCES="include/*.hpp          \
         include/*/*.hpp        \
         include/*/*.hpp.in     \
         src/*/*.cpp            \
         prog/*.hpp             \
         prog/*.cpp             \
         test/*.hpp             \
         test/*.cpp             \
         test/catch2/catch2-*   \
         tutorial/*.cpp"

clang-format --verbose -i $SOURCES
