[![DOI](https://zenodo.org/badge/4569/aeantipov/pomerol.svg)](
http://dx.doi.org/10.5281/zenodo.17900
)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-red)](
https://aeantipov.github.io/pomerol/)

**pomerol** is an exact diagonalization (full ED) code written in C++ aimed at
solving condensed matter second-quantized models of interacting fermions
and bosons on finite size lattices at finite temperatures.
It is designed to compute thermal expectation values of observables, single- and
two-particle Green's functions as well as susceptibilities.

##  Features

  * High performance exact calculation of Green's functions, two-particle
    Green's functions and susceptibilities in Matsubara domain.
  * Many-body Hamiltonians can be specified in a natural mathematical form using
    [libcommute's](https://krivenko.github.io/libcommute/) Domain-Specific
    Language. Hamiltonian presets for commonly used lattice models are also
    available.
  * Automatic symmetry analysis of the many-body Hamiltonians drastically
    reduces computational costs.
  * [Eigen 3](http://eigen.tuxfamily.org) template library is used for numerical
    linear algebra.
  * [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface) +
    [OpenMP](https://en.wikipedia.org/wiki/OpenMP) support.

## Installation
### From source

  - Check the *dependencies*:

    * A C++11 conformant compiler
    * CMake >= 3.1.0
    * Boost >= 1.54.0 (only headers are required)
    * Eigen >= 3.1.0
    * [libcommute >= 0.7](https://github.com/krivenko/libcommute)
    * An MPI 3.0 implementation
    * Git to fetch the sources

  - Download the latest sources:

    ```
    git clone https://github.com/aeantipov/pomerol.git
    ```

  - Create a (temporary) build directory and change to it:

    ```
    mkdir build && cd build
    ```

  - In this build directory, run

    ```
    cmake <path_to_pomerol_sources> -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<installation_path>
    ```
    * If CMake configuration file of libcommute cannot be found, add
      `-Dlibcommute_DIR=<libcommute_installation_path>/lib/cmake` to the command line.
    * Add `-DTesting=OFF` to disable compilation of unit tests (not recommended).
    * Add `-DProgs=ON` to compile provided executables (from `progs`
      directory). Some of the executables depend on the
      [gftools](https://github.com/aeantipov/gftools) library, which will be
      automatically downloaded in case it cannot be found by CMake (use
      `-Dgftools_DIR` to specify its installation path). gftools supports saving
      to HDF5 through [ALPSCore](http://alpscore.org).
    * Add `-DDocumentation=OFF` to disable generation of reference
      documentation.
    * Add `-DUSE_OPENMP=OFF` to disable OpenMP optimization for two-particle GF
      calculation.
    * Add `-DBUILD_SHARED_LIBS=OFF` to compile static instead of shared libraries.
  - `make`
  - `make test` (if unit tests are compiled)
  - `make install`
  - `make doc` generates the Doxygen reference documentation in the `doc/html`
    subdirectory.

The library, _libpomerol_ is built. It can be used for linking with executables.
Some working executables are given in `prog` subdirectory.

## Interfacing with your own code and other libraries

Check the `tutorial` directory for an example of a pomerol-based code that is
linked to external libraries.

The interface to [TRIQS library](https://triqs.github.io/triqs/latest/) is
readily available: <https://github.com/krivenko/pomerol2triqs>.

# Documentation
Check <https://aeantipov.github.io/pomerol/html/> or type `make doc` during
compilation stage to build the reference documentation.

[V2-NEWS.md](V2-NEWS.md) lists main changes introduced in pomerol 2.0
compared to the 1.x branch.

# License
This Source Code Form is subject to the terms of the Mozilla Public License,
v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain
one at <http://mozilla.org/MPL/2.0/>.

# Academic usage

Please, attribute this work by a citation to
<http://dx.doi.org/10.5281/zenodo.17900>.

# Authors & Contributors
  * Andrey Antipov <Andrey.E.Antipov\at\gmail.com>
  * Igor Krivenko <igor.s.krivenko\at\gmail.com>
  * Mikhail Alejnikov
  * Alexey Rubtsov
  * Christoph Jung
  * Aljoscha Wilhelm
  * Junya Otsuki
  * Sergei Iskakov
  * Hiroshi Shinaoka
  * Nils Wentzell
  * Hugo U.R. Strand

We acknowledge [NRC Kurchatov Institute](http://eng.nrcki.ru/) for providing
computational resources.

# Development/Help
Please, feel free to contact us and to contribute!
