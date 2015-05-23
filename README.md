[![DOI](https://zenodo.org/badge/4569/aeantipov/pomerol.svg)](http://dx.doi.org/10.5281/zenodo.17900)

**pomerol** is an exact diagonalization (full-ED) code written in C++ aimed at solving condensed matter second-quantized models of interacting fermions on finite size lattices at finite temperatures. It is designed to produce single and two-particle Greens functions.

##### Documentation
Check http://pomerol.sourceforge.net or type `make doc` during compilation stage for the reference documentation.

The library, _libpomerol_ is built. It then can be used to linking with executables. The example of the latter is given in example section and some working executables are given in prog subdirectory.
Documentation can be compiled with a `make doc` command.
#####  Features
  * High performance exact calculation of a Green's function and a two-particle Green's function in Matsubara domain.
  * Written in C++: iterators are used to avoid zero matrix elements and vanishing combinations. 
  * Symmetry analysis. The commutation relations between operators are taken into account.
  * Fermionic operators algebra to diagonalize any fermionic Hamiltonian.
  * [Eigen3](http://eigen.tuxfamily.org) template library for linear algebra is used (mostly its Sparse module).
  * [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface) support. 
  * [CMake](http://www.cmake.org) is used for the installation.

##### Installation
  Check the *dependencies*: c++ compiler, CMake, Eigen3, Boost (with Boost::mpi and serialization), mpi and git to fetch the sources. For compiling the binaries from prog you'll need a tclap header-only library. 
  - Checkout the latest sources `git clone https://github.com/aeantipov/pomerol.git`
  - In a directory run ` cmake <path_to_pomerol> -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<path>` 
    * add `-DTesting=ON` for compiling tests. Default = ON.
    * add `-DProgs=ON` for compiling provided binaries (from progs directory). These include diagonalizations of the Anderson impurity and Hubbard 2d cluster. Default = OFF.
    * add `-DExamples=ON` for compiling provided binary examples (from examples directory). Default = OFF.
    * add `-Duse_complex=ON` for allowing complex matrix elements in the Hamiltonian. Default = OFF.
  - ` make`
  - ` make test` (if tests are compiled)
  - ` make install`
    * Shared library _libpomerol_ will be in `<path>/lib`.
  - ` make doc` generates the documentation in the `doc` subfolder.

##### License 
The software is released under GPLv2 license. 
Academic usage : please attribute this work by a citation to http://dx.doi.org/10.5281/zenodo.17900.
##### Authors
  * Andrey Antipov <Andrey.E.Antipov\at\gmail.com>
  * Igor Krivenko <igor\at\shg.ru>

##### Development/Help 
Please feel free to contact and contribute!
