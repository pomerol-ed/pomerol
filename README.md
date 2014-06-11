*pomerol* is an exact diagonalization (full-ED) code written in C++ aimed at solving condensed matter second-quantized models of interacting fermions on finite size lattices at finite temperatures. It is mostly developed to produce single and two-particle Greens functions and corresponding vertex functions of a Hubbard model on a cluster.

_[10.01.2014]_ Updated pomerol to 2.0 version, that supports MPI. Removed obsolete hdf5 bindings.  
_[10.06.2014]_ Moved repo to github. 

##### Documentation
Check http://pomerol.sourceforge.net for the reference documentation.

The library, _libpomerol_ is built. It then can be used to linking with executables. The example of the latter is given in example section and some working executables are given in prog subdirectory.
Documentation can be compiled with a `make doc` command.
#####  Features
  * Written in C++: iterators are used to avoid zero matrix elements and vanishing combinations. 
  * Symmetry analysis. The commutation relations between operators are taken into account.
  * [http://eigen.tuxfamily.org Eigen3] template library for linear algebra is used (mostly its Sparse module).
  * High performance exact calculation of a Green's function and a two-particle Green's function in Matsubara domain.
  * [http://en.wikipedia.org/wiki/Message_Passing_Interface MPI] support. 
  * [http://www.cmake.org CMake] used for the installation.

##### Installation
  # Check the *dependencies*: c++11 - compatible compiler (clang++-3.2 / intel c++ - 14.0 / g++ - 4.8), CMake, Eigen3, Boost (with Boost::mpi and serialization), mpi and git to fetch the sources. For compiling the binaries from prog you'll need a tclap header-only library. 
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
The software is released under MIT license. You may modify and redistribute it according to the license. We however consider that it contains scientific value and kindly ask for an acknowledgment in case of publishing results obtained with help of *pomerol* code.
##### Authors
  * Andrey Antipov <Andrey.E.Antipov\at\gmail.com>
  * Igor Krivenko <igor\at\shg.ru>

##### Development/Help 
Please feel free to contact and contribute!
