# Changelog

All notable changes to this project will be documented in this file.

## [2.1] - Unreleased

- New classes to compute the 3-point fermion-boson susceptibilities,
  `ThreePointSusceptibility`, `ThreePointSusceptibilityPart` and
  `ThreePointSusceptibilityContainer`.
- `QuadraticOperator` can now be a product of two creators or two annihilators.
- Renamed type aliases `FreqTuple` -> `FreqTuple3` and `FreqVec` -> `FreqVec3`.
  The old names are still usable but marked as deprecated.
- Bump required libcommute version to 0.7.2.
- Install CMake configuration files into
  ``${CMAKE_INSTALL_PREFIX}/lib/cmake/pomerol``, which is the recommended
  location.

## [2.0] - 2021-11-30

- The Mozilla Public License Version 2.0 has been adopted.

- Pomerol 2.0 requires a C++11 compatible compiler to build.

- Dependence on Boost.MPI and Boost.Serialization has been dropped. Pomerol 2.0
  still depends on a few header-only Boost libraries and an MPI-3.0
  implementation that provides a working `<mpi.h>`.

- Changed extension of all header files from `.h` to `.hpp`.

- The input layer and diagonalization routines have been rewritten to benefit
  from facilities provided by
  [libcommute](https://github.com/krivenko/libcommute). Hamiltonians to be
  solved are now specified as libcommute's expressions with arbitrary static
  types of indices carried by creation/annihilation operators.

- Thanks to the use of libcommute's expressions, it is now possible to solve
  models, whose Hamiltonians involve bosonic degrees of freedom.

- The `Lattice` class has been retired. The `LatticePresets` class has been
  turned into a namespace with preset functions returning their respective
  expressions of Hamiltonian terms. For the sake of backward compatibility,
  operators in the expressions returned by these functions carry the
  traditional (site label, orbital index, spin projection) index triplets.

- The `spin` enumeration type is now declared in the `LatticePresets` namespace.
  The enumeration has also been extended with an extra `undef` value, which is
  meant to be used on bosonic creation/annihilation operators.

- The `POMEROL_COMPLEX_MATRIX_ELEMENTS` CMake option has been removed.
  A proper matrix storage format is selected at runtime depending on the types
  of input expressions for Hamiltonians and operators of physical observables.

- API of class `IndexClassification` has been made more generic (it is now
  templated on the operator index types to accommodate the flexibility of
  libcommute's expressions).

- Class `IndexHamiltonian` has been removed.

- Functionality of the `Symmetrizer` and `StatesClassification` classes has been
  redistributed between `StatesClassification` and a new class `HilbertSpace`.
  The notion of quantum numbers has been abandoned since partition of a
  Hilbert space into sectors is now performed by libcommute's `space_partition`
  algorithm.

- Class `FieldOperator` has been generalized and renamed into
  `MonomialOperator`. It can now compute and store matrices of arbitrary
  monomial operators, i.e. operators that are products of creation/annihilation
  operators, possibly with a real or complex prefactor.

- It is now possible to use `EnsembleAverage` and `Susceptibility` to compute
  averages/dynamical fluctuations of monomial operators instead of just
  quadratic operators.

- Renamed `EnsembleAverage::getResult()` to `EnsembleAverage::operator()()`.

- The outdated `ENABLE_SAVE_PLAINTEXT` macro has been removed.

- The `pomerol/first_include.h.in` header has been renamed into
  `pomerol/Version.hpp.in`.

- API reference documentation and the tutorial have been updated and cleaned up.

- Unit tests have been ported from Google Test to Catch2. The header file of
  Catch2 is bundled to the source code.
