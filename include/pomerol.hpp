//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol.hpp
/// \brief Main "include-all" header of the library.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_POMEROL_HPP
#define POMEROL_INCLUDE_POMEROL_HPP

#include "mpi_dispatcher/mpi_dispatcher.hpp"
#include "mpi_dispatcher/mpi_skel.hpp"

#include "pomerol/DensityMatrix.hpp"
#include "pomerol/EnsembleAverage.hpp"
#include "pomerol/FieldOperatorContainer.hpp"
#include "pomerol/GFContainer.hpp"
#include "pomerol/Hamiltonian.hpp"
#include "pomerol/Index.hpp"
#include "pomerol/IndexClassification.hpp"
#include "pomerol/LatticePresets.hpp"
#include "pomerol/Misc.hpp"
#include "pomerol/MonomialOperator.hpp"
#include "pomerol/Operators.hpp"
#include "pomerol/StatesClassification.hpp"
#include "pomerol/Susceptibility.hpp"
#include "pomerol/ThreePointSusceptibility.hpp"
#include "pomerol/ThreePointSusceptibilityContainer.hpp"
#include "pomerol/TwoParticleGF.hpp"
#include "pomerol/TwoParticleGFContainer.hpp"

namespace Pomerol {

/// \mainpage
/// The source code and installation instructions are located in
/// <a href="https://github.com/pomerol-ed/pomerol">project's GitHub repository</a>.
///
/// \section ref_API libpomerol API
/// Most of pomerol's API is defined in \ref Pomerol namespace with an exception of the MPI tools
/// that are defined in \ref pMPI.
///
/// The entire library is available via a single include file \ref pomerol.hpp.
///
/// A pomerol-based exact diagonalization program would normally proceed as follows.
/// \li Define an expression corresponding to the Hamiltonian of the system either in terms of
///     <a href="https://krivenko.github.io/libcommute/expression/factories.html">
///     libcommute's creation/annihilation operators</a> (these are injected into the
///     \ref Pomerol namespace in \ref pomerol/Operators.hpp) or by using \ref LatticePresets "lattice presets".
/// \li Construct an \ref IndexClassification object that carries information about all
///     fermionic single-particle states of the system (typically, by calling \ref MakeIndexClassification()).
/// \li Construct system's Hilbert space by calling \ref MakeHilbertSpace().
/// \li Employ a \ref StatesClassification object to reveal invariant subspaces (sectors) of the Hamiltonian.
/// \li Define, prepare and diagonalize the \ref Hamiltonian matrix of the system.
/// \li Define, prepare and compute the many-body \ref DensityMatrix of the system.
/// \li Use a \ref FieldOperatorContainer to allocate matrices of fermionic creation/annihilation operators
///     (\ref CreationOperator / \ref AnnihilationOperator).
///     and transform them into the eigenbasis of the Hamiltonian.
/// \li Finally, compute some of the following physically relevant quantities.
///     - Gibbs ensemble averages of \ref MonomialOperator "operators" of physical observables (\ref EnsembleAverage).
///     - Single-particle fermionic Green's functions (\ref GreensFunction, \ref GFContainer).
///     - Dynamical susceptibilities -- correlators of two \ref MonomialOperator's (\ref Susceptibility).
///     - 3-point susceptibilities -- correlators of two fermionic and one bosonic quadratic operator
///       (\ref ThreePointSusceptibility, \ref ThreePointSusceptibilityContainer).
///     - Two-particle fermionic Green's functions and irreducible vertices (\ref TwoParticleGF,
///       \ref TwoParticleGFContainer, \ref Vertex4).
///
/// All available classes and functions are grouped into a few modules, so you can check out
/// the <a href="modules.html">Modules</a> page to get an overview of libpomerol's API.

}

#endif // #ifndef POMEROL_INCLUDE_POMEROL_HPP
