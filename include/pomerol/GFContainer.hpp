//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/GFContainer.hpp
/// \brief Storage for multiple fermionic single-particle Matsubara Green's functions.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_GFCONTAINER_HPP
#define POMEROL_INCLUDE_POMEROL_GFCONTAINER_HPP

#include "DensityMatrix.hpp"
#include "FieldOperatorContainer.hpp"
#include "GreensFunction.hpp"
#include "Hamiltonian.hpp"
#include "Index.hpp"
#include "IndexClassification.hpp"
#include "IndexContainer2.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

#include <memory>
#include <set>

namespace Pomerol {

/// \addtogroup GF
///@{

/// \brief Container for instances of \ref GreensFunction.
///
/// This class stores elements of a matrix-valued fermionic single-particle Matsubara Green's function
/// \f[
///  G_{ij}(i\omega_n) = -\int_0^\beta d\tau e^{i\omega_n\tau} Tr[\mathcal{T}_\tau \hat\rho c_i(\tau) c_j^\dagger(0)].
/// \f]
class GFContainer : public IndexContainer2<GreensFunction, GFContainer>, public Thermal {

public:
    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by the creation and annihilation operators.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] DM Many-body density matrix \f$\hat\rho\f$.
    /// \param[in] Ops A set of creation/annihilation operators \f$c^\dagger_j\f$/\f$c_i\f$.
    template <typename... IndexTypes>
    GFContainer(IndexClassification<IndexTypes...> const& IndexInfo,
                StatesClassification const& S,
                Hamiltonian const& H,
                DensityMatrix const& DM,
                FieldOperatorContainer const& Ops)
        : IndexContainer2<GreensFunction, GFContainer>(*this, IndexInfo),
          Thermal(DM),
          S(S),
          H(H),
          DM(DM),
          Operators(Ops) {}

    /// Prepare a set of matrix elements \f$G_{ij}\f$.
    /// \param[in] Indices Set of index combinations of the elements \f$G_{ij}\f$ to be prepared.
    ///            An empty set results in creation of elements for all possible index combinations.
    void prepareAll(std::set<IndexCombination2> const& Indices = {});
    /// Compute all prepared matrix elements \f$G_{ij}\f$.
    /// \pre \ref prepareAll() has been called.
    void computeAll();

protected:
    friend class IndexContainer2<GreensFunction, GFContainer>;

    /// Create a single matrix element \f$G_{ij}\f$.
    /// \param[in] Indices Index combination \f$(i,j)\f$.
    std::shared_ptr<GreensFunction> createElement(IndexCombination2 const& Indices) const;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;

    /// The Hamiltonian.
    Hamiltonian const& H;
    /// Many-body density matrix \f$\hat\rho\f$.
    DensityMatrix const& DM;
    /// A set of creation/annihilation operators \f$c^\dagger_j\f$/\f$c_i\f$.
    FieldOperatorContainer const& Operators;
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_GFCONTAINER_HPP
