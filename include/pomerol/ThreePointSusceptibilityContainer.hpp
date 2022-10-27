//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/ThreePointSusceptibilityContainer.hpp
/// \brief Storage for multiple 3-point susceptibilities in the Matsubara representation.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITYCONTAINER_HPP
#define POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITYCONTAINER_HPP

#include "DensityMatrix.hpp"
#include "FieldOperatorContainer.hpp"
#include "Hamiltonian.hpp"
#include "Index.hpp"
#include "IndexClassification.hpp"
#include "IndexContainer2.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"
#include "ThreePointSusceptibility.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace Pomerol {

/// \addtogroup 3PSusc
///@{

/// \brief Container for instances of \ref ThreePointSusceptibility.
///
/// This class stores multiple \f$(i, j, \hat B)\f$-elements of a 3-point susceptibility.
/// A quadratic operator \f$\hat B\f$ determines the channel the susceptibility is defined in,
/// and is shared by all stored elements.
class ThreePointSusceptibilityContainer
    : public IndexContainer2<ThreePointSusceptibility, ThreePointSusceptibilityContainer>,
      public Thermal {
public:
    /// A difference in energies with magnitude below this value is treated as zero.
    RealType ReduceResonanceTolerance = 1e-8;
    /// Minimal magnitude of the coefficient of a term for it to be taken into account.
    RealType CoefficientTolerance = 1e-16;

    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by the creation and annihilation operators.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] DM Many-body density matrix \f$\hat\rho\f$.
    /// \param[in] Ops A set of creation/annihilation operators \f$c^\dagger\f$/\f$c\f$.
    /// \param[in] B A quadratic operator \f$\hat B\f$. See \ref ThreePointSusceptibility for more details.
    template <typename... IndexTypes>
    ThreePointSusceptibilityContainer(IndexClassification<IndexTypes...> const& IndexInfo,
                                      StatesClassification const& S,
                                      Hamiltonian const& H,
                                      DensityMatrix const& DM,
                                      FieldOperatorContainer const& Ops,
                                      QuadraticOperator const& B)
        : IndexContainer2<ThreePointSusceptibility, ThreePointSusceptibilityContainer>(*this, IndexInfo),
          Thermal(DM),
          S(S),
          H(H),
          DM(DM),
          Operators(Ops),
          B(B),
          Channel(ThreePointSusceptibility::selectChannel(B)) {}

    /// Prepare a set of elements \f$\chi^{(3)}_{ij,\hat B}\f$.
    /// \param[in] Indices Set of index combinations of the elements \f$\chi^{(3)}_{ij,\hat B}\f$ to be prepared.
    ///            An empty set results in creation of elements for all possible index combinations \f$(i,j)\f$.
    void prepareAll(std::set<IndexCombination2> const& Indices = {});
    /// Compute all prepared elements \f$\chi^{(3)}_{ij,\hat B}\f$.
    /// \param[in] clear If true, computed \ref ThreePointSusceptibilityPart's of all elements will be destroyed
    ///                  immediately after filling the precomputed value cache.
    /// \param[in] freqs List of frequency duplets \f$(\omega_{n_1},\omega_{n_2})\f$ for value pre-computation.
    /// \pre \ref prepareAll() has been called.
    std::map<IndexCombination2, std::vector<ComplexType>>
    computeAll(bool clearTerms = false, FreqVec2 const& freqs = {}, MPI_Comm const& comm = MPI_COMM_WORLD);

protected:
    friend class IndexContainer2<ThreePointSusceptibility, ThreePointSusceptibilityContainer>;

    /// Create a single element \f$\chi^{(3)}_{ij,\hat B}\f$.
    /// \param[in] Indices Index combination \f$(i,j)\f$.
    std::shared_ptr<ThreePointSusceptibility> createElement(IndexCombination2 const& Indices) const;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;

    /// The Hamiltonian.
    Hamiltonian const& H;
    /// Many-body density matrix \f$\hat\rho\f$.
    DensityMatrix const& DM;
    /// A set of creation/annihilation operators \f$c^\dagger\f$/\f$c\f$.
    FieldOperatorContainer const& Operators;
    /// Quadratic operator B.
    QuadraticOperator const& B;
    /// Channel of \f$\chi^{(3)}\f$.
    ThreePointSusceptibility::Channel Channel;
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITYCONTAINER_HPP
