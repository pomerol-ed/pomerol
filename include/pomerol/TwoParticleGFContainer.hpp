//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/TwoParticleGFContainer.hpp
/// \brief Storage for multiple fermionic two-particle Matsubara Green's functions.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_TWOPARTICLEGFCONTAINER_HPP
#define POMEROL_INCLUDE_TWOPARTICLEGFCONTAINER_HPP

#include "DensityMatrix.hpp"
#include "FieldOperatorContainer.hpp"
#include "Hamiltonian.hpp"
#include "Index.hpp"
#include "IndexClassification.hpp"
#include "IndexContainer4.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"
#include "TwoParticleGF.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace Pomerol {

/// \addtogroup 2PGF
///@{

/// \brief Container for instances of \ref TwoParticleGF.
///
/// This class stores multiple \f$(i,j,k,l)\f$-elements of a fermionic two-particle Matsubara Green's function
/// \f[ \chi_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_1}+\omega_{n_2}-\omega_{n_3}) =
///   \int_0^\beta
///   Tr[\mathcal{T}_\tau \hat\rho c_i(\tau_1)c_j(\tau_2)c^\dagger_k(\tau_3)c^\dagger_l(0)]
///   e^{i\omega_{n_1}\tau_1+i\omega_{n_2}\tau_2-i\omega_{n_3}\tau_3}
///   d\tau_1 d\tau_2 d\tau_3.
/// \f]
class TwoParticleGFContainer : public IndexContainer4<TwoParticleGF, TwoParticleGFContainer>, public Thermal {
public:
    /// A difference in energies with magnitude below this value is treated as zero.
    RealType ReduceResonanceTolerance = 1e-8;
    /// Minimal magnitude of the coefficient of a term for it to be taken into account.
    RealType CoefficientTolerance = 1e-16;
    /// Minimal magnitude of the coefficient of a term for it to be taken into account with respect to
    /// the amount of terms.
    RealType MultiTermCoefficientTolerance = 1e-5;

    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by the creation and annihilation operators.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] DM Many-body density matrix \f$\hat\rho\f$.
    /// \param[in] Ops A set of creation/annihilation operators \f$c^\dagger\f$/\f$c\f$.
    template <typename... IndexTypes>
    TwoParticleGFContainer(IndexClassification<IndexTypes...> const& IndexInfo,
                           StatesClassification const& S,
                           Hamiltonian const& H,
                           DensityMatrix const& DM,
                           FieldOperatorContainer const& Ops)
        : IndexContainer4<TwoParticleGF, TwoParticleGFContainer>(*this, IndexInfo),
          Thermal(DM),
          S(S),
          H(H),
          DM(DM),
          Operators(Ops) {}

    /// Prepare a set of elements \f$\chi_{ijkl}\f$.
    /// \param[in] Indices Set of index combinations of the elements \f$\chi_{ijkl}\f$ to be prepared.
    ///            An empty set results in creation of elements for all possible index combinations \f$(i,j,k,l)\f$.
    void prepareAll(std::set<IndexCombination4> const& Indices = {});
    /// Compute all prepared elements \f$\chi_{ijkl}\f$.
    /// \param[in] clearTerms If true, computed \ref TwoParticleGFPart's of all elements will be destroyed
    ///                       immediately after filling the precomputed value cache.
    /// \param[in] freqs List of frequency triplets \f$(\omega_{n_1},\omega_{n_2},\omega_{n_3})\f$
    ///                  for value pre-computation.
    /// \param[in] comm MPI communicator used to parallelize the computation.
    /// \param[in] split Enable MPI parallelization.
    /// \pre \ref prepareAll() has been called.
    std::map<IndexCombination4, std::vector<ComplexType>> computeAll(bool clearTerms = false,
                                                                     FreqVec const& freqs = {},
                                                                     MPI_Comm const& comm = MPI_COMM_WORLD,
                                                                     bool split = true);

protected:
    friend class IndexContainer4<TwoParticleGF, TwoParticleGFContainer>;

    /// Create a single element \f$\chi_{ijkl}\f$.
    /// \param[in] Indices Index combination \f$(i,j,k,l)\f$.
    std::shared_ptr<TwoParticleGF> createElement(IndexCombination4 const& Indices) const;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;

    /// The Hamiltonian.
    Hamiltonian const& H;
    /// Many-body density matrix \f$\hat\rho\f$.
    DensityMatrix const& DM;
    /// A set of creation/annihilation operators \f$c^\dagger\f$/\f$c\f$.
    FieldOperatorContainer const& Operators;

private:
    // Implementation details.
    std::map<IndexCombination4, std::vector<ComplexType>>
    computeAll_nosplit(bool clearTerms, FreqVec const& freqs = {}, MPI_Comm const& comm = MPI_COMM_WORLD);
    std::map<IndexCombination4, std::vector<ComplexType>>
    computeAll_split(bool clearTerms, FreqVec const& freqs = {}, MPI_Comm const& comm = MPI_COMM_WORLD);
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_TWOPARTICLEGFCONTAINER_HPP
