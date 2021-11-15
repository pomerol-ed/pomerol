//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

class TwoParticleGFContainer : public IndexContainer4<TwoParticleGF, TwoParticleGFContainer>, public Thermal {
public:
    /** A difference in energies with magnitude less than this value is treated as zero. default = 1e-8. */
    RealType ReduceResonanceTolerance = 1e-8;
    /** Minimal magnitude of the coefficient of a term to take it into account. default = 1e-16. */
    RealType CoefficientTolerance = 1e-16;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. default = 1e-5. */
    RealType MultiTermCoefficientTolerance = 1e-5;

    template <typename... IndexTypes>
    TwoParticleGFContainer(IndexClassification<IndexTypes...> const& IndexInfo,
                           StatesClassification const& S,
                           Hamiltonian const& H,
                           DensityMatrix const& DM,
                           FieldOperatorContainer const& Operators)
        : IndexContainer4<TwoParticleGF, TwoParticleGFContainer>(*this, IndexInfo),
          Thermal(DM),
          S(S),
          H(H),
          DM(DM),
          Operators(Operators) {}

    void prepareAll(std::set<IndexCombination4> const& InitialIndices = {});
    std::map<IndexCombination4, std::vector<ComplexType>> computeAll(bool clearTerms = false,
                                                                     FreqVec const& freqs = {},
                                                                     MPI_Comm const& comm = MPI_COMM_WORLD,
                                                                     bool split = true);
    std::map<IndexCombination4, std::vector<ComplexType>>
    computeAll_nosplit(bool clearTerms, FreqVec const& freqs = {}, MPI_Comm const& comm = MPI_COMM_WORLD);
    std::map<IndexCombination4, std::vector<ComplexType>>
    computeAll_split(bool clearTerms, FreqVec const& freqs = {}, MPI_Comm const& comm = MPI_COMM_WORLD);

protected:
    friend class IndexContainer4<TwoParticleGF, TwoParticleGFContainer>;
    std::shared_ptr<TwoParticleGF> createElement(IndexCombination4 const& Indices) const;

    StatesClassification const& S;

    Hamiltonian const& H;
    DensityMatrix const& DM;
    FieldOperatorContainer const& Operators;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_TWOPARTICLEGFCONTAINER_HPP
