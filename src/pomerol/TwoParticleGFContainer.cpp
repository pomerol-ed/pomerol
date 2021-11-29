//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/TwoParticleGFContainer.cpp
/// \brief Storage for multiple fermionic two-particle Matsubara Green's functions (implementation).
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#include "pomerol/TwoParticleGFContainer.hpp"

#include <cmath>
#include <cstddef>

namespace Pomerol {

void TwoParticleGFContainer::prepareAll(std::set<IndexCombination4> const& InitialIndices) {
    fill(InitialIndices);
    for(auto& el : ElementsMap) {
        auto& g = static_cast<TwoParticleGF&>(el.second);
        g.ReduceResonanceTolerance = ReduceResonanceTolerance;
        g.CoefficientTolerance = CoefficientTolerance;
        g.MultiTermCoefficientTolerance = MultiTermCoefficientTolerance;
        g.prepare();
    }
}

std::map<IndexCombination4, std::vector<ComplexType>>
TwoParticleGFContainer::computeAll(bool clearTerms, FreqVec const& freqs, MPI_Comm const& comm, bool split) {
    if(split)
        return computeAll_split(clearTerms, freqs, comm);
    else
        return computeAll_nosplit(clearTerms, freqs, comm);
}

std::map<IndexCombination4, std::vector<ComplexType>>
TwoParticleGFContainer::computeAll_nosplit(bool clearTerms, FreqVec const& freqs, MPI_Comm const& comm) {
    std::map<IndexCombination4, std::vector<ComplexType>> out;
    for(auto& el : ElementsMap) {
        INFO("Computing 2PGF for " << el.first);
        out.emplace(el.first, static_cast<TwoParticleGF&>(el.second).compute(clearTerms, freqs, comm));
    }
    return out;
}

std::map<IndexCombination4, std::vector<ComplexType>>
TwoParticleGFContainer::computeAll_split(bool clearTerms, FreqVec const& freqs, MPI_Comm const& comm) {
    std::map<IndexCombination4, std::vector<ComplexType>> out;
    std::map<IndexCombination4, std::vector<ComplexType>> storage;

    int comm_size = pMPI::size(comm);
    int comm_rank = pMPI::rank(comm);

    // split communicator
    int ncomponents = static_cast<int>(NonTrivialElements.size());
    int ncolors = std::min(comm_size, int(NonTrivialElements.size()));
    RealType color_size = 1.0 * comm_size / static_cast<RealType>(ncolors);
    std::map<int, int> proc_colors;
    std::map<int, int> elem_colors;
    std::map<int, int> color_roots;
    for(int p = 0; p < comm_size; ++p) {
        int color = int(1.0 * p / color_size);
        proc_colors[p] = color;
        color_roots[color] = p;
    }
    for(int i = 0; i < ncomponents; ++i) {
        int color = i * ncolors / ncomponents;
        elem_colors[i] = color;
    }

    if(!comm_rank) {
        INFO("Splitting " << ncomponents << " components in " << ncolors << " communicators");
        for(std::size_t i = 0; i < ncomponents; ++i)
            INFO("2pgf " << i << " color: " << elem_colors[i] << " color_root: " << color_roots[elem_colors[i]]);
    }
    MPI_Barrier(comm);
    int comp = 0;

    MPI_Comm comm_split = nullptr;
    MPI_Comm_split(comm, proc_colors[comm_rank], comm_rank, &comm_split);

    for(auto iter = NonTrivialElements.begin(); iter != NonTrivialElements.end(); iter++, comp++) {
        if(elem_colors[comp] == proc_colors[comm_rank]) {
            INFO("C" << elem_colors[comp] << "p" << comm_rank << ": computing 2PGF for " << iter->first);
            storage[iter->first] = static_cast<TwoParticleGF&>(*(iter->second)).compute(clearTerms, freqs, comm_split);
        }
    }
    MPI_Barrier(comm);
    // distribute data
    if(!comm_rank)
        INFO_NONEWLINE("Distributing 2PGF container...");
    comp = 0;
    for(auto iter = NonTrivialElements.begin(); iter != NonTrivialElements.end(); iter++, comp++) {
        int sender = color_roots[elem_colors[comp]];
        TwoParticleGF& chi = *((iter)->second);
        for(std::size_t p = 0; p < chi.parts.size(); p++) {
            chi.parts[p].NonResonantTerms.broadcast(comm, sender);
            chi.parts[p].ResonantTerms.broadcast(comm, sender);
            std::vector<ComplexType> freq_data;
            int freq_data_size = {};
            if(comm_rank == sender) {
                freq_data = storage[iter->first];
                freq_data_size = static_cast<int>(freq_data.size());
                MPI_Bcast(&freq_data_size, 1, MPI_LONG, sender, comm);
                MPI_Bcast(freq_data.data(), freq_data_size, MPI_CXX_DOUBLE_COMPLEX, sender, comm);
            } else {
                MPI_Bcast(&freq_data_size, 1, MPI_LONG, sender, comm);
                freq_data.resize(freq_data_size);
                MPI_Bcast(freq_data.data(), freq_data_size, MPI_CXX_DOUBLE_COMPLEX, sender, comm);
            }
            out[iter->first] = freq_data;

            if(comm_rank != sender) {
                chi.setStatus(TwoParticleGF::Computed);
            }
        }
    }
    MPI_Barrier(comm);
    if(!comm_rank)
        INFO("done.");
    return out;
}

std::shared_ptr<TwoParticleGF> TwoParticleGFContainer::createElement(IndexCombination4 const& Indices) const {
    AnnihilationOperator const& C1 = Operators.getAnnihilationOperator(Indices.Index1);
    AnnihilationOperator const& C2 = Operators.getAnnihilationOperator(Indices.Index2);
    CreationOperator const& CX3 = Operators.getCreationOperator(Indices.Index3);
    CreationOperator const& CX4 = Operators.getCreationOperator(Indices.Index4);

    return std::make_shared<TwoParticleGF>(S, H, C1, C2, CX3, CX4, DM);
}

} // namespace Pomerol
