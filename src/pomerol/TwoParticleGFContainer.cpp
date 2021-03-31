//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


#include "pomerol/TwoParticleGFContainer.h"

#include "mpi_dispatcher/misc.hpp"

namespace Pomerol{

template<bool Complex>
TwoParticleGFContainer<Complex>::TwoParticleGFContainer(const IndexClassification<Complex>& IndexInfo,
                                                        const StatesClassification<Complex> &S,
                                                        const Hamiltonian<Complex> &H,
                                                        const DensityMatrix<Complex> &DM,
                                                        const FieldOperatorContainer<Complex>& Operators) :
    ContainerBase(this, IndexInfo), Thermal(DM),
    S(S),H(H),DM(DM), Operators(Operators),
    ReduceResonanceTolerance (1e-8),//1e-16),
    CoefficientTolerance (1e-16),//1e-16),
    MultiTermCoefficientTolerance (1e-5)//1e-5),
{}

template<bool Complex>
void TwoParticleGFContainer<Complex>::prepareAll(const std::set<IndexCombination4>& InitialIndices)
{
    ContainerBase::fill(InitialIndices);
    for(auto iter = ContainerBase::ElementsMap.begin();
        iter != ContainerBase::ContainerBaseElementsMap.end(); iter++) {
        static_cast<ElementT&>(iter->second).ReduceResonanceTolerance = ReduceResonanceTolerance;
        static_cast<ElementT&>(iter->second).CoefficientTolerance = CoefficientTolerance;
        static_cast<ElementT&>(iter->second).MultiTermCoefficientTolerance = MultiTermCoefficientTolerance;
        static_cast<ElementT&>(iter->second).prepare();
    };
}

template<bool Complex>
std::map<IndexCombination4, std::vector<ComplexType> > TwoParticleGFContainer<Complex>::computeAll(bool clearTerms, std::vector<std::tuple<ComplexType, ComplexType, ComplexType> > const& freqs, const MPI_Comm& comm, bool split)
{
    if (split)
        return computeAll_split(clearTerms, freqs, comm);
    else
        return computeAll_nosplit(clearTerms, freqs, comm);
}

template<bool Complex>
std::map<IndexCombination4,std::vector<ComplexType> > TwoParticleGFContainer<Complex>::computeAll_nosplit(bool clearTerms, std::vector<std::tuple<ComplexType, ComplexType, ComplexType> > const& freqs, const MPI_Comm& comm)
{
    std::map<IndexCombination4,std::vector<ComplexType> > out;
    for(auto iter = ContainerBase::ElementsMap.begin();
        iter != ContainerBase::ElementsMap.end(); iter++) {
        INFO("Computing 2PGF for " << iter->first);
        out.insert(std::make_pair(iter->first, static_cast<ElementT&>(iter->second).compute(clearTerms, freqs, comm)));
    };
    return out;
}

template<bool Complex>
std::map<IndexCombination4,std::vector<ComplexType> > TwoParticleGFContainer<Complex>::computeAll_split(bool clearTerms, std::vector<std::tuple<ComplexType, ComplexType, ComplexType> > const& freqs, const MPI_Comm& comm)
{
    std::map<IndexCombination4,std::vector<ComplexType> > out;
    std::map<IndexCombination4,std::vector<ComplexType> > storage;

    int comm_size = pMPI::size(comm);
    int rank = pMPI::rank(comm);

    // split communicator
    size_t ncomponents = ContainerBase::NonTrivialElements.size();
    size_t ncolors = std::min(comm_size, int(ContainerBase::NonTrivialElements.size()));
    RealType color_size = 1.0*comm_size/ncolors;
    std::map<int,int> proc_colors;
    std::map<int,int> elem_colors;
    std::map<int,int> color_roots;
    bool calc = false;
    for (size_t p=0; p<comm_size; p++) {
        int color = int (1.0*p / color_size);
        proc_colors[p] = color;
        color_roots[color]=p;
    }
    for (size_t i=0; i<ncomponents; i++) {
        int color = i*ncolors/ncomponents;
        elem_colors[i] = color;
    };

    if (!rank) {
        INFO("Splitting " << ncomponents << " components in " << ncolors << " communicators");
        for (size_t i=0; i<ncomponents; i++)
            INFO("2pgf " << i << " color: " << elem_colors[i] << " color_root: "
                         << color_roots[elem_colors[i]]);
    }
    MPI_Barrier(comm);
    int comp = 0;

    MPI_Comm comm_split;
    MPI_Comm_split(comm, proc_colors[rank], rank, &comm_split);

    for(auto iter = ContainerBase::NonTrivialElements.begin(); iter != ContainerBase::NonTrivialElements.end(); iter++, comp++) {
        bool calc = (elem_colors[comp] == proc_colors[rank]);
        if (calc) {
            INFO("C" << elem_colors[comp] << "p" << rank << ": computing 2PGF for " << iter->first);
            if (calc)
                storage[iter->first] = static_cast<ElementT&>(*(iter->second)).compute(clearTerms, freqs, comm_split);
        }
    }
    MPI_Barrier(comm);
    // distribute data
    if (!rank)
        INFO_NONEWLINE("Distributing 2PGF container...");
    comp = 0;
    for(auto iter = ContainerBase::NonTrivialElements.begin(); iter != ContainerBase::NonTrivialElements.end(); iter++, comp++) {
        int sender = color_roots[elem_colors[comp]];
        ElementT& chi = *((iter)->second);
        for (size_t p = 0; p < chi.parts.size(); p++) {
            chi.parts[p]->NonResonantTerms.broadcast(comm, sender);
            chi.parts[p]->ResonantTerms.broadcast(comm, sender);
            std::vector<ComplexType> freq_data;
            long freq_data_size;
            if (rank == sender) {
                freq_data = storage[iter->first];
                freq_data_size = freq_data.size();
                MPI_Bcast(&freq_data_size, 1, MPI_LONG, sender, comm);
                MPI_Bcast(freq_data.data(), freq_data_size, MPI_CXX_DOUBLE_COMPLEX, sender, comm);
            } else {
                MPI_Bcast(&freq_data_size, 1, MPI_LONG, sender, comm);
                freq_data.resize(freq_data_size);
                MPI_Bcast(freq_data.data(), freq_data_size, MPI_CXX_DOUBLE_COMPLEX, sender, comm);
            }
            out[iter->first] = freq_data;

            if (rank != sender) {
                chi.setStatus(ElementT::Computed);
            }
        }
    }
    MPI_Barrier(comm);
    if (!rank) INFO("done.");
    return out;
}

template<bool Complex>
auto TwoParticleGFContainer<Complex>::createElement(const IndexCombination4& Indices) const -> ElementT*
{
    const AnnihilationOperator<Complex> &C1 = Operators.getAnnihilationOperator(Indices.Index1);
    const AnnihilationOperator<Complex> &C2 = Operators.getAnnihilationOperator(Indices.Index2);
    const CreationOperator<Complex> &CX3 = Operators.getCreationOperator   (Indices.Index3);
    const CreationOperator<Complex> &CX4 = Operators.getCreationOperator   (Indices.Index4);

    return new ElementT(S,H,C1,C2,CX3,CX4,DM);
}

} // end of namespace Pomerol
