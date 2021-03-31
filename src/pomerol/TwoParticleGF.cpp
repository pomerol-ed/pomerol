#include "pomerol/TwoParticleGF.h"

#include "mpi_dispatcher/mpi_skel.hpp"

#include <vector>
#include <tuple>

namespace Pomerol{

template<bool Complex>
TwoParticleGF<Complex>::TwoParticleGF(const StatesClassification<Complex>& S, const Hamiltonian<Complex>& H,
                const AnnihilationOperator<Complex>& C1, const AnnihilationOperator<Complex>& C2,
                const CreationOperator<Complex>& CX3, const CreationOperator<Complex>& CX4,
                const DensityMatrix<Complex>& DM) :
    Thermal(DM.beta), ComputableObject(),
    S(S), H(H), C1(C1), C2(C2), CX3(CX3), CX4(CX4), DM(DM),
    parts(0), Vanishing(true),
    ReduceResonanceTolerance (1e-8),
    CoefficientTolerance (1e-16),
    MultiTermCoefficientTolerance (1e-5)
{
}

template<bool Complex>
TwoParticleGF<Complex>::~TwoParticleGF()
{
      for(auto iter = parts.begin(); iter != parts.end(); iter++)
          delete *iter;
}

template<bool Complex>
BlockNumber TwoParticleGF<Complex>::getLeftIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber RightIndex) const
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getLeftIndex(RightIndex);
        case 1: return C2.getLeftIndex(RightIndex);
        case 2: return CX3.getLeftIndex(RightIndex);
        default: return ERROR_BLOCK_NUMBER;
    }
}

template<bool Complex>
BlockNumber TwoParticleGF<Complex>::getRightIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex) const
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getRightIndex(LeftIndex);
        case 1: return C2.getRightIndex(LeftIndex);
        case 2: return CX3.getRightIndex(LeftIndex);
        default: return ERROR_BLOCK_NUMBER;
    }
}

template<bool Complex>
auto TwoParticleGF<Complex>::OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex) const -> const FieldOperatorPart<Complex>&
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getPartFromLeftIndex(LeftIndex);
        case 1: return C2.getPartFromLeftIndex(LeftIndex);
        case 2: return CX3.getPartFromLeftIndex(LeftIndex);
        default: assert(0);
    }
    throw std::logic_error("TwoParticleGF : could not find operator part");
}

template<bool Complex>
void TwoParticleGF<Complex>::prepare()
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of CX4.
    typename FieldOperator<Complex>::BlocksBimap const& CX4NontrivialBlocks = CX4.getBlockMapping();
    for(auto outer_iter = CX4NontrivialBlocks.right.begin();
        outer_iter != CX4NontrivialBlocks.right.end(); outer_iter++){ // Iterate over the outermost index.
            for(size_t p=0; p<6; ++p){ // Choose a permutation
                  BlockNumber LeftIndices[4];
                  LeftIndices[0] = outer_iter->first;
                  LeftIndices[3] = outer_iter->second;
                  LeftIndices[2] = getLeftIndex(p,2,LeftIndices[3]);
                  LeftIndices[1] = getRightIndex(p,0,LeftIndices[0]);
                  // < LeftIndices[0] | O_1 | LeftIndices[1] >
                  // < LeftIndices[1] | O_2 | getRightIndex(p,1,LeftIndices[1]) >
                  // < LeftIndices[2]| O_3 | LeftIndices[3] >
                  // < LeftIndices[3] | CX4 | LeftIndices[0] >
                  // Select a relevant 'world stripe' (sequence of blocks).
                  if(getRightIndex(p,1,LeftIndices[1]) == LeftIndices[2] && LeftIndices[1].isCorrect() && LeftIndices[2].isCorrect()){
                      { // check if retained blocks are included. If not, do not push.
                          bool include_block_retained=false;
                          for(int k=0; k<4; k++)
                              if (DM.isRetained(LeftIndices[k]))  include_block_retained=true;
                          if(!include_block_retained)  continue;
                      }
                      // DEBUG
                      /*DEBUG("new part: "  << S.getBlockInfo(LeftIndices[0]) << " "
                                          << S.getBlockInfo(LeftIndices[1]) << " "
                                          << S.getBlockInfo(LeftIndices[2]) << " "
                                          << S.getBlockInfo(LeftIndices[3]) << " "
                      <<"BlockNumbers part: "  << LeftIndices[0] << " " << LeftIndices[1] << " " << LeftIndices[2] << " " << LeftIndices[3]);
                      */
                      parts.push_back(new PartT(
                            OperatorPartAtPosition(p,0,LeftIndices[0]),
                            OperatorPartAtPosition(p,1,LeftIndices[1]),
                            OperatorPartAtPosition(p,2,LeftIndices[2]),
                            (CreationOperatorPart<Complex>&)CX4.getPartFromLeftIndex(LeftIndices[3]),
                            H.getPart(LeftIndices[0]), H.getPart(LeftIndices[1]), H.getPart(LeftIndices[2]), H.getPart(LeftIndices[3]),
                            DM.getPart(LeftIndices[0]), DM.getPart(LeftIndices[1]), DM.getPart(LeftIndices[2]), DM.getPart(LeftIndices[3]),
                      permutations3[p]));

                      (*parts.rbegin())->ReduceResonanceTolerance = ReduceResonanceTolerance;
                      (*parts.rbegin())->CoefficientTolerance = CoefficientTolerance;
                      (*parts.rbegin())->MultiTermCoefficientTolerance = MultiTermCoefficientTolerance;
                      }
            }
    }
    if ( parts.size() > 0 ) {
        Vanishing = false;
        INFO("TwoParticleGF(" << getIndex(0) << getIndex(1) << getIndex(2) << getIndex(3) << "): " << parts.size() << " parts will be calculated");
        }
    Status = Prepared;
}

template<bool Complex>
bool TwoParticleGF<Complex>::isVanishing(void) const
{
    return Vanishing;
}


// An mpi adapter to 1) compute 2pgf terms; 2) convert them to a Matsubara Container; 3) purge terms
typedef std::tuple<ComplexType, ComplexType, ComplexType> freq_tuple;
typedef std::vector<freq_tuple> freq_vec_t;
template<bool Complex = false>
struct ComputeAndClearWrap
{
    void run(){
        p->compute();
        if (fill_) {
            int wsize = freqs_->size();
            #ifdef POMEROL_USE_OPENMP
            #pragma omp parallel for
            #endif
            for (int w = 0; w < wsize; ++w) {
                (*data_)[w] += (*p)(std::get<0>((*freqs_)[w]), std::get<1>((*freqs_)[w]), std::get<2>((*freqs_)[w]));
                }
            #ifdef POMEROL_USE_OPENMP
            #pragma omp barrier
            #endif
            }
        if (clear_) p->clear();
    };
    ComputeAndClearWrap(freq_vec_t const* freqs, std::vector<ComplexType> *data,  TwoParticleGFPart<Complex> *p, bool clear, bool fill, int complexity = 1):
        freqs_(freqs), data_(data), p(p),clear_(clear), fill_(fill), complexity(complexity){};
    int complexity;
protected:
    freq_vec_t const* freqs_;
    std::vector<ComplexType>* data_;
    TwoParticleGFPart<Complex> *p;
    bool clear_;
    bool fill_;
};

template<bool Complex>
std::vector<ComplexType> TwoParticleGF<Complex>::compute(bool clear, std::vector<std::tuple<ComplexType, ComplexType, ComplexType> > const& freqs, const MPI_Comm& comm)
{
    std::vector<ComplexType> m_data;
    if (Status < Prepared) throw (exStatusMismatch());
    if (Status >= Computed) return m_data;
    if (!Vanishing) {
        // Create a "skeleton" class with pointers to part that can call a compute method
        pMPI::mpi_skel<ComputeAndClearWrap<Complex>> skel;
        bool fill_container = freqs.size() > 0;
        skel.parts.reserve(parts.size());
        m_data.resize(freqs.size(), 0.0);
        for (size_t i=0; i<parts.size(); i++) {
            skel.parts.push_back(ComputeAndClearWrap<Complex>(&freqs, &m_data, parts[i], clear, fill_container, 1));
        }
        std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true); // actual running - very costly

        // Start distributing data
        MPI_Barrier(comm);

        MPI_Allreduce(MPI_IN_PLACE,
                      m_data.data(),
                      m_data.size(),
                      MPI_CXX_DOUBLE_COMPLEX,
                      MPI_SUM,
                      comm);

        // Optionally distribute terms to other processes
        if (!clear) {
            for (size_t p = 0; p < parts.size(); p++) {
                parts[p]->NonResonantTerms.broadcast(comm, job_map[p]);
                parts[p]->ResonantTerms.broadcast(comm, job_map[p]);
                parts[p]->Status = TwoParticleGFPart<Complex>::Computed;
            }
            MPI_Barrier(comm);
        }
    }
    Status = Computed;
    return m_data;
}

template<bool Complex>
ParticleIndex TwoParticleGF<Complex>::getIndex(size_t Position) const
{
    switch(Position){
        case 0: return C1.getIndex();
        case 1: return C2.getIndex();
        case 2: return CX3.getIndex();
        case 3: return CX4.getIndex();
        default: assert(0);
    }
    throw std::logic_error("TwoParticleGF : could not get operator index");
}

template<bool Complex>
unsigned short TwoParticleGF<Complex>::getPermutationNumber ( const Permutation3& in )
{
    for (unsigned short i=0; i<6; ++i) if (in == permutations3[i]) return i;
    ERROR("TwoParticleGF: Permutation " << in << " not found in all permutations3");
    return 0;
}

template class TwoParticleGF<false>;
template class TwoParticleGF<true>;

} // end of namespace Pomerol
