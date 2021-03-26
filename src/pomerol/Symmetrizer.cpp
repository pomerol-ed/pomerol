#include "pomerol/Symmetrizer.h"
#include "pomerol/OperatorPresets.h"


namespace Pomerol {

//
//Symmetrizer::IndexPermutation
//

template<bool Complex>
Symmetrizer<Complex>::IndexPermutation::IndexPermutation(const DynamicIndexCombination &in):N(in.getNumberOfIndices())
{
    if ( checkConsistency(in) && checkIrreducibility(in) ) {
        Combinations.push_back(new DynamicIndexCombination(in));
        CycleLength=0;
        calculateCycleLength();
    }
    else throw ( DynamicIndexCombination::exWrongIndices()) ;
}

template<bool Complex>
bool Symmetrizer<Complex>::IndexPermutation::checkConsistency(const DynamicIndexCombination &in)
{

    for (ParticleIndex i=0; i<N; ++i) {
        if (in.getIndex(i)>=N) {
            ERROR("Indices in IndexPermutation should belong to the interval 0..N-1");
            return false;
            };
        for (ParticleIndex j=i+1; j<N; ++j)
            if (in.getIndex(i)==in.getIndex(j)) {
            ERROR("Found equal indices in given combination");
            return false;
            }
        };
    return true;
}

template<bool Complex>
bool Symmetrizer<Complex>::IndexPermutation::checkIrreducibility(const DynamicIndexCombination &in)
{
    std::map<ParticleIndex, unsigned int> nontrivial_indices; // Here collected indices, which are changed
    std::vector<ParticleIndex> trivial_indices;
    ParticleIndex current_index=0;

    while ((nontrivial_indices.size() + trivial_indices.size())!=N) {
        //DEBUG("Current index : " << current_index << "-->" << in.getIndex(current_index));
        if (current_index == in.getIndex(current_index)) { // Index is a trivial index - no loop for it is done
            trivial_indices.push_back(current_index);
            current_index++;
            }
        else if (nontrivial_indices.find(current_index)!=nontrivial_indices.end()) // Check that this index was found during previous iterations.
            current_index++;
            else if ( nontrivial_indices.size()==0 ) { // This is a first nontrivial index found - start a loop then.
                nontrivial_indices[current_index]=1;
                ParticleIndex result=in.getIndex(current_index);
                //DEBUG("Begin with" << current_index);
                //DEBUG("-->" << result);
                while ( result!=current_index) { // make a small loop in indices to determine the length of found cycle
                    nontrivial_indices[result]=1;
                    result=in.getIndex(result);
                    //DEBUG("-->" << result);
                    };
                current_index++;
            }
        else {
            ERROR("Permutation " << in << " is reducible");
            return false; // Once an index which is not a part of previously checked loop is found means a permutation is reducible.
            };
    };
    if ( trivial_indices.size() == N ) { ERROR("Identity permutation " << in << " is rejected."); return false; } // reject trivial identity permutation.
    return true;
}

template<bool Complex>
void Symmetrizer<Complex>::IndexPermutation::calculateCycleLength()
{
    DynamicIndexCombination initial(**(Combinations.begin()));
    DynamicIndexCombination current(**(Combinations.begin()));
    DynamicIndexCombination trivial (Symmetrizer::generateTrivialCombination(N)); // #warning think of better static implementation
    DynamicIndexCombination next(N);
    for (ParticleIndex i=0; i<N; ++i) trivial[i] = i;
    bool exit_loop=false;
    while (!exit_loop) {
        for (ParticleIndex i=0; i<N; ++i) next[i]=current[current[i]];
        CycleLength++;
        exit_loop = ( next == initial || next == trivial );
        current = next;
        }
}

template<bool Complex>
const DynamicIndexCombination& Symmetrizer<Complex>::IndexPermutation::getIndices( unsigned int cycle_number ) const
{
    return *Combinations[cycle_number];
}

template<bool Complex>
const unsigned int Symmetrizer<Complex>::IndexPermutation::getCycleLength() const
{
    return CycleLength;
}

template<bool Complex>
const char* Symmetrizer<Complex>::IndexPermutation::exEqualIndices::what() const throw(){
    return "Cannot have equal indices in the Symmetrizer index combination";
};

//
// Symmetrizer::QuantumNumbers
//

template<bool Complex>
Symmetrizer<Complex>::QuantumNumbers::QuantumNumbers(int amount):
  amount(amount),
  numbers( std::vector<MelemType<Complex>>(amount)),
  NumbersHash(numbers_hash_generator(numbers))
{
};

template<bool Complex>
bool Symmetrizer<Complex>::QuantumNumbers::set ( int pos, MelemType<Complex> val )
{
    if (pos<amount) {
        numbers[pos] = val;
        NumbersHash = numbers_hash_generator(numbers);
        }
    else {
        ERROR("Tried to insert element " << val << " to wrong position " << pos << " in " << __PRETTY_FUNCTION__ );
        return false;
        };
    return true;
}

template<bool Complex>
bool Symmetrizer<Complex>::QuantumNumbers::operator< (const Symmetrizer::QuantumNumbers& rhs) const
{
    return (NumbersHash<rhs.NumbersHash);
}

template<bool Complex>
bool Symmetrizer<Complex>::QuantumNumbers::operator== (const Symmetrizer::QuantumNumbers& rhs) const
{
    return (NumbersHash==rhs.NumbersHash);
}

template<bool Complex>
bool Symmetrizer<Complex>::QuantumNumbers::operator!= (const Symmetrizer::QuantumNumbers& rhs) const
{
    return (NumbersHash!=rhs.NumbersHash);
}

//
// Symmetrizer
//

template<bool Complex>
Symmetrizer<Complex>::Symmetrizer(const IndexClassification<Complex> &IndexInfo, const IndexHamiltonian<Complex> &Storage):
    ComputableObject(),
    IndexInfo(IndexInfo),
    Storage(Storage),
    NSymmetries(0)
{
}

template<bool Complex>
const DynamicIndexCombination& Symmetrizer<Complex>::generateTrivialCombination(ParticleIndex N)
{
    static DynamicIndexCombination trivial(N);
    for (ParticleIndex i=0; i<N; ++i) trivial[i] = i;
    return trivial;
}

template<bool Complex>
const std::vector<std::shared_ptr<Operator<Complex>> >& Symmetrizer<Complex>::getOperations() const
{
    return Operations;
}

template<bool Complex>
bool Symmetrizer<Complex>::checkSymmetry(const Operator<Complex> &in)
{
    std::shared_ptr<Operator<Complex>> OP1 ( new Operator<Complex>(in));
    // Check that OP1 is an integrals of motion
    if (!Storage.commutes(*OP1)) return false;

    // Check that all Fock states are eigenstates of OP1
    // Otherwise, it's unsuitable for Hilbert space partitioning
    for(ParticleIndex i = 0; i < IndexSize; ++i) {
        if (!OperatorPresets::n<Complex>(i).commutes(*OP1)) return false;
    }

    Operations.push_back(OP1);
    NSymmetries++;
    return true;
}

template<bool Complex>
void Symmetrizer<Complex>::compute(const std::vector<Operator<Complex>>& integrals_of_motion)
{
    if (Status>=Computed) return;
    IndexSize = IndexInfo.getIndexSize();

    for(int i = 0; i < integrals_of_motion.size(); ++i) {
        const Operator<Complex>& in = integrals_of_motion[i];
        if (checkSymmetry(in)) INFO("[ H ," << in << " ]=0");
    }

    Status = Computed;
}

template<bool Complex>
void Symmetrizer<Complex>::compute(bool ignore_symmetries)
{
    if (Status>=Computed) return;
    IndexSize = IndexInfo.getIndexSize();
    if (!ignore_symmetries) {
        // Check particle number conservation
        Operator<Complex> op_n = Pomerol::OperatorPresets::N<Complex>(IndexSize);
        if (this->checkSymmetry(op_n)) INFO("[ H ," << op_n << " ]=0");

        // Check Sz conservation
        bool valid_sz = true;
        for (ParticleIndex i=0; i<IndexSize && valid_sz; ++i)
            valid_sz = valid_sz && (IndexInfo.getInfo(i).Spin == up || IndexInfo.getInfo(i).Spin == down);
        if (valid_sz) {
            std::vector<ParticleIndex> SpinUpIndices;
            for (ParticleIndex i=0; i<IndexSize; ++i) {
                unsigned short Spin = IndexInfo.getInfo(i).Spin;
                if ( Spin == up ) SpinUpIndices.push_back(i);
            }
            Operator<Complex> op_sz = Pomerol::OperatorPresets::Sz<Complex>(IndexSize, SpinUpIndices);
            if (this->checkSymmetry(op_sz)) INFO("[ H ," << op_sz << " ]=0");
        };
    };

    Status = Computed;
}

template<bool Complex>
typename Symmetrizer<Complex>::QuantumNumbers Symmetrizer<Complex>::getQuantumNumbers() const
{
    return Symmetrizer::QuantumNumbers(NSymmetries);
}

template<bool Complex>
std::ostream& operator<<(std::ostream& output, const typename Symmetrizer<Complex>::QuantumNumbers& out)
{
    output << "[";
    for (int i=0 ;i<out.amount-1; ++i) output << out.numbers[i] << ",";
    if ( out.amount ) output << out.numbers[out.amount-1];
    output << "]";
    return output;
}

} // end of namespace Pomerol
