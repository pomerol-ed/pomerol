#include "pomerol/OperatorPresets.h"

namespace Pomerol {
namespace OperatorPresets {

//
// Operator N
//

template<bool Complex>
N<Complex>::N(ParticleIndex Nmodes):Operator<Complex>(),Nmodes(Nmodes)
{
    for (ParticleIndex index=0; index<Nmodes; ++index) {
        (*this)+=n<Complex>(index);
    };
};

template<bool Complex>
std::map<FockState,MelemType<Complex>> N<Complex>::actRight(const FockState &ket) const
{
    std::map<FockState,MelemType<Complex>> output;
    output[ket]=this->getMatrixElement(ket);
    return output;
}

template<bool Complex>
MelemType<Complex> N<Complex>::getMatrixElement(const FockState &bra, const FockState &ket) const
{
    return (bra!=ket)?0:getMatrixElement(ket);
}

template<bool Complex>
MelemType<Complex> N<Complex>::getMatrixElement(const FockState &ket) const
{
    return ket.count();
}

//
// Operator Sz
//

template<bool Complex>
Sz<Complex>::Sz(ParticleIndex Nmodes, const std::vector<ParticleIndex> & SpinUpIndices)
    : Operator<Complex>(),Nmodes(Nmodes),SpinUpIndices(SpinUpIndices)
{
    SpinDownIndices.reserve(Nmodes/2);
    for (ParticleIndex i=0; i<Nmodes; i++) {
        if (std::find(SpinUpIndices.begin(), SpinUpIndices.end(), i)==SpinUpIndices.end()) SpinDownIndices.push_back(i);
        };
    if (SpinUpIndices.size() != SpinDownIndices.size() ) { throw (typename Pomerol::Operator<Complex>::exWrongLabel() ); ERROR("Sz operator requires even number of indices"); };
    generateTerms();
}

template<bool Complex>
Sz<Complex>::Sz(const std::vector<ParticleIndex> & SpinUpIndices, const std::vector<ParticleIndex> & SpinDownIndices)
    : Operator<Complex>(),Nmodes(SpinUpIndices.size() + SpinDownIndices.size()),SpinUpIndices(SpinUpIndices), SpinDownIndices(SpinDownIndices)
{
    if (SpinUpIndices.size() != SpinDownIndices.size() ) { throw (typename Pomerol::Operator<Complex>::exWrongLabel() ); ERROR("Sz operator requires even number of indices"); };
    generateTerms();
}

template<bool Complex>
void Sz<Complex>::generateTerms()
{
    for (ParticleIndex i=0; i<SpinUpIndices.size(); ++i) {
        (*this)+=n<Complex>(SpinUpIndices[i])*0.5;
        (*this)-=n<Complex>(SpinDownIndices[i])*0.5;
    }
}

template<bool Complex>
MelemType<Complex> Sz<Complex>::getMatrixElement(const FockState &ket) const
{
    int up_value=0, down_value=0;
    std::vector<ParticleIndex>::const_iterator it_up=SpinUpIndices.begin(), it_down=SpinDownIndices.begin();
    for (;it_up!=SpinUpIndices.end();it_up++) up_value+=ket.test(*it_up);
    for (;it_down!=SpinDownIndices.end();it_down++) down_value+=ket.test(*it_down);
    return 0.5*(up_value-down_value);
}

template<bool Complex>
MelemType<Complex> Sz<Complex>::getMatrixElement(const FockState &bra, const FockState &ket) const
{
    return (bra!=ket)?0:getMatrixElement(ket);
}

template<bool Complex>
std::map<FockState,MelemType<Complex>> Sz<Complex>::actRight(const FockState &ket) const
{
    std::map<FockState,MelemType<Complex>> output;
    output[ket]=this->getMatrixElement(ket);
    return output;
}

// Explicit instantiations: Real case

template class N<false>;
template class Sz<false>;
template class Cdag<false>;
template class C<false>;
template class N_offdiag<false>;

// Explicit instantiations: Complex case

template class N<true>;
template class Sz<true>;
template class Cdag<true>;
template class C<true>;
template class N_offdiag<true>;

} // end of namespace OperatorPresets
} // end of namespace Pomerol
