/** \file OperatorPresets.h
**  \brief Declarations of the OperatorPresets class. This is some workaround to generate easy classes
**
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_OPERATOR_PRESETS_H__
#define __INCLUDE_OPERATOR_PRESETS_H__

#include "Misc.h"
#include "Operator.h"

namespace Pomerol {
namespace OperatorPresets {

template<bool Complex = false>
class N : public Operator<Complex> {
private:
    const ParticleIndex Nmodes;
public:

    using MelemT = MelemType<Complex>;

    N(ParticleIndex Nmodes);
    std::map <FockState, MelemT> actRight(const FockState &ket) const;
    MelemT getMatrixElement(const FockState &bra, const FockState &ket) const;
    MelemT getMatrixElement(const FockState &ket) const;
};

template<bool Complex = false>
class Sz : public Operator<Complex> {
private:
    const int Nmodes;
    std::vector<ParticleIndex> SpinUpIndices;
    std::vector<ParticleIndex> SpinDownIndices;
    void generateTerms();
public:

    using MelemT = MelemType<Complex>;

    Sz(ParticleIndex Nmodes, const std::vector<ParticleIndex> & SpinUpIndices);
    Sz(const std::vector<ParticleIndex> & SpinUpIndices, const std::vector<ParticleIndex> & SpinDownIndices);
    std::map <FockState, MelemT> actRight(const FockState &ket) const;
    MelemT getMatrixElement(const FockState &bra, const FockState &ket) const;
    MelemT getMatrixElement(const FockState &ket) const;
};

template<bool Complex = false>
class Cdag : public Operator<Complex> {
private:
    ParticleIndex index;
public:
    Cdag(ParticleIndex index) : Operator<Complex>(c_dag<Complex>(index)){};
};

template<bool Complex = false>
class C : public Operator<Complex> {
private:
    ParticleIndex index;
public:
    C(ParticleIndex index):Operator<Complex>(c<Complex>(index)){};
};

template<bool Complex = false>
class N_offdiag : public Operator<Complex> {
private:
    ParticleIndex index1, index2;
public:
    N_offdiag(ParticleIndex index1, ParticleIndex index2):
      Operator<Complex>(n_offdiag<Complex>(index1, index2)){};
};

// External templates: Real case

extern template class N<false>;
extern template class Sz<false>;
extern template class Cdag<false>;
extern template class C<false>;
extern template class N_offdiag<false>;

// External templates: Complex case

extern template class N<true>;
extern template class Sz<true>;
extern template class Cdag<true>;
extern template class C<true>;
extern template class N_offdiag<true>;

} // end of namespace OperatorPresets
} // end of namespace Pomerol

#endif // endif :: #ifndef __INCLUDE_OPERATOR_PRESETS_H__
