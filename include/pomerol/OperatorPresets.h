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

class N : public Operator {
private:
    const ParticleIndex Nmodes;
public:
    N(ParticleIndex Nmodes);
    std::map <FockState,MelemType> actRight(const FockState &ket) const;
    MelemType getMatrixElement(const FockState &bra, const FockState &ket) const;
    MelemType getMatrixElement(const FockState &ket) const;
};

class Sz : public Operator {
private:
    const int Nmodes;
    std::vector<ParticleIndex> SpinUpIndices; 
    std::vector<ParticleIndex> SpinDownIndices; 
    void generateTerms();
public:
    Sz(ParticleIndex Nmodes, const std::vector<ParticleIndex> & SpinUpIndices); 
    Sz(const std::vector<ParticleIndex> & SpinUpIndices, const std::vector<ParticleIndex> & SpinDownIndices);
    std::map <FockState,MelemType> actRight(const FockState &ket) const;
    MelemType getMatrixElement(const FockState &bra, const FockState &ket) const;
    MelemType getMatrixElement(const FockState &ket) const;
};

 
class Cdag : public Operator {
private:
    ParticleIndex index;
public:
    Cdag(ParticleIndex index):Operator(c_dag(index)){};
};

class C : public Operator {
private:
    ParticleIndex index;
public:
    C(ParticleIndex index):Operator(c(index)){};
};

class N_offdiag : public Operator {
private:
    ParticleIndex index1, index2;
public:
    N_offdiag(ParticleIndex index1, ParticleIndex index2):Operator(n_offdiag(index1, index2)){};
};



} // end of namespace OperatorPresets
} // end of namespace Pomerol

#endif // endif :: #ifndef __INCLUDE_OPERATOR_PRESETS_H__
