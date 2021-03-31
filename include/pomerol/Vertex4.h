/** \file include/pomerol/Vertex4.h
** \brief Irreducible two-particle vertex in the Matsubara representation.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_VERTEX4_H
#define __INCLUDE_VERTEX4_H

#include"Misc.h"
#include"GFContainer.h"
#include"TwoParticleGF.h"
#include"MatsubaraContainers.h"

namespace Pomerol{

template<bool Complex = false>
class Vertex4 : public Thermal, public ComputableObject {

    TwoParticleGF<Complex> &Chi4;
    GreensFunction<Complex> &G13;
    GreensFunction<Complex> &G24;
    GreensFunction<Complex> &G14;
    GreensFunction<Complex> &G23;

    /** Storage for precomputed values. */
    mutable MatsubaraContainer4<Vertex4> Storage;
    friend class MatsubaraContainer4<Vertex4>;


public:

    Vertex4(TwoParticleGF<Complex>& Chi4,
            GreensFunction<Complex>& G13, GreensFunction<Complex>& G24,
            GreensFunction<Complex>& G14, GreensFunction<Complex>& G23);

    void compute(long NumberOfMatsubaras = 0);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
    ComplexType value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    bool isVanishing(void) const;
};

extern template class Vertex4<false>;
extern template class Vertex4<true>;

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_VERTEX4_H
