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
#include"TwoParticleGFContainer.h"
#include"MatsubaraContainers.h"

namespace Pomerol{

class Vertex4 : public Thermal, public ComputableObject {

    TwoParticleGF &Chi4;
    GreensFunction &G13;
    GreensFunction &G24;
    GreensFunction &G14;
    GreensFunction &G23;

    /** Storage for precomputed values. */
    mutable MatsubaraContainer4<Vertex4> Storage;
    friend class MatsubaraContainer4<Vertex4>;


public:

    Vertex4(TwoParticleGF& Chi4,
            GreensFunction& G13, GreensFunction& G24,
            GreensFunction& G14, GreensFunction& G23);

    void compute(long NumberOfMatsubaras = 0);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
    ComplexType value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    bool isVanishing(void) const;
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_VERTEX4_H
