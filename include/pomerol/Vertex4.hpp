/** \file include/pomerol/Vertex4.h
** \brief Irreducible two-particle vertex in the Matsubara representation.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_VERTEX4_H
#define POMEROL_INCLUDE_VERTEX4_H

#include "ComputableObject.hpp"
#include "GreensFunction.hpp"
#include "MatsubaraContainers.hpp"
#include "Misc.hpp"
#include "Thermal.hpp"
#include "TwoParticleGF.hpp"

namespace Pomerol {

class Vertex4 : public Thermal, public ComputableObject {

    TwoParticleGF const& Chi4;
    GreensFunction const& G13;
    GreensFunction const& G24;
    GreensFunction const& G14;
    GreensFunction const& G23;

    /** Storage for precomputed values. */
    mutable MatsubaraContainer4<Vertex4> Storage;
    friend class MatsubaraContainer4<Vertex4>;

public:
    Vertex4(TwoParticleGF const& Chi4,
            GreensFunction const& G13,
            GreensFunction const& G24,
            GreensFunction const& G14,
            GreensFunction const& G23);

    void compute(long NumberOfMatsubaras = 0);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
    ComplexType value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_VERTEX4_H
