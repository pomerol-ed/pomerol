/** \file include/pomerol/GFContainer.h
** \brief Storage of GF for multiple indices (obsolete, remove)
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/


#ifndef __INCLUDE_GFCONTAINER_H
#define __INCLUDE_GFCONTAINER_H

#include"Misc.h"
#include"GreensFunction.h"
#include"FieldOperatorContainer.h"
#include"IndexContainer2.h"

#include <memory>

namespace Pomerol{

template<bool Complex = false>
class GFContainer: public IndexContainer2<GreensFunction<Complex>, GFContainer<Complex>>, public Thermal
{

public:

    using ElementT = GreensFunction<Complex>;
    using ContainerBase = IndexContainer2<ElementT, GFContainer<Complex>>;

    GFContainer(const IndexClassification<Complex>& IndexInfo,
                const StatesClassification<Complex> &S,
                const Hamiltonian<Complex> &H,
                const DensityMatrix<Complex> &DM,
                const FieldOperatorContainer<Complex>& Operators);


    void prepareAll(const std::set<IndexCombination2>& InitialIndices = std::set<IndexCombination2>());
    void computeAll();

protected:

    friend ContainerBase;
    ElementT* createElement(const IndexCombination2& Indices) const;

    const StatesClassification<Complex> &S;

    const Hamiltonian<Complex> &H;
    const DensityMatrix<Complex> &DM;
    const FieldOperatorContainer<Complex> &Operators;
};

extern template class GFContainer<false>;
extern template class GFContainer<true>;

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GFCONTAINER_H
