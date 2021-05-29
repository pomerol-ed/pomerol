/** \file include/pomerol/GFContainer.h
** \brief Storage of GF for multiple indices (obsolete, remove)
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/


#ifndef __INCLUDE_GFCONTAINER_H
#define __INCLUDE_GFCONTAINER_H

#include"Misc.h"
#include"Index.h"
#include"IndexClassification.h"
#include"GreensFunction.h"
#include"FieldOperatorContainer.h"
#include"IndexContainer2.h"

#include <memory>

namespace Pomerol {

class GFContainer: public IndexContainer2<GreensFunction,GFContainer>, public Thermal
{

public:
    template<typename... IndexTypes>
    GFContainer(const IndexClassification<IndexTypes...>& IndexInfo,
                const StatesClassification &S,
                const Hamiltonian &H,
                const DensityMatrix &DM,
                const FieldOperatorContainer& Operators) :
        IndexContainer2<GreensFunction, GFContainer>(this, IndexInfo),
        Thermal(DM), S(S), H(H), DM(DM), Operators(Operators)
    {}

    void prepareAll(const std::set<IndexCombination2>& InitialIndices = {});
    void computeAll();

protected:

    friend class IndexContainer2<GreensFunction,GFContainer>;
    std::shared_ptr<GreensFunction> createElement(const IndexCombination2& Indices) const;

    const StatesClassification &S;

    const Hamiltonian &H;
    const DensityMatrix &DM;
    const FieldOperatorContainer &Operators;
};

} // namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GFCONTAINER_H
