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

typedef std::shared_ptr<GreensFunction> GFPointer;

class GFContainer: public IndexContainer2<GreensFunction,GFContainer>, public Thermal
{

public:
    GFContainer(const IndexClassification& IndexInfo,
                const StatesClassification &S,
                const Hamiltonian &H, const DensityMatrix &DM, const FieldOperatorContainer& Operators);


    void prepareAll(const std::set<IndexCombination2>& InitialIndices = std::set<IndexCombination2>());
    void computeAll();

protected:

    friend class IndexContainer2<GreensFunction,GFContainer>;
    GreensFunction* createElement(const IndexCombination2& Indices) const;

    const StatesClassification &S;

    const Hamiltonian &H;
    const DensityMatrix &DM;
    const FieldOperatorContainer &Operators;
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GFCONTAINER_H
