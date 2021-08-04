/** \file include/pomerol/GFContainer.h
** \brief Storage of GF for multiple indices (obsolete, remove)
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_GFCONTAINER_H
#define POMEROL_INCLUDE_POMEROL_GFCONTAINER_H

#include "DensityMatrix.hpp"
#include "FieldOperatorContainer.hpp"
#include "GreensFunction.hpp"
#include "Hamiltonian.hpp"
#include "Index.hpp"
#include "IndexClassification.hpp"
#include "IndexContainer2.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

#include <memory>
#include <set>

namespace Pomerol {

class GFContainer: public IndexContainer2<GreensFunction,GFContainer>, public Thermal
{

public:
    template<typename... IndexTypes>
    GFContainer(IndexClassification<IndexTypes...> const& IndexInfo,
                StatesClassification const& S,
                Hamiltonian const& H,
                DensityMatrix const& DM,
                FieldOperatorContainer const& Operators) :
        IndexContainer2<GreensFunction, GFContainer>(*this, IndexInfo),
        Thermal(DM), S(S), H(H), DM(DM), Operators(Operators)
    {}

    void prepareAll(std::set<IndexCombination2> const& InitialIndices = {});
    void computeAll();

protected:

    friend class IndexContainer2<GreensFunction,GFContainer>;
    std::shared_ptr<GreensFunction> createElement(IndexCombination2 const& Indices) const;

    StatesClassification const& S;

    Hamiltonian const& H;
    DensityMatrix const& DM;
    FieldOperatorContainer const& Operators;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_GFCONTAINER_H
