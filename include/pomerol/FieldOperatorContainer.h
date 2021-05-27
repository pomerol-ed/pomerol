/** \file include/pomerol/FieldOperatorContainer.h
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_FIELDOPERATORCONTAINER_H
#define __INCLUDE_FIELDOPERATORCONTAINER_H

#include"Misc.h"
#include"IndexClassification.h"
#include"MonomialOperator.h"
#include"Hamiltonian.h"
#include"StatesClassification.h"

#include <set>
#include <unordered_map>

namespace Pomerol {

/** This class represents a container to store and retrieve FieldOperators ( CreationOperator or AnnihilationOperator
 * rotated to eigenvector basis of Hamiltonian H ) for a given Index.
 * If no field operator is yet initialized then calculation of the field operator is done.
 */
class FieldOperatorContainer
{

    /** A map which gives a link to the CreationOperator for a given index */
    std::unordered_map<ParticleIndex, CreationOperator> mapCreationOperators;
    /** A map which gives a link to the AnnihilationOperator for a given index */
    std::unordered_map<ParticleIndex, AnnihilationOperator> mapAnnihilationOperators;

public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] IndexInfo A reference to a IndexClassification
     */
    template<typename... IndexTypes>
    FieldOperatorContainer(const IndexClassification<IndexTypes...> &IndexInfo,
                           const HilbertSpace<IndexTypes...> &HS,
                           StatesClassification &S,
                           const Hamiltonian &H,
                           std::set<ParticleIndex> in = {})
    {
        if(in.empty()) {
            for (ParticleIndex p = 0; p < IndexInfo.getIndexSize(); ++p) {
                in.insert(p);
            }
        }
        for(auto p : in) {
            mapCreationOperators.emplace(p, CreationOperator(IndexInfo, HS, S, H, p));
            mapAnnihilationOperators.emplace(p, AnnihilationOperator(IndexInfo, HS, S, H, p));
        }
    }

    template<typename... IndexTypes>
    void prepareAll(const HilbertSpace<IndexTypes...> &HS)
    {
        for(auto & CX : mapCreationOperators)
            CX.second.prepare(HS);
        for(auto & C : mapAnnihilationOperators)
            C.second.prepare(HS);
    }
    void computeAll();

    /** Returns the CreationOperator for a given Index. Makes on-demand computation. */
    const CreationOperator& getCreationOperator(ParticleIndex in) const;
    /** Returns the AnnihilationOperator for a given Index. Makes on-demand computation */
    const AnnihilationOperator& getAnnihilationOperator(ParticleIndex in) const;
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_FIELDOPERATORCONTAINER_H
