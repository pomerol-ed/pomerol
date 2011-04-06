/** \file src/FieldOperatorContainer.h
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __FIELD_OPERATOR_CONTAINER__
#define __FIELD_OPERATOR_CONTAINER__

#include "config.h"
#include "FieldOperator.h"
#include "Hamiltonian.h"
#include "StatesClassification.h"

/** This class represents a container to store and retrieve FieldOperators ( CreationOperator or AnnihilationOperator 
 * rotated to eigenvector basis of Hamiltonian H ) for a given Index.
 * If no field operator is yet initialized then calculation of the field operator is done.
 */
class FieldOperatorContainer
{
private:
    /** A reference to a states classification object. */
    StatesClassification &S;
    /** A reference to a Hamiltonian. */
    Hamiltonian &H;
    /** A reference to a BitClassification object in order to check the input indices. */
    BitClassification &IndexInfo;
    /** A map which gives a link to the CreationOperator for a given index */
    std::map <ParticleIndex, CreationOperator*>     mapCreationOperators;
    /** A map which gives a link to the AnnihilationOperator for a given index */
    std::map <ParticleIndex, AnnihilationOperator*> mapAnnihilationOperators;
public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] IndexInfo A reference to a BitClassification
     */
    FieldOperatorContainer(StatesClassification &S, 
        Hamiltonian &H, BitClassification &IndexInfo);

    /** Returns the CreationOperator for a given Index */
    CreationOperator& getCreationOperator(ParticleIndex in);
    /** Returns the AnnihilationOperator for a given Index */
    AnnihilationOperator& getAnnihilationOperator(ParticleIndex in);
};

#endif // endif :: #ifndef __FIELD_OPERATOR_CONTAINER__
