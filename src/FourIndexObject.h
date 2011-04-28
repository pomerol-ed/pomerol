/** \file src/FourIndexObject.h
** \brief A prototype class for all objects, depending on four indices
**
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_FOURINDEXOBJECT_H 
#define __INCLUDE_FOURINDEXOBJECT_H 

#include "Misc.h"
#include "TwoParticleGFPart.h"

/** This class is a prototype for every object which is defined by a combination of 4 arbitrary Particle Indices
 *  It defines main subobjects, which are used by the 4 index-dependent quantites: IndexCombination and MatsubaraContainer
 */
class FourIndexObject {
public:
    /** A combination of four indices. The notation is ccc^*c^* */
    struct IndexCombination;
    /** A storage of an object over Matsubara frequencies */
    class MatsubaraContainer;
};

/** A prototype class for objects which are constructed for a given set of 4 indices, i.e Chi_{ijkl} */
class FourIndexSingleObject : public FourIndexObject {
public:
    /** Return the combination of indices, for which the object is constructed */
    IndexCombination getIndices();

    /** Return the value for a given Matsubara numbers (not frequencies themselves)
     * \param[in] MatsubaraNumber1 An index of the 1st Matsubara frequency.
     * \param[in] MatsubaraNumber2 An index of the 2nd Matsubara frequency.
     * \param[in] MatsubaraNumber3 An index of the 3rd Matsubara frequency.
     */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
    };

/** A prototype container class of Four Indices Objects - stores all values for all possible Particle Indices 
 * The notation is ccc^*c^*
 */ 
class FourIndexContainerObject : public FourIndexObject {
protected:
    /** A set of trivial permutations, which just lead to a change of sign */
    static const Permutation4 TrivialOperatorPermutations[];
public:
    /** Returns the value for a given set of indices and Matsubara numbers (not frequencies themselves)
     * \param[in] In An IndexCombination at which the value should be get.
     * \param[in] MatsubaraNumber1 An index of the 1st Matsubara frequency.
     * \param[in] MatsubaraNumber2 An index of the 2nd Matsubara frequency.
     * \param[in] MatsubaraNumber3 An index of the 3rd Matsubara frequency.
     */
    ComplexType operator()(const IndexCombination& In, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);
};

/** A structure to handle a combination of 4 Particle Index. The notation is ccc^*c^* */
struct FourIndexObject::IndexCombination
{
    /** Actual indices */
    ParticleIndex Indices[4];
    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& output, const FourIndexObject::IndexCombination& out);
    /** Operator < - comparison method for IndexCombination */
    bool operator< (const FourIndexObject::IndexCombination& rhs) const ;
    /** Operator == */
    bool operator==(const FourIndexObject::IndexCombination& rhs) const ;
    /** Operator != */
    bool operator!=(const FourIndexObject::IndexCombination& rhs) const ;
    /** Constructor
     * \param[in] cindex1 - Index of a 1st operator ( annihilation )
     * \param[in] cindex2 - Index of a 2nd operator ( annihilation )
     * \param[in] cdagindex3 - Index of a 3rd operator ( creation )
     * \param[in] cdagindex4 - Index of a 4th operator ( creation )
     */
    IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4);
};


/**
* A subclass of FourIndexObject - it stores values of a 4 FourIndexObject  over Matsubara frequencies. 
* The data is stored in a 1 bosonic and 2 fermionic degrees flavour : 
* (OMEGA,nu,nu'), 
* where OMEGA=w1+w2 - bosonic frequency, nu=w1, nu'=w4
*/
class FourIndexObject::MatsubaraContainer 
{
    /** This is a i\frac{\pi}{\beta} - an interval between 2 adjacent matsubaras. 
     * It is more convenient to store Inverse Temperature with a spacing, rather than to invert it afterwards 
     */
    ComplexType MatsubaraSpacing;

    /* An amount of positive Matsubaras on which to store the values. The range of values [-MatsubaraSpacing;MatsubaraSpacing-1] will be available */
    long NumberOfMatsubaras;

    /* The storage - an array of Matrices for nu,nu' space, stored as a vector which is dependent on a bosonic frequency index */
    std::vector<MatrixType> Data;
    /* A first non-vanishing fermionic index - used for correspondence between number of element in a matrix and real Matsubara numbers */
    std::vector<long> FermionicFirstIndex;
public:
    /** Constructor
     * \param[in] beta An inverse temperature.
     */
    MatsubaraContainer(RealType beta);

    /** Allocate memory for a storage
     * \param[in] NumberOfMatsubaras An amount of positive Matsubara frequencies that will be held in a MatsubaraContainer.
     */
    void prepare(long NumberOfMatsubaras);

    /** Returns the value for a given Matsubara numbers (not frequencies themselves)
     * \param[in] MatsubaraNumber1 An index of the 1st Matsubara frequency.
     * \param[in] MatsubaraNumber2 An index of the 2nd Matsubara frequency.
     * \param[in] MatsubaraNumber3 An index of the 3rd Matsubara frequency.
     */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /** Fill container from a list of terms 
     * \param[in] NonResonantTerms A list of NonResonant Terms.
     * \param[in] ResonantTerms A list of Resonant Terms.
     * \param[in] Permutation A permutation of input Matsubara frequencies
     */
    void fill(const std::list<TwoParticleGFPart::NonResonantTerm>& NonResonantTerms, const std::list<TwoParticleGFPart::ResonantTerm>& ResonantTerms, Permutation3 Permutation);

    /** Operator+= : adds to a current MatsubaraContainer another one
     * \param[in] rhs Right hand side of the equation Matsubara Container to add.
     */
    MatsubaraContainer& operator+=(const MatsubaraContainer& rhs);

    /** Empty memory */
    void clear();
};

#endif // endif :: #ifndef __INCLUDE_FOURINDEXOBJECT_H 
