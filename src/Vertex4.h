/** \file src/Vertex4.h
** \brief Irreducible two-particle vertex in the Matsubara representation.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef ____VERTEX4____
#define ____VERTEX4____

#include "GFContainer.h"
#include "TwoParticleGFContainer.h"

/** Objects of this class just transforms a two-particle Green's function into
 * an irreducible vertex part or into an amputated irreducible vertex.
 */
class Vertex4 {

    /** A reference to a two-particle Green's function. */
    TwoParticleGFContainer &Chi;
    /** A reference to a Green's function container */
    GFContainer &g;
    /** A reference to a bit classification object */
    const BitClassification &IndexInfo;
    
    /** Precomputed inverted Green's function matrices calculated at different Matsubara frequencies. */
    std::vector<MatrixType> InvertedGFs;

    /** Amount of computed matsubaras in TwoParticleGF */
    long NumberOfMatsubaras;

    /** A storage for unamputated values */
    std::map<TwoParticleGFContainer::IndexCombination,TwoParticleGFPart::MatsubaraContainer*> mapUnAmputatedValues;

    /** A vector of all nontrivial combinations to compute */
    std::vector<TwoParticleGFContainer::IndexCombination*> NonTrivialAmputatedCombinations;

    /** A storage for amputated values */
    std::map<TwoParticleGFContainer::IndexCombination,TwoParticleGFPart::MatsubaraContainer*> mapAmputatedValues;

public:
    /** Constructor.
     * \param[in] IndexInfo A reference to a bit classification object.
     * \param[in] Chi A reference to a two-particle Green's function.
     * \param[in] g1 A reference to a Green's function container.
     */
    Vertex4(const BitClassification &IndexInfo, TwoParticleGFContainer &Chi, GFContainer &g);

    //============================= UnAmputated methods ==============================//

    /** Do some preparation procedures : prepare storage */
    void prepareUnAmputated();

    /** Compute unamputated values
     */
    void computeUnAmputated();

     /** Returns the value of the unamputated irreducible vertex calculated at given frequencies.
     * \param[in] Requested indices.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType getUnAmputatedValue(const TwoParticleGFContainer::IndexCombination& in,
                             long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);


    //============================= Amputated methods ==============================//
     /** Do some preparation procedures : calculate inverted GF's, prepare storage */
    void prepareAmputated(std::vector<TwoParticleGFContainer::IndexCombination*>&);

     /** Returns the value of the amputated irreducible vertex calculated at given frequencies.
     * \param[in] Requested indices.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType getAmputatedValue(const TwoParticleGFContainer::IndexCombination& in,
                             long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);

    void computeAmputated();

    //==============================================================================//
    /** Returns the value of the irreducible vertex calculated at given frequencies.
     * \param[in] in Requested indices.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType operator()(const TwoParticleGFContainer::IndexCombination& in, 
                           long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);


private:
    bool vanishes(const TwoParticleGFContainer::IndexCombination& in); 

    /** Compute unamputated values
     * \param[in] Requested indices.
     */
    void computeUnAmputated(const TwoParticleGFContainer::IndexCombination& in);

    /** Compute amputated values
     * \param[in] Requested indices.
     */
    void computeAmputated(const TwoParticleGFContainer::IndexCombination& in);
};

#endif // endif :: #ifndef ____VERTEX4____
