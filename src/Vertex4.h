/** \file src/Vertex4.h
** \brief Irreducible two-particle vertex in the Matsubara representation.
**
** \author Igor Krivenko (igor@shg.ru)
*/
#ifndef ____VERTEX4____
#define ____VERTEX4____

#include "GreensFunction.h"
#include "TwoParticleGF.h"

/** Objects of this class just transforms a two-particle Green's function into
 * an irreducible vertex part or into an amputated irreducible vertex.
 */
class Vertex4 {

    /** A reference to a two-particle Green's function. */
    TwoParticleGF &Chi;
    /** A reference to a Green's function corresponding to the first index of Chi. */
    GreensFunction &g1;
    /** A reference to a Green's function corresponding to the second index of Chi. */
    GreensFunction &g2;
    /** A reference to a Green's function corresponding to the third index of Chi. */
    GreensFunction &g3;
    /** A reference to a Green's function corresponding to the fourth index of Chi. */
    GreensFunction &g4;

    /** 'Bits' (indices) of the two-particle Green's functions. */
    unsigned short ChiBit1, ChiBit2, ChiBit3, ChiBit4;

public:
    /** Constructor.
     * \param[in] Chi A reference to a two-particle Green's function.
     * \param[in] g1 A reference to a Green's function corresponding to the first index of Chi.
     * \param[in] g2 A reference to a Green's function corresponding to the second index of Chi.
     * \param[in] g3 A reference to a Green's function corresponding to the third index of Chi.
     * \param[in] g4 A reference to a Green's function corresponding to the fourth index of Chi.
     */
    Vertex4(TwoParticleGF &Chi, GreensFunction &g1, GreensFunction &g2, GreensFunction &g3, GreensFunction &g4);

    /** Returns the value of the irreducible vertex calculated at given frequencies.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);
    /** Returns the value of the amputated irreducible vertex calculated at given frequencies.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType getAmputated(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);

    //void dumpMatsubara(unsigned short points);
    //void dumpAmputatedMatsubara(unsigned short points);
};

#endif // endif :: #ifndef ____VERTEX4____
