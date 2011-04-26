/** \file src/GreensFunction.h
** \brief Thermal Green's function.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_GREENSFUNCTION_H
#define __INCLUDE_GREENSFUNCTION_H

#include <sstream>
#include <fstream>
//#include <iomanip>
//#include <cmath>

#include"Misc.h"
#include"ComputableObject.h"
#include"output.h"
#include"StatesClassification.h"
#include"FieldOperator.h"
#include"DensityMatrix.h"
#include"GreensFunctionPart.h"

/** This class represents a thermal Green's function in the Matsubara representation.
 *
 * Exact definition:
 * 
 * \f[
 *      G(\omega_n) = -\int_0^\beta \langle\mathbf{T}c_i(\tau)c^+_j(0)\rangle e^{i\omega_n\tau} d\tau
 * \f]
 * 
 * It is actually a container class for a collection of parts (most of real calculations
 * take place inside the parts). A pair of parts, one part of an annihilation operator and
 * another from a creation operator, corresponds to a part of the Green's function.
 */
class GreensFunction : public ComputableObject {
    /** A reference to a states classification object. */
    StatesClassification& S;
    /** A reference to a Hamiltonian. */
    Hamiltonian& H;
    /** A reference to an annihilation operator. */
    AnnihilationOperator& C;
    /** A reference to a creation operator. */
    CreationOperator& CX;
    /** A reference to a density matrix. */
    DensityMatrix& DM;
    /** A flag to represent if Greens function vanishes, i.e. identical to 0 */
    bool vanish;

    /** A list of pointers to parts (every part corresponds to a part of the annihilation operator
     * and a part of the creation operator).
     */
    std::list<GreensFunctionPart*> parts;

public:
     /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] C A reference to an annihilation operator.
     * \param[in] CX A reference to a creation operator.
     * \param[in] DM A reference to a density matrix.
     */
    GreensFunction(StatesClassification& S, Hamiltonian& H,
                   AnnihilationOperator& C, CreationOperator& CX, DensityMatrix& DM);
    /** Destructor. */
    ~GreensFunction();

    /** Chooses relevant parts of C and CX and allocates resources for the parts of the Green's function. */
    void prepare(void);
    /** Actually computes the parts. */
    void compute(void);

    /** Returns the 'bit' (index) of the operator C or CX.
     * \param[in] Position Use C for Position==0 and CX for Position==1.
     */
    unsigned short getIndex(size_t Position) const;

     /** Returns the value of the Green's function calculated at a given frequency.
     * \param[in] MatsubaraNum Number of the Matsubara frequency (\f$ \omega_n = \pi(2n+1)/\beta \f$).
     */
    ComplexType operator()(long MatsubaraNum);

    void dumpToPlainText(long wn);
    /** Returns the path of the output directory associated with this Green's function. */
    std::string getPath();
    /** Dumps the Green's function for a range of the Matsubara frequencies.
     * \param[in] point The number of points in the range.
     */
    void dumpMatsubara(unsigned short points);
    /** Returns true if current Greens function is identical to zero */
    bool vanishes();
};

#endif // endif :: #ifndef __INCLUDE_GREENSFUNCTION_H

