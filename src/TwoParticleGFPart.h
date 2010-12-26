#ifndef ____DEFINE_2PGF_PART____
#define ____DEFINE_2PGF_PART____

#include <list>

#include "config.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "FieldOperator.h"
#include "DensityMatrixPart.h"

class TwoParticleGFPart {

public:

    enum ComputationMethod {ChasingIndices1, ChasingIndices2};

    struct TwoParticleGFTerm{
        size_t z1, z2, z3;

        ComplexType CoeffZ2;
        ComplexType CoeffZ4;
        ComplexType CoeffZ1Z2Res;
        ComplexType CoeffZ1Z2NonRes;
        ComplexType CoeffZ2Z3Res;
        ComplexType CoeffZ2Z3NonRes;

        ComplexType Poles[3];

        TwoParticleGFTerm(ComplexType Coeff, RealType beta,
                         RealType Ei, RealType Ej, RealType Ek, RealType El,
                         RealType Wi, RealType Wj, RealType Wk, RealType Wl,
                         Permutation3& Permutation);

        // \frac{1}{(z1-Poles[1])(z3-Poles[3])}*
        //      (\frac{CoeffZ4}{z1+z2+z3-Poles[1]-Poles[2]-Poles[3]} - \frac{CoeffZ2}{z2-Poles[2]}
        //      + CoeffZ1Z2Res*IsZ1Z2Resonant*\delta(z1+z2)
        //      + CoeffZ1Z2NonRes*\frac{1 - IsZ1Z2Resonant*\delta(z1+z2)}{z1+z2-Poles[1]-Poles[2]}
        //      + CoeffZ2Z3Res*IsZ2Z3Resonant*\delta(z2+z3)
        //      + CoeffZ2Z3NonRes*\frac{1 - IsZ2Z3Resonant*\delta(z2+z3)}{z2+z3-Poles[2]-Poles[3]})
        ComplexType operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3) const;

        static const RealType MatrixElementTolerance = 1e-10;
        static const RealType ResonanceTolerance = 1e-16;
        static bool IsRelevant(const ComplexType &MatrixElement);
    };

    /**
     * A miniclass to store value of Chi over Matsubara frequencies. Stores data in a (OMEGA,nu,nu'), where OMEGA=w1+w2 - bosonic frequency
     * and nu=w1, nu'=w4
     */
    class MatsubaraContainer{

        ComplexType MatsubaraSpacing;
        long NumberOfMatsubaras;
        std::vector<MatrixType> Data;
        std::vector<long> FermionicFirstIndex;
    public:
        MatsubaraContainer(RealType beta);
        void prepare(long NumberOfMatsubaras);
        ComplexType& operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);
        void fill(std::list<TwoParticleGFTerm> &Terms);
        MatsubaraContainer& operator+=(const MatsubaraContainer& rhs);
        void clear();
    };

private:

    FieldOperatorPart& O1;
    FieldOperatorPart& O2;
    FieldOperatorPart& O3;
    CreationOperatorPart& CX4;

    HamiltonianPart& Hpart1;
    HamiltonianPart& Hpart2;
    HamiltonianPart& Hpart3;
    HamiltonianPart& Hpart4;

    DensityMatrixPart& DMpart1; 
    DensityMatrixPart& DMpart2;
    DensityMatrixPart& DMpart3;
    DensityMatrixPart& DMpart4;

    Permutation3 Permutation;

    std::list<TwoParticleGFTerm> Terms;

    MatsubaraContainer *Storage;

    void computeChasing1(long NumberOfMatsubaras);
    void computeChasing2(long NumberOfMatsubaras);

public:
    TwoParticleGFPart(FieldOperatorPart& O1, FieldOperatorPart& O2, FieldOperatorPart& O3, CreationOperatorPart& CX4,
                HamiltonianPart& Hpart1, HamiltonianPart& Hpart2, HamiltonianPart& Hpart3, HamiltonianPart& Hpart4,
                DensityMatrixPart& DMpart1, DensityMatrixPart& DMpart2, DensityMatrixPart& DMpart3, DensityMatrixPart& DMpart4,
                Permutation3 Permutation);

    void compute(long NumberOfMatsubaras, ComputationMethod method = ChasingIndices2);
    void clear();
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);
    const MatsubaraContainer& getMatsubaraContainer();
};

#endif // endif :: #ifndef ____DEFINE_2PGF_PART____
