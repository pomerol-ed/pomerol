#include "Vertex4Part.h"

Vertex4Part::Vertex4Part(AnnihilationOperatorPart& C0, AnnihilationOperatorPart& C1,
                         CreationOperatorPart& CX2, CreationOperatorPart& CX3,
                         HamiltonianPart& HpartOuter, HamiltonianPart& Hpart23,
                         HamiltonianPart& Hpart12, HamiltonianPart& Hpart01,
                         DensityMatrixPart& DMpartOuter, DensityMatrixPart& DMpart23,
                         DensityMatrixPart& DMpart12, DensityMatrixPart& DMpart01,
                         size_t PermutationNumber) :
C0(C0), C1(C1), CX2(CX2), CX3(CX3), 
HpartOuter(HpartOuter), Hpart23(Hpart23), Hpart12(Hpart12), Hpart01(Hpart01),
DMpartOuter(DMpartOuter), DMpart23(DMpart23), DMpart12(DMpart12), DMpart01(DMpart01),
PermutationNumber(PermutationNumber)
{
    // TODO
}

void Vertex4Part::compute(void)
{
    // TODO
}

ComplexType Vertex4Part::operator()(ComplexType Frequency0, ComplexType Frequency1, ComplexType Frequency2)
{
    // TODO
    return 0;
}