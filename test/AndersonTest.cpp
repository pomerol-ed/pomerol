#include "Misc.hpp"
#include "LatticePresets.hpp"
#include "Operators.hpp"
#include "IndexClassification.hpp"
#include "HilbertSpace.hpp"
#include "StatesClassification.hpp"
#include "HamiltonianPart.hpp"
#include "Hamiltonian.hpp"
#include "DensityMatrix.hpp"
#include "FieldOperatorContainer.hpp"
#include "GFContainer.hpp"

#include "./Utility.hpp"

#include <vector>

using namespace Pomerol;

RealType beta = 20;
RealType U = 3.7;
RealType mu = U/2*0;
RealType h = 0.0;
RealType V = 1.0;
RealType epsilon = 2.3;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    using namespace LatticePresets;

    auto HExpr = CoulombS("C", U, -mu);
    HExpr += Magnetization("C", 2*h);
    HExpr += Level("0", -epsilon);
    HExpr += Level("1", epsilon);
    HExpr += Hopping("C", "0", V);
    HExpr += Hopping("C", "1", V);

    auto IndexInfo = MakeIndexClassification(HExpr);
    print_section("Indices");
    std::cout << IndexInfo << std::endl;

    INFO("Hamiltonian");
    INFO(HExpr);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    for (BlockNumber i=0; i < S.getNumberOfBlocks(); i++) {
        std::vector<QuantumState> st = S.getFockStates(i);
        for (int i=0; i<st.size(); ++i) INFO(st[i]);
        INFO("");
    };

    H.compute(MPI_COMM_WORLD);

    DensityMatrix rho(S, H, beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex down_index = IndexInfo.getIndex("C", 0, down);
    ParticleIndex up_index = IndexInfo.getIndex("C", 0, up);

    DEBUG(down_index);
    DEBUG(up_index);

    std::cout << IndexInfo << std::endl;
    for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); i++) {
        INFO("C^+_"<<i);
    }

    GreensFunction GF_down(S, H,
        Operators.getAnnihilationOperator(down_index),
        Operators.getCreationOperator(down_index),
    rho);

    GreensFunction GF_up(S, H,
        Operators.getAnnihilationOperator(up_index),
        Operators.getCreationOperator(up_index),
    rho);

    GF_down.prepare(); DEBUG(""); GF_up.prepare();
    GF_down.compute(); DEBUG(""); GF_up.compute();

    for (size_t i=0; i<10; i++) {
        INFO(GF_down(i) << " " << GF_up(i));
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
