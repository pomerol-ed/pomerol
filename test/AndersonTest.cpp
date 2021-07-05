#include "Misc.hpp"
#include "LatticePresets.hpp"
#include "IndexClassification.hpp"
#include "HilbertSpace.hpp"
#include "StatesClassification.hpp"
#include "Hamiltonian.hpp"
#include "DensityMatrix.hpp"
#include "FieldOperatorContainer.hpp"
#include "GreensFunction.hpp"

#undef INFO // Catch2 has its own INFO() macro
#include "catch2/catch.hpp"

using namespace Pomerol;

TEST_CASE("Anderson model with 2 bath sites", "[Anderson]") {
    RealType U = 3.7;
    RealType mu = 0.6 * U;
    RealType h = 0.1;
    RealType V = 1.0;
    RealType epsilon = 2.3;
    RealType beta = 20;

    // Reference Green's functions
    ComplexVectorType G_ref_up(10);
    G_ref_up << -0.7545439 -0.14723373*I,
                -0.59517353-0.34478922*I,
                -0.42646689-0.4031622*I,
                -0.30758605-0.39519013*I,
                -0.23068661-0.36811109*I,
                -0.18013326-0.33885065*I,
                -0.14541077-0.31208195*I,
                -0.12043127-0.28860297*I,
                -0.10171394-0.26815136*I,
                -0.08721336-0.25025992*I;
    ComplexVectorType G_ref_down(10);
    G_ref_down << 0.49196891-0.07241433*I,
                  0.44396903-0.18681652*I,
                  0.37248532-0.24764566*I,
                  0.30425235-0.26969548*I,
                  0.24921656-0.27135953*I,
                  0.20705484-0.26399883*I,
                  0.1748209 -0.25319165*I,
                  0.14976473-0.2414213*I,
                  0.12986709-0.22973732*I,
                  0.11373954-0.21855949*I;

    using namespace LatticePresets;

    auto HExpr = CoulombS("C", U, -mu);
    HExpr += Magnetization("C", h);
    HExpr += Level("0", -epsilon);
    HExpr += Level("1", epsilon);
    HExpr += Hopping("C", "0", V);
    HExpr += Hopping("C", "1", V);
    INFO("Hamiltonian\n" << HExpr);

    auto IndexInfo = MakeIndexClassification(HExpr);
    INFO("Indices\n" << IndexInfo);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);
    INFO("Energy levels " << H.getEigenValues());
    INFO("The value of ground energy is " << H.getGroundEnergy());

    DensityMatrix rho(S, H, beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex C_down_index = IndexInfo.getIndex("C", 0, down);
    ParticleIndex C_up_index = IndexInfo.getIndex("C", 0, up);

    GreensFunction GF_down(S, H,
        Operators.getAnnihilationOperator(C_down_index),
        Operators.getCreationOperator(C_down_index),
    rho);

    GreensFunction GF_up(S, H,
        Operators.getAnnihilationOperator(C_up_index),
        Operators.getCreationOperator(C_up_index),
    rho);

    GF_down.prepare(); GF_up.prepare();
    GF_down.compute(); GF_up.compute();

    for(int n = 0; n < G_ref_up.size(); ++n) {
        auto result = GF_up(n);
        auto ref = G_ref_up[n];
        REQUIRE(result.real() == Approx(ref.real()).margin(1e-12));
        REQUIRE(result.imag() == Approx(ref.imag()).margin(1e-12));
    }
    for(int n = 0; n < G_ref_down.size(); ++n) {
        auto result = GF_down(n);
        auto ref = G_ref_down[n];
        REQUIRE(result.real() == Approx(ref.real()).margin(1e-12));
        REQUIRE(result.imag() == Approx(ref.imag()).margin(1e-12));
    }
}
