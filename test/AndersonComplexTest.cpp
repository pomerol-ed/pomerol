// Include the pomerol library
#include <pomerol.h>

#include "./Utility.h"

using namespace Pomerol;

// Parameters
double mu = 1.2;
double U = 2.0;
double beta = 10;

struct GF_ref {
    std::vector<double> E;
    std::vector<double> w;
    double Z;
    double phi;
    GF_ref(ComplexType J) : phi(std::arg(J)), Z(0) {
        E.push_back(0);
        E.push_back(-mu-std::abs(J));
        E.push_back(-mu+std::abs(J));
        E.push_back(-2*mu+U);

        for(int i = 0; i < E.size(); ++i) {
            w.push_back(std::exp(-beta*E[i]));
            Z += w.back();
        }
        for(int i = 0; i < E.size(); ++i) w[i] /= Z;
    }

    ComplexType operator()(spin s1, spin s2, int n) const {
        ComplexType iw = ComplexType(0,M_PI*(2*n+1)/beta);

        ComplexType res = 0;
        if(s1 == up && s2 == up) {
            res += 0.5*(w[1] + w[0])/(iw - (E[1] - E[0]));
            res += 0.5*(w[2] + w[0])/(iw - (E[2] - E[0]));
            res += 0.5*(w[3] + w[1])/(iw - (E[3] - E[1]));
            res += 0.5*(w[3] + w[2])/(iw - (E[3] - E[2]));
        } else if(s1 == up && s2 == down) {
            ComplexType u = std::exp(ComplexType(0,phi));
            res += -u*0.5*(w[1] + w[0])/(iw - (E[1] - E[0]));
            res += u*0.5*(w[2] + w[0])/(iw - (E[2] - E[0]));
            res += u*0.5*(w[3] + w[1])/(iw - (E[3] - E[1]));
            res += -u*0.5*(w[3] + w[2])/(iw - (E[3] - E[2]));
        } else if(s1 == down && s2 == up) {
            ComplexType u = std::exp(ComplexType(0,-phi));
            res += -u*0.5*(w[1] + w[0])/(iw - (E[1] - E[0]));
            res += u*0.5*(w[2] + w[0])/(iw - (E[2] - E[0]));
            res += u*0.5*(w[3] + w[1])/(iw - (E[3] - E[1]));
            res += -u*0.5*(w[3] + w[2])/(iw - (E[3] - E[2]));
        } else {
            res += 0.5*(w[1] + w[0])/(iw - (E[1] - E[0]));
            res += 0.5*(w[2] + w[0])/(iw - (E[2] - E[0]));
            res += 0.5*(w[3] + w[1])/(iw - (E[3] - E[1]));
            res += 0.5*(w[3] + w[2])/(iw - (E[3] - E[2]));
        }
        return res;
    }
};

bool run_test(ComplexType J /* spin-flip amplitude */) {

    INFO("J = " << J);

    using namespace LatticePresets;

    // Add a site with a name "C", that has 1 orbital and 2 spins.

    ComplexExpr HExpr = CoulombS("C", U, -mu);
    HExpr += Hopping("C", "C", J, 0, 0, up, down);

    int rank = pMPI::rank(MPI_COMM_WORLD);

    auto IndexInfo = MakeIndexClassification(HExpr);
    if (!rank) {
        print_section("Indices");
        IndexInfo.printIndices();
    }

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    // Reference
    GF_ref ref(J);

    for(spin s1 : {down, up}){
    for(spin s2 : {down, up}){
        INFO("s1 = " << s1 << " s2 = " << s2);
        ParticleIndex index1 = IndexInfo.getIndex("C",0,s1);
        ParticleIndex index2 = IndexInfo.getIndex("C",0,s2);
        GreensFunction GF(S,H,
            Operators.getAnnihilationOperator(index1),
            Operators.getCreationOperator(index2),
        rho);
        GF.prepare();
        GF.compute();

        for (int n=0; n<10; ++n)
            if(std::abs(GF(n) - ref(s1,s2,n)) > 1e-10) return false;
    }}
    return true;
}

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    bool result =
     run_test(ComplexType( 0.1,0)) &&
     run_test(ComplexType(-0.1,0)) &&
     run_test(ComplexType(0, 0.1)) &&
     run_test(ComplexType(0,-0.1)) &&
     run_test(ComplexType(0.1, 0.1)) &&
     run_test(ComplexType(0.1,-0.1)) &&
     run_test(ComplexType(-0.1,0.1)) &&
     run_test(ComplexType(-0.1,-0.1));

     MPI_Finalize();

     return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
