// Include the pomerol library
#include <pomerol.h>

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

    Lattice L;
    // Add a site with a name "C", that has 1 orbital and 2 spins.
    L.addSite(new Lattice::Site("C",1,2));

    LatticePresets::addCoulombS(&L, "C", U, -mu);
    LatticePresets::addHopping(&L, "C", "C", J, 0, 0, up, down);

    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare();
    IndexInfo.printIndices();

    IndexHamiltonian HStorage(&L,IndexInfo);
    HStorage.prepare();

    Symmetrizer Symm(IndexInfo, HStorage);
    Symm.compute(true);

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    Hamiltonian H(IndexInfo, HStorage, S);
    H.prepare();
    H.compute(MPI_COMM_WORLD);

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, S, H);
    Operators.prepareAll();
    Operators.computeAll();

    IndexInfo.printIndices();

    std::vector<spin> all_spins;
    all_spins.push_back(down);
    all_spins.push_back(up);

    // Reference
    GF_ref ref(J);

    BOOST_FOREACH(spin s1, all_spins){
    BOOST_FOREACH(spin s2, all_spins){
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

    return (
     run_test(ComplexType( 0.1,0)) &&
     run_test(ComplexType(-0.1,0)) &&
     run_test(ComplexType(0, 0.1)) &&
     run_test(ComplexType(0,-0.1)) &&
     run_test(ComplexType(0.1, 0.1)) &&
     run_test(ComplexType(0.1,-0.1)) &&
     run_test(ComplexType(-0.1,0.1)) &&
     run_test(ComplexType(-0.1,-0.1))
     ) ? EXIT_SUCCESS : EXIT_FAILURE;
}
