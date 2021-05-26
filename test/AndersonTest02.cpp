// Hubbard 2x2
// Antipov, 2013

#pragma clang diagnostic ignored "-Wc++11-extensions"
#pragma clang diagnostic ignored "-Wgnu"

#include <string>
#include <tuple>
#include <iostream>
#include <algorithm>

#include <cstdlib>
#include <fstream>

#include <pomerol.h>

#include "./Utility.h"

using namespace Pomerol;

/* Auxiliary routines - implemented in the bottom. */
template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    print_section("Anderson ");

    int wn;
    RealType U=0.5, mu = 0.25, beta = 26, reduce_tol = 1e-5, coeff_tol = 1e-8;
    bool calc_gf = true, calc_2pgf = true;
    size_t L=2;

    std::vector<double> levels(4);
    levels[0] = 1.02036910873357;
    levels[1] = -1.02036910873357;
    levels[2] = 0.140037222821207;
    levels[3] =  -0.140037222821207;
    std::vector<double> hoppings(4);
    hoppings[0] = 0.296439333614347;
    hoppings[1] = 0.296439333614347;
    hoppings[2] =  0.229348742868022;
    hoppings[3] = 0.229348742868022;
    L = std::min(L,levels.size());

    INFO("Diagonalization of 1+" << L << " sites");

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu);

    std::vector<std::string> names(L);
    for (size_t i=0; i<L; i++)
        {
            std::stringstream s; s << i;
            names[i] = "b" + s.str();
            HExpr += Hopping("A", names[i], hoppings[i]);
            HExpr += Level(names[i], levels[i]);
        };

    int rank = pMPI::rank(MPI_COMM_WORLD);
    if (!rank) {
        INFO("Hamiltonian");
        INFO(HExpr);
    }

    auto IndexInfo = MakeIndexClassification(HExpr);
    if (!rank) {
        print_section("Indices");
        std::cout << IndexInfo << std::endl;
    }

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    RealVectorType evals (H.getEigenValues());
    std::sort(evals.data(), evals.data() + H.getEigenValues().size());

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    if (calc_gf) {
        INFO("1-particle Green's functions calc");
        std::set<ParticleIndex> f;
        std::set<IndexCombination2> indices2;
        ParticleIndex d0 = IndexInfo.getIndex("A",0,down);
        ParticleIndex u0 = IndexInfo.getIndex("A",0,up);
        f.insert(u0);
        f.insert(d0);

        indices2.insert(IndexCombination2(d0,d0));

        FieldOperatorContainer Operators(IndexInfo, HS, S, H, f);
        Operators.prepareAll(HS);
        Operators.computeAll();
        GFContainer G(IndexInfo,S,H,rho,Operators);

        G.prepareAll(indices2);
        G.computeAll();

        if(calc_2pgf) {
            print_section("2-Particle Green's function calc");
            std::set<IndexCombination4> indices4;
            indices4.insert(IndexCombination4(u0,u0,u0,u0));
            indices4.insert(IndexCombination4(u0,d0,u0,d0));
            indices4.insert(IndexCombination4(d0,d0,d0,d0));
            TwoParticleGFContainer Chi4(IndexInfo,S,H,rho,Operators);
            /** A difference in energies with magnitude less than this value is treated as zero. */
            Chi4.ReduceResonanceTolerance = reduce_tol;
            /** Minimal magnitude of the coefficient of a term to take it into account. */
            Chi4.CoefficientTolerance = coeff_tol;
            /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
            Chi4.MultiTermCoefficientTolerance = 1e-6;
            Chi4.prepareAll(indices4);
            MPI_Barrier(MPI_COMM_WORLD);
            std::vector<std::tuple<ComplexType, ComplexType, ComplexType> > freqs;
            Chi4.computeAll(false, freqs, MPI_COMM_WORLD, true);

            std::vector<ComplexType> chi_uuuu_vals = {-2.342841271771e+01,
                                                       0.000000000000e+00,
                                                       6.932231165814e-03,
                                                       2.037522082872e-03,
                                                      -2.150424835716e-03,
                                                      -4.384848776411e-03,
                                                      -5.253420668000e-03,
                                                      -5.370700986029e-03,
                                                      -5.126175681822e-03,
                                                      -4.732777836189e-03};

            const TwoParticleGF& chi_uuuu = Chi4(IndexCombination4(u0,u0,u0,u0));
            const TwoParticleGF& chi_udud = Chi4(IndexCombination4(u0,d0,u0,d0));
            const TwoParticleGF& chi_dddd = Chi4(IndexCombination4(d0,d0,d0,d0));

            ComplexType Omega = I*2.*M_PI/beta;
            ComplexType omega = I*M_PI/beta;

            for (size_t v = 0; v < chi_uuuu_vals.size(); v++) {
                ComplexType w_p = I*(2.*v+1.)*M_PI/beta;
                ComplexType chi_uuuu_val =  chi_uuuu(omega+Omega, w_p, omega);
                ComplexType chi_dddd_val =  chi_dddd(omega+Omega, w_p, omega);

                std::cout << "(" << imag(omega+Omega) << "," << imag(w_p) << "," << imag(omega) << "): " << std::flush;
                bool success = is_equal(chi_uuuu_val, chi_uuuu_vals[v], 1e-6);
                std::cout << "uuuu: " << chi_uuuu_val << " == (exact_uuuu) " << chi_uuuu_vals[v] << " == " << std::boolalpha << success << std::endl;
                if (!success) {
                    MPI_Finalize();
                    return EXIT_FAILURE;
                }
            }
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
