// Hubbard 2x2
// Antipov, 2013

#pragma clang diagnostic ignored "-Wc++11-extensions"
#pragma clang diagnostic ignored "-Wgnu"

#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#include <string>
#include <iostream>
#include <algorithm>

#include<cstdlib>
#include <fstream>

#include <pomerol.h>

using namespace Pomerol;

extern boost::mpi::environment env;
boost::mpi::communicator comm;

/* Auxiliary routines - implemented in the bottom. */
bool compare(ComplexType a, ComplexType b);
void print_section (const std::string& str);
template <typename F1, typename F2>
    bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7);
template <typename T1> void savetxt(std::string fname, T1 in);
struct my_logic_error;

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator comm;

    Lattice Lat;
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

    /* Add sites */
    Lat.addSite(new Lattice::Site("A",1,2));
    LatticePresets::addCoulombS(&Lat, "A", U, -mu);

    std::vector<std::string> names(L);
    for (size_t i=0; i<L; i++)
        {
            std::stringstream s; s << i; 
            names[i]="b"+s.str();
            Lat.addSite(new Lattice::Site(names[i],1,2));
            LatticePresets::addHopping(&Lat, "A", names[i], hoppings[i]);
            LatticePresets::addLevel(&Lat, names[i], levels[i]);
        };
    
    INFO("Sites");
    Lat.printSites();

    int rank = comm.rank();
    if (!rank) {
        INFO("Terms with 2 operators");
        Lat.printTerms(2);

        INFO("Terms with 4 operators");
        Lat.printTerms(4);
        };

    IndexClassification IndexInfo(Lat.getSiteMap());
    IndexInfo.prepare(false);
    if (!rank) { print_section("Indices"); IndexInfo.printIndices(); };
    int index_size = IndexInfo.getIndexSize();

    print_section("Matrix element storage");
    IndexHamiltonian Storage(&Lat,IndexInfo);
    Storage.prepare();
    print_section("Terms");
    if (!rank) INFO(Storage);

    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute();

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    Hamiltonian H(IndexInfo, Storage, S);
    H.prepare(comm);
    H.compute(comm);

    RealVectorType evals (H.getEigenValues());
    std::sort(evals.data(), evals.data() + H.getEigenValues().size());

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    INFO("<N> = " << rho.getAverageOccupancy());
    
    FieldOperatorContainer Operators(IndexInfo, S, H);

    if (calc_gf) {
        INFO("1-particle Green's functions calc");
        std::set<ParticleIndex> f; 
        std::set<IndexCombination2> indices2;
        ParticleIndex d0 = IndexInfo.getIndex("A",0,down); 
        ParticleIndex u0 = IndexInfo.getIndex("A",0,up);
        f.insert(u0);
        f.insert(d0);
        
        indices2.insert(IndexCombination2(d0,d0));
        Operators.prepareAll(f);
        Operators.computeAll();
        GFContainer G(IndexInfo,S,H,rho,Operators);

        G.prepareAll(indices2);
        G.computeAll();

        if (calc_2pgf) {   
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
            /** Knob that controls the caching frequency. */
            Chi4.ReduceInvocationThreshold = 1e5;
            /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
            Chi4.MultiTermCoefficientTolerance = 1e-6;
            Chi4.prepareAll(indices4);
            comm.barrier();
            std::vector<boost::tuple<ComplexType, ComplexType, ComplexType> > freqs;
            Chi4.computeAll(false, freqs, comm, true);


            std::vector<ComplexType> chi_uuuu_vals(10);
            chi_uuuu_vals[0] = -2.342841271771e+01;
            chi_uuuu_vals[1] = 0.000000000000e+00;
            chi_uuuu_vals[2] = 6.932231165814e-03;
            chi_uuuu_vals[3] = 2.037522082872e-03;
            chi_uuuu_vals[4] = -2.150424835716e-03;
            chi_uuuu_vals[5] = -4.384848776411e-03;
            chi_uuuu_vals[6] = -5.253420668000e-03;
            chi_uuuu_vals[7] = -5.370700986029e-03;
            chi_uuuu_vals[8] = -5.126175681822e-03;
            chi_uuuu_vals[9] = -4.732777836189e-03;

            const TwoParticleGF& chi_uuuu = Chi4(IndexCombination4(u0,u0,u0,u0));
            const TwoParticleGF& chi_udud = Chi4(IndexCombination4(u0,u0,u0,u0));
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
                if (!success) return EXIT_FAILURE;
                
                };
            };
    };
}


void print_section (const std::string& str)
{
    if (!comm.rank()) { 
        std::cout << std::string(str.size(),'=') << std::endl;
        std::cout << str << std::endl;
        std::cout << std::string(str.size(),'=') << std::endl;
        };
}

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance)
{
    return (std::abs(x-y)<tolerance);
}

template <typename T1> void savetxt(std::string fname, T1 in){std::ofstream out(fname.c_str()); out << in << std::endl; out.close();};

struct my_logic_error : public std::logic_error { my_logic_error (const std::string& what_arg):logic_error(what_arg){}; };



