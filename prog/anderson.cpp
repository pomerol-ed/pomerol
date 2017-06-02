//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2014 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2014 Igor Krivenko <Igor.S.Krivenko@gmail.com>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


/** \file prog/anderson.cpp
** \brief Diagonalization of the Anderson impurity model (1 impurity coupled to a set of non-interacting bath sites)
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/lexical_cast.hpp>

#include <string>
#include <iostream>
#include <algorithm>

#include<cstdlib>
#include <fstream>

#include <pomerol.h>
#include "mpi_dispatcher/mpi_dispatcher.hpp"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace Pomerol;

extern boost::mpi::environment env;
boost::mpi::communicator comm;
#define mpi_cout if(!comm.rank()) std::cout

/* Auxiliary routines - implemented in the bottom. */
bool compare(ComplexType a, ComplexType b);
void print_section (const std::string& str);
template <typename F1, typename F2>
    bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7);
template <typename T1> void savetxt(std::string fname, T1 in);
struct my_logic_error;
double FMatsubara(int n, double beta){return M_PI/beta*(2.*n+1);}
double BMatsubara(int n, double beta){return M_PI/beta*(2.*n);}

template <typename T>
po::options_description define(po::options_description& opts, std::string name, T def_val, std::string desc) {
    opts.add_options()(name.c_str(), po::value<T>()->default_value(T(def_val)), desc.c_str()); return opts; }

template <typename T>
po::options_description define_vec(po::options_description& opts, std::string name, T def_val, std::string desc) {
    opts.add_options()(name.c_str(), po::value<T>()->multitoken()->default_value(T(def_val),""), desc.c_str()); return opts; }

// cmdline parser
po::variables_map cmdline_params(int argc, char* argv[])
{
    po::options_description p("Full-ED of the Anderson model");
    //po::variables_map p(argc, (const char**)argv);
    define_vec <std::vector<double> > (p, "levels", std::vector<double>(), "energy levels of the bath sites");
    define_vec <std::vector<double> > (p, "hoppings", std::vector<double>(), "hopping to the bath sites");

    define<double> (p, "U", 10.0, "Value of U");
    define<double> (p, "beta", 1, "Value of inverse temperature");
    define<double> (p, "ed", 0.0, "Value of energy level of the impurity");

    define<int>(p, "calc_gf", false, "Calculate Green's functions");
    define<int>(p, "calc_2pgf", false, "Calculate 2-particle Green's functions");
    define<int>(p, "wf_min", -20, "Minimum fermionic Matsubara freq");
    define<int>(p, "wf_max", 20, "Maximum fermionic Matsubara freq (4x for GF)");
    define<int>(p, "wb_min", 0, "Minimum bosonic Matsubara freq");
    define<int>(p, "wb_max", 0, "Maximum bosonic Matsubara freq");

    define<double>(p, "gf.eta", 0.05, "GF: Offset from the real axis for Green's function calculation");
    define<double>(p, "gf.step", 0.01, "GF: step of the real-freq grid");
    define<double>(p, "gf.D", 6, "GF: length of the real-freq grid");
    define<int>(p, "gf.ntau", 100, "GF: amount of points on the imag-time grid");

    define<double>(p, "2pgf.reduce_tol", 1e-5, "Energy resonance resolution in 2pgf");
    define<double>(p, "2pgf.coeff_tol",  1e-12, "Tolerance on nominators in 2pgf");
    define<double>(p, "2pgf.multiterm_tol", 1e-6, "How often to reduce terms in 2pgf");

    std::vector<size_t> default_inds(4,0);
    define_vec<std::vector<size_t> >(p, "2pgf.indices", default_inds, "2pgf index combination");

    p.add_options()("help","help");

    po::variables_map vm;
    //po::store(po::parse_command_line(argc, argv, p), vm);
    po::store(po::command_line_parser(argc, argv).options(p).style(
        po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);

    po::notify(vm);

    if (vm.count("help")) { std::cerr << p << "\n"; MPI_Finalize(); exit(0); }

    return vm;
}


int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator comm;

    // all the params of the model
    RealType e0, U, beta;
    bool calc_gf, calc_2pgf;
    size_t L;
    std::vector<double> levels;
    std::vector<double> hoppings;

    po::variables_map p = cmdline_params(argc, argv);

    U = p["U"].as<double>();
    e0 = p["ed"].as<double>();
    boost::tie(beta, calc_gf, calc_2pgf) = boost::make_tuple(
        p["beta"].as<double>(), p["calc_gf"].as<int>(), p["calc_2pgf"].as<int>());
    calc_gf = calc_gf || calc_2pgf;


    if (p.count("levels")) {
        levels = p["levels"].as<std::vector<double> >();
        hoppings = p["hoppings"].as<std::vector<double> >();
    }

    if (levels.size() != hoppings.size()) {MPI_Finalize(); throw (std::logic_error("number of levels != number of hoppings")); }

    L = levels.size();
    mpi_cout << "Diagonalization of 1+" << L << " sites" << std::endl;

    /* Add sites */
    Lattice Lat;
    Lat.addSite(new Lattice::Site("A",1,2));
    LatticePresets::addCoulombS(&Lat, "A", U, e0);

    std::vector<std::string> names(L);
    for (size_t i=0; i<L; i++)
        {
            std::stringstream s; s << i;
            names[i]="b"+s.str();
            Lat.addSite(new Lattice::Site(names[i],1,2));
            LatticePresets::addHopping(&Lat, "A", names[i], hoppings[i]);
            LatticePresets::addLevel(&Lat, names[i], levels[i]);
        };

    mpi_cout << "Sites" << std::endl;
    if (!comm.rank()) Lat.printSites();

    int rank = comm.rank();
    if (!rank) {
        mpi_cout << "Terms with 2 operators" << std::endl;
        Lat.printTerms(2);

        mpi_cout << "Terms with 4 operators" << std::endl;
        Lat.printTerms(4);
        };

    IndexClassification IndexInfo(Lat.getSiteMap());
    IndexInfo.prepare(false); // Create index space
    if (!rank) { print_section("Indices"); IndexInfo.printIndices(); };
    int index_size = IndexInfo.getIndexSize();

    print_section("Matrix element storage");
    IndexHamiltonian Storage(&Lat,IndexInfo);
    Storage.prepare(); // Write down the Hamiltonian as a symbolic formula
    print_section("Terms");
    mpi_cout << Storage << std::endl;

    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute(); // Find symmetries of the problem

    StatesClassification S(IndexInfo,Symm); // Introduce Fock space and classify states to blocks
    S.compute();

    Hamiltonian H(IndexInfo, Storage, S); // Hamiltonian in the basis of Fock Space
    H.prepare(); // enter the Hamiltonian matrices
    H.compute(); // compute eigenvalues and eigenvectors

    RealVectorType evals (H.getEigenValues());
    std::sort(evals.data(), evals.data() + H.getEigenValues().size());
    savetxt("spectrum.dat", evals); // dump eigenvalues

    DensityMatrix rho(S,H,beta); // create Density Matrix
    rho.prepare();
    rho.compute(); // evaluate thermal weights with respect to ground energy, i.e exp(-beta(e-e_0))/Z

    mpi_cout << "<N> = " << rho.getAverageOccupancy() << std::endl; // get average total particle number
    mpi_cout << "<H> = " << rho.getAverageEnergy() << std::endl; // get average energy
    ParticleIndex d0 = IndexInfo.getIndex("A",0,down); // find the indices of the impurity, i.e. spin up index
    ParticleIndex u0 = IndexInfo.getIndex("A",0,up);
    mpi_cout << "<N_{" << IndexInfo.getInfo(u0) << "}N_{"<< IndexInfo.getInfo(u0) << "}> = " << rho.getAverageDoubleOccupancy(u0,d0) << std::endl; // get double occupancy

    for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); i++) {
        mpi_cout << "<N_{" << IndexInfo.getInfo(i) << "[" << i <<"]}> = " << rho.getAverageOccupancy(i) << std::endl; // get average total particle number
        }

    if (!comm.rank()) savetxt("N_T.dat",rho.getAverageOccupancy());

    // Green's function calculation starts here

    FieldOperatorContainer Operators(IndexInfo, S, H); // Create a container for c and c^+ in the eigenstate basis

    if (calc_gf) {
        print_section("1-particle Green's functions calc");
        std::set<ParticleIndex> f; // a set of indices to evaluate c and c^+
        std::set<IndexCombination2> indices2; // a set of pairs of indices to evaluate Green's function

        int ntau; double eta, step, hbw;
        boost::tie(ntau, eta, step, hbw) = boost::make_tuple(p["gf.ntau"].as<int>(), p["gf.eta"].as<double>(), p["gf.step"].as<double>(), p["gf.D"].as<double>());
        int wf_min, wf_max, wb_min, wb_max;
        boost::tie(wf_min, wf_max, wb_min, wb_max) = boost::make_tuple(p["wf_min"].as<int>(), p["wf_max"].as<int>(), p["wb_min"].as<int>(), p["wb_max"].as<int>());

        // Take only impurity spin up and spin down indices
        f.insert(u0);
        f.insert(d0);
        indices2.insert(IndexCombination2(d0,d0)); // evaluate only G_{\down \down}

        Operators.prepareAll(f);
        Operators.computeAll(); // evaluate c, c^+ for chosen indices

        GFContainer G(IndexInfo,S,H,rho,Operators);

        G.prepareAll(indices2); // identify all non-vanishing block connections in the Green's function
        G.computeAll(); // Evaluate all GF terms, i.e. resonances and weights of expressions in Lehmans representation of the Green's function

        if (!comm.rank()) // dump gf into a file
        // loops over all components (pairs of indices) of the Green's function
        for (std::set<IndexCombination2>::const_iterator it = indices2.begin(); it != indices2.end(); ++it) {
            IndexCombination2 ind2 = *it;
            // Save Matsubara GF from pi/beta to pi/beta*(4*wf_max + 1)
            std::cout << "Saving imfreq G" << ind2 << " on "<< 4*wf_max << " Matsubara freqs. " << std::endl;
            std::ofstream gw_im(("gw_imag"+ boost::lexical_cast<std::string>(ind2.Index1)+ boost::lexical_cast<std::string>(ind2.Index2)+".dat").c_str());
            const GreensFunction & GF = G(ind2);
            for (int wn = 0; wn < wf_max*4; wn++) {
                ComplexType val = GF(I*FMatsubara(wn,beta)); // this comes from Pomerol - see GreensFunction::operator()
                gw_im << std::scientific << std::setprecision(12) << FMatsubara(wn,beta) << "   " << real(val) << " " << imag(val) << std::endl;
            };
            gw_im.close();
            // Save Retarded GF on the real axis
            std::ofstream gw_re(("gw_real"+boost::lexical_cast<std::string>(ind2.Index1)+boost::lexical_cast<std::string>(ind2.Index2)+".dat").c_str());
            std::cout << "Saving real-freq GF " << ind2 << " in energy space [" << e0-hbw << ":" << e0+hbw << ":" << step << "] + I*" << eta << "." << std::endl;
            for (double w = e0-hbw; w < e0+hbw; w+=step) {
                ComplexType val = GF(ComplexType(w) + I*eta);
                gw_re << std::scientific << std::setprecision(12) << w << "   " << real(val) << " " << imag(val) << std::endl;
            };
            gw_re.close();
        }

        // Start Two-particle GF calculation

        if (calc_2pgf) {
            print_section("2-Particle Green's function calc");

            std::vector<size_t> indices_2pgf = p["2pgf.indices"].as<std::vector<size_t> >();
            if (indices_2pgf.size() != 4) throw std::logic_error("Need 4 indices for 2pgf");

            // a set of four indices to evaluate the 2pgf
            IndexCombination4 index_comb(indices_2pgf[0], indices_2pgf[1], indices_2pgf[2], indices_2pgf[3]);

            std::set<IndexCombination4> indices4;
            // 2PGF = <T c c c^+ c^+>
            indices4.insert(index_comb);
            std::string ind_str = boost::lexical_cast<std::string>(index_comb.Index1)
                                + boost::lexical_cast<std::string>(index_comb.Index2)
                                + boost::lexical_cast<std::string>(index_comb.Index3)
                                + boost::lexical_cast<std::string>(index_comb.Index4);

            AnnihilationOperator const& C1 = Operators.getAnnihilationOperator(index_comb.Index1);
            AnnihilationOperator const& C2 = Operators.getAnnihilationOperator(index_comb.Index2);
            CreationOperator const&    CX3 = Operators.getCreationOperator(index_comb.Index3);
            CreationOperator const&    CX4 = Operators.getCreationOperator(index_comb.Index4);
            TwoParticleGF G4(S, H, C1, C2, CX3, CX4, rho);

            /* Some knobs to make calc faster - the larger the values of tolerances, the faster is calc, but rounding errors may show up. */
            /** A difference in energies with magnitude less than this value is treated as zero - resolution of energy resonances. */
            G4.ReduceResonanceTolerance = p["2pgf.reduce_tol"].as<double>();
            /** Minimal magnitude of the coefficient of a term to take it into account - resolution of thermal weight. */
            G4.CoefficientTolerance = p["2pgf.coeff_tol"].as<double>();
            /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
            G4.MultiTermCoefficientTolerance = p["2pgf.multiterm_tol"].as<double>();

            G4.prepare();
            comm.barrier(); // MPI::BARRIER

            // Fill a vector of tuples of fermionic Matsubara frequencies - these will be evaluated in-place
            std::vector<boost::tuple<ComplexType, ComplexType, ComplexType> > freqs_2pgf;
            for (int W_index = -wb_min; W_index <= wb_max; W_index++) { // loop over bosonic freq
                ComplexType W = I*BMatsubara(W_index, beta);
                for (int w3_index = -wf_max; w3_index<wf_max; w3_index++) { // loop over first fermionic frequency
                    ComplexType w3 = I*FMatsubara(w3_index, beta);
                    for (int w2_index = -wf_max; w2_index<wf_max; w2_index++) { // loop over second fermionic
                        ComplexType w2 = I*FMatsubara(w2_index, beta);
                        ComplexType w1 = W+w3;
                        freqs_2pgf.push_back(boost::make_tuple(w1,w2,w3));
                    }
                }
            }
            mpi_cout << "2PGF : " << freqs_2pgf.size() << " freqs to evaluate" << std::endl;

            std::vector<ComplexType> chi_freq_data = G4.compute(true, freqs_2pgf, comm); // mdata[ind];

            // Save terms of two particle GF
            std::ofstream term_res_stream(("terms_res"+ind_str+".pom").c_str());
            std::ofstream term_nonres_stream(("terms_nonres"+ind_str+".pom").c_str());
            boost::archive::text_oarchive oa_res(term_res_stream);
            boost::archive::text_oarchive oa_nonres(term_nonres_stream);
            for(std::vector<TwoParticleGFPart*>::const_iterator iter = G4.parts.begin(); iter != G4.parts.end(); iter++) {
                oa_nonres << ((*iter)->getNonResonantTerms());
                oa_res << ((*iter)->getResonantTerms());
                };

            if (!rank) {
                size_t w = 0;
                for (int W_index = -wb_min; W_index <= wb_max; W_index++) { // loop over bosonic freq
                    ComplexType W = I*BMatsubara(W_index, beta);
                    std::ofstream chi_stream (("chi"+ind_str+"_W"+boost::lexical_cast<std::string>(std::imag(W))+".dat").c_str());
                    for (int w3_index = -wf_max; w3_index<wf_max; w3_index++) { // loop over first fermionic frequency
                        ComplexType w3 = I*FMatsubara(w3_index, beta);
                        for (int w2_index = -wf_max; w2_index<wf_max; w2_index++) { // loop over second fermionic
                            ComplexType w2 = I*FMatsubara(w2_index, beta);
                            ComplexType w1 = W+w3;
                            if (std::abs(w1 - boost::get<0>(freqs_2pgf[w])) > 1e-8) throw std::logic_error("2pgf freq mismatch");

                            std::complex<double> val = chi_freq_data[w];

                            chi_stream << std::scientific << std::setprecision(12)
                               << w3.real() << " " << w3.imag() << " " << w2.real() << " " << w2.imag() << "   " << std::real(val) << " " << std::imag(val) << std::endl;
                            ++w;
                    }
                }
                chi_stream << std::endl;
                chi_stream.close();
                }
            }
        }
    }
}


bool compare(ComplexType a, ComplexType b)
{
    return abs(a-b) < 1e-5;
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



