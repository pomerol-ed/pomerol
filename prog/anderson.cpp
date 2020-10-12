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

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/mpi.hpp>

#include <string>
#include <iostream>
#include <algorithm>
#include <tuple>

#include <pomerol.h>

#undef DEBUG
#include <gftools.hpp>
using gftools::tools::is_float_equal;
using gftools::grid_object;
using gftools::fmatsubara_grid;
using gftools::bmatsubara_grid;
using gftools::real_grid;

namespace po = boost::program_options;

void print_section (const std::string& str, boost::mpi::communicator comm = boost::mpi::communicator());

//boost::mpi::environment env;
using namespace Pomerol;
#define mpi_cout if(!comm.rank()) std::cout

template <typename T>
po::options_description define(po::options_description& opts, std::string name, T def_val, std::string desc) {
    opts.add_options()(name.c_str(), po::value<T>()->default_value(T(def_val)), desc.c_str()); return opts; }

template <typename T>
po::options_description define_vec(po::options_description& opts, std::string name, T&& def_val, std::string desc) {
    opts.add_options()(name.c_str(), po::value<T>()->multitoken()->default_value(T(def_val),""), desc.c_str()); return opts; }

// cmdline parser
po::variables_map cmdline_params(int argc, char* argv[])
{
    po::options_description p("Full-ED of the Anderson model");
    //po::variables_map p(argc, (const char**)argv);
    define_vec <std::vector<double> > (p, "levels", { }, "energy levels of the bath sites");
    define_vec <std::vector<double> > (p, "hoppings", { }, "hopping to the bath sites");

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
    define_vec<std::vector<size_t> >(p, "2pgf.indices", std::vector<size_t>({0, 0, 0, 0 }), "2pgf index combination");

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

    print_section("Anderson model ED");

    // all the params of the model
    RealType e0, U, beta;
    bool calc_gf, calc_2pgf;
    size_t L;
    std::vector<double> levels;
    std::vector<double> hoppings;

    po::variables_map p = cmdline_params(argc, argv);

    U = p["U"].as<double>();
    e0 = p["ed"].as<double>();
    std::tie(beta, calc_gf, calc_2pgf) = std::make_tuple(
        p["beta"].as<double>(), p["calc_gf"].as<int>(), p["calc_2pgf"].as<int>());
    calc_gf = calc_gf || calc_2pgf;


    if (p.count("levels")) {
        levels = p["levels"].as<std::vector<double>>();
        hoppings = p["hoppings"].as<std::vector<double>>();
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

    int rank = comm.rank();
    mpi_cout << "Sites" << std::endl;
    if (!rank) Lat.printSites();

    if (!rank) {
        mpi_cout << "Terms with 2 operators" << std::endl;
        Lat.printTerms(2);

        mpi_cout << "Terms with 4 operators" << std::endl;
        Lat.printTerms(4);
        };

    IndexClassification IndexInfo(Lat.getSiteMap());
    // Create index space
    IndexInfo.prepare(false);
    if (!rank) { print_section("Indices"); IndexInfo.printIndices(); };
    int index_size = IndexInfo.getIndexSize();

    print_section("Matrix element storage");
    IndexHamiltonian Storage(&Lat,IndexInfo);
    // Write down the Hamiltonian as a symbolic formula
    Storage.prepare();
    print_section("Terms");
    mpi_cout << Storage << std::endl;

    Symmetrizer Symm(IndexInfo, Storage);
    // Find symmetries of the problem
    Symm.compute();

    // Introduce Fock space and classify states to blocks
    StatesClassification S(IndexInfo,Symm);
    S.compute();

    // Hamiltonian in the basis of Fock Space
    Hamiltonian H(IndexInfo, Storage, S);
    // enter the Hamiltonian matrices
    H.prepare();
    // compute eigenvalues and eigenvectors
    H.compute();

    // Save spectrum
    if (!rank) {
        gftools::grid_object<double, gftools::enum_grid> evals1(gftools::enum_grid(0, S.getNumberOfStates()));
        RealVectorType evals (H.getEigenValues());
        std::sort(evals.data(), evals.data() + H.getEigenValues().size());
        std::copy(evals.data(), evals.data() + S.getNumberOfStates(), evals1.data().data());
        evals1.savetxt("spectrum.dat");
    }
    //savetxt("spectrum.dat", evals); // dump eigenvalues

    DensityMatrix rho(S,H,beta); // create Density Matrix
    rho.prepare();
    rho.compute(); // evaluate thermal weights with respect to ground energy, i.e exp(-beta(e-e_0))/Z

    ParticleIndex d0 = IndexInfo.getIndex("A",0,down); // find the indices of the impurity, i.e. spin up index
    ParticleIndex u0 = IndexInfo.getIndex("A",0,up);

    if (!rank) {
        // get average total particle number
        mpi_cout << "<N> = " << rho.getAverageOccupancy() << std::endl;
        // get average energy
        mpi_cout << "<H> = " << rho.getAverageEnergy() << std::endl;
        // get double occupancy
        mpi_cout << "<N_{" << IndexInfo.getInfo(u0) << "}N_{"<< IndexInfo.getInfo(u0) << "}> = " << rho.getAverageDoubleOccupancy(u0,d0) << std::endl;
        // get average total particle number per index
        for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); i++) {
            std::cout << "<N_{" << IndexInfo.getInfo(i) << "[" << i <<"]}> = " << rho.getAverageOccupancy(i) << std::endl;
            }
        double n_av = rho.getAverageOccupancy();
        gftools::num_io<double>(n_av).savetxt("N_T.dat");
        }

    // Green's function calculation starts here

    FieldOperatorContainer Operators(IndexInfo, S, H); // Create a container for c and c^+ in the eigenstate basis

    if (calc_gf) {
        int ntau; double eta, step, hbw;
        std::tie(ntau, eta, step, hbw) = std::make_tuple(p["gf.ntau"].as<int>(), p["gf.eta"].as<double>(), p["gf.step"].as<double>(), p["gf.D"].as<double>());
        int wf_min, wf_max, wb_min, wb_max;
        std::tie(wf_min, wf_max, wb_min, wb_max) = std::make_tuple(p["wf_min"].as<int>(), p["wf_max"].as<int>(), p["wb_min"].as<int>(), p["wb_max"].as<int>());

        mpi_cout << "1-particle Green's functions calc" << std::endl;
        std::set<ParticleIndex> f; // a set of indices to evaluate c and c^+
        std::set<IndexCombination2> indices2; // a set of pairs of indices to evaluate Green's function

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
            const GreensFunction & GF = G(ind2);

            mpi_cout << "Saving imfreq G" << ind2 << " on "<< 4*wf_max << " Matsubara freqs. " << std::endl;
            grid_object<std::complex<double>, fmatsubara_grid> gf_imfreq (fmatsubara_grid(wf_min, wf_max*4, beta));
            std::string ind_str = boost::lexical_cast<std::string>(ind2.Index1)+ boost::lexical_cast<std::string>(ind2.Index2);
            for (auto p : gf_imfreq.grid().points()) { gf_imfreq[p] = GF(p.value()); }
            gf_imfreq.savetxt("gw_imfreq_"+ ind_str +".dat");

            real_grid freq_grid(-hbw, hbw, 2*hbw/step+1, true);
            grid_object<std::complex<double>, real_grid> gf_refreq(freq_grid);
            for (auto p : freq_grid.points()) {
                ComplexType val = GF(ComplexType(p.value()) + I*eta);
                gf_refreq[p] = val;
            };
            gf_refreq.savetxt("gw_refreq_"+ ind_str +".dat");
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

            std::vector<std::tuple<ComplexType, ComplexType, ComplexType> > freqs_2pgf;
            fmatsubara_grid fgrid(wf_min, wf_max, beta);
            bmatsubara_grid bgrid(wb_min, wb_max + 1, beta);
            freqs_2pgf.reserve(fgrid.size() * fgrid.size() * bgrid.size());
            for (auto W : bgrid.values()) {
                for (auto w3 : fgrid.values()) {
                    for (auto w2 : fgrid.values()) {
                        ComplexType w1 = W+w3;
                        freqs_2pgf.push_back(boost::make_tuple(w1,w2,w3));
                        }
                    }
                }
            mpi_cout << "2PGF : " << freqs_2pgf.size() << " freqs to evaluate" << std::endl;

            // ! The most important routine - actually calculate the 2PGF
            auto chi_freq_data = G4.compute(true, freqs_2pgf, comm);

            // dump 2PGF into files - loop through 2pgf components
            if (!comm.rank()) {
                mpi_cout << "Saving 2PGF " << index_comb << std::endl;
                grid_object<std::complex<double>, bmatsubara_grid, fmatsubara_grid, fmatsubara_grid> full_vertex(std::forward_as_tuple(bgrid, fgrid, fgrid));
                grid_object<std::complex<double>, fmatsubara_grid, fmatsubara_grid> full_vertex_1freq(std::forward_as_tuple(fgrid, fgrid));
                size_t w_ind = 0;
                for (auto W : bgrid.points()) {
                    for (auto w3 : fgrid.points()) {
                        for (auto w2 : fgrid.points()) {
                            std::complex<double> val = chi_freq_data[w_ind];
                            full_vertex[W][w3.index()][w2.index()] = val;
                            full_vertex_1freq[w3.index()][w2.index()] = val;
                            if (!is_float_equal(boost::get<0>(freqs_2pgf[w_ind]), W.value()+w3.value())) throw std::logic_error("2pgf freq mismatch");
                            ++w_ind;
                            }
                        }
                    std::string fv1_name = "chi"+ind_str+"_W"+boost::lexical_cast<std::string>(std::imag(W.value()))+".dat";
                    full_vertex_1freq.savetxt(fv1_name);
                    }
                }
            // with the help of alpscore grid_object could be dumped to hdf5 - see save_grid_object(alps::hdf5::archive&, grid_object const&, group)
            }
    }
}

void print_section (const std::string& str, boost::mpi::communicator comm)
{
    mpi_cout << std::string(str.size(),'=') << std::endl;
    mpi_cout << str << std::endl;
    mpi_cout << std::string(str.size(),'=') << std::endl;
}

