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

#pragma clang diagnostic ignored "-Wc++11-extensions"
#pragma clang diagnostic ignored "-Wgnu"

#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/lexical_cast.hpp>


#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>

#include<cstdlib>
#include <fstream>

#include <pomerol.h>
#include "mpi_dispatcher/mpi_dispatcher.hpp"

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
double FMatsubara(int n, double beta){return M_PI/beta*(2.*n+1);}
double BMatsubara(int n, double beta){return M_PI/beta*(2.*n);}

template <typename T>
ComplexType chi_bfreq_f(T const& chi, double W, double w1, double w2) { 
    return chi(I*(W+w1), I*w2, I*w1); // this comes from Pomerol - see TwoParticleGF::operator()
};
int job_to_bfreq_index(int job, int wbmax) { return -wbmax + job+1; }

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator comm;

    print_section("Hubbard nxn");
    
    int wf_max, wb_max;
    RealType e0, U, beta, reduce_tol, coeff_tol;
    bool calc_gf, calc_2pgf;
    size_t L;

    double eta, hbw, step; // for evaluation of GF on real axis 

    std::vector<double> levels;
    std::vector<double> hoppings;

    try { // command line parser
        TCLAP::CmdLine cmd("Hubbard nxn diag", ' ', "");
        TCLAP::ValueArg<RealType> U_arg("U","U","Value of U",true,10.0,"RealType",cmd);
        TCLAP::ValueArg<RealType> beta_arg("b","beta","Inverse temperature",true,100.,"RealType");
        TCLAP::ValueArg<RealType> T_arg("T","T","Temperature",true,0.01,"RealType");
        cmd.xorAdd(beta_arg,T_arg);

        TCLAP::MultiArg<double> level_args("l", "level", "level on auxiliary site", false,"RealType", cmd );
        TCLAP::MultiArg<double> hopping_args("t", "hopping", "hopping to an auxiliary site", false,"RealType", cmd );

        TCLAP::ValueArg<size_t> wn_arg("","wf","Number of positive fermionic Matsubara Freqs",false,64,"int",cmd);
        TCLAP::ValueArg<size_t> wb_arg("","wb","Number of positive bosonic Matsubara Freqs",false,1,"int",cmd);
        TCLAP::SwitchArg gf_arg("","calcgf","Calculate Green's functions",cmd, false);
        TCLAP::SwitchArg twopgf_arg("","calc2pgf","Calculate 2-particle Green's functions",cmd, false);
        TCLAP::ValueArg<RealType> reduce_tol_arg("","reducetol","Energy resonance resolution in 2pgf",false,1e-5,"RealType",cmd);
        TCLAP::ValueArg<RealType> coeff_tol_arg("","coefftol","Total weight tolerance",false,1e-12,"RealType",cmd);
        TCLAP::ValueArg<RealType> e0_arg("e","e0","Energy level of the impurity",false,0.0,"RealType",cmd);

        TCLAP::ValueArg<RealType> eta_arg("","eta","Offset from the real axis for Green's function calculation",false,0.05,"RealType",cmd);
        TCLAP::ValueArg<RealType> hbw_arg("D","hbw","Half-bandwidth. Default = U",false,0.0,"RealType",cmd);
        TCLAP::ValueArg<RealType> step_arg("","step","Step on a real axis. Default : 0.01",false,0.01,"RealType",cmd);

        cmd.parse( argc, argv ); // parse arguments
 
        U = U_arg.getValue();
        e0 = (e0_arg.isSet()?e0_arg.getValue():-U/2.0);
        boost::tie(beta, calc_gf, calc_2pgf, reduce_tol, coeff_tol) = boost::make_tuple( beta_arg.getValue(), 
            gf_arg.getValue(), twopgf_arg.getValue(), reduce_tol_arg.getValue(), coeff_tol_arg.getValue());
        boost::tie(wf_max, wb_max) = boost::make_tuple(wn_arg.getValue(), wb_arg.getValue());
        boost::tie(eta, hbw, step) = boost::make_tuple(eta_arg.getValue(), (hbw_arg.isSet()?hbw_arg.getValue():2.*U), step_arg.getValue());
        calc_gf = calc_gf || calc_2pgf;

        levels = level_args.getValue(); 
        hoppings = hopping_args.getValue();

        if (levels.size() != hoppings.size()) throw (std::logic_error("number of levels != number of hoppings"));
        }
    catch (TCLAP::ArgException &e)  // catch parsing exceptions
        { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; exit(1);}
    catch (std::exception &e)  // catch standard exceptions
        { std::cerr << "error: " << e.what() << std::endl; exit(1);}


    L = levels.size();
    INFO("Diagonalization of 1+" << L << " sites");

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
    IndexInfo.prepare(false); // Create index space
    if (!rank) { print_section("Indices"); IndexInfo.printIndices(); };
    int index_size = IndexInfo.getIndexSize();

    print_section("Matrix element storage");
    IndexHamiltonian Storage(&Lat,IndexInfo); 
    Storage.prepare(); // Write down the Hamiltonian as a symbolic formula
    print_section("Terms");
    if (!rank) INFO(Storage);

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

    INFO("<N> = " << rho.getAverageOccupancy()); // get average total particle number
    INFO("<H> = " << rho.getAverageEnergy()); // get average energy
    ParticleIndex d0 = IndexInfo.getIndex("A",0,down); // find the indices of the impurity, i.e. spin up index
    ParticleIndex u0 = IndexInfo.getIndex("A",0,up);
    INFO("<N_{" << IndexInfo.getInfo(u0) << "}N_{"<< IndexInfo.getInfo(u0) << "}> = " << rho.getAverageDoubleOccupancy(u0,d0)); // get double occupancy

    for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); i++) {  
        INFO("<N_{" << IndexInfo.getInfo(i) << "[" << i <<"]}> = " << rho.getAverageOccupancy(i)); // get average total particle number
        }
        
    savetxt("N_T.dat",rho.getAverageOccupancy());

    // Green's function calculation starts here
    
    FieldOperatorContainer Operators(IndexInfo, S, H); // Create a container for c and c^+ in the eigenstate basis

    if (calc_gf) {
        INFO("1-particle Green's functions calc");
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
            std::set<IndexCombination4> indices4; // a set of four indices to evaluate the 2pgf
            // 2PGF = <T c c c^+ c^+>
            indices4.insert(IndexCombination4(u0,u0,u0,u0)); // register up-up-up-up component for evaluation
            indices4.insert(IndexCombination4(u0,d0,u0,d0)); // register up-down-up-down component for evaluation
            //indices4.insert(IndexCombination4(d0,d0,d0,d0));

            TwoParticleGFContainer Chi4(IndexInfo,S,H,rho,Operators);
            /* Some knobs to make calc faster - the larger the values of tolerances, the faster is calc, but rounding errors may show up. */
            /** A difference in energies with magnitude less than this value is treated as zero - resolution of energy resonances. */
            Chi4.ReduceResonanceTolerance = reduce_tol;
            /** Minimal magnitude of the coefficient of a term to take it into account - resolution of thermal weight. */
            Chi4.CoefficientTolerance = coeff_tol;
            /** Knob that controls the caching frequency. */
            Chi4.ReduceInvocationThreshold = 1e5;
            /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
            Chi4.MultiTermCoefficientTolerance = 1e-6;
            
            Chi4.prepareAll(indices4); // find all non-vanishing block connections inside 2pgf
            comm.barrier(); // MPI::BARRIER

            // ! The most important routine - actually calculate the 2PGF
            Chi4.computeAll(comm, true); 

            // dump 2PGF into files - loop through 2pgf components
            for (std::set<IndexCombination4>::const_iterator it = indices4.begin(); it != indices4.end(); ++it) { 
                IndexCombination4 ind = *it;
                if (!comm.rank()) std::cout << "Saving 2PGF " << ind << std::endl;
                std::string ind_str = boost::lexical_cast<std::string>(ind.Index1) + boost::lexical_cast<std::string>(ind.Index2) +boost::lexical_cast<std::string>(ind.Index3) +boost::lexical_cast<std::string>(ind.Index4);
                const TwoParticleGF &chi = Chi4(ind);

                // Save terms of two particle GF
                std::ofstream term_res_stream(("terms_res"+ind_str+".pom").c_str());
                std::ofstream term_nonres_stream(("terms_nonres"+ind_str+".pom").c_str());
                boost::archive::text_oarchive oa_res(term_res_stream);
                boost::archive::text_oarchive oa_nonres(term_nonres_stream);
                for(std::vector<TwoParticleGFPart*>::const_iterator iter = chi.parts.begin(); iter != chi.parts.end(); iter++) {
                    oa_nonres << ((*iter)->getNonResonantTerms());
                    oa_res << ((*iter)->getResonantTerms());
                    };

                // start output of vertex

                // wb_max -  number of non-negative bosonic freqs
                // wf_max -  number of positive fermionic freqs

                // give 2pgf in bosonic-fermionic-fermionic freq notation

                 { // dispatch and save two-particle GF data - MPI parallelization in bosonic freqs

                    // Master-slave scheme to distribute the bosonic frequencies on different processes
                    int ROOT = 0;
                    int ntasks = std::max(2*wb_max-1,0);

                    boost::scoped_ptr<pMPI::MPIMaster> disp;

                    if (comm.rank() == ROOT) { 
                        DEBUG("Master at " << comm.rank());
                        disp.reset(new pMPI::MPIMaster(comm,ntasks,true));
                    };
                    comm.barrier();
                    
                
                    for (pMPI::MPIWorker worker(comm,ROOT);!worker.is_finished();) {
                        if (rank == ROOT) disp->order(); 
                        worker.receive_order(); 
                        //DEBUG((worker.Status == WorkerTag::Pending));
                        if (worker.is_working()) { 
                            // this is what every process executes
                            pMPI::JobId p = worker.current_job();

                            double W = BMatsubara(job_to_bfreq_index(p, wb_max), beta); // get current bosonic frequency
                            std::cout << "["<<p+1<<"/" << ntasks << "] p" << comm.rank() << " Omega = " << W << std::endl;

                            std::ofstream chi_stream (("chi"+ind_str+"_W"+boost::lexical_cast<std::string>(W)+".dat").c_str());

                            // Most important part 2 - loop over fermionic frequencies. Consider parallelizing one of the loop
                            for (int w1_index = -wf_max; w1_index<wf_max; w1_index++) { // loop over first fermionic frequency 
                                double w1 = FMatsubara(w1_index, beta);        
                                for (int w2_index = -wf_max; w2_index<wf_max; w2_index++) { // loop over second fermionic
                                    double w2 = FMatsubara(w2_index, beta);        

                                    ComplexType val = chi_bfreq_f(chi, W, w1, w2);

                                    chi_stream << std::scientific << std::setprecision(12)  
                                               << w1 << " " << w2 << "   " << std::real(val) << " " << std::imag(val) << std::endl;
                                }
                                chi_stream << std::endl;
                            };
                            chi_stream.close();
                            worker.report_job_done(); 
                            };
                        if (rank == ROOT) disp->check_workers(); // check if there are free workers 
                        };
                    comm.barrier();
                    //if (comm.rank() == ROOT) { disp.release(); DEBUG("Released master"); };
                    }
                }
            };
        };
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



