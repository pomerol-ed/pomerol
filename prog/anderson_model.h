//
// Created by iskakoff on 05/12/16.
//

#ifndef POMEROL_ANDERSON_MODEL_H
#define POMEROL_ANDERSON_MODEL_H

#include "quantum_model.h"

/**
 * @brief anderson_model class
 *
 * @author iskakoff
 */
class anderson_model: public quantum_model {

public:
  anderson_model(int argc, char ** argv): quantum_model(argc, argv){
    p = cmdline_params(argc, argv);
    init_parameters();
    init_lattice();
  };

  virtual void init_parameters() {
    quantum_model::init_parameters();
    if (p.count("levels")) {
      levels = p["levels"].as<std::vector<double> >();
      hoppings = p["hoppings"].as<std::vector<double> >();
    }

    if (levels.size() != hoppings.size()) {MPI_Finalize(); throw (std::logic_error("number of levels != number of hoppings")); }

    L = levels.size();
    int rank = comm.rank();
    mpi_cout << "Diagonalization of 1+" << L << " sites" << std::endl;
  }


  virtual po::variables_map cmdline_params(int argc, char* argv[]) {
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
    define<size_t>(p, "2pgf.reduce_freq", 1e5, "How often to reduce terms in 2pgf");
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

  virtual void init_lattice() {
    /* Add sites */
    Lat.addSite(new Lattice::Site("A",1,2));
    LatticePresets::addCoulombS(&Lat, "A", U(), e0());

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
    if (!rank) Lat.printSites();

    if (!rank) {
      mpi_cout << "Terms with 2 operators" << std::endl;
      Lat.printTerms(2);

      mpi_cout << "Terms with 4 operators" << std::endl;
      Lat.printTerms(4);
    };
  }
  size_t L;
  std::vector<double> levels;
  std::vector<double> hoppings;
};


#endif //POMEROL_ANDERSON_MODEL_H
