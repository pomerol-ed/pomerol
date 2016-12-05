//
// Created by iskakoff on 05/12/16.
//

#ifndef POMEROL_HUBBARD2D_MODEL_H
#define POMEROL_HUBBARD2D_MODEL_H


#include "quantum_model.h"

/**
 * @brief hubbard2d_model class
 *
 * @author iskakoff
 */
class hubbard2d_model: public quantum_model {

public:
  hubbard2d_model(int argc, char ** argv): quantum_model(argc, argv){
    p = cmdline_params(argc, argv);
    init_parameters();
    init_lattice();
  };

  virtual void init_parameters() {
    quantum_model::init_parameters();
    _mu = p.count("mu") ? p["mu"].as<double>() :U()/2;
    _t = p["t"].as<double>();
    boost::tie(size_x, size_y) = boost::make_tuple(p["x"].as<int>(), p["y"].as<int>());
  }

  virtual void init_lattice() {
    int L = size_x*size_y;
    INFO("Diagonalization of " << L << "=" << size_x << "*" << size_y << " sites");
    /* Add sites */
    std::vector<std::string> names(L);
    //auto SiteIndexF = [size_x](size_t x, size_t y){return y*size_x + x;};
    for (size_t y=0; y<size_y; y++) {
      for (size_t x = 0; x < size_x; x++) {
        auto i = SiteIndexF(x, y);
        std::stringstream s;
        s << i;
        names[i] = "S" + s.str();
        Lat.addSite(new Lattice::Site(names[i], 1, 2));
      };
    }

    INFO("Sites");
    Lat.printSites();

    /* Add interaction on each site*/
    for (size_t i=0; i<L; i++) LatticePresets::addCoulombS(&Lat, names[i], U(), -_mu);

    /* Add hopping */
    for (size_t y=0; y<size_y; y++) {
      for (size_t x=0; x<size_x; x++) {
        auto pos = SiteIndexF(x,y);
        auto pos_right = SiteIndexF((x+1)%size_x,y); /*if (x == size_x - 1) pos_right = SiteIndexF(0,y); */
        auto pos_up = SiteIndexF(x,(y+1)%size_y);
        mpi_cout<<pos<<"->"<<pos_right<<"\n|\n"<<pos_up<<std::endl;
        if (size_x > 1) LatticePresets::addHopping(&Lat, std::min(names[pos], names[pos_right]), std::max(names[pos], names[pos_right]), -_t);
        if (size_y > 1) LatticePresets::addHopping(&Lat, std::min(names[pos], names[pos_up]), std::max(names[pos], names[pos_up]), -_t);
      };
    };

    auto rank = comm.rank();
    if (!rank) {
      INFO("Terms with 2 operators");
      Lat.printTerms(2);
      INFO("Terms with 4 operators");
      Lat.printTerms(4);
    };
  }

  virtual po::variables_map cmdline_params(int argc, char**argv) {
    po::options_description p("Hubbard nxn diag");
    //po::variables_map p(argc, (const char**)argv);
    define_vec <std::vector<double> > (p, "levels", std::vector<double>(), "energy levels of the bath sites");
    define_vec <std::vector<double> > (p, "hoppings", std::vector<double>(), "hopping to the bath sites");

    define<double> (p, "U", 10.0, "Value of U");
    define<double> (p, "mu", 0.0, "Global chemical potential");
    define<double> (p, "t", 1.0, "Value of t");
    define<double> (p, "beta", 100, "Value of inverse temperature");
    define<double> (p, "T", .01, "Value of temperature");
    define<double> (p, "ed", 0.0, "Value of energy level of the impurity");
    define<int> (p,"x", 2,"Size over x");
    define<int> (p,"y", 2,"Size over y");

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

private:
  int size_x;
  int size_y;
  RealType _mu;
  RealType _t;
  size_t SiteIndexF(size_t x, size_t y) {
    return y*size_x + x;
  }
};


#endif //POMEROL_HUBBARD2D_MODEL_H
