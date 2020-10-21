//
// Created by iskakoff on 05/12/16.
//

#ifndef POMEROL_QUANTUM_MODEL_H
#define POMEROL_QUANTUM_MODEL_H

#include <iostream>
#include <string>
#include <algorithm>
#include <tuple>

#include <cstdlib>
#include <fstream>

#include <pomerol.h>
#include "mpi_dispatcher/mpi_dispatcher.hpp"
#include <boost/program_options.hpp>

#include <set>

#define mpi_cout if(!pMPI::rank(comm)) std::cout

namespace po = boost::program_options;

using namespace Pomerol;

#undef DEBUG
#include <gftools.hpp>
using gftools::tools::is_float_equal;
using gftools::grid_object;
using gftools::fmatsubara_grid;
using gftools::bmatsubara_grid;
using gftools::real_grid;

/**
 * @brief Base class for quantum model full-ED calculation
 *
 * @author iskakoff
 */
class quantum_model {
private:
  struct my_logic_error : public std::logic_error { my_logic_error (const std::string& what_arg):logic_error(what_arg){}; };

public:

  quantum_model(int argc, char ** argv) : comm(MPI_COMM_WORLD) {
      MPI_Init(&argc, &argv);
  }

  ~quantum_model() {
      MPI_Finalize();
  }

  virtual void init_parameters() {
    _U = p["U"].as<double>();
    _e0 = p["ed"].as<double>();
    std::tie(beta, calc_gf, calc_2pgf) = std::make_tuple(p["beta"].as<double>(), p["calc_gf"].as<int>(), p["calc_2pgf"].as<int>());
    calc_gf = calc_gf || calc_2pgf;
    rank = pMPI::rank(comm);
  }

  virtual void init_lattice() = 0;

  void compute();

  virtual std::pair<ParticleIndex, ParticleIndex> get_node(const IndexClassification &IndexInfo) = 0;

  double FMatsubara(int n, double beta){return M_PI/beta*(2.*n+1);}
  double BMatsubara(int n, double beta){return M_PI/beta*(2.*n);}

  template <typename T>
  po::options_description define(po::options_description& opts, std::string name, T def_val, std::string desc) {
    opts.add_options()(name.c_str(), po::value<T>()->default_value(T(def_val)), desc.c_str()); return opts; }

  template <typename T>
  po::options_description define_vec(po::options_description& opts, std::string name, T def_val, std::string desc) {
    opts.add_options()(name.c_str(), po::value<T>()->multitoken()->default_value(T(def_val),""), desc.c_str()); return opts; }

// cmdline parser
  virtual po::variables_map cmdline_params(int argc, char* argv[]) = 0;

  RealType U() {return _U;};
  RealType e0() {return _e0;};

  virtual void prepare_indices(ParticleIndex d0, ParticleIndex u0, std::set < IndexCombination2 > &indices2, std::set<ParticleIndex>& f, const IndexClassification &IndexInfo) = 0;
private:
  // all the params of the model
  RealType _e0;
  RealType _U;
  RealType beta;
  bool calc_gf;
  bool calc_2pgf;

protected:
  MPI_Comm comm;
  int rank;
  po::variables_map p;
  Lattice Lat;


  bool compare(ComplexType a, ComplexType b)
  {
    return abs(a-b) < 1e-5;
  }

  void print_section (const std::string& str)
  {
    if (!rank) {
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
};


#endif //POMEROL_QUANTUM_MODEL_H
