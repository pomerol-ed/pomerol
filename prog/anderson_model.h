//
// Created by iskakoff on 05/12/16.
//

#ifndef POMEROL_PROG_ANDERSON_MODEL_H
#define POMEROL_PROG_ANDERSON_MODEL_H

#include "quantum_model.h"

using namespace Pomerol;

/**
 * @brief anderson_model class
 *
 * @author iskakoff
 */
class anderson_model: public quantum_model {

public:

  anderson_model(int argc, char* argv[]) :
    quantum_model(argc, argv, "Full-ED of the Anderson model"),
    args_options{{args_parser, "U", "Interaction constant U", {"U"}, 10.0},
                 {args_parser, "ed", "Energy level of the impurity", {"ed"}, 0},
                 {args_parser, "levels", "Energy levels of the bath sites", {"levels"}, {}},
                 {args_parser, "hoppings", "Hopping to the bath sites", {"hoppings"}, {}}}
  {
    parse_args(argc, argv);

    levels = args::get(args_options.levels);
    hoppings = args::get(args_options.hoppings);

    if(levels.size() != hoppings.size()) {
      MPI_Finalize();
      std::cerr << "Number of levels != number of hoppings" << std::endl;
      exit(2);
    }

    L = levels.size();
    mpi_cout << "Diagonalization of 1+" << L << " sites" << std::endl;

    init_hamiltonian();
  }

  virtual std::pair<ParticleIndex, ParticleIndex>
  get_node(const IndexInfoType & IndexInfo) override {
    ParticleIndex d0 = IndexInfo.getIndex("A",0,down);
    ParticleIndex u0 = IndexInfo.getIndex("A",0,up);
    return std::make_pair(d0, u0);
  }

  virtual void prepare_indices(ParticleIndex d0,
                               ParticleIndex u0,
                               std::set<IndexCombination2> &indices2,
                               std::set<ParticleIndex>& f,
                               const IndexInfoType &IndexInfo) override {
    indices2.insert(IndexCombination2(d0, d0)); // evaluate only G_{\down \down}
  }

  virtual void init_hamiltonian() override {
    HExpr += LatticePresets::CoulombS("A",
                                      args::get(args_options.U),
                                      args::get(args_options.ed));

    std::vector<std::string> names(L);
    for (size_t i=0; i<L; i++) {
      names[i] = "b" + std::to_string(i);
      HExpr += LatticePresets::Hopping("A", names[i], hoppings[i]);
      HExpr += LatticePresets::Level(names[i], levels[i]);
    }

    if (!rank) mpi_cout << "Hamiltonian:\n" << HExpr << std::endl;
  }

  struct {
    args::ValueFlag<double> U;
    args::ValueFlag<double> ed;
    args::ValueFlag<std::vector<double>, VectorReader> levels;
    args::ValueFlag<std::vector<double>, VectorReader> hoppings;
  } args_options;

  size_t L;
  std::vector<double> levels;
  std::vector<double> hoppings;
};

#endif // #ifndef POMEROL_PROG_ANDERSON_MODEL_H
