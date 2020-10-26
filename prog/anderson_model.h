//
// Created by iskakoff on 05/12/16.
//

#ifndef POMEROL_ANDERSON_MODEL_H
#define POMEROL_ANDERSON_MODEL_H

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

    init_lattice();
  }

  virtual std::pair<ParticleIndex, ParticleIndex>
  get_node(const IndexClassification & IndexInfo) override {
    ParticleIndex d0 = IndexInfo.getIndex("A",0,down);
    ParticleIndex u0 = IndexInfo.getIndex("A",0,up);
    return std::make_pair(d0, u0);
  }

  virtual void prepare_indices(ParticleIndex d0,
                               ParticleIndex u0,
                               std::set<IndexCombination2> &indices2,
                               std::set<ParticleIndex>& f,
                               const IndexClassification &IndexInfo) override {
    indices2.insert(IndexCombination2(d0, d0)); // evaluate only G_{\down \down}
  }

  virtual void init_lattice() override {
    /* Add sites */
    Lat.addSite(new Lattice::Site("A",1,2));
    LatticePresets::addCoulombS(&Lat,
                                "A",
                                args::get(args_options.U),
                                args::get(args_options.ed));

    std::vector<std::string> names(L);
    for (size_t i=0; i<L; i++) {
      names[i] = "b" + std::to_string(i);
      Lat.addSite(new Lattice::Site(names[i],1,2));
      LatticePresets::addHopping(&Lat, "A", names[i], hoppings[i]);
      LatticePresets::addLevel(&Lat, names[i], levels[i]);
    }

    mpi_cout << "Sites" << std::endl;
    if (!rank) Lat.printSites();

    if (!rank) {
      mpi_cout << "Terms with 2 operators" << std::endl;
      Lat.printTerms(2);

      mpi_cout << "Terms with 4 operators" << std::endl;
      Lat.printTerms(4);
    }
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

#endif //POMEROL_ANDERSON_MODEL_H
