//
// Created by iskakoff on 05/12/16.
//

#ifndef POMEROL_PROG_HUBBARD2D_MODEL_H
#define POMEROL_PROG_HUBBARD2D_MODEL_H

#include "quantum_model.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>
#include <utility>

using namespace Pomerol;

/**
 * @brief hubbard2d_model class
 *
 * @author iskakoff
 */
class hubbard2d_model: public quantum_model {

public:

  hubbard2d_model(int argc, char* argv[]) :
    quantum_model(argc, argv, "Full-ED of the MxM Hubbard cluster"),
    args_options_hubbard2d{{args_parser, "U", "Hubbard constant U", {"U"}, 10.0},
                           {args_parser, "mu", "Chemical potential", {"mu"}, std::numeric_limits<double>::quiet_NaN()},
                           {args_parser, "t", "NN hopping constant t", {"t"}, 1.0},
                           {args_parser, "tp", "NNN hopping constant t'", {"tp"}, 0.0},
                           {args_parser, "x", "Size over x", {"x"}, 2},
                           {args_parser, "y", "Size over y", {"y"}, 2}}
  {
    args_options_hubbard2d.mu.HelpDefault("U/2");
    parse_args(argc, argv);

    size_x = args::get(args_options_hubbard2d.x);
    size_y = args::get(args_options_hubbard2d.y);

    init_hamiltonian();
  }

  virtual std::pair<ParticleIndex, ParticleIndex>
  get_node(IndexInfoType const& IndexInfo) override {
    ParticleIndex d0 = IndexInfo.getIndex("S0",0,LatticePresets::down);
    ParticleIndex u0 = IndexInfo.getIndex("S0",0,LatticePresets::up);
    return std::make_pair(d0, u0);
  };

  virtual void init_hamiltonian() override {
    int L = size_x*size_y;
    INFO("Diagonalization of " << L << "=" << size_x << "*" << size_y << " sites");

    /* Add sites */
    names.resize(L);
    for (std::size_t y=0; y<size_y; y++) {
      for (std::size_t x = 0; x < size_x; x++) {
        auto i = SiteIndexF(x, y);
        names[i] = "S" + std::to_string(i);
      }
    }

    /* Add interaction on each site*/
    double U = args::get(args_options_hubbard2d.U);
    double mu = args::get(args_options_hubbard2d.mu);
    if(std::isnan(mu)) mu = U / 2;

    for (std::size_t i=0; i<L; i++)
      HExpr += LatticePresets::CoulombS(names[i], U, -mu);

    /* Add hopping */
    double t = args::get(args_options_hubbard2d.t);
    double tp = args::get(args_options_hubbard2d.tp);

    for (std::size_t y=0; y<size_y; y++) {
      for (std::size_t x=0; x<size_x; x++) {
        auto pos = SiteIndexF(x,y);
        auto pos_right = SiteIndexF((x+1)%size_x,y); /*if (x == size_x - 1) pos_right = SiteIndexF(0,y); */
        auto pos_up = SiteIndexF(x,(y+1)%size_y);
        auto pos_dia_right = SiteIndexF((x+1)%size_x,(y+1)%size_y);
        auto pos_dia_left  = SiteIndexF((x+1)%size_x,(y-1)%size_y);
        std::cout<<pos_dia_right<<" "<<pos_dia_left<<std::endl;
        if (size_x > 1)
          HExpr += LatticePresets::Hopping(std::min(names[pos], names[pos_right]), std::max(names[pos], names[pos_right]), -t);
        if (size_y > 1)
          HExpr += LatticePresets::Hopping(std::min(names[pos], names[pos_up]), std::max(names[pos], names[pos_up]), -t);
        if (std::abs(tp)>1e-10 && size_x > 1 && size_y>1) {
          HExpr += LatticePresets::Hopping(std::min(names[pos], names[pos_dia_right]), std::max(names[pos], names[pos_dia_right]), tp);
          HExpr += LatticePresets::Hopping(std::min(names[pos], names[pos_dia_left]), std::max(names[pos], names[pos_dia_left]), tp);
        }
      }
    }

    auto rank = pMPI::rank(comm);
    if (!rank) mpi_cout << "Hamiltonian:\n" << HExpr << std::endl;
  }

  virtual void prepare_indices(ParticleIndex d0,
                               ParticleIndex u0,
                               std::set<IndexCombination2> &indices2,
                               std::set<ParticleIndex>& f,
                               IndexInfoType const& IndexInfo) override {
    for (std::size_t x=0; x<size_x; x++) {
      ParticleIndex ind = IndexInfo.getIndex(names[SiteIndexF(x,0)],0,LatticePresets::down);
      f.insert(ind);
      indices2.insert(IndexCombination2(d0,ind));
    }
  }

private:

  struct {
    args::ValueFlag<double> U;
    args::ValueFlag<double> mu;
    args::ValueFlag<double> t;
    args::ValueFlag<double> tp;
    args::ValueFlag<int> x;
    args::ValueFlag<int> y;
  } args_options_hubbard2d;

  int size_x;
  int size_y;
  std::vector<std::string> names;

  std::size_t SiteIndexF(std::size_t x, std::size_t y) { return y * size_x + x; }
};

#endif // #ifndef POMEROL_PROG_HUBBARD2D_MODEL_H
