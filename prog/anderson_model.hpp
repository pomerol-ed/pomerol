//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file prog/anderson_model.hpp
/// \brief ED calculations for the Anderson impurity model.
/// \author Sergei Iskakov (sir.iskakoff@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_PROG_ANDERSON_MODEL_HPP
#define POMEROL_PROG_ANDERSON_MODEL_HPP

#include "quantum_model.hpp"

#include <cstddef>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace Pomerol;

/// A full-ED calculation on the Anderson impurity model.
class anderson_model : public quantum_model {

public:
    /// Constructor.
    /// \param[in] argc The number of command line arguments.
    /// \param[in] argv The command line arguments.
    anderson_model(int argc, char* argv[])
        : quantum_model(argc, argv, "Full-ED of the Anderson model"),
          // clang-format off
    args_options_anderson{{args_parser, "U", "Interaction constant U", {"U"}, 10.0},
                          {args_parser, "ed", "Energy level of the impurity", {"ed"}, 0},
                          {args_parser, "levels", "Energy levels of the bath sites", {"levels"}, {}},
                          {args_parser, "hoppings", "Hopping to the bath sites", {"hoppings"}, {}}}
    // clang-format on
    {
        parse_args(argc, argv);

        levels = args::get(args_options_anderson.levels);
        hoppings = args::get(args_options_anderson.hoppings);

        if(levels.size() != hoppings.size()) {
            MPI_Finalize();
            std::cerr << "Number of levels != number of hoppings" << std::endl;
            exit(2);
        }

        L = levels.size();
        mpi_cout << "Diagonalization of 1+" << L << " sites" << std::endl;

        init_hamiltonian();
    }

    /// Return the (spin down, spin up) pair of indices of the correlated atom.
    /// \param[in] IndexInfo Classification of operator index tuples.
    virtual std::pair<ParticleIndex, ParticleIndex> get_node(IndexInfoType const& IndexInfo) override {
        ParticleIndex d0 = IndexInfo.getIndex("A", 0, LatticePresets::down);
        ParticleIndex u0 = IndexInfo.getIndex("A", 0, LatticePresets::up);
        return std::make_pair(d0, u0);
    }

    /// Prepare a set of indices to evaluate annihilation/creation operators and
    /// a set of pairs of indices to evaluate single-particle Green's functions.
    /// \param[in] d0 The spin-down single particle index of the correlated atom.
    /// \param[in] u0 The spin-up single particle index of the correlated atom.
    /// \param[out] indices2 The set of index pairs for Green's function computation.
    /// \param[out] f The set of indices for annihilation/creation operator computation.
    /// \param[in] IndexInfo Classification of operator index tuples.
    virtual void prepare_indices(ParticleIndex d0,
                                 ParticleIndex u0,
                                 std::set<IndexCombination2>& indices2,
                                 std::set<ParticleIndex>& f,
                                 IndexInfoType const& IndexInfo) override {
        indices2.insert(IndexCombination2(d0, d0)); // evaluate only G_{\down \down}
    }

    /// \brief Construct Hamiltonian of the single impurity Anderson model.
    void init_hamiltonian() {
        HExpr += LatticePresets::CoulombS("A", args::get(args_options_anderson.U), args::get(args_options_anderson.ed));

        std::vector<std::string> names(L);
        for(std::size_t i = 0; i < L; i++) {
            names[i] = "b" + std::to_string(i);
            HExpr += LatticePresets::Hopping("A", names[i], hoppings[i]);
            HExpr += LatticePresets::Level(names[i], levels[i]);
        }

        if(!rank)
            mpi_cout << "Hamiltonian:\n" << HExpr << std::endl;
    }

private:
    /// Description of command line options.
    struct {
        /// Coulomb repulsion constant.
        args::ValueFlag<double> U;
        /// Energy level localized on the correlated atom.
        args::ValueFlag<double> ed;
        /// Bath levels.
        args::ValueFlag<std::vector<double>, VectorReader> levels;
        /// Hopping constants between the correlated atom and the bath sites.
        args::ValueFlag<std::vector<double>, VectorReader> hoppings;
    } args_options_anderson;

    /// The number of bath sites.
    std::size_t L;
    /// Bath levels.
    std::vector<double> levels;
    /// Hopping constants between the correlated atom and the bath sites.
    std::vector<double> hoppings;
};

#endif // #ifndef POMEROL_PROG_ANDERSON_MODEL_HPP
