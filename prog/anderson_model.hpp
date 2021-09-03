//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//
// Created by iskakoff on 05/12/16.
//

#ifndef POMEROL_PROG_ANDERSON_MODEL_H
#define POMEROL_PROG_ANDERSON_MODEL_H

#include "quantum_model.hpp"

#include <cstddef>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace Pomerol;

/**
 * @brief anderson_model class
 *
 * @author iskakoff
 */
class anderson_model : public quantum_model {

public:
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

    virtual std::pair<ParticleIndex, ParticleIndex> get_node(IndexInfoType const& IndexInfo) override {
        ParticleIndex d0 = IndexInfo.getIndex("A", 0, LatticePresets::down);
        ParticleIndex u0 = IndexInfo.getIndex("A", 0, LatticePresets::up);
        return std::make_pair(d0, u0);
    }

    virtual void prepare_indices(ParticleIndex d0,
                                 ParticleIndex u0,
                                 std::set<IndexCombination2>& indices2,
                                 std::set<ParticleIndex>& f,
                                 IndexInfoType const& IndexInfo) override {
        indices2.insert(IndexCombination2(d0, d0)); // evaluate only G_{\down \down}
    }

    virtual void init_hamiltonian() override {
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

    struct {
        args::ValueFlag<double> U;
        args::ValueFlag<double> ed;
        args::ValueFlag<std::vector<double>, VectorReader> levels;
        args::ValueFlag<std::vector<double>, VectorReader> hoppings;
    } args_options_anderson;

    std::size_t L;
    std::vector<double> levels;
    std::vector<double> hoppings;
};

#endif // #ifndef POMEROL_PROG_ANDERSON_MODEL_H
