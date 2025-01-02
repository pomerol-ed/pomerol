//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file prog/hubbard2d_model.hpp
/// \brief ED calculations for the two-dimensional Hubbard model on an MxM cluster.
/// \author Sergei Iskakov (sir.iskakoff@gmail.com)
/// \author Igor Krivenko

#ifndef POMEROL_PROG_HUBBARD2D_MODEL_HPP
#define POMEROL_PROG_HUBBARD2D_MODEL_HPP

#include "quantum_model.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace Pomerol;

/// A full-ED calculation for the two-dimensional Hubbard model on an \f$N_x \times N_y\f$ cluster.
class hubbard2d_model : public quantum_model {

public:
    /// Constructor.
    /// \param[in] argc The number of command line arguments.
    /// \param[in] argv The command line arguments.
    hubbard2d_model(int argc, char* argv[])
        : quantum_model(argc, argv, "Full-ED of the MxM Hubbard cluster"),
          // clang-format off
    args_options_hubbard2d{{args_parser, "U", "Hubbard constant U", {"U"}, 10.0},
                           {args_parser, "mu", "Chemical potential", {"mu"}, std::numeric_limits<double>::quiet_NaN()},
                           {args_parser, "t", "NN hopping constant t", {"t"}, 1.0},
                           {args_parser, "tp", "NNN hopping constant t'", {"tp"}, 0.0},
                           {args_parser, "x", "Size over x", {"x"}, 2},
                           {args_parser, "y", "Size over y", {"y"}, 2}} // clang-format on
    {
        args_options_hubbard2d.mu.HelpDefault("U/2");
        parse_args(argc, argv);

        size_x = args::get(args_options_hubbard2d.x);
        size_y = args::get(args_options_hubbard2d.y);

        init_hamiltonian();
    }

    /// Return the (spin down, spin up) pair of indices belonging to one of cluster's sites.
    /// \param[in] IndexInfo Classification of operator index tuples.
    virtual std::pair<ParticleIndex, ParticleIndex> get_node(IndexInfoType const& IndexInfo) override {
        ParticleIndex d0 = IndexInfo.getIndex("S0", 0, LatticePresets::down);
        ParticleIndex u0 = IndexInfo.getIndex("S0", 0, LatticePresets::up);
        return std::make_pair(d0, u0);
    };

    /// \brief Construct Hamiltonian of the two-dimensional Hubbard model on an \f$N_x \times N_y\f$ cluster.
    void init_hamiltonian() {
        int L = size_x * size_y;
        INFO("Diagonalization of " << L << "=" << size_x << "*" << size_y << " sites");

        /* Add sites */
        names.resize(L);
        for(std::size_t y = 0; y < size_y; y++) {
            for(std::size_t x = 0; x < size_x; x++) {
                auto i = SiteIndexF(x, y);
                names[i] = "S" + std::to_string(i);
            }
        }

        /* Add interaction on each site*/
        double U = args::get(args_options_hubbard2d.U);
        double mu = args::get(args_options_hubbard2d.mu);
        if(std::isnan(mu))
            mu = U / 2;

        for(std::size_t i = 0; i < L; i++)
            HExpr += LatticePresets::CoulombS(names[i], U, -mu);

        /* Add hopping */
        double t = args::get(args_options_hubbard2d.t);
        double tp = args::get(args_options_hubbard2d.tp);

        for(std::size_t y = 0; y < size_y; y++) {
            for(std::size_t x = 0; x < size_x; x++) {
                auto pos = SiteIndexF(x, y);
                auto pos_right = SiteIndexF((x + 1) % size_x, y); /*if (x == size_x - 1) pos_right = SiteIndexF(0,y); */
                auto pos_up = SiteIndexF(x, (y + 1) % size_y);
                auto pos_dia_right = SiteIndexF((x + 1) % size_x, (y + 1) % size_y);
                auto pos_dia_left = SiteIndexF((x + 1) % size_x, (y - 1) % size_y);
                std::cout << pos_dia_right << " " << pos_dia_left << '\n';
                if(size_x > 1)
                    HExpr += LatticePresets::Hopping(std::min(names[pos], names[pos_right]),
                                                     std::max(names[pos], names[pos_right]),
                                                     -t);
                if(size_y > 1)
                    HExpr += LatticePresets::Hopping(std::min(names[pos], names[pos_up]),
                                                     std::max(names[pos], names[pos_up]),
                                                     -t);
                if(std::abs(tp) > 1e-10 && size_x > 1 && size_y > 1) {
                    HExpr += LatticePresets::Hopping(std::min(names[pos], names[pos_dia_right]),
                                                     std::max(names[pos], names[pos_dia_right]),
                                                     tp);
                    HExpr += LatticePresets::Hopping(std::min(names[pos], names[pos_dia_left]),
                                                     std::max(names[pos], names[pos_dia_left]),
                                                     tp);
                }
            }
        }

        auto rank = pMPI::rank(comm);
        if(!rank)
            mpi_cout << "Hamiltonian:\n" << HExpr << '\n';
    }

    /// Prepare a set of indices to evaluate annihilation/creation operators and
    /// a set of pairs of indices to evaluate single-particle Green's functions.
    /// \param[in] d0 The spin-down single particle index belonging to one of cluster's sites.
    /// \param[in] u0 The spin-up single particle index belonging to one of cluster's sites.
    /// \param[out] indices2 The set of index pairs for Green's function computation.
    /// \param[out] f The set of indices for annihilation/creation operator computation.
    /// \param[in] IndexInfo Classification of operator index tuples.
    virtual void prepare_indices(ParticleIndex d0,
                                 ParticleIndex u0,
                                 std::set<IndexCombination2>& indices2,
                                 std::set<ParticleIndex>& f,
                                 IndexInfoType const& IndexInfo) override {
        for(std::size_t x = 0; x < size_x; x++) {
            ParticleIndex ind = IndexInfo.getIndex(names[SiteIndexF(x, 0)], 0, LatticePresets::down);
            f.insert(ind);
            indices2.insert(IndexCombination2(d0, ind));
        }
    }

private:
    /// Description of command line options.
    struct {
        /// Hubbard interaction constant.
        args::ValueFlag<double> U;
        /// Chemical potential.
        args::ValueFlag<double> mu;
        /// Nearest neighbor hopping constant.
        args::ValueFlag<double> t;
        /// Next-nearest neighbor hopping constant.
        args::ValueFlag<double> tp;
        /// Linear size of the cluster along the \f$x\f$-direction, \f$N_x\f$.
        args::ValueFlag<int> x;
        /// Linear size of the cluster along the \f$y\f$-direction, \f$N_y\f$.
        args::ValueFlag<int> y;
    } args_options_hubbard2d;

    /// Linear size of the cluster along the \f$x\f$-direction, \f$N_x\f$.
    int size_x;
    /// Linear size of the cluster along the \f$y\f$-direction, \f$N_y\f$.
    int size_y;
    /// Names of cluster's sites.
    std::vector<std::string> names;

    /// Translate coordinates of cluster's sites into a linear index.
    /// \param[in] x Site coordinate along the \f$x\f$-direction.
    /// \param[in] y Site coordinate along the \f$y\f$-direction.
    /// \return Flattened index.
    std::size_t SiteIndexF(std::size_t x, std::size_t y) { return y * size_x + x; }
};

#endif // #ifndef POMEROL_PROG_HUBBARD2D_MODEL_HPP
