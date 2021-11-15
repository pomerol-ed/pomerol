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

#ifndef POMEROL_PROG_QUANTUM_MODEL_HPP
#define POMEROL_PROG_QUANTUM_MODEL_HPP

#include <pomerol.hpp>

#include "args.hxx"

#include <cstddef>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#define mpi_cout                                                                                                       \
    if(!pMPI::rank(comm))                                                                                              \
    std::cout

// Vector value reader for args
struct VectorReader {
    template <typename T> void operator()(std::string const&, std::string const& value, std::vector<T>& destination) {
        std::stringstream ss(value);
        std::string token;
        T element;
        while(std::getline(ss, token, ',')) {
            std::stringstream(token) >> element;
            destination.emplace_back(element);
        }
    }
};

/**
 * @brief Base class for quantum model full-ED calculation
 *
 * @author iskakoff
 */
class quantum_model {

public:
    using IndexInfoType = Pomerol::IndexClassification<std::string, unsigned short, Pomerol::LatticePresets::spin>;

    quantum_model(int argc, char* argv[], std::string const& prog_desc);
    ~quantum_model();

    void parse_args(int argc, char* argv[]);

    virtual void init_hamiltonian() = 0;

    void compute();

    virtual std::pair<Pomerol::ParticleIndex, Pomerol::ParticleIndex> get_node(IndexInfoType const& IndexInfo) = 0;

    double FMatsubara(int n, double beta) { return M_PI / beta * (2. * n + 1); }
    double BMatsubara(int n, double beta) { return M_PI / beta * (2. * n); }

    virtual void prepare_indices(Pomerol::ParticleIndex d0,
                                 Pomerol::ParticleIndex u0,
                                 std::set<Pomerol::IndexCombination2>& indices2,
                                 std::set<Pomerol::ParticleIndex>& f,
                                 IndexInfoType const& IndexInfo) = 0;

private:
    // Simulation parameters
    Pomerol::RealType beta = {};
    bool calc_gf = false;
    bool calc_2pgf = false;

protected:
    // Command line parsing
    args::ArgumentParser args_parser;
    struct {
        args::HelpFlag help;
        args::ValueFlag<double> beta;
        args::Flag calc_gf;
        args::Flag calc_2pgf;
        args::ValueFlag<double> gf_eta;
        args::ValueFlag<double> gf_step;
        args::ValueFlag<double> gf_D;
        args::ValueFlag<int> wf_min;
        args::ValueFlag<int> wf_max;
        args::ValueFlag<int> wb_min;
        args::ValueFlag<int> wb_max;
        args::ValueFlag<std::vector<std::size_t>, VectorReader> _2pgf_indices;
        args::ValueFlag<double> _2pgf_reduce_tol;
        args::ValueFlag<double> _2pgf_coeff_tol;
        args::ValueFlag<double> _2pgf_multiterm_tol;
    } args_options;

    MPI_Comm comm;
    int rank;

    Pomerol::LatticePresets::RealExpr HExpr;

    void print_section(std::string const& str) {
        if(!rank) {
            std::cout << std::string(str.size(), '=') << std::endl;
            std::cout << str << std::endl;
            std::cout << std::string(str.size(), '=') << std::endl;
        }
    }
};

#endif // #ifndef POMEROL_PROG_QUANTUM_MODEL_HPP
