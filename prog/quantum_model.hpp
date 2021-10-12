//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file prog/quantum_model.hpp
/// \brief Base class for ED calculations of finite quantum many-body models.
/// \author Sergei Iskakov (sir.iskakoff@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

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

/// Standard output stream that discards all output on all MPI ranks except for rank 0.
#define mpi_cout                                                                                                       \
    if(!pMPI::rank(comm))                                                                                              \
    std::cout

/// std::vector<T> value reader for args.hxx
struct VectorReader {
    /// Read a value of type std::vector<T> from a string.
    /// \tparam T Element type of the vector to be read.
    /// \param[in] value Input string.
    /// \param[out] destination Output vector.
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

/// Base class for a full-ED calculation on a finite quantum many-body model.
class quantum_model {

public:
    /// IndexClassification type for (site label, orbital index, spin projection) index tuples.
    using IndexInfoType = Pomerol::IndexClassification<std::string, unsigned short, Pomerol::LatticePresets::spin>;

    /// Constructor.
    /// \param[in] argc The number of command line arguments.
    /// \param[in] argv The command line arguments.
    /// \param[in] prog_desc Brief description of the ED calculation.
    quantum_model(int argc, char* argv[], std::string const& prog_desc);
    /// Destructor.
    ~quantum_model();

    /// Parse the command line arguments common to all ED calculations.
    /// \param[in] argc The number of command line arguments.
    /// \param[in] argv The command line arguments.
    void parse_args(int argc, char* argv[]);

    /// \brief Construct system's Hamiltonian.
    ///
    /// To be implemented by calculations for specific quantum models.
    virtual void init_hamiltonian() = 0;

    /// Diagonalize the model and compute Green's functions.
    void compute();

    /// Return a model-dependent pair of indices for the single-particle Green's function calculation.
    /// \param[in] IndexInfo Model-dependent classification of operator index tuples.
    virtual std::pair<Pomerol::ParticleIndex, Pomerol::ParticleIndex> get_node(IndexInfoType const& IndexInfo) = 0;

    /// Return fermionic Matsubara frequency \f$\omega_n = \pi(2n+1)/\beta\f$.
    /// \param[in] n Index of the Matsubara frequency \f$n\f$.
    /// \param[in] beta Inverse temperature \f$\beta\f$.
    inline double FMatsubara(int n, double beta) { return M_PI / beta * (2. * n + 1); }
    /// Return bosonic Matsubara frequency \f$\nu_n = \pi(2n)/\beta\f$.
    /// \param[in] n Index of the Matsubara frequency \f$n\f$.
    /// \param[in] beta Inverse temperature \f$\beta\f$.
    inline double BMatsubara(int n, double beta) { return M_PI / beta * (2. * n); }

    /// Prepare a set of indices to evaluate annihilation/creation operators and
    /// a set of pairs of indices to evaluate single-particle Green's functions.
    /// \param[in] d0 The spin-down single particle index.
    /// \param[in] u0 The spin-up single particle index.
    /// \param[out] indices2 The set of index pairs for Green's function computation.
    /// \param[out] f The set of indices for annihilation/creation operator computation.
    /// \param[in] IndexInfo Model-dependent classification of operator index tuples.
    virtual void prepare_indices(Pomerol::ParticleIndex d0,
                                 Pomerol::ParticleIndex u0,
                                 std::set<Pomerol::IndexCombination2>& indices2,
                                 std::set<Pomerol::ParticleIndex>& f,
                                 IndexInfoType const& IndexInfo) = 0;

private:
    /// Inverse temperature.
    Pomerol::RealType beta = {};
    /// Whether to compute the single-particle Matsubara Green's function.
    bool calc_gf = false;
    /// Whether to compute the two-particle Matsubara Green's function.
    bool calc_2pgf = false;

protected:
    /// Parser for command line arguments.
    args::ArgumentParser args_parser;
    /// Description of command line options.
    struct {
        args::HelpFlag help;             ///< Request a usage message.
        args::ValueFlag<double> beta;    ///< Inverse temperature.
        args::Flag calc_gf;              ///< Whether to compute the single-particle Matsubara Green's function.
        args::Flag calc_2pgf;            ///< Whether to compute the two-particle Matsubara Green's function.
        args::ValueFlag<double> gf_eta;  ///< GF: Offset from the real axis for Green's function calculation.
        args::ValueFlag<double> gf_step; ///< GF: step of the real frequency grid.
        args::ValueFlag<double> gf_D;    ///< GF: length of the real frequency grid.
        args::ValueFlag<int> wf_min;     ///< Minimum fermionic Matsubara frequency.
        args::ValueFlag<int> wf_max;     ///< Maximum fermionic Matsubara frequency.
        args::ValueFlag<int> wb_min;     ///< Minimum bosonic Matsubara frequency.
        args::ValueFlag<int> wb_max;     ///< Maximum bosonic Matsubara frequency.
        /// Index combination of the two-particle Green's function.
        args::ValueFlag<std::vector<std::size_t>, VectorReader> _2pgf_indices;
        /// 2PGF: Energy resonance resolution.
        args::ValueFlag<double> _2pgf_reduce_tol;
        /// 2PGF: Tolerance on numerators.
        args::ValueFlag<double> _2pgf_coeff_tol;
        /// 2PGF: How often to reduce terms in 2PGF.
        args::ValueFlag<double> _2pgf_multiterm_tol;
    } args_options;

    /// MPI communicator for the calculation.
    MPI_Comm comm;
    /// MPI rank of the process owning this object.
    int rank;

    /// Expression of system's Hamiltonian.
    Pomerol::LatticePresets::RealExpr HExpr;

    /// Print a line preceded and forwarded by two horizontal rules.
    void print_section(std::string const& str) {
        if(!rank) {
            std::cout << std::string(str.size(), '=') << std::endl;
            std::cout << str << std::endl;
            std::cout << std::string(str.size(), '=') << std::endl;
        }
    }
};

#endif // #ifndef POMEROL_PROG_QUANTUM_MODEL_HPP
