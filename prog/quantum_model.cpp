//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file prog/quantum_model.hpp
/// \brief Base class for ED calculations of finite quantum many-body models (implementation).
/// \author Sergei Iskakov (sir.iskakoff@gmail.com)
/// \author Igor Krivenko

#include "quantum_model.hpp"

#undef DEBUG
#include <gftools.hpp>

#include <algorithm>
#include <complex>
#include <stdexcept>
#include <tuple>

using namespace Pomerol;

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
quantum_model::quantum_model(int argc, char* argv[], std::string const& prog_desc)
    : args_parser(prog_desc),
      // clang-format off
  args_options{{args_parser, "help", "Display this help menu", {'h', "help"}},
               {args_parser, "beta", "Inverse temperature", {"beta"}, 1.0},
               {args_parser, "calc_gf", "Calculate Green's functions", {"calc_gf"}},
               {args_parser, "calc_2gf", "Calculate 2-particle Green's functions", {"calc_2gf"}},
               {args_parser, "eta", "GF: Offset from the real axis for Green's function calculation", {"gf.eta"}, 0.05},
               {args_parser, "step", "GF: step of the real frequency grid", {"gf.step"}, 0.01},
               {args_parser, "D", "GF: length of the real frequency grid", {"gf.D"}, 6.0},
               {args_parser, "wf_min", "Minimum fermionic Matsubara frequency", {"wf_min"}, -20},
               {args_parser, "wf_max", "Maximum fermionic Matsubara frequency (4x for GF)", {"wf_max"}, 20},
               {args_parser, "wb_min", "Minimum bosonic Matsubara frequency", {"wb_min"}, 0},
               {args_parser, "wb_max", "Maximum bosonic Matsubara frequency", {"wb_max"}, 0},
               {args_parser, "indices", "2PGF index combination", {"2pgf.indices"}, {0, 0, 0, 0}},
               {args_parser, "tol", "Energy resonance resolution in 2PGF", {"2pgf.reduce_tol"}, 1e-5},
               {args_parser, "tol", "Tolerance on numerators in 2PGF", {"2pgf.coeff_tol"}, 1e-12}
  },
      // clang-format on
      comm(MPI_COMM_WORLD) {
    MPI_Init(&argc, &argv);
    // NOLINTNEXTLINE(cppcoreguidelines-prefer-member-initializer)
    rank = pMPI::rank(comm);
    // Print defaults in the usage message
    args_parser.helpParams.addDefault = true;
}

quantum_model::~quantum_model() {
    MPI_Finalize();
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
void quantum_model::parse_args(int argc, char const* const argv[]) {
    try {
        args_parser.ParseCLI(argc, argv);
    } catch(args::Help) {
        std::cout << args_parser;
        exit(0);
    } catch(args::ParseError& e) {
        std::cerr << e.what() << '\n';
        std::cerr << args_parser;
        exit(1);
    } catch(args::ValidationError& e) {
        std::cerr << e.what() << '\n';
        std::cerr << args_parser;
        exit(1);
    }

    beta = args::get(args_options.beta);
    calc_gf = args::get(args_options.calc_gf);
    calc_2pgf = args::get(args_options.calc_2pgf);
    calc_gf = calc_gf || calc_2pgf;
}

void quantum_model::compute() {

    using gftools::tools::is_float_equal;
    using gftools::grid_object;
    using gftools::fmatsubara_grid;
    using gftools::bmatsubara_grid;
    using gftools::real_grid;

    IndexInfoType IndexInfo = MakeIndexClassification(HExpr);
    if(!rank) {
        print_section("Indices");
        std::cout << IndexInfo << '\n';
    };

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();

    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    if(!rank) {
        gftools::grid_object<double, gftools::enum_grid> evals1(
            gftools::enum_grid(0, static_cast<int>(S.getNumberOfStates())));
        RealVectorType evals(H.getEigenValues());
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        std::sort(evals.data(), evals.data() + H.getEigenValues().size());
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        std::copy(evals.data(), evals.data() + S.getNumberOfStates(), evals1.data().data());
        evals1.savetxt("spectrum.dat");
    }
    DensityMatrix rho(S, H, beta); // Create density matrix.
    rho.prepare();
    rho.compute(); // Evaluate thermal weights with respect to ground energy, i.e exp(-beta(e-e_0))/Z.

    std::pair<ParticleIndex, ParticleIndex> pair = get_node(IndexInfo);
    // cppcheck-suppress variableScope
    ParticleIndex d0 = pair.first;
    // cppcheck-suppress variableScope
    ParticleIndex u0 = pair.second;

    // Green's function calculation starts here.

    if(calc_gf) {
        print_section("1-particle Green's functions calc");
        std::set<ParticleIndex> f;            // A set of indices to evaluate c and c^+.
        std::set<IndexCombination2> indices2; // A set of pairs of indices to evaluate Green's function.

        double eta = args::get(args_options.gf_eta);
        double step = args::get(args_options.gf_step);
        double hbw = args::get(args_options.gf_D);

        int wf_min = args::get(args_options.wf_min);
        int wf_max = args::get(args_options.wf_max);
        int wb_min = args::get(args_options.wb_min);
        int wb_max = args::get(args_options.wb_max);

        // Take only impurity spin up and spin down indices.
        f.insert(u0);
        f.insert(d0);
        prepare_indices(d0, u0, indices2, f, IndexInfo);

        // Create a container for c and c^+ in the eigenstate basis.
        FieldOperatorContainer Operators(IndexInfo, HS, S, H, f);
        Operators.prepareAll(HS);
        Operators.computeAll(); // evaluate c, c^+ for chosen indices.

        GFContainer G(IndexInfo, S, H, rho, Operators);

        // Identify all non-vanishing block connections in the Green's function.
        G.prepareAll(indices2);
        // Evaluate all GF terms, i.e. resonances and weights of expressions in Lehmann's representation
        // of the Green's function.
        G.computeAll();

        if(!rank) // Dump gf into a file
            // Loops over all components (pairs of indices) of the Green's function.
            for(auto const& ind2 : indices2) {
                GreensFunction const& GF = G(ind2);
                // Save Matsubara GF from pi/beta to pi/beta*(4*wf_max + 1)
                std::cout << "Saving imfreq G" << ind2 << " on " << 4 * wf_max << " Matsubara freqs.\n";
                grid_object<std::complex<double>, fmatsubara_grid> gf_imfreq(
                    fmatsubara_grid(wf_min, wf_max * 4, beta, true));
                std::string ind_str = std::to_string(ind2.Index1) + std::to_string(ind2.Index2);
                for(auto p : gf_imfreq.grid().points()) {
                    gf_imfreq[p] = GF(p.value());
                }
                gf_imfreq.savetxt("gw_imfreq_" + ind_str + ".dat");

                real_grid freq_grid(-hbw, hbw, 2 * static_cast<std::size_t>(hbw / step) + 1, true);
                grid_object<std::complex<double>, real_grid> gf_refreq(freq_grid);
                for(auto p : freq_grid.points()) {
                    ComplexType val = GF(ComplexType(p.value()) + I * eta);
                    gf_refreq[p] = val;
                };
                gf_refreq.savetxt("gw_refreq_" + ind_str + ".dat");
            }

        // Start Two-particle GF calculation.

        if(calc_2pgf) {
            print_section("2-Particle Green's function calc");

            std::vector<std::size_t> indices_2pgf = args::get(args_options._2pgf_indices);

            if(indices_2pgf.size() != 4)
                throw std::runtime_error("Need 4 indices for 2PGF");

            // a set of four indices to evaluate the 2pgf
            IndexCombination4 index_comb(indices_2pgf[0], indices_2pgf[1], indices_2pgf[2], indices_2pgf[3]);

            std::set<IndexCombination4> indices4;
            // 2PGF = <T c c c^+ c^+>
            indices4.insert(index_comb);
            std::string ind_str = std::to_string(index_comb.Index1) + std::to_string(index_comb.Index2) +
                                  std::to_string(index_comb.Index3) + std::to_string(index_comb.Index4);

            AnnihilationOperator const& C1 = Operators.getAnnihilationOperator(index_comb.Index1);
            AnnihilationOperator const& C2 = Operators.getAnnihilationOperator(index_comb.Index2);
            CreationOperator const& CX3 = Operators.getCreationOperator(index_comb.Index3);
            CreationOperator const& CX4 = Operators.getCreationOperator(index_comb.Index4);
            TwoParticleGF G4(S, H, C1, C2, CX3, CX4, rho);

            // Some knobs to make calc faster - the larger the values of tolerances, the faster is calc,
            // but rounding errors may show up.

            // A difference in energies with magnitude less than this value is treated as zero
            // - resolution of energy resonances.
            G4.ReduceResonanceTolerance = args::get(args_options._2pgf_reduce_tol);
            // Minimal magnitude of the coefficient of a term to take it into account - resolution of thermal weight.
            G4.CoefficientTolerance = args::get(args_options._2pgf_coeff_tol);

            G4.prepare();
            MPI_Barrier(comm);
            std::vector<std::tuple<ComplexType, ComplexType, ComplexType>> freqs_2pgf;
            fmatsubara_grid fgrid(wf_min, wf_max, beta, true);
            bmatsubara_grid bgrid(wb_min, wb_max, beta, true);
            freqs_2pgf.reserve(fgrid.size() * fgrid.size() * bgrid.size());
            for(auto W : bgrid.values()) {
                for(auto w3 : fgrid.values()) {
                    for(auto w2 : fgrid.values()) {
                        ComplexType w1 = W + w3;
                        freqs_2pgf.emplace_back(w1, w2, w3);
                    }
                }
            }
            mpi_cout << "2PGF : " << freqs_2pgf.size() << " freqs to evaluate\n";

            std::vector<ComplexType> chi_freq_data = G4.compute(true, freqs_2pgf, comm);

            // dump 2PGF into files - loop through 2pgf components
            if(!rank) {
                mpi_cout << "Saving 2PGF " << index_comb << '\n';
                grid_object<std::complex<double>, bmatsubara_grid, fmatsubara_grid, fmatsubara_grid> full_vertex(
                    std::forward_as_tuple(bgrid, fgrid, fgrid));
                grid_object<std::complex<double>, fmatsubara_grid, fmatsubara_grid> full_vertex_1freq(
                    std::forward_as_tuple(fgrid, fgrid));
                std::size_t w_ind = 0;
                for(auto W : bgrid.points()) {
                    for(auto w3 : fgrid.points()) {
                        for(auto w2 : fgrid.points()) {
                            std::complex<double> val = chi_freq_data[w_ind];
                            full_vertex[W][w3.index()][w2.index()] = val;
                            full_vertex_1freq[w3.index()][w2.index()] = val;
                            if(!is_float_equal(std::get<0>(freqs_2pgf[w_ind]), W.value() + w3.value()))
                                throw std::logic_error("2PGF freq mismatch");
                            ++w_ind;
                        }
                    }
                    std::string fv1_name = "chi" + ind_str + "_W" + std::to_string(std::imag(W.value())) + ".dat";
                    full_vertex_1freq.savetxt(fv1_name);
                }
                std::string fv_name = "chi" + ind_str + ".dat";
                full_vertex.savetxt(fv_name);
            }
        }
    }
}
