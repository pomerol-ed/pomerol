//
// Created by iskakoff on 05/12/16.
//

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <tuple>

#include "quantum_model.h"

#undef DEBUG
#include <gftools.hpp>

using namespace Pomerol;

quantum_model::quantum_model(int argc, char* argv[], const std::string &prog_desc) :
  args_parser(prog_desc),
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
               {args_parser, "tol", "Tolerance on nominators in 2PGF", {"2pgf.coeff_tol"}, 1e-12},
               {args_parser, "tol", "How often to reduce terms in 2PGF", {"2pgf.multiterm_tol"}, 1e-6}
  },
  comm(MPI_COMM_WORLD) {
  MPI_Init(&argc, &argv);
  rank = pMPI::rank(comm);
  // Print defaults in the usage message
  args_parser.helpParams.addDefault = true;
}

quantum_model::~quantum_model() {
  MPI_Finalize();
}

void quantum_model::parse_args(int argc, char* argv[]) {
  try {
    args_parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << args_parser;
    exit(0);
  }
  catch (args::ParseError& e)
  {
    std::cerr << e.what() << std::endl;
    std::cerr << args_parser;
    exit(1);
  }
  catch (args::ValidationError& e)
  {
    std::cerr << e.what() << std::endl;
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
  if (!rank) {
    print_section("Indices");
    std::cout << IndexInfo << std::endl;
  };

  auto HS = MakeHilbertSpace(IndexInfo, HExpr);
  HS.compute();

  StatesClassification S;
  S.compute(HS);

  Hamiltonian H(S);
  H.prepare(HExpr, HS, MPI_COMM_WORLD);
  H.compute(MPI_COMM_WORLD);

  if(!rank) {
    gftools::grid_object<double, gftools::enum_grid> evals1(gftools::enum_grid(0, S.getNumberOfStates()));
    RealVectorType evals (H.getEigenValues());
    std::sort(evals.data(), evals.data() + H.getEigenValues().size());
    std::copy(evals.data(), evals.data() + S.getNumberOfStates(), evals1.data().data());
    evals1.savetxt("spectrum.dat");
  }
  DensityMatrix rho(S,H,beta); // create Density Matrix
  rho.prepare();
  rho.compute(); // evaluate thermal weights with respect to ground energy, i.e exp(-beta(e-e_0))/Z

  std::pair<ParticleIndex, ParticleIndex> pair = get_node(IndexInfo);
  ParticleIndex d0 = pair.first;//IndexInfo.getIndex("A",0,down); // find the indices of the impurity, i.e. spin up index
  ParticleIndex u0 = pair.second;//IndexInfo.getIndex("A",0,up);

  // Green's function calculation starts here

  if (calc_gf) {
    print_section("1-particle Green's functions calc");
    std::set<ParticleIndex> f; // a set of indices to evaluate c and c^+
    std::set<IndexCombination2> indices2; // a set of pairs of indices to evaluate Green's function

    double eta = args::get(args_options.gf_eta);
    double step = args::get(args_options.gf_step);
    double hbw = args::get(args_options.gf_D);

    int wf_min = args::get(args_options.wf_min);
    int wf_max = args::get(args_options.wf_max);
    int wb_min = args::get(args_options.wb_min);
    int wb_max = args::get(args_options.wb_max);

    // Take only impurity spin up and spin down indices
    f.insert(u0);
    f.insert(d0);
    prepare_indices(d0, u0, indices2, f, IndexInfo);

    // Create a container for c and c^+ in the eigenstate basis
    FieldOperatorContainer Operators(IndexInfo, HS, S, H, f);
    Operators.prepareAll(HS);
    Operators.computeAll(); // evaluate c, c^+ for chosen indices

    GFContainer G(IndexInfo,S,H,rho,Operators);

    G.prepareAll(indices2); // identify all non-vanishing block connections in the Green's function
    G.computeAll(); // Evaluate all GF terms, i.e. resonances and weights of expressions in Lehmans representation of the Green's function

    if (!rank) // dump gf into a file
      // loops over all components (pairs of indices) of the Green's function
      for (std::set<IndexCombination2>::const_iterator it = indices2.begin(); it != indices2.end(); ++it) {
        IndexCombination2 ind2 = *it;
        const GreensFunction & GF = G(ind2);
        // Save Matsubara GF from pi/beta to pi/beta*(4*wf_max + 1)
        std::cout << "Saving imfreq G" << ind2 << " on " << 4 * wf_max << " Matsubara freqs. " << std::endl;
        grid_object<std::complex<double>, fmatsubara_grid> gf_imfreq (fmatsubara_grid(wf_min, wf_max*4, beta, true));
        std::string ind_str = std::to_string(ind2.Index1)+ std::to_string(ind2.Index2);
        for (auto p : gf_imfreq.grid().points()) { gf_imfreq[p] = GF(p.value()); }
        gf_imfreq.savetxt("gw_imfreq_"+ ind_str +".dat");

        real_grid freq_grid(-hbw, hbw, 2*hbw/step+1, true);
        grid_object<std::complex<double>, real_grid> gf_refreq(freq_grid);
        for (auto p : freq_grid.points()) {
          ComplexType val = GF(ComplexType(p.value()) + I*eta);
          gf_refreq[p] = val;
        };
        gf_refreq.savetxt("gw_refreq_"+ ind_str +".dat");
      }

    // Start Two-particle GF calculation

    if (calc_2pgf) {
      print_section("2-Particle Green's function calc");

      std::vector<size_t> indices_2pgf = args::get(args_options._2pgf_indices);

      if (indices_2pgf.size() != 4)
        throw std::logic_error("Need 4 indices for 2PGF");

      // a set of four indices to evaluate the 2pgf
      IndexCombination4 index_comb(indices_2pgf[0], indices_2pgf[1], indices_2pgf[2], indices_2pgf[3]);

      std::set<IndexCombination4> indices4;
      // 2PGF = <T c c c^+ c^+>
      indices4.insert(index_comb);
      std::string ind_str = std::to_string(index_comb.Index1)
                          + std::to_string(index_comb.Index2)
                          + std::to_string(index_comb.Index3)
                          + std::to_string(index_comb.Index4);

      AnnihilationOperator const& C1 = Operators.getAnnihilationOperator(index_comb.Index1);
      AnnihilationOperator const& C2 = Operators.getAnnihilationOperator(index_comb.Index2);
      CreationOperator const&    CX3 = Operators.getCreationOperator(index_comb.Index3);
      CreationOperator const&    CX4 = Operators.getCreationOperator(index_comb.Index4);
      TwoParticleGF G4(S, H, C1, C2, CX3, CX4, rho);

      /* Some knobs to make calc faster - the larger the values of tolerances, the faster is calc, but rounding errors may show up. */
      /** A difference in energies with magnitude less than this value is treated as zero - resolution of energy resonances. */
      G4.ReduceResonanceTolerance = args::get(args_options._2pgf_reduce_tol);
      /** Minimal magnitude of the coefficient of a term to take it into account - resolution of thermal weight. */
      G4.CoefficientTolerance = args::get(args_options._2pgf_coeff_tol);
      /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
      G4.MultiTermCoefficientTolerance = args::get(args_options._2pgf_multiterm_tol);

      G4.prepare();
      MPI_Barrier(comm);
      std::vector<std::tuple<ComplexType, ComplexType, ComplexType>> freqs_2pgf;
      fmatsubara_grid fgrid(wf_min, wf_max, beta, true);
      bmatsubara_grid bgrid(wb_min, wb_max, beta, true);
      freqs_2pgf.reserve(fgrid.size() * fgrid.size() * bgrid.size());
      for (auto W : bgrid.values()) {
        for (auto w3 : fgrid.values()) {
          for (auto w2 : fgrid.values()) {
            ComplexType w1 = W+w3;
            freqs_2pgf.push_back(std::make_tuple(w1,w2,w3));
          }
        }
      }
      mpi_cout << "2PGF : " << freqs_2pgf.size() << " freqs to evaluate" << std::endl;

      std::vector<ComplexType> chi_freq_data = G4.compute(true, freqs_2pgf, comm);

      // dump 2PGF into files - loop through 2pgf components
      if (!rank) {
        mpi_cout << "Saving 2PGF " << index_comb << std::endl;
        grid_object<std::complex<double>, bmatsubara_grid, fmatsubara_grid, fmatsubara_grid> full_vertex(std::forward_as_tuple(bgrid, fgrid, fgrid));
        grid_object<std::complex<double>, fmatsubara_grid, fmatsubara_grid> full_vertex_1freq(std::forward_as_tuple(fgrid, fgrid));
        size_t w_ind = 0;
        for (auto W : bgrid.points()) {
          for (auto w3 : fgrid.points()) {
            for (auto w2 : fgrid.points()) {
              std::complex<double> val = chi_freq_data[w_ind];
              full_vertex[W][w3.index()][w2.index()] = val;
              full_vertex_1freq[w3.index()][w2.index()] = val;
              if (!is_float_equal(std::get<0>(freqs_2pgf[w_ind]), W.value()+w3.value()))
                throw std::logic_error("2PGF freq mismatch");
              ++w_ind;
            }
          }
          std::string fv1_name = "chi"+ind_str+"_W"+std::to_string(std::imag(W.value()))+".dat";
          full_vertex_1freq.savetxt(fv1_name);
        }
      }
    }
  }
}
