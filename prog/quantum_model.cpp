//
// Created by iskakoff on 05/12/16.
//

#include "quantum_model.h"

void quantum_model::compute() {
  IndexClassification IndexInfo(Lat.getSiteMap());
  IndexInfo.prepare(false); // Create index space
  if (!rank) { print_section("Indices"); IndexInfo.printIndices(); };
  int index_size = IndexInfo.getIndexSize();

  print_section("Matrix element storage");
  IndexHamiltonian Storage(&Lat,IndexInfo);
  Storage.prepare(); // Write down the Hamiltonian as a symbolic formula
  print_section("Terms");
  mpi_cout << Storage << std::endl;

  Symmetrizer Symm(IndexInfo, Storage);
  Symm.compute(); // Find symmetries of the problem

  StatesClassification S(IndexInfo,Symm); // Introduce Fock space and classify states to blocks
  S.compute();

  Hamiltonian H(IndexInfo, Storage, S); // Hamiltonian in the basis of Fock Space
  H.prepare(); // enter the Hamiltonian matrices
  H.compute(); // compute eigenvalues and eigenvectors

  if(!rank) {
#ifdef POMEROL_CXX11
    gftools::grid_object<double, gftools::enum_grid> evals1(gftools::enum_grid(0, S.getNumberOfStates()));
    RealVectorType evals (H.getEigenValues());
    std::sort(evals.data(), evals.data() + H.getEigenValues().size());
    std::copy(evals.data(), evals.data() + S.getNumberOfStates(), evals1.data().data());
    evals1.savetxt("spectrum.dat");
#else
    RealVectorType evals(H.getEigenValues());
    std::sort(evals.data(), evals.data() + H.getEigenValues().size());
    savetxt("spectrum.dat", evals); // dump eigenvalues
#endif
  }
  DensityMatrix rho(S,H,beta); // create Density Matrix
  rho.prepare();
  rho.compute(); // evaluate thermal weights with respect to ground energy, i.e exp(-beta(e-e_0))/Z

  mpi_cout << "<N> = " << rho.getAverageOccupancy() << std::endl; // get average total particle number
  mpi_cout << "<H> = " << rho.getAverageEnergy() << std::endl; // get average energy
  std::pair<ParticleIndex, ParticleIndex> pair = get_node(IndexInfo);
  ParticleIndex d0 = pair.first;//IndexInfo.getIndex("A",0,down); // find the indices of the impurity, i.e. spin up index
  ParticleIndex u0 = pair.second;//IndexInfo.getIndex("A",0,up);
  mpi_cout << "<N_{" << IndexInfo.getInfo(u0) << "}N_{"<< IndexInfo.getInfo(u0) << "}> = " << rho.getAverageDoubleOccupancy(u0,d0) << std::endl; // get double occupancy

  for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); i++) {
    mpi_cout << "<N_{" << IndexInfo.getInfo(i) << "[" << i <<"]}> = " << rho.getAverageOccupancy(i) << std::endl; // get average total particle number
  }

  if (!comm.rank()) {
#ifdef POMEROL_CXX11
    double n_av = rho.getAverageOccupancy();
    gftools::num_io<double>(n_av).savetxt("N_T.dat");
#else
    savetxt("N_T.dat",rho.getAverageOccupancy());
#endif
  }

  // Green's function calculation starts here

  FieldOperatorContainer Operators(IndexInfo, S, H); // Create a container for c and c^+ in the eigenstate basis

  if (calc_gf) {
    print_section("1-particle Green's functions calc");
    std::set<ParticleIndex> f; // a set of indices to evaluate c and c^+
    std::set<IndexCombination2> indices2; // a set of pairs of indices to evaluate Green's function

    int ntau; double eta, step, hbw;
    boost::tuples::tie(ntau, eta, step, hbw) = boost::tuples::make_tuple(p["gf.ntau"].as<int>(), p["gf.eta"].as<double>(), p["gf.step"].as<double>(), p["gf.D"].as<double>());
    int wf_min, wf_max, wb_min, wb_max;
    boost::tuples::tie(wf_min, wf_max, wb_min, wb_max) = boost::tuples::make_tuple(p["wf_min"].as<int>(), p["wf_max"].as<int>(), p["wb_min"].as<int>(), p["wb_max"].as<int>());

    // Take only impurity spin up and spin down indices
    f.insert(u0);
    f.insert(d0);
    prepare_indices(d0, u0, indices2, f, IndexInfo);

    Operators.prepareAll(f);
    Operators.computeAll(); // evaluate c, c^+ for chosen indices

    GFContainer G(IndexInfo,S,H,rho,Operators);

    G.prepareAll(indices2); // identify all non-vanishing block connections in the Green's function
    G.computeAll(); // Evaluate all GF terms, i.e. resonances and weights of expressions in Lehmans representation of the Green's function

    if (!comm.rank()) // dump gf into a file
      // loops over all components (pairs of indices) of the Green's function
      for (std::set<IndexCombination2>::const_iterator it = indices2.begin(); it != indices2.end(); ++it) {
        IndexCombination2 ind2 = *it;
        const GreensFunction & GF = G(ind2);
        // Save Matsubara GF from pi/beta to pi/beta*(4*wf_max + 1)
        std::cout << "Saving imfreq G" << ind2 << " on " << 4 * wf_max << " Matsubara freqs. " << std::endl;
#ifdef POMEROL_CXX11
        grid_object<std::complex<double>, fmatsubara_grid> gf_imfreq (fmatsubara_grid(wf_min, wf_max*4, beta, true));
        std::string ind_str = boost::lexical_cast<std::string>(ind2.Index1)+ boost::lexical_cast<std::string>(ind2.Index2);
        for (auto p : gf_imfreq.grid().points()) { gf_imfreq[p] = GF(p.value()); }
        gf_imfreq.savetxt("gw_imfreq_"+ ind_str +".dat");

        real_grid freq_grid(-hbw, hbw, 2*hbw/step+1, true);
        grid_object<std::complex<double>, real_grid> gf_refreq(freq_grid);
        for (auto p : freq_grid.points()) {
          ComplexType val = GF(ComplexType(p.value()) + I*eta);
          gf_refreq[p] = val;
        };
        gf_refreq.savetxt("gw_refreq_"+ ind_str +".dat");
#else
        std::ofstream gw_im(("gw_imag" + boost::lexical_cast< std::string>(ind2.Index1) + boost::lexical_cast< std::string>(ind2.Index2) + ".dat").c_str());
        for (int wn = 0; wn <= wf_max*4; wn++) {
          ComplexType val = GF(I*FMatsubara(wn,beta)); // this comes from Pomerol - see GreensFunction::operator()
          gw_im << std::scientific << std::setprecision(12) << FMatsubara(wn, beta) << "   " << real(val) << " " << imag(val) << std::endl;
        };
        gw_im.close();
        // Save Retarded GF on the real axis
        std::ofstream gw_re(("gw_real" + boost::lexical_cast< std::string>(ind2.Index1) + boost::lexical_cast< std::string>(ind2.Index2) + ".dat").c_str());
        std::cout << "Saving real-freq GF " << ind2 << " in energy space [" << _e0 - hbw << ":" << _e0 + hbw << ":" << step << "] + I*" << eta << "." << std::endl;
        for (double w = _e0-hbw; w < _e0+hbw; w+=step) {
          ComplexType val = GF(ComplexType(w) + I*eta);
          gw_re << std::scientific << std::setprecision(12) << w << "   " << real(val) << " " << imag(val) << std::endl;
        };
        gw_re.close();
#endif
      }

    // Start Two-particle GF calculation

    if (calc_2pgf) {
      print_section("2-Particle Green's function calc");

      std::vector<size_t> indices_2pgf = p["2pgf.indices"].as< std::vector<size_t> >();
      if (indices_2pgf.size() != 4) throw std::logic_error("Need 4 indices for 2pgf");

      // a set of four indices to evaluate the 2pgf
      IndexCombination4 index_comb(indices_2pgf[0], indices_2pgf[1], indices_2pgf[2], indices_2pgf[3]);

      std::set<IndexCombination4> indices4;
      // 2PGF = <T c c c^+ c^+>
      indices4.insert(index_comb);
      std::string ind_str = boost::lexical_cast< std::string>(index_comb.Index1)
                            + boost::lexical_cast< std::string>(index_comb.Index2)
                            + boost::lexical_cast< std::string>(index_comb.Index3)
                            + boost::lexical_cast< std::string>(index_comb.Index4);

      AnnihilationOperator const& C1 = Operators.getAnnihilationOperator(index_comb.Index1);
      AnnihilationOperator const& C2 = Operators.getAnnihilationOperator(index_comb.Index2);
      CreationOperator const&    CX3 = Operators.getCreationOperator(index_comb.Index3);
      CreationOperator const&    CX4 = Operators.getCreationOperator(index_comb.Index4);
      TwoParticleGF G4(S, H, C1, C2, CX3, CX4, rho);

      /* Some knobs to make calc faster - the larger the values of tolerances, the faster is calc, but rounding errors may show up. */
      /** A difference in energies with magnitude less than this value is treated as zero - resolution of energy resonances. */
      G4.ReduceResonanceTolerance = p["2pgf.reduce_tol"].as<double>();
      /** Minimal magnitude of the coefficient of a term to take it into account - resolution of thermal weight. */
      G4.CoefficientTolerance = p["2pgf.coeff_tol"].as<double>();
      /** Knob that controls the caching frequency. */
      G4.ReduceInvocationThreshold = p["2pgf.reduce_freq"].as<size_t>();
      /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
      G4.MultiTermCoefficientTolerance = p["2pgf.multiterm_tol"].as<double>();

      G4.prepare();
      comm.barrier(); // MPI::BARRIER
#ifdef POMEROL_CXX11
      std::vector<boost::tuple<ComplexType, ComplexType, ComplexType> > freqs_2pgf;
      fmatsubara_grid fgrid(wf_min, wf_max, beta, true);
      bmatsubara_grid bgrid(wb_min, wb_max, beta, true);
      freqs_2pgf.reserve(fgrid.size() * fgrid.size() * bgrid.size());
      for (auto W : bgrid.values()) {
        for (auto w3 : fgrid.values()) {
          for (auto w2 : fgrid.values()) {
            ComplexType w1 = W+w3;
            freqs_2pgf.push_back(boost::make_tuple(w1,w2,w3));
          }
        }
      }
#else
      // Fill a vector of tuples of fermionic Matsubara frequencies - these will be evaluated in-place
      std::vector< boost::tuples::tuple<ComplexType, ComplexType, ComplexType> > freqs_2pgf;
      for (int W_index = -wb_min; W_index <= wb_max; W_index++) { // loop over bosonic freq
        ComplexType W = I*BMatsubara(W_index, beta);
        for (int w3_index = -wf_max; w3_index<=wf_max; w3_index++) { // loop over first fermionic frequency
          ComplexType w3 = I*FMatsubara(w3_index, beta);
          for (int w2_index = -wf_max; w2_index<=wf_max; w2_index++) { // loop over second fermionic
            ComplexType w2 = I*FMatsubara(w2_index, beta);
            ComplexType w1 = W+w3;
            freqs_2pgf.push_back(boost::tuples::make_tuple(w1, w2, w3));
          }
        }
      }
#endif
      mpi_cout << "2PGF : " << freqs_2pgf.size() << " freqs to evaluate" << std::endl;

      std::vector<ComplexType> chi_freq_data = G4.compute(true, freqs_2pgf, comm); // mdata[ind];

#ifdef POMEROL_CXX11
      // dump 2PGF into files - loop through 2pgf components
      if (!comm.rank()) {
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
              if (!is_float_equal(boost::get<0>(freqs_2pgf[w_ind]), W.value()+w3.value())) throw std::logic_error("2pgf freq mismatch");
              ++w_ind;
            }
          }
          std::string fv1_name = "chi"+ind_str+"_W"+boost::lexical_cast<std::string>(std::imag(W.value()))+".dat";
          full_vertex_1freq.savetxt(fv1_name);
        }
      }
#else
      // with the help of alpscore grid_object could be dumped to hdf5 - see save_grid_object(alps::hdf5::archive&, grid_object const&, group)
      // Save terms of two particle GF
      std::ofstream term_res_stream(("terms_res" + ind_str + ".pom").c_str());
      std::ofstream term_nonres_stream(("terms_nonres" + ind_str + ".pom").c_str());
      boost::archive::text_oarchive oa_res(term_res_stream);
      boost::archive::text_oarchive oa_nonres(term_nonres_stream);
      for(std::vector<TwoParticleGFPart*>::const_iterator iter = G4.parts.begin(); iter != G4.parts.end(); iter++) {
        oa_nonres << ((*iter)->getNonResonantTerms());
        oa_res << ((*iter)->getResonantTerms());
      };

      if (!rank) {
        size_t w = 0;
        for (int W_index = -wb_min; W_index <= wb_max; W_index++) { // loop over bosonic freq
          ComplexType W = I*BMatsubara(W_index, beta);
          std::ofstream chi_stream (("chi" + ind_str + "_W" + boost::lexical_cast< std::string>(imag(W)) + ".dat").c_str());
          for (int w3_index = -wf_max; w3_index<=wf_max; w3_index++) { // loop over first fermionic frequency
            ComplexType w3 = I*FMatsubara(w3_index, beta);
            for (int w2_index = -wf_max; w2_index<=wf_max; w2_index++) { // loop over second fermionic
              ComplexType w2 = I*FMatsubara(w2_index, beta);
              ComplexType w1 = W+w3;
              if (::std::abs(w1 - boost::tuples::get<0>(freqs_2pgf[w])) > 1e-8) throw std::logic_error("2pgf freq mismatch");

              std::complex<double> val = chi_freq_data[w];

              chi_stream << std::scientific << std::setprecision(12)
                         << w3.real() << " " << w3.imag() << " " << w2.real() << " " << w2.imag() << "   " << real(val) << " " << imag(val) << std::endl;
              ++w;
            }
          }
          chi_stream << std::endl;
          chi_stream.close();
        }
      }
#endif
    }
  }
}