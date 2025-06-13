//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file tutorial/example2site.cpp
/// \brief An example of how to use the pomerol library.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

// In this file we provide a tutorial example of how to actually use the pomerol library.

#include <iostream>
#include <string>
#include <tuple>
#include <vector>

// Include the pomerol library
#include <pomerol.hpp>

// Use Pomerol's main namespace.
// Some MPI-related classes and functions are declared in the pMPI namespace.
using namespace Pomerol;

// Small routine to make fancy screen output for text.
void print_section(std::string const& str);

// Generic tips:
// The calculation is done by computing a set of objects in the following order:
// Hamiltonian expression -> IndexClassification -> HilbertSpace ->
// -> StatesClassification -> Hamiltonian -> CreationOperator/AnnihilationOperator.
//
// (for thermal objects, such as GFs in Matsubara domain)
// -> DensityMatrix -> EnsembleAverage
//                  -> GreensFunction
//                  -> TwoParticleGF -> Vertex4
//                  -> Susceptibility
// The detailed explanation of each class is given below.

int main(int argc, char* argv[]) {
    // Initializing MPI.
    MPI_Init(&argc, &argv);

    // Let us construct a two-site lattice with sites labeled "A" and "B".
    // The sites are connected by a hopping term with matrix element -1.
    RealType t = 1.0;

    // Expression of system's Hamiltonian.
    LatticePresets::RealExpr HExpr = LatticePresets::Hopping("A", "B", -t);

    // Now add interaction. We will use Hubbard-type n_{up} n_{down} interaction
    // on each site. For this and some other common interactions, such as SzSz
    // or SS couplings, shortcuts are provided in the Pomerol::LatticePresets
    // namespace.
    RealType U = 2.0;  // Hubbard interaction constant.
    RealType mu = 1.0; // Chemical potential.

    // LatticePresets::CoulombS adds U n_{up}n_{down} - mu(n_{up} + n_{down})
    // for 1 orbital and 2 spins.
    HExpr += LatticePresets::CoulombS("A", U, -mu);
    HExpr += LatticePresets::CoulombS("B", U, -mu);

    // It is possible to add arbitrary custom terms to 'HExpr' by building them
    // out of fermionic and bosonic creation/annihilation operators. Such operators
    // are returned by functions c(), c_dag(), n(), a() and a_dag() declared in
    // the Pomerol::Operators namespace.
    // More information about expression objects similar to 'HExpr' and means to
    // build and manipulate them can be found in libcommute's documentation under
    // https://krivenko.github.io/libcommute/expression/expression.html.

    // Let us now print 'HExpr'.
    if(pMPI::rank(MPI_COMM_WORLD) == 0)
        std::cout << "HExpr = " << HExpr << '\n';

    // In order to go further, we need to introduce the single-particle index space.
    // A single-particle index is an integer that uniquely identifies a combination of
    // indices carried by a creation/annihilation operator in 'HExpr'. When functions
    // from LatticePresets are used to construct the Hamiltonian expression,
    // each operator in 'HExpr' carries a combination of three indices
    // (site label, orbital index, spin projection). Note, however, that combinations
    // of indices with arbitrary types can be used when 'HExpr' is manually built
    // out of libcommute's creation/annihilation operators.
    // The object that takes care of handling single-particle indices is called
    // IndexClassification.

    // Construct IndexClassification.
    auto IndexInfo = MakeIndexClassification(HExpr);
    // Print which indices we have.
    if(pMPI::rank(MPI_COMM_WORLD) == 0)
        std::cout << "Indices:\n" << IndexInfo << '\n';

    // Let us make a test that our Hamiltonian expression commutes with an operator
    // that represents the total number of particles in the system.
    auto NExpr = Operators::n("A", (unsigned short)0, LatticePresets::up) +
                 Operators::n("A", (unsigned short)0, LatticePresets::down) +
                 Operators::n("B", (unsigned short)0, LatticePresets::up) +
                 Operators::n("B", (unsigned short)0, LatticePresets::down);
    if(pMPI::rank(MPI_COMM_WORLD) == 0) {
        std::cout << "NExpr = " << NExpr << '\n';
        std::cout << "[HExpr, NExpr] = " << (HExpr * NExpr - NExpr * HExpr) << '\n';
    }

    // Having created the Hamiltonian expression and the IndexClassification object
    // we can now introduce system's Hilbert space.
    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();

    // Important remark 1!
    //
    // Many of the objects defined within Pomerol library have the following semantics.
    // They can be constructed, prepared and computed.
    // This means
    //   - constructed: No operations are done except from initializing references to
    //                  other objects that current class depends on.
    //   - prepared: Typically, this is when all memory allocation takes place.
    //   - computed: The actual computation. This is the most costly operation.

    // The Hamiltonian has a set of symmetries. These symmetries allow to partition the
    // Hilbert space into invariant subspaces (sectors) of the Hamiltonian and to
    // effectively reduce its matrix to a block-diagonal form.
    // The StatesClassification object uses a special algorithm to reveal the sectors.

    StatesClassification S;
    S.compute(HS); // Find the invariant subspaces.

    // We shall proceed now with obtaining the spectrum of the Hamiltonian.
    // The Hamiltonian class converts an expression into its block-diagonal matrix
    // representation.

    Hamiltonian H(S);
    // Allocate all diagonal blocks of the Hamiltonian.
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    // Diagonalize the blocks.
    H.compute(MPI_COMM_WORLD);

    // Get ground energy energy.
    std::cout << "The value of ground energy is " << H.getGroundEnergy() << '\n';

    // Important remark 2!
    //
    // All further calculations done in the pomerol code take into account the
    // block structure of the Hamiltonian. All objects that handle matrices
    // and all thermal objects, such as Green's functions are in fact a
    // set of pieces (called "parts") that operate on a certain block or a set of
    // blocks. As such all actual computations are done within these parts and
    // their encompassing objects like Green's functions or Hamiltonian basically
    // just loop over the parts and tell them to call prepare() or compute() methods.

    // At this stage the Hamiltonian is diagonalized and its spectrum and
    // eigenvectors can be directly accessed to calculate some observables.
    //
    // We shall now proceed to the calculations of thermal quantities, e.i.
    // assume that our finite-size system was adiabatically connected to a thermal
    // reservoir that sets certain temperature (in fact, inverse temperature
    // \beta). This means that expectation values of the observables in the system
    // should be calculated with a Gibbs density matrix exp(-\beta H) / Z, rather than
    // by averaging with the ground state. In the eigenbasis of the Hamiltonian the
    // calculation of a density matrix is straightforward - it is just
    // \exp(-\beta (E_i - E_0)) / Z, where E_i is an energy of an excited state,
    // E_0 is the ground state energy, and Z is the partition function.
    // The procedure is done as following:

    // Define inverse temperature
    RealType beta = 10.0;

    // Create the density matrix.
    DensityMatrix rho(S, H, beta);
    // Allocate all internal parts of the density matrix.
    rho.prepare();
    // Actually compute the density matrix.
    rho.compute();
    // Truncate blocks that have only negligibly small contributions.
    rho.truncateBlocks(1e-15);

    // Lehmann representation of the Green's function requires matrices of creation and
    // annihilation operators calculated in the eigenbasis of the Hamiltonian.
    // CreationOperator/AnnihilationOperator are the classes that compute the matrices.

    // Let us create c^\dagger_{"A",up}, c^\dagger_{"A",down} and their conjugates
    ParticleIndex up_index = IndexInfo.getIndex("A", 0, LatticePresets::up);
    ParticleIndex dn_index = IndexInfo.getIndex("A", 0, LatticePresets::down);

    CreationOperator CX_up(IndexInfo, HS, S, H, up_index), CX_dn(IndexInfo, HS, S, H, dn_index);
    CX_up.prepare(HS);
    CX_up.compute();
    CX_dn.prepare(HS);
    CX_dn.compute();

    AnnihilationOperator C_up(IndexInfo, HS, S, H, up_index), C_dn(IndexInfo, HS, S, H, dn_index);
    C_up.prepare(HS);
    C_up.compute();
    C_dn.prepare(HS);
    C_dn.compute();

    print_section("Single-particle Green's function");

    // The local Green's function in the Matsubara domain G_{"A",up}(i\omega_n)
    GreensFunction GF(S, H, C_up, CX_up, rho);
    // Allocate GF parts.
    GF.prepare();
    // Calculate the GF.
    GF.compute();

    for(int n = 0; n < 10; ++n) {
        std::cout << n << " | " << GF(n) << "\n";
    }

    print_section("Two-particle Green's function");

    // The two-particle GF is constructed in analogy to the single-particle GF,
    // it requires 4 operators to be provided though.
    TwoParticleGF Chi(S, H, C_up, C_up, CX_up, CX_up, rho);

    // Some knobs to make calculation faster; The larger the values of tolerances,
    // the faster is the calculation, but rounding errors may show.
    // Here are some settings that give very high precision. If you want to make
    // things faster, and when many values for different frequencies are required,
    // change ReduceResonanceTolerance to something like 10^{-4}.
    //
    // A difference in energies with magnitude below this value is treated as zero.
    Chi.ReduceResonanceTolerance = 1e-8;
    // Minimal magnitude of the coefficient of a term for it to be taken into account.
    Chi.CoefficientTolerance = 1e-16;
    // Minimal magnitude of the coefficient of a term for it to be taken into account
    // with respect to the amount of terms.
    Chi.MultiTermCoefficientTolerance = 1e-6;

    Chi.prepare();
    std::vector<std::tuple<ComplexType, ComplexType, ComplexType>> freqs_2pgf;
    Chi.compute(false, freqs_2pgf, MPI_COMM_WORLD);

    if(pMPI::rank(MPI_COMM_WORLD) == 0) {
        int nm = 2;
        for(int n1 = -nm; n1 < nm; ++n1) {
            for(int n2 = -nm; n2 < nm; ++n2) {
                for(int n3 = -nm; n3 < nm; ++n3) {
                    std::cout << n1 << " " << n2 << " " << n3 << "|" << Chi(n1, n2, n3) << "\n";
                }
            }
        }
    }

    print_section("Quadratic operator");

    // We define a quadratic operator O_{ij} = c^+_i c_j to compute its ensemble average
    // and its fluctuations (dynamical susceptibility).
    // QuadraticOperator is the class that computes and stores the matrix of O_{ij}.

    // Define a quadratic operator O = c^+_{up} c_{up}.
    QuadraticOperator N_up(IndexInfo, HS, S, H, up_index, up_index);
    N_up.prepare(HS);
    N_up.compute();

    print_section("Ensemble average");
    // Compute an ensemble average, <O>

    EnsembleAverage EA(N_up, rho);
    EA.compute();
    RealType occup_up = real(EA());
    if(pMPI::rank(MPI_COMM_WORLD) == 0)
        std::cout << "Occupation number of up spin is " << occup_up << '\n';

    print_section("Quartic operator");

    // It is also possible to compute ensemble average of a quartic operator
    // O_{ijkl} = c^+_i c^+_j c_k c_l.
    // QuarticOperator is the class that computes and stores the matrix of O_{ijkl}.

    // Define a quartic operator O = c^+_{up} c^+_{dn} c_{dn} c_{up}.
    QuarticOperator N_up_N_dn(IndexInfo, HS, S, H, up_index, dn_index, dn_index, up_index);
    N_up_N_dn.prepare(HS);
    N_up_N_dn.compute();

    EnsembleAverage EA2(N_up_N_dn, rho);
    EA2.compute();
    RealType double_occ = real(EA2());
    if(pMPI::rank(MPI_COMM_WORLD) == 0)
        std::cout << "Double occupancy is " << double_occ << '\n';

    print_section("Dynamical susceptibility");

    // The dynamical susceptibility is computed by Susceptibility class.
    // One can obtain either F[ <A(\tau)B> ] or F[ <A(\tau)B> - <A><B> ],
    // where F denotes Fourier transform from \tau to Matsubara frequency.
    // To choose the latter quantity, call subtractDisconnected() method.
    //
    // There are 3 variants of subtractDisconnected() method.
    // 1. <A> and <B> are computed in Susceptibility class;
    // 2. Use precomputed <A> and <B>;
    // 3. Use predefined EnsembleAverage instances for A and B.

    Susceptibility Sus(S, H, N_up, N_up, rho);
    Sus.prepare();
    Sus.compute();
    // Subtract <n_up><n_uo>
    Sus.subtractDisconnected(); // 1
    //Sus.subtractDisconnected(occup_up, occup_up);  // 2
    //Sus.subtractDisconnected(EA, EA);  // 3
    if(pMPI::rank(MPI_COMM_WORLD) == 0) {
        for(int n = 0; n < 10; n++) {
            std::cout << n << "|" << Sus(n) << "\n";
        }
    }

    print_section("3-points susceptibility");
    // The 3-point susceptibility is computed by ThreePointSusceptibility class.
    //
    // It can be defined in one of the following three channels.
    // 1. Particle-particle channel:
    //   \chi^{(3)}_{pp}(\omega_{n_1},\omega_{n_2}) =
    //   \int_0^\beta d\tau_1 d\tau_2 e^{-i\omega_{n_1}\tau_1} e^{-i\omega_{n_2}\tau_2}
    //   Tr[\mathcal{T}_\tau \hat\rho c^+_1(\tau_1) c_2(0^+) c^+_3(\tau_2) c_4(0)]
    //
    // 2. Particle-hole channel:
    //   \chi^{(3)}_{ph}(\omega_{n_1},\omega_{n_2}) =
    //   \int_0^\beta d\tau_1 d\tau_2 e^{-i\omega_{n_1}\tau_1} e^{i\omega_{n_2}\tau_2}
    //   Tr[\mathcal{T}_\tau \hat\rho c^+_1(\tau_1) c_2(\tau_2) c^+_3(0^+) c_4(0)]
    //
    // 3. Crossed particle-hole channel:
    //   \chi^{(3)}_{\bar{ph}}(\omega_{n_1},\omega_{n_2}) =
    //   \int_0^\beta d\tau_1 d\tau_2 e^{-i\omega_{n_1}\tau_1} e^{i\omega_{n_2}\tau_2}
    //   Tr[\mathcal{T}_\tau \hat\rho c^+_1(\tau_1) c_2(0) c^+_3(0^+) c_4(\tau_2)]

    // Particle-particle channel
    // The PP-susceptibility object with indices (up, up, down, down)
    ThreePointSusceptibility chi3pp(PP, S, H, CX_up, C_up, CX_dn, C_dn, rho);
    chi3pp.prepare();
    chi3pp.compute();
    if(pMPI::rank(MPI_COMM_WORLD) == 0) {
        for(int n1 = 0; n1 < 3; n1++) {
            for(int n2 = 0; n2 < 3; n2++) {
                std::cout << n1 << "," << n2 << "|" << chi3pp(n1, n2) << "\n";
            }
        }
    }

    // Particle-hole channel
    // The PH-susceptibility object with indices (up, up, down, down)
    ThreePointSusceptibility chi3ph(PH, S, H, CX_up, C_up, CX_dn, C_dn, rho);
    chi3ph.prepare();
    chi3ph.compute();
    if(pMPI::rank(MPI_COMM_WORLD) == 0) {
        for(int n1 = 0; n1 < 3; n1++) {
            for(int n2 = 0; n2 < 3; n2++) {
                std::cout << n1 << "," << n2 << "|" << chi3ph(n1, n2) << "\n";
            }
        }
    }

    // Crossed particle-hole channel
    // The xPH-susceptibility object with indices (up, up, down, down)
    ThreePointSusceptibility chi3xph(xPH, S, H, CX_up, C_up, CX_dn, C_dn, rho);
    chi3xph.prepare();
    chi3xph.compute();
    if(pMPI::rank(MPI_COMM_WORLD) == 0) {
        for(int n1 = 0; n1 < 3; n1++) {
            for(int n2 = 0; n2 < 3; n2++) {
                std::cout << n1 << "," << n2 << "|" << chi3xph(n1, n2) << "\n";
            }
        }
    }

    // Shut down MPI.
    MPI_Finalize();
}

void print_section(std::string const& str) {
    if(!pMPI::rank(MPI_COMM_WORLD)) {
        std::cout << std::string(str.size(), '=') << '\n';
        std::cout << str << '\n';
        std::cout << std::string(str.size(), '=') << '\n';
    }
}
