//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/LatticePresets.hpp
/// \brief Factory functions for terms commonly used to construct various lattice Hamiltonians.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_POMEROL_LATTICEPRESETS_HPP
#define POMEROL_INCLUDE_POMEROL_LATTICEPRESETS_HPP

#include "Misc.hpp"
#include "Operators.hpp"

#include <ostream>
#include <string>

namespace Pomerol {

/// This namespace encloses factory functions for various terms most commonly used to write a lattice Hamiltonian.
namespace LatticePresets {

/// Possible values of spin-1/2 z-projection.
enum spin : short {
    undef = -1, ///< Undefined (useful for bosonic degrees of freedom)
    down = 0,   ///< Spin down
    up = 1      ///< Spin up
};

/// Output stream insertion operator for the spin projection values.
/// \param[out] os Output stream.
/// \param[in] s Spin projection to be inserted.
/// \return Reference to the output stream.
std::ostream& operator<<(std::ostream& os, spin s);

/// Real-valued expression built out of lattice creation/annihilation operators.
/// Each operator in the expression carries a site name label (a string index),
/// an integer orbital index and a spin index.
using RealExpr = Operators::expression<RealType, std::string, unsigned short, spin>;
/// Complex-valued expression built out of lattice creation/annihilation operators.
/// Each operator in the expression carries a site name label (a string index),
/// an integer orbital index and a spin index.
using ComplexExpr = Operators::expression<ComplexType, std::string, unsigned short, spin>;

//
// Overloads of Level()
//

/// \defgroup LatticePresets Factory functions to construct widely used lattice models

/// \defgroup Level Factory functions for single fermion level terms
/// \ingroup LatticePresets
///@{

/// Make a single energy level term \f$\varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma}\f$
/// for a fermion on a given site for a given spin and orbital.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] Eps Real energy level \f$\varepsilon\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
/// \param[in] Spin Spin component \f$\sigma\f$.
RealExpr Level(std::string const& Label, RealType Eps, unsigned short Orbital, spin Spin);
/// Make a single energy level term \f$\varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma}\f$
/// for a fermion on a given site for a given spin and orbital.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] Eps Complex energy level \f$\varepsilon\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
/// \param[in] Spin Spin component \f$\sigma\f$.
ComplexExpr Level(std::string const& Label, ComplexType Eps, unsigned short Orbital, spin Spin);

/// Make a sum of real energy fermionic terms
/// \f$ \sum\limits_{\alpha, \sigma} \varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma}\f$.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] Eps Energy level \f$\varepsilon\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
RealExpr Level(std::string const& Label, RealType Eps, unsigned short NOrbitals = 1);
/// Make a sum of complex energy fermionic terms
/// \f$ \sum\limits_{\alpha, \sigma} \varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma}\f$.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] Eps Energy level \f$\varepsilon\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
ComplexExpr Level(std::string const& Label, ComplexType Eps, unsigned short NOrbitals = 1);

///@}

//
// Overloads of Hopping()
//

/// \defgroup Hopping Factory functions for fermionic hopping terms
/// \ingroup LatticePresets
///@{

/// Make a fermionic hopping term \f$t c^\dagger_{i\alpha_1\sigma_1}c_{j\alpha_2\sigma_2} + h.c.\f$
/// between two lattice sites \f$i \neq j\f$ with a real hopping matrix element.
/// \param[in] Label1 The first lattice site \f$i\f$ connected by this term.
/// \param[in] Label2 The second lattice site \f$j\f$ connected by this term.
/// \param[in] t Hopping matrix element \f$t\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$ connected by this term.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$ connected by this term.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$ on site \f$i\f$ connected by this term.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$ on site \f$j\f$ connected by this term.
RealExpr Hopping(std::string const& Label1,
                 std::string const& Label2,
                 RealType t,
                 unsigned short Orbital1,
                 unsigned short Orbital2,
                 spin Spin1,
                 spin Spin2);
/// Make a fermionic hopping term \f$t c^\dagger_{i\alpha_1\sigma_1}c_{j\alpha_2\sigma_2} + h.c.\f$
/// between two lattice sites \f$i \neq j\f$ with a complex hopping matrix element.
/// \param[in] Label1 The first lattice site \f$i\f$ connected by this term.
/// \param[in] Label2 The second lattice site \f$j\f$ connected by this term.
/// \param[in] t Hopping matrix element \f$t\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$ connected by this term.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$ connected by this term.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$ on site \f$i\f$ connected by this term.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$ on site \f$j\f$ connected by this term.
ComplexExpr Hopping(std::string const& Label1,
                    std::string const& Label2,
                    ComplexType t,
                    unsigned short Orbital1,
                    unsigned short Orbital2,
                    spin Spin1,
                    spin Spin2);

/// Make a fermionic hopping term \f$t c^\dagger_{i\alpha\sigma}c_{j\alpha\sigma} + h.c.\f$
/// between two lattice sites \f$i \neq j\f$ with a real hopping matrix element.
/// \param[in] Label1 The first lattice site \f$i\f$ connected by this term.
/// \param[in] Label2 The second lattice site \f$j\f$ connected by this term.
/// \param[in] t Hopping matrix element \f$t\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
/// \param[in] Spin Spin component \f$\sigma\f$.
RealExpr Hopping(std::string const& Label1, std::string const& Label2, RealType t, unsigned short Orbital, spin Spin);
/// Make a fermionic hopping term \f$t c^\dagger_{i\alpha\sigma}c_{j\alpha\sigma} + h.c.\f$
/// between two lattice sites \f$i \neq j\f$ with a complex hopping matrix element.
/// \param[in] Label1 The first lattice site \f$i\f$ connected by this term.
/// \param[in] Label2 The second lattice site \f$j\f$ connected by this term.
/// \param[in] t Hopping matrix element \f$t\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
/// \param[in] Spin Spin component \f$\sigma\f$.
ComplexExpr
Hopping(std::string const& Label1, std::string const& Label2, ComplexType t, unsigned short Orbital, spin Spin);

/// Make a fermionic hopping term \f$t \sum_\sigma c^\dagger_{i\alpha_1\sigma}c_{j\alpha_2\sigma} + h.c.\f$
/// between two lattice sites \f$i \neq j\f$ with a real hopping matrix element.
/// \param[in] Label1 The first lattice site \f$i\f$ connected by this term.
/// \param[in] Label2 The second lattice site \f$j\f$ connected by this term.
/// \param[in] t Hopping matrix element \f$t\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$ connected by this term.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$ connected by this term.
RealExpr Hopping(std::string const& Label1,
                 std::string const& Label2,
                 RealType t,
                 unsigned short Orbital1,
                 unsigned short Orbital2);
/// Make a fermionic hopping term \f$t \sum_\sigma c^\dagger_{i\alpha_1\sigma}c_{j\alpha_2\sigma} + h.c.\f$
/// between two lattice sites \f$i \neq j\f$ with a complex hopping matrix element.
/// \param[in] Label1 The first lattice site \f$i\f$ connected by this term.
/// \param[in] Label2 The second lattice site \f$j\f$ connected by this term.
/// \param[in] t Hopping matrix element \f$t\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$ connected by this term.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$ connected by this term.
ComplexExpr Hopping(std::string const& Label1,
                    std::string const& Label2,
                    ComplexType t,
                    unsigned short Orbital1,
                    unsigned short Orbital2);

/// Make a fermionic hopping term \f$t \sum_{\alpha\sigma} c^\dagger_{i\alpha\sigma}c_{j\alpha\sigma} + h.c.\f$
/// between two lattice sites \f$i \neq j\f$ with a real hopping matrix element.
/// \param[in] Label1 The first lattice site \f$i\f$ connected by this term.
/// \param[in] Label2 The second lattice site \f$j\f$ connected by this term.
/// \param[in] t Hopping matrix element \f$t\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
RealExpr Hopping(std::string const& Label1, std::string const& Label2, RealType t, unsigned short NOrbitals = 1);
/// Make a fermionic hopping term \f$t \sum_{\alpha\sigma} c^\dagger_{i\alpha\sigma}c_{j\alpha\sigma} + h.c.\f$
/// between two lattice sites \f$i \neq j\f$ with a complex hopping matrix element.
/// \param[in] Label1 The first lattice site \f$i\f$ connected by this term.
/// \param[in] Label2 The second lattice site \f$j\f$ connected by this term.
/// \param[in] t Hopping matrix element \f$t\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
ComplexExpr Hopping(std::string const& Label1, std::string const& Label2, ComplexType t, unsigned short NOrbitals = 1);

///@}

//
// Overloads of Magnetization()
//

/// \defgroup Magnetization Factory functions for magnetization terms
/// \ingroup LatticePresets
///@{

/// Make a magnetic splitting term \f$H \sum_\alpha (n_{i\alpha\uparrow} - n_{i\alpha\downarrow})\f$
/// with a real magnetization constant.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] H Magnetization constant \f$H\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
RealExpr Magnetization(std::string const& Label, RealType H, unsigned short NOrbitals = 1);
/// Make a magnetic splitting term \f$H \sum_\alpha (n_{i\alpha\uparrow} - n_{i\alpha\downarrow})\f$
/// with a complex magnetization constant.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] H Magnetization constant \f$H\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
ComplexExpr Magnetization(std::string const& Label, ComplexType H, unsigned short NOrbitals = 1);

///@}

//
// Overloads of Pairing()
//

/// \defgroup Pairing Factory functions for pairing terms
/// \ingroup LatticePresets
///@{

/// Make a pairing term \f$\Delta c^\dagger_{i\alpha_1\sigma_1}c^\dagger_{j\alpha_2\sigma_2} + h.c.\f$
/// with a real pairing amplitude \f$\Delta\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] Delta Pairing amplitude \f$\Delta\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$ on site \f$i\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$ on site \f$j\f$.
RealExpr Pairing(std::string const& Label1,
                 std::string const& Label2,
                 RealType Delta,
                 unsigned short Orbital1,
                 unsigned short Orbital2,
                 spin Spin1,
                 spin Spin2);
/// Make a pairing term \f$\Delta c^\dagger_{i\alpha_1\sigma_1}c^\dagger_{j\alpha_2\sigma_2} + h.c.\f$
/// with a complex pairing amplitude \f$\Delta\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] Delta Pairing amplitude \f$\Delta\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$ on site \f$i\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$ on site \f$j\f$.
ComplexExpr Pairing(std::string const& Label1,
                    std::string const& Label2,
                    ComplexType Delta,
                    unsigned short Orbital1,
                    unsigned short Orbital2,
                    spin Spin1,
                    spin Spin2);

/// Make a pairing term \f$\Delta c^\dagger_{i\alpha_1\uparrow}c^\dagger_{j\alpha_2\downarrow} + h.c.\f$
/// with a real pairing amplitude \f$\Delta\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] Delta Pairing amplitude \f$\Delta\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$.
RealExpr Pairing(std::string const& Label1,
                 std::string const& Label2,
                 RealType Delta,
                 unsigned short Orbital1,
                 unsigned short Orbital2);
/// Make a pairing term \f$\Delta c^\dagger_{i\alpha_1\uparrow}c^\dagger_{j\alpha_2\downarrow} + h.c.\f$
/// with a complex pairing amplitude \f$\Delta\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] Delta Pairing amplitude \f$\Delta\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$.
ComplexExpr Pairing(std::string const& Label1,
                    std::string const& Label2,
                    ComplexType Delta,
                    unsigned short Orbital1,
                    unsigned short Orbital2);

/// Make a local pairing term \f$\Delta \sum_\alpha c^\dagger_{i\alpha\uparrow}c^\dagger_{i\alpha\downarrow} + h.c.\f$
/// with a real pairing amplitude \f$\Delta\f$.
/// \param[in] Label The lattice site \f$i\f$.
/// \param[in] Delta Pairing amplitude \f$\Delta\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
RealExpr Pairing(std::string const& Label, RealType Delta, unsigned short NOrbitals = 1);
/// Make a local pairing term \f$\Delta \sum_\alpha c^\dagger_{i\alpha\uparrow}c^\dagger_{i\alpha\downarrow} + h.c.\f$
/// with a complex pairing amplitude \f$\Delta\f$.
/// \param[in] Label The lattice site \f$i\f$.
/// \param[in] Delta Pairing amplitude \f$\Delta\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
ComplexExpr Pairing(std::string const& Label, ComplexType Delta, unsigned short NOrbitals = 1);

///@}

//
// Overloads of NupNdown()
//

/// \defgroup NupNdown Factory functions for fermionic density-density interaction terms
/// \ingroup LatticePresets
///@{

/// Make a fermionic density-density interaction term \f$ U n_{i\alpha_1\sigma_1}n_{j\alpha_2\sigma_2}\f$
/// with a real interaction strength \f$U\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$ on site \f$i\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$ on site \f$j\f$.
RealExpr NupNdown(std::string const& Label1,
                  std::string const& Label2,
                  RealType U,
                  unsigned short Orbital1,
                  unsigned short Orbital2,
                  spin Spin1,
                  spin Spin2);
/// Make a fermionic density-density interaction term \f$ U n_{i\alpha_1\sigma_1}n_{j\alpha_2\sigma_2}\f$
/// with a complex interaction strength \f$U\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$ on site \f$i\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$ on site \f$j\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$ on site \f$i\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$ on site \f$j\f$.
ComplexExpr NupNdown(std::string const& Label1,
                     std::string const& Label2,
                     ComplexType U,
                     unsigned short Orbital1,
                     unsigned short Orbital2,
                     spin Spin1,
                     spin Spin2);

/// Make a fermionic density-density interaction term \f$ U n_{i\alpha_1\sigma_1}n_{i\alpha_2\sigma_2}\f$
/// with a real interaction strength \f$U\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$.
RealExpr NupNdown(std::string const& Label,
                  RealType U,
                  unsigned short Orbital1,
                  unsigned short Orbital2,
                  spin Spin1,
                  spin Spin2);
/// Make a fermionic density-density interaction term \f$ U n_{i\alpha_1\sigma_1}n_{i\alpha_2\sigma_2}\f$
/// with a complex interaction strength \f$U\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$.
ComplexExpr NupNdown(std::string const& Label,
                     ComplexType U,
                     unsigned short Orbital1,
                     unsigned short Orbital2,
                     spin Spin1,
                     spin Spin2);

/// Make a fermionic density-density interaction term \f$ U n_{i\alpha_1\uparrow}n_{i\alpha_2\downarrow}\f$
/// with a real interaction strength \f$U\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$.
RealExpr NupNdown(std::string const& Label, RealType U, unsigned short Orbital1, unsigned short Orbital2);
/// Make a fermionic density-density interaction term \f$ U n_{i\alpha_1\uparrow}n_{i\alpha_2\downarrow}\f$
/// with a complex interaction strength \f$U\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$.
ComplexExpr NupNdown(std::string const& Label, ComplexType U, unsigned short Orbital1, unsigned short Orbital2);

/// Make a fermionic density-density interaction term \f$ U n_{i\alpha\sigma_1}n_{i\alpha\sigma_2}\f$
/// with a real interaction strength \f$U\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$.
RealExpr NupNdown(std::string const& Label, RealType U, unsigned short Orbital, spin Spin1 = up, spin Spin2 = down);
/// Make a fermionic density-density interaction term \f$ U n_{i\alpha\sigma_1}n_{i\alpha\sigma_2}\f$
/// with a complex interaction strength \f$U\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$.
ComplexExpr
NupNdown(std::string const& Label, ComplexType U, unsigned short Orbital, spin Spin1 = up, spin Spin2 = down);

///@}

/// \defgroup SpinflipPairHopping Factory functions for spin-flip and pair-hopping terms
/// \ingroup LatticePresets
///@{

//
// Overloads of Spinflip()
//

/// Make a spin-flip term
/// \f$J c^\dagger_{i\alpha_1\sigma_1}c^\dagger_{i\alpha_2\sigma_2}c_{i\alpha_2\sigma_1}c_{i\alpha_1\sigma_2}\f$,
/// with \f$\alpha_1 \neq \alpha_2, \sigma_1 \neq \sigma_2\f$ and a real exchange constant \f$J\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$.
RealExpr Spinflip(std::string const& Label,
                  RealType J,
                  unsigned short Orbital1,
                  unsigned short Orbital2,
                  spin Spin1 = up,
                  spin Spin2 = down);
/// Make a spin-flip term
/// \f$J c^\dagger_{i\alpha_1\sigma_1}c^\dagger_{i\alpha_2\sigma_2}c_{i\alpha_2\sigma_1}c_{i\alpha_1\sigma_2}\f$,
/// with \f$\alpha_1 \neq \alpha_2, \sigma_1 \neq \sigma_2\f$ and a complex exchange constant \f$J\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$.
ComplexExpr Spinflip(std::string const& Label,
                     ComplexType J,
                     unsigned short Orbital1,
                     unsigned short Orbital2,
                     spin Spin1 = up,
                     spin Spin2 = down);

//
// Overloads of PairHopping()
//

/// Make a pair-hopping term
/// \f$J c^\dagger_{i\alpha_1\sigma_1}c^\dagger_{i\alpha_1\sigma_2}c_{i\alpha_2\sigma_1}c_{i\alpha_2\sigma_2}\f$,
/// with \f$\alpha_1 \neq \alpha_2, \sigma_1 \neq \sigma_2\f$ and a real exchange constant \f$J\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$.
RealExpr PairHopping(std::string const& Label,
                     RealType J,
                     unsigned short Orbital1,
                     unsigned short Orbital2,
                     spin Spin1 = up,
                     spin Spin2 = down);
/// Make a pair-hopping term
/// \f$J c^\dagger_{i\alpha_1\sigma_1}c^\dagger_{i\alpha_1\sigma_2}c_{i\alpha_2\sigma_1}c_{i\alpha_2\sigma_2}\f$,
/// with \f$\alpha_1 \neq \alpha_2, \sigma_1 \neq \sigma_2\f$ and a complex exchange constant \f$J\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] Orbital1 Orbital index \f$\alpha_1\f$.
/// \param[in] Orbital2 Orbital index \f$\alpha_2\f$.
/// \param[in] Spin1 Spin component \f$\sigma_1\f$.
/// \param[in] Spin2 Spin component \f$\sigma_2\f$.
ComplexExpr PairHopping(std::string const& Label,
                        ComplexType J,
                        unsigned short Orbital1,
                        unsigned short Orbital2,
                        spin Spin1 = up,
                        spin Spin2 = down);

///@}

/// \defgroup SS Factory functions for spin coupling terms
/// \ingroup LatticePresets
///@{

//
// Overloads of SplusSminus()
//

/// Make a fermionic \f$S_+ S_-\f$-coupling term
/// \f$J S_{+,i\alpha} S_{-,j\alpha}  = J c^\dagger_{i\alpha\uparrow} c_{i\alpha\downarrow}
/// c^\dagger_{j\alpha\downarrow} c_{j\alpha\uparrow}\f$ with a real exchange constant \f$J\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
RealExpr SplusSminus(std::string const& Label1, std::string const& Label2, RealType J, unsigned short Orbital);
/// Make a fermionic \f$S_+ S_-\f$-coupling term
/// \f$J S_{+,i\alpha} S_{-,j\alpha}  = J c^\dagger_{i\alpha\uparrow} c_{i\alpha\downarrow}
/// c^\dagger_{j\alpha\downarrow} c_{j\alpha\uparrow}\f$ with a complex exchange constant \f$J\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
ComplexExpr SplusSminus(std::string const& Label1, std::string const& Label2, ComplexType J, unsigned short Orbital);

//
// Overloads of SminusSplus()
//

/// Make a fermionic \f$S_- S_+\f$-coupling term
/// \f$J S_{-,i\alpha} S_{+,j\alpha} = J c^\dagger_{i\alpha\downarrow} c_{i\alpha\uparrow}
/// c^\dagger_{j\alpha\uparrow} c_{j\alpha\downarrow}\f$ with a real exchange constant \f$J\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
RealExpr SminusSplus(std::string const& Label1, std::string const& Label2, RealType J, unsigned short Orbital);
/// Make a fermionic \f$S_- S_+\f$-coupling term
/// \f$J S_{-,i\alpha} S_{+,j\alpha} = J c^\dagger_{i\alpha\downarrow} c_{i\alpha\uparrow}
/// c^\dagger_{j\alpha\uparrow} c_{j\alpha\downarrow}\f$ with a complex exchange constant \f$J\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
ComplexExpr SminusSplus(std::string const& Label1, std::string const& Label2, ComplexType J, unsigned short Orbital);

//
// Overloads of SzSz()
//

/// Make a fermionic \f$S_z S_z\f$-coupling term
/// \f[
///  J S_{z,i} S_{z,j} = \frac{J}{4} \sum_\alpha
///    (n_{i\alpha\uparrow} - n_{i\alpha\downarrow}) (n_{j\alpha\uparrow} - n_{j\alpha\downarrow})
/// \f]
/// with a real exchange constant \f$J\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
RealExpr SzSz(std::string const& Label1, std::string const& Label2, RealType J, unsigned short NOrbitals = 1);
/// Make a fermionic \f$S_z S_z\f$-coupling term
/// \f[
///  J S_{z,i} S_{z,j} = \frac{J}{4} \sum_\alpha
///    (n_{i\alpha\uparrow} - n_{i\alpha\downarrow}) (n_{j\alpha\uparrow} - n_{j\alpha\downarrow})
/// \f]
/// with a complex exchange constant \f$J\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
ComplexExpr SzSz(std::string const& Label1, std::string const& Label2, ComplexType J, unsigned short NOrbitals = 1);

//
// Overloads of SS()
//

/// Make a fermionic \f$\mathbf{S S}\f$-coupling term
/// \f[
///  J \mathbf{S}_{i} \mathbf{S}_{j} = J \sum_\alpha \left[ S_{z,i\alpha} S_{z,j\alpha} +
///    \frac{1}{2} S_{+,i\alpha} S_{-,j\alpha} + \frac{1}{2} S_{-,i\alpha} S_{+,j\alpha}\right]
/// \f]
/// with a real exchange constant \f$J\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
RealExpr SS(std::string const& Label1, std::string const& Label2, RealType J, unsigned short NOrbitals = 1);
/// Make a fermionic \f$\mathbf{S S}\f$-coupling term
/// \f[
///  J \mathbf{S}_{i} \mathbf{S}_{j} = J \sum_\alpha \left[ S_{z,i\alpha} S_{z,j\alpha} +
///    \frac{1}{2} S_{+,i\alpha} S_{-,j\alpha} + \frac{1}{2} S_{-,i\alpha} S_{+,j\alpha}\right]
/// \f]
/// with a complex exchange constant \f$J\f$.
/// \param[in] Label1 The first lattice site \f$i\f$.
/// \param[in] Label2 The second lattice site \f$j\f$.
/// \param[in] J Exchange constant \f$J\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
ComplexExpr SS(std::string const& Label1, std::string const& Label2, ComplexType J, unsigned short NOrbitals = 1);

///@}

/// \defgroup Coulomb Factory functions for Coulomb interaction terms
/// \ingroup LatticePresets
///@{

//
// Overloads of CoulombS()
//

/// Make a Coulomb interaction term of the following form,
/// \f[
///   U \sum_\alpha n_{i\alpha\uparrow} n_{i\alpha\downarrow} +
///   \varepsilon \sum_{\alpha,\sigma} n_{i\alpha\sigma}
/// \f]
/// with a real interaction constant \f$U\f$ and a real energy level \f$\varepsilon\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Eps Energy level \f$\varepsilon\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
RealExpr CoulombS(std::string const& Label, RealType U, RealType Eps, unsigned short NOrbitals = 1);
/// Make a Coulomb interaction term of the following form,
/// \f[
///   U \sum_\alpha n_{i\alpha\uparrow} n_{i\alpha\downarrow} +
///   \varepsilon \sum_{\alpha,\sigma} n_{i\alpha\sigma}
/// \f]
/// with a complex interaction constant \f$U\f$ and a complex energy level \f$\varepsilon\f$.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] Eps Energy level \f$\varepsilon\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha\f$ to sum over.
ComplexExpr CoulombS(std::string const& Label, ComplexType U, ComplexType Eps, unsigned short NOrbitals = 1);

//
// Overloads of CoulombP()
//

/// Make a Hubbard-Kanamori interaction term of the following form,
/// \f[
/// U \sum_{\alpha, \sigma > \sigma'} n_{i\alpha\sigma}n_{i\alpha\sigma'} +
/// U' \sum_{\alpha\neq\alpha',\sigma > \sigma'} n_{i\alpha\sigma} n_{i\alpha'\sigma'} +
/// \frac{U'-J}{2} \sum_{\alpha\neq\alpha',\sigma} n_{i\alpha\sigma} n_{i\alpha'\sigma} -
/// J \sum_{\alpha\neq\alpha',\sigma > \sigma'} (c^\dagger_{i\alpha \sigma}
///   c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'} +
///   c^\dagger_{i\alpha'\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha\sigma}c_{i\alpha\sigma'}) +
/// \varepsilon \sum_{\alpha,\sigma} n_{i\alpha\sigma}
/// \f]
/// with real interaction constants \f$U\f$, \f$U_p\f$, \f$J\f$ and a real energy level \f$\varepsilon\f$.
/// The number of orbitals must be at least 2.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Hubbard-Kanamori interaction constant \f$U\f$.
/// \param[in] U_p Hubbard-Kanamori interaction constant \f$U_p\f$.
/// \param[in] J Hund's coupling \f$J\f$.
/// \param[in] Eps Energy level \f$\varepsilon\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha, \alpha'\f$ to sum over.
RealExpr
CoulombP(std::string const& Label, RealType U, RealType U_p, RealType J, RealType Eps, unsigned short NOrbitals = 3);
/// Make a Hubbard-Kanamori interaction term of the following form,
/// \f[
/// U \sum_{\alpha, \sigma > \sigma'} n_{i\alpha\sigma}n_{i\alpha\sigma'} +
/// U' \sum_{\alpha\neq\alpha',\sigma > \sigma'} n_{i\alpha\sigma} n_{i\alpha'\sigma'} +
/// \frac{U'-J}{2} \sum_{\alpha\neq\alpha',\sigma} n_{i\alpha\sigma} n_{i\alpha'\sigma} -
/// J \sum_{\alpha\neq\alpha',\sigma > \sigma'} (c^\dagger_{i\alpha \sigma}
///   c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'} +
///   c^\dagger_{i\alpha'\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha\sigma}c_{i\alpha\sigma'}) +
/// \varepsilon \sum_{\alpha,\sigma} n_{i\alpha\sigma}
/// \f]
/// with complex interaction constants \f$U\f$, \f$U_p\f$, \f$J\f$ and a complex energy level \f$\varepsilon\f$.
/// The number of orbitals must be at least 2.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Hubbard-Kanamori interaction constant \f$U\f$.
/// \param[in] U_p Hubbard-Kanamori interaction constant \f$U_p\f$.
/// \param[in] J Hund's coupling \f$J\f$.
/// \param[in] Eps Energy level \f$\varepsilon\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha, \alpha'\f$ to sum over.
ComplexExpr CoulombP(std::string const& Label,
                     ComplexType U,
                     ComplexType U_p,
                     ComplexType J,
                     ComplexType Eps,
                     unsigned short NOrbitals = 3);

/// Make a Hubbard-Kanamori interaction term of the following form,
/// \f[
/// U \sum_{\alpha, \sigma > \sigma'} n_{i\alpha\sigma}n_{i\alpha\sigma'} +
/// (U-2J) \sum_{\alpha\neq\alpha',\sigma > \sigma'} n_{i\alpha\sigma} n_{i\alpha'\sigma'} +
/// \frac{U-3J}{2} \sum_{\alpha\neq\alpha',\sigma} n_{i\alpha\sigma} n_{i\alpha'\sigma} -
/// J \sum_{\alpha\neq\alpha',\sigma > \sigma'} (c^\dagger_{i\alpha \sigma}
///   c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'} +
///   c^\dagger_{i\alpha'\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha\sigma}c_{i\alpha\sigma'}) +
/// \varepsilon \sum_{\alpha,\sigma} n_{i\alpha\sigma}
/// \f]
/// with real interaction constants \f$U\f$, \f$J\f$ and a real energy level \f$\varepsilon\f$.
/// The number of orbitals must be at least 2.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Hubbard-Kanamori interaction constant \f$U\f$.
/// \param[in] J Hund's coupling \f$J\f$.
/// \param[in] Eps Energy level \f$\varepsilon\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha, \alpha'\f$ to sum over.
RealExpr CoulombP(std::string const& Label, RealType U, RealType J, RealType Eps, unsigned short NOrbitals = 3);
/// Make a Hubbard-Kanamori interaction term of the following form,
/// \f[
/// U \sum_{\alpha, \sigma > \sigma'} n_{i\alpha\sigma}n_{i\alpha\sigma'} +
/// (U-2J) \sum_{\alpha\neq\alpha',\sigma > \sigma'} n_{i\alpha\sigma} n_{i\alpha'\sigma'} +
/// \frac{U-3J}{2} \sum_{\alpha\neq\alpha',\sigma} n_{i\alpha\sigma} n_{i\alpha'\sigma} -
/// J \sum_{\alpha\neq\alpha',\sigma > \sigma'} (c^\dagger_{i\alpha \sigma}
///   c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'} +
///   c^\dagger_{i\alpha'\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha\sigma}c_{i\alpha\sigma'}) +
/// \varepsilon \sum_{\alpha,\sigma} n_{i\alpha\sigma}
/// \f]
/// with complex interaction constants \f$U\f$, \f$J\f$ and a complex energy level \f$\varepsilon\f$.
/// The number of orbitals must be at least 2.
/// \param[in] Label Lattice site \f$i\f$.
/// \param[in] U Hubbard-Kanamori interaction constant \f$U\f$.
/// \param[in] J Hund's coupling \f$J\f$.
/// \param[in] Eps Energy level \f$\varepsilon\f$.
/// \param[in] NOrbitals Number of orbitals \f$\alpha, \alpha'\f$ to sum over.
ComplexExpr
CoulombP(std::string const& Label, ComplexType U, ComplexType J, ComplexType Eps, unsigned short NOrbitals = 3);

/// @}

//
// Bosons
//

/// \defgroup Boson Factory functions for terms with bosonic degrees of freedom
/// \ingroup LatticePresets
///@{

/// Make a single energy level term \f$\varepsilon a^{\dagger}_{i\alpha} a_{i\alpha}\f$
/// for a boson on a given site with a given additional index.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] Eps Real energy level \f$\varepsilon\f$.
/// \param[in] ExtraIndex Additional index \f$\alpha\f$.
RealExpr BosonLevel(std::string const& Label, RealType Eps, unsigned short ExtraIndex);
/// Make a single energy level term \f$\varepsilon a^{\dagger}_{i\alpha} a_{i\alpha}\f$
/// for a boson on a given site with a given additional index.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] Eps Complex energy level \f$\varepsilon\f$.
/// \param[in] ExtraIndex Additional index \f$\alpha\f$.
ComplexExpr BosonLevel(std::string const& Label, ComplexType Eps, unsigned short ExtraIndex);

/// Make a bosonic interaction term
/// \f$\frac{U}{2}a^{\dagger}_{i\alpha} a_{i\alpha} (a^{\dagger}_{i\alpha} a_{i\alpha} - 1)\f$
/// with a real interaction strength \f$U\f$.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] ExtraIndex Additional index \f$\alpha\f$.
RealExpr BosonInteraction(std::string const& Label, RealType U, unsigned short ExtraIndex);
/// Make a bosonic interaction term
/// \f$\frac{U}{2}a^{\dagger}_{i\alpha} a_{i\alpha} (a^{\dagger}_{i\alpha} a_{i\alpha} - 1)\f$
/// with a complex interaction strength \f$U\f$.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] U Interaction constant \f$U\f$.
/// \param[in] ExtraIndex Additional index \f$\alpha\f$.
ComplexExpr BosonInteraction(std::string const& Label, ComplexType U, unsigned short ExtraIndex);

/// Make a Holstein fermion-boson coupling term of the following form,
/// \f[
/// \lambda (n_{i\alpha\uparrow} + n_{i\alpha\downarrow}) (a^\dagger_{i\beta} + a_{i\beta})
/// \f]
/// with a real coupling constant \f$\lambda\f$.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] Lambda Coupling constant \f$\lambda\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
/// \param[in] BosonExtraIndex Additional bosonic index \f$\beta\f$.
RealExpr
HolsteinInteraction(std::string const& Label, RealType Lambda, unsigned short Orbital, unsigned short BosonExtraIndex);
/// Make a Holstein fermion-boson coupling term of the following form,
/// \f[
/// \lambda (n_{i\alpha\uparrow} + n_{i\alpha\downarrow}) (a^\dagger_{i\beta} + a_{i\beta})
/// \f]
/// with a complex coupling constant \f$\lambda\f$.
/// \param[in] Label Site label \f$i\f$.
/// \param[in] Lambda Coupling constant \f$\lambda\f$.
/// \param[in] Orbital Orbital index \f$\alpha\f$.
/// \param[in] BosonExtraIndex Additional bosonic index \f$\beta\f$.
ComplexExpr HolsteinInteraction(std::string const& Label,
                                ComplexType Lambda,
                                unsigned short Orbital,
                                unsigned short BosonExtraIndex);

///@}

} // namespace LatticePresets
} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_LATTICEPRESETS_HPP
