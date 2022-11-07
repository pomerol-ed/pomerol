//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/LatticePresets.cpp
/// \brief Factory functions for terms commonly used to construct various lattice Hamiltonians (implementation).
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#include "pomerol/LatticePresets.hpp"

#include <stdexcept>

namespace Pomerol {
namespace LatticePresets {

using Operators::c_dag;
using Operators::c;
using Operators::n;
using Operators::a_dag;
using Operators::a;

//
// spin
//

std::ostream& operator<<(std::ostream& os, spin s) {
    switch(s) {
    case undef: return os;
    case up: return os << "up";
    case down: return os << "dn";
    default: return os;
    }
}

//
// Overloads of Level()
//

RealExpr Level(std::string const& Label, RealType Eps, unsigned short Orbital, spin Spin) {
    return Eps * n(Label, Orbital, Spin);
}
ComplexExpr Level(std::string const& Label, ComplexType Eps, unsigned short Orbital, spin Spin) {
    return Eps * n(Label, Orbital, Spin);
}

RealExpr Level(std::string const& Label, RealType Eps, unsigned short NOrbitals) {
    RealExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::Level(Label, Eps, Orbital, up);
        res += LatticePresets::Level(Label, Eps, Orbital, down);
    }
    return res;
}

ComplexExpr Level(std::string const& Label, ComplexType Eps, unsigned short NOrbitals) {
    ComplexExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::Level(Label, Eps, Orbital, up);
        res += LatticePresets::Level(Label, Eps, Orbital, down);
    }
    return res;
}

//
// Overloads of Hopping()
//

RealExpr Hopping(std::string const& Label1,
                 std::string const& Label2,
                 RealType t,
                 unsigned short Orbital1,
                 unsigned short Orbital2,
                 spin Spin1,
                 spin Spin2) {
    return (t * c_dag(Label1, Orbital1, Spin1) * c(Label2, Orbital2, Spin2)) + Operators::hc;
}
ComplexExpr Hopping(std::string const& Label1,
                    std::string const& Label2,
                    ComplexType t,
                    unsigned short Orbital1,
                    unsigned short Orbital2,
                    spin Spin1,
                    spin Spin2) {
    return (t * c_dag(Label1, Orbital1, Spin1) * c(Label2, Orbital2, Spin2)) + Operators::hc;
}

RealExpr Hopping(std::string const& Label1, std::string const& Label2, RealType t, unsigned short Orbital, spin Spin) {
    return Hopping(Label1, Label2, t, Orbital, Orbital, Spin, Spin);
}
ComplexExpr
Hopping(std::string const& Label1, std::string const& Label2, ComplexType t, unsigned short Orbital, spin Spin) {
    return Hopping(Label1, Label2, t, Orbital, Orbital, Spin, Spin);
}

RealExpr Hopping(std::string const& Label1,
                 std::string const& Label2,
                 RealType t,
                 unsigned short Orbital1,
                 unsigned short Orbital2) {
    return Hopping(Label1, Label2, t, Orbital1, Orbital2, up, up) +
           Hopping(Label1, Label2, t, Orbital1, Orbital2, down, down);
}
ComplexExpr Hopping(std::string const& Label1,
                    std::string const& Label2,
                    ComplexType t,
                    unsigned short Orbital1,
                    unsigned short Orbital2) {
    return Hopping(Label1, Label2, t, Orbital1, Orbital2, up, up) +
           Hopping(Label1, Label2, t, Orbital1, Orbital2, down, down);
}

RealExpr Hopping(std::string const& Label1, std::string const& Label2, RealType t, unsigned short NOrbitals) {
    RealExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::Hopping(Label1, Label2, t, Orbital, up);
        res += LatticePresets::Hopping(Label1, Label2, t, Orbital, down);
    }
    return res;
}
ComplexExpr Hopping(std::string const& Label1, std::string const& Label2, ComplexType t, unsigned short NOrbitals) {
    ComplexExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::Hopping(Label1, Label2, t, Orbital, up);
        res += LatticePresets::Hopping(Label1, Label2, t, Orbital, down);
    }
    return res;
}

//
// Overloads of Magnetization()
//

RealExpr Magnetization(std::string const& Label, RealType H, unsigned short NOrbitals) {
    RealExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::Level(Label, H, Orbital, up);
        res += LatticePresets::Level(Label, -H, Orbital, down);
    }
    return res;
}

ComplexExpr Magnetization(std::string const& Label, ComplexType H, unsigned short NOrbitals) {
    ComplexExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::Level(Label, H, Orbital, up);
        res += LatticePresets::Level(Label, -H, Orbital, down);
    }
    return res;
}

//
// Overloads of NupNdown()
//

RealExpr NupNdown(std::string const& Label1,
                  std::string const& Label2,
                  RealType U,
                  unsigned short Orbital1,
                  unsigned short Orbital2,
                  spin Spin1,
                  spin Spin2) {
    return U * n(Label1, Orbital1, Spin1) * n(Label2, Orbital2, Spin2);
}
ComplexExpr NupNdown(std::string const& Label1,
                     std::string const& Label2,
                     ComplexType U,
                     unsigned short Orbital1,
                     unsigned short Orbital2,
                     spin Spin1,
                     spin Spin2) {
    return U * n(Label1, Orbital1, Spin1) * n(Label2, Orbital2, Spin2);
}

RealExpr NupNdown(std::string const& Label,
                  RealType U,
                  unsigned short Orbital1,
                  unsigned short Orbital2,
                  spin Spin1,
                  spin Spin2) {
    return NupNdown(Label, Label, U, Orbital1, Orbital2, Spin1, Spin2);
}
ComplexExpr NupNdown(std::string const& Label,
                     ComplexType U,
                     unsigned short Orbital1,
                     unsigned short Orbital2,
                     spin Spin1,
                     spin Spin2) {
    return NupNdown(Label, Label, U, Orbital1, Orbital2, Spin1, Spin2);
}

RealExpr NupNdown(std::string const& Label, RealType U, unsigned short Orbital1, unsigned short Orbital2) {
    return NupNdown(Label, Label, U, Orbital1, Orbital2, up, down);
}
ComplexExpr NupNdown(std::string const& Label, ComplexType U, unsigned short Orbital1, unsigned short Orbital2) {
    return NupNdown(Label, Label, U, Orbital1, Orbital2, up, down);
}

RealExpr NupNdown(std::string const& Label, RealType U, unsigned short Orbital, spin Spin1, spin Spin2) {
    return NupNdown(Label, Label, U, Orbital, Orbital, Spin1, Spin2);
}
ComplexExpr NupNdown(std::string const& Label, ComplexType U, unsigned short Orbital, spin Spin1, spin Spin2) {
    return NupNdown(Label, Label, U, Orbital, Orbital, Spin1, Spin2);
}

//
// Overloads of Spinflip()
//

RealExpr Spinflip(std::string const& Label,
                  RealType J,
                  unsigned short Orbital1,
                  unsigned short Orbital2,
                  spin Spin1,
                  spin Spin2) {
    return J * c_dag(Label, Orbital1, Spin1) * c_dag(Label, Orbital2, Spin2) * c(Label, Orbital2, Spin1) *
           c(Label, Orbital1, Spin2);
}
ComplexExpr Spinflip(std::string const& Label,
                     ComplexType J,
                     unsigned short Orbital1,
                     unsigned short Orbital2,
                     spin Spin1,
                     spin Spin2) {
    return J * c_dag(Label, Orbital1, Spin1) * c_dag(Label, Orbital2, Spin2) * c(Label, Orbital2, Spin1) *
           c(Label, Orbital1, Spin2);
}

//
// Overloads of PairHopping()
//

RealExpr PairHopping(std::string const& Label,
                     RealType J,
                     unsigned short Orbital1,
                     unsigned short Orbital2,
                     spin Spin1,
                     spin Spin2) {
    return J * c_dag(Label, Orbital1, Spin1) * c_dag(Label, Orbital1, Spin2) * c(Label, Orbital2, Spin1) *
           c(Label, Orbital2, Spin2);
}
ComplexExpr PairHopping(std::string const& Label,
                        ComplexType J,
                        unsigned short Orbital1,
                        unsigned short Orbital2,
                        spin Spin1,
                        spin Spin2) {
    return J * c_dag(Label, Orbital1, Spin1) * c_dag(Label, Orbital1, Spin2) * c(Label, Orbital2, Spin1) *
           c(Label, Orbital2, Spin2);
}

//
// Overloads of SplusSminus()
//

RealExpr SplusSminus(std::string const& Label1, std::string const& Label2, RealType J, unsigned short Orbital) {
    return J * c_dag(Label1, Orbital, up) * c(Label1, Orbital, down) * c_dag(Label2, Orbital, down) *
           c(Label2, Orbital, up);
}
ComplexExpr SplusSminus(std::string const& Label1, std::string const& Label2, ComplexType J, unsigned short Orbital) {
    return J * c_dag(Label1, Orbital, up) * c(Label1, Orbital, down) * c_dag(Label2, Orbital, down) *
           c(Label2, Orbital, up);
}

//
// Overloads of SminusSplus()
//

RealExpr SminusSplus(std::string const& Label1, std::string const& Label2, RealType J, unsigned short Orbital) {
    return J * c_dag(Label1, Orbital, down) * c(Label1, Orbital, up) * c_dag(Label2, Orbital, up) *
           c(Label2, Orbital, down);
}
ComplexExpr SminusSplus(std::string const& Label1, std::string const& Label2, ComplexType J, unsigned short Orbital) {
    return J * c_dag(Label1, Orbital, down) * c(Label1, Orbital, up) * c_dag(Label2, Orbital, up) *
           c(Label2, Orbital, down);
}

//
// Overloads of SzSz()
//

RealExpr SzSz(std::string const& Label1, std::string const& Label2, RealType J, unsigned short NOrbitals) {
    RealExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        LatticePresets::NupNdown(Label1, Label2, -J / 4., Orbital, Orbital, up, down);
        LatticePresets::NupNdown(Label1, Label2, -J / 4., Orbital, Orbital, down, up);
        if(Label1 != Label2) {
            res += LatticePresets::NupNdown(Label1, Label2, J / 4., Orbital, Orbital, up, up);
            res += LatticePresets::NupNdown(Label1, Label2, J / 4., Orbital, Orbital, down, down);
        } else {
            res += LatticePresets::Level(Label1, J / 4., Orbital, up);
            res += LatticePresets::Level(Label1, J / 4., Orbital, down);
        }
    }
    return res;
}

ComplexExpr SzSz(std::string const& Label1, std::string const& Label2, ComplexType J, unsigned short NOrbitals) {
    ComplexExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        LatticePresets::NupNdown(Label1, Label2, -J / 4., Orbital, Orbital, up, down);
        LatticePresets::NupNdown(Label1, Label2, -J / 4., Orbital, Orbital, down, up);
        if(Label1 != Label2) {
            res += LatticePresets::NupNdown(Label1, Label2, J / 4., Orbital, Orbital, up, up);
            res += LatticePresets::NupNdown(Label1, Label2, J / 4., Orbital, Orbital, down, down);
        } else {
            res += LatticePresets::Level(Label1, J / 4., Orbital, up);
            res += LatticePresets::Level(Label1, J / 4., Orbital, down);
        }
    }
    return res;
}

//
// Overloads of SS()
//

RealExpr SS(std::string const& Label1, std::string const& Label2, RealType J, unsigned short NOrbitals) {
    RealExpr res = SzSz(Label1, Label2, J, NOrbitals);
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::SplusSminus(Label1, Label2, J / 2., Orbital);
        res += LatticePresets::SminusSplus(Label1, Label2, J / 2., Orbital);
    }
    return res;
}

ComplexExpr SS(std::string const& Label1, std::string const& Label2, ComplexType J, unsigned short NOrbitals) {
    ComplexExpr res = SzSz(Label1, Label2, J, NOrbitals);
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::SplusSminus(Label1, Label2, J / 2., Orbital);
        res += LatticePresets::SminusSplus(Label1, Label2, J / 2., Orbital);
    }
    return res;
}

//
// Overloads of CoulombS()
//

RealExpr CoulombS(std::string const& Label, RealType U, RealType Eps, unsigned short NOrbitals) {
    RealExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::Level(Label, Eps, Orbital, up) + LatticePresets::Level(Label, Eps, Orbital, down) +
               LatticePresets::NupNdown(Label, U, Orbital, Orbital, up, down);
    }
    return res;
}
ComplexExpr CoulombS(std::string const& Label, ComplexType U, ComplexType Eps, unsigned short NOrbitals) {
    ComplexExpr res;
    for(unsigned short Orbital = 0; Orbital < NOrbitals; ++Orbital) {
        res += LatticePresets::Level(Label, Eps, Orbital, up) + LatticePresets::Level(Label, Eps, Orbital, down) +
               LatticePresets::NupNdown(Label, U, Orbital, Orbital, up, down);
    }
    return res;
}

//
// Overloads of CoulombP()
//

RealExpr
CoulombP(std::string const& Label, RealType U, RealType U_p, RealType J, RealType Eps, unsigned short NOrbitals) {
    if(NOrbitals < 2)
        throw std::runtime_error("Cannot add multiorbital interaction to a site with 1 orbital");
    RealExpr res;
    for(unsigned short Orbital1 = 0; Orbital1 < NOrbitals; ++Orbital1) {
        for(spin s1 : {up, down}) {
            res += LatticePresets::Level(Label, Eps, Orbital1, s1);
            for(unsigned short Orbital2 = 0; Orbital2 < NOrbitals; ++Orbital2)
                if(Orbital1 != Orbital2)
                    res += LatticePresets::NupNdown(Label, (U_p - J) / 2.0, Orbital1, Orbital2, s1, s1);
            for(spin s2 : {up, down}) {
                if(s2 >= s1)
                    continue;
                res += LatticePresets::NupNdown(Label, U, Orbital1, Orbital1, s1, s2);
                for(unsigned short Orbital2 = 0; Orbital2 < NOrbitals; ++Orbital2) {
                    if(Orbital1 != Orbital2) {
                        res += LatticePresets::NupNdown(Label, U_p, Orbital1, Orbital2, s1, s2);
                        res += LatticePresets::Spinflip(Label, -J, Orbital1, Orbital2, s1, s2);
                        res += LatticePresets::PairHopping(Label, -J, Orbital1, Orbital2, s1, s2);
                    }
                }
            }
        }
    }
    return res;
}
ComplexExpr CoulombP(std::string const& Label,
                     ComplexType U,
                     ComplexType U_p,
                     ComplexType J,
                     ComplexType Eps,
                     unsigned short NOrbitals) {
    if(NOrbitals < 2)
        throw std::runtime_error("Cannot add multiorbital interaction to a site with 1 orbital");
    ComplexExpr res;
    for(unsigned short Orbital1 = 0; Orbital1 < NOrbitals; ++Orbital1) {
        for(spin s1 : {up, down}) {
            res += LatticePresets::Level(Label, Eps, Orbital1, s1);
            for(unsigned short Orbital2 = 0; Orbital2 < NOrbitals; ++Orbital2)
                if(Orbital1 != Orbital2)
                    res += LatticePresets::NupNdown(Label, (U_p - J) / 2.0, Orbital1, Orbital2, s1, s1);
            for(spin s2 : {up, down}) {
                if(s2 >= s1)
                    continue;
                res += LatticePresets::NupNdown(Label, U, Orbital1, Orbital1, s1, s2);
                for(unsigned short Orbital2 = 0; Orbital2 < NOrbitals; ++Orbital2) {
                    if(Orbital1 != Orbital2) {
                        res += LatticePresets::NupNdown(Label, U_p, Orbital1, Orbital2, s1, s2);
                        res += LatticePresets::Spinflip(Label, -J, Orbital1, Orbital2, s1, s2);
                        res += LatticePresets::PairHopping(Label, -J, Orbital1, Orbital2, s1, s2);
                    }
                }
            }
        }
    }
    return res;
}

RealExpr CoulombP(std::string const& Label, RealType U, RealType J, RealType Eps, unsigned short NOrbitals) {
    return CoulombP(Label, U, U - 2.0 * J, J, Eps, NOrbitals);
}
ComplexExpr
addCoulombP(std::string const& Label, ComplexType U, ComplexType J, ComplexType Eps, unsigned short NOrbitals) {
    return CoulombP(Label, U, U - 2.0 * J, J, Eps, NOrbitals);
}

//
// Bosons
//

RealExpr BosonLevel(std::string const& Label, RealType Eps, unsigned short ExtraIndex) {
    return Eps * a_dag(Label, ExtraIndex, undef) * a(Label, ExtraIndex, undef);
}

ComplexExpr BosonLevel(std::string const& Label, ComplexType Eps, unsigned short ExtraIndex) {
    return Eps * a_dag(Label, ExtraIndex, undef) * a(Label, ExtraIndex, undef);
}

RealExpr BosonInteraction(std::string const& Label, RealType U, unsigned short ExtraIndex) {
    auto nb = a_dag(Label, ExtraIndex, undef) * a(Label, ExtraIndex, undef);
    return 0.5 * U * nb * (nb - 1.0);
}

ComplexExpr BosonInteraction(std::string const& Label, ComplexType U, unsigned short ExtraIndex) {
    auto nb = a_dag(Label, ExtraIndex, undef) * a(Label, ExtraIndex, undef);
    return 0.5 * U * nb * (nb - 1.0);
}

RealExpr
HolsteinInteraction(std::string const& Label, RealType Lambda, unsigned short Orbital, unsigned short BosonExtraIndex) {
    auto N = n(Label, Orbital, up) + n(Label, Orbital, down);
    return Lambda * N * (a_dag(Label, BosonExtraIndex, undef) + a(Label, BosonExtraIndex, undef));
}

ComplexExpr HolsteinInteraction(std::string const& Label,
                                ComplexType Lambda,
                                unsigned short Orbital,
                                unsigned short BosonExtraIndex) {
    auto N = n(Label, Orbital, up) + n(Label, Orbital, down);
    return Lambda * N * (a_dag(Label, BosonExtraIndex, undef) + a(Label, BosonExtraIndex, undef));
}

} // namespace LatticePresets
} // namespace Pomerol
