//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.

/** \file src/Lattice.h
** \brief A lattice handler. Reads and stores the lattice from a JSON file. Designed to be extended to other formats.
** 
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_LATTICE_SITES_PRESETS_H
#define __INCLUDE_LATTICE_SITES_PRESETS_H

#include "Misc.h"
#include "Logger.h"
#include "Lattice.h"

namespace Pomerol{

/** Some generic presets for spin-1/2 models */
const Lattice::Site sSite (std::string(""), 1, 2); 
const Lattice::Site pSite (std::string(""), 3, 2); 

class Lattice::Term::Presets{
private:
public:
    static Term* Hopping     ( const std::string& Label1, const std::string& Label2, RealType Value, unsigned short orbital, unsigned short spin);
    static Term* Level       ( const std::string& Label, RealType Value, unsigned short orbital, unsigned short spin);

    static Term* NupNdown    ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2);
    static Term* NupNdown    ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2);
    static Term* NupNdown    ( const std::string& Label, RealType Value, unsigned short orbital, unsigned short spin1 = up, unsigned short spin2 = down);

    static Term* Spinflip ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1 = up, unsigned short spin2 = down);
    static Term* PairHopping (const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1 = up, unsigned short spin2 = down);
};

class Lattice::Presets {
public:
    static void addSSite(Lattice *L, const std::string& label, RealType U, RealType Level, unsigned short Orbitals=1, unsigned short Spins=2);

    static void addPSite(Lattice *L, const std::string& label, RealType U, RealType U_p, RealType J, RealType Level, unsigned short Orbitals, unsigned short Spins);
    static void addPSite(Lattice *L, const std::string& label, RealType U, RealType J, RealType Level, unsigned short Orbitals, unsigned short Spins);
    static void addPSite(Lattice *L, const std::string& label, RealType U, RealType J, RealType Level, unsigned Orbitals);

    static void addMagnetization( Lattice *L, const std::string& label, RealType Magnetization, unsigned short Orbitals, unsigned short Spins=2);
};


}; // end of namespace Pomerol

#endif // endif : ifndef __INCLUDE_LATTICE_SITES_PRESETS_H
