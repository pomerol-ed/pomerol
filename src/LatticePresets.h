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

/** Some generic presets for spin-1/2 models. */
const Lattice::Site sSite (std::string(""), 1, 2); 
const Lattice::Site t2gSite (std::string(""), 3, 2); 
const Lattice::Site egSite (std::string(""), 2, 2); 

/** This is a set of presets of different Terms, most commonly used while writing a Hamiltonian. */
class Lattice::Term::Presets{
private:
public:
    /** Generates a hopping term \f$ t c^{\dagger}_{i\alpha\sigma}c_{j\alpha\sigma}, j \neq i \f$ between two sites. 
     * \param[in] Label1 \f$i\f$ - the first site which is connected by this term.
     * \param[in] Label2 \f$j\f$ - the second site which is connected by this term.
     * \param[in] Value \f$t\f$ - matrix element of a term.
     * \param[in] orbital \f$m\f$ - orbitals of sites, which are connected by this term.
     * \param[in] spin \f$\sigma\f$ - spins of sites, which are connected by this term.
     */
    static Term* Hopping     ( const std::string& Label1, const std::string& Label2, RealType Value, unsigned short orbital, unsigned short spin);

    /** Generates a single energy level term \f$\varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma} \f$ on a local site for a given spin and orbital. 
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$\varepsilon\f$ - the energy level. 
     * \param[in] orbital \f$m\f$ - affected orbital of the site.
     * \param[in] spin \f$\sigma\f$ - affected spin component.
     */
    static Term* Level       ( const std::string& Label, RealType Value, unsigned short orbital, unsigned short spin);

    /** Generates a local density-density 4-point term \f$ U n_{i\alpha\sigma}n_{i\alpha'\sigma'} \f$.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$U\f$ - matrix element of the term.
     * \param[in] orbital1 \f$m\f$ - the orbital affected by the first density operator.
     * \param[in] orbital2 \f$m'\f$ - the orbital affected by the second density operator.
     * \param[in] spin1 \f$\sigma\f$ - the spin component affected by the first density operator.
     * \param[in] spin2 \f$\sigma'\f$ - the spin component affected by the second density operator.
     */
    static Term* NupNdown    ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2);
    /** A shortcut to Pomerol::Lattice::Term::Presets::NupNdown \f$ U n_{i\alpha\uparrow}n_{i\alpha'\downarrow'} \f$ term for spin1 = \f$\uparrow\f$, spin2 = \f$\downarrow\f$. */
    static Term* NupNdown    ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2);
    /** A shortcut to Pomerol::Lattice::Term::Presets::NupNdown \f$ U n_{i\alpha\uparrow}n_{i\alpha\downarrow'} \f$ term for the same orbital \f$m=m'\f$ and default parameters spin1 = \f$\uparrow\f$, spin2 = \f$\downarrow\f$. */
    static Term* NupNdown    ( const std::string& Label, RealType Value, unsigned short orbital, unsigned short spin1 = up, unsigned short spin2 = down);


    /** Generates a spinflip \f$ J c^\dagger_{i\alpha\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'}, \alpha \neq \alpha', \sigma \neq \sigma' \f$ term.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$J\f$ - matrix element of the term.
     * \param[in] orbital \f$\alpha\f$ - first orbital affected by this term.
     * \param[in] orbital \f$\alpha'\f$ - second orbital affected by this term.
     * \param[in] spin1 \f$\sigma\f$ - first affected spin component. By default set to \f$\uparrow\f$.
     * \param[in] spin2 \f$\sigma'\f$ - second affected spin component. By default set to \f$\downarrow\f$.
     */
    static Term* Spinflip ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1 = up, unsigned short spin2 = down);


    /** Generates a pair-hopping \f$ J c^\dagger_{i\alpha\sigma}c^\dagger_{i\alpha\sigma'}c_{i\alpha'\sigma}c_{i\alpha'\sigma'}, \alpha \neq \alpha', \sigma \neq \sigma' \f$ term.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$J\f$ - matrix element of the term.
     * \param[in] orbital \f$\alpha\f$ - first orbital affected by this term.
     * \param[in] orbital \f$\alpha'\f$ - second orbital affected by this term.
     * \param[in] spin1 \f$\sigma\f$ - first affected spin component. By default set to \f$\uparrow\f$.
     * \param[in] spin2 \f$\sigma'\f$ - second affected spin component. By default set to \f$\downarrow\f$.
     */
    static Term* PairHopping (const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1 = up, unsigned short spin2 = down);
};

class Lattice::Presets {
public:
    /** Adds a site with a hamiltonian \f[ \sum\limits_{\alpha, \sigma > \sigma'} Un_{i\alpha\sigma}Un_{i\alpha\sigma'} + \sum\limits_{\alpha,\sigma} \varepsilon n_{i\alpha\sigma}. \f] 
     * \param[in] L A pointer to the Lattice to add the site.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] U \f$U\f$ - value of the onsite Coulomb interaction.
     * \param[in] Level \f$\varepsilon\f$ - the local energy level on the site.
     * \param[in] Orbitals Total amount of orbitals on the site. By default equal to 1.
     * \param[in] Spins Total amount of spin components on the site. By default equal to 2.
     */
    static void addSSite(Lattice *L, const std::string& label, RealType U, RealType Level, unsigned short Orbitals=1, unsigned short Spins=2);

    /** Adds a site with a hamiltonian \f[ U \sum_{\alpha, \sigma > \sigma'} n_{i\alpha\sigma}n_{i\alpha\sigma'} + U' \sum_{\alpha\neq\alpha',\sigma > \sigma'} n_{i\alpha\sigma} n_{i\alpha'\sigma'} + \frac{U'-J}{2} \sum_{\alpha\neq\alpha',\sigma} n_{i\alpha\sigma} n_{i\alpha'\sigma} - J \sum_{\alpha\neq\alpha',\sigma > \sigma'} (c^\dagger_{i\alpha \sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'} + c^\dagger_{i\alpha'\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha\sigma}c_{i\alpha\sigma'}). \f] 
     * \param[in] L A pointer to the Lattice to add the site.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] U \f$U\f$ - value of the onsite Coulomb interaction.
     * \param[in] U_p \f$U'\f$ - value of Kanamori parameter. 
     * \param[in] J \f$J\f$ - value of the Hund's coupling.
     * \param[in] Level \f$\varepsilon\f$ - the local energy level on the site.
     * \param[in] Orbitals Total amount of orbitals on the site. By default equal to 1.
     * \param[in] Spins Total amount of spin components on the site. By default equal to 2.
     */
    static void addPSite(Lattice *L, const std::string& label, RealType U, RealType U_p, RealType J, RealType Level, unsigned short Orbitals, unsigned short Spins);
    /** A shortcut to Lattice::Presets::addPSite with \f$U'=U-2J\f$, i.e. U_p = U - 2.0* J */
    static void addPSite(Lattice *L, const std::string& label, RealType U, RealType J, RealType Level, unsigned short Orbitals, unsigned short Spins);
    /** A shortcut to Lattice::Presets::addPSite with \f$U'=U-2J\f$, i.e. U_p = U - 2.0* J, and 2 spins */
    static void addPSite(Lattice *L, const std::string& label, RealType U, RealType J, RealType Level, unsigned Orbitals);

    /** Adds a magnetic \f$ \sum\limits_\alpha mH \frac{1}{2} (n_{i\alpha\uparrow} - n_{i\alpha\downarrow}) \f$ splitting to a given site.
     * \param[in] L A pointer to the Lattice to add the site.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] Magnetization \f$mH\f$ - magnetization to add.
     * \param[in] Orbitals Total amount of orbitals on the site. By default equal to 1.
     * \param[in] Spins Total amount of spin components on the site. By default equal to 2. Works only for 2 spins.
     */
    static void addMagnetization( Lattice *L, const std::string& label, RealType Magnetization, unsigned short Orbitals, unsigned short Spins=2);
};


}; // end of namespace Pomerol

#endif // endif : ifndef __INCLUDE_LATTICE_SITES_PRESETS_H
