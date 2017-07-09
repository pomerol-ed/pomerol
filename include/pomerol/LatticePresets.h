/** \file src/LatticePresets.h
** \brief A set of preset methods to simplify Pomerol::Lattice entering. 
** 
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_LATTICE_SITES_PRESETS_H
#define __INCLUDE_LATTICE_SITES_PRESETS_H

#include "Misc.h"
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
    /** Generates a hopping term \f$ t c^{\dagger}_{i\alpha\sigma}c_{j\alpha'\sigma}, j \neq i \f$ between two sites. 
     * \param[in] Label1 \f$i\f$ - the first site which is connected by this term.
     * \param[in] Label2 \f$j\f$ - the second site which is connected by this term.
     * \param[in] Value \f$t\f$ - matrix element of a term.
     * \param[in] orbital \f$\alpha\f$ - orbital of site \f$i\f$, which is connected by this term.
     * \param[in] orbital \f$\alpha'\f$ - orbital of site \f$j\f$, which is connected by this term.
     * \param[in] spin \f$\sigma\f$ - spins of sites, which are connected by this term.
     */
    static Term* Hopping     ( const std::string& Label1, const std::string& Label2, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1 , unsigned short spin2);

    /** A shortcut to hopping Lattice::Term \f$ t c^{\dagger}_{i\alpha\sigma}c_{j\alpha\sigma}, j \neq i \f$ */
    static Term* Hopping     ( const std::string& Label1, const std::string& Label2, MelemType Value, unsigned short orbital, unsigned short spin);

    /** Generates a single energy level term \f$\varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma} \f$ on a local site for a given spin and orbital. 
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$\varepsilon\f$ - the energy level. 
     * \param[in] orbital \f$\alpha\f$ - affected orbital of the site.
     * \param[in] spin \f$\sigma\f$ - affected spin component.
     */
    static Term* Level       ( const std::string& Label, MelemType Value, unsigned short orbital, unsigned short spin);

    /** Generates a local density-density 4-point term \f$ U n_{i\alpha\sigma}n_{j\alpha'\sigma'} \f$.
     * \param[in] Label1 \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Label2 \f$j\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$U\f$ - matrix element of the term.
     * \param[in] orbital1 \f$\alpha\f$ - the orbital affected by the first density operator.
     * \param[in] orbital2 \f$\alpha'\f$ - the orbital affected by the second density operator.
     * \param[in] spin1 \f$\sigma\f$ - the spin component affected by the first density operator.
     * \param[in] spin2 \f$\sigma'\f$ - the spin component affected by the second density operator.
     */
    static Term* NupNdown    ( const std::string& Label1, const std::string& Label2, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2);

    /** A shortcut to Pomerol::Lattice::Term::Presets::NupNdown \f$ U n_{i\alpha\uparrow}n_{j\alpha'\downarrow'} \f$ term for \f$i=j\f$. */
    static Term* NupNdown    ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2);
    /** A shortcut to Pomerol::Lattice::Term::Presets::NupNdown \f$ U n_{i\alpha\uparrow}n_{i\alpha'\downarrow'} \f$ term for spin1 = \f$\uparrow\f$, spin2 = \f$\downarrow\f$. */
    static Term* NupNdown    ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2);
    /** A shortcut to Pomerol::Lattice::Term::Presets::NupNdown \f$ U n_{i\alpha\uparrow}n_{i\alpha\downarrow'} \f$ term for the same orbital \f$m=m'\f$ and default parameters spin1 = \f$\uparrow\f$, spin2 = \f$\downarrow\f$. */
    static Term* NupNdown    ( const std::string& Label, MelemType Value, unsigned short orbital, unsigned short spin1 = up, unsigned short spin2 = down);


    /** Generates a spinflip \f$ J c^\dagger_{i\alpha\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'}, \alpha \neq \alpha', \sigma \neq \sigma' \f$ term.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$J\f$ - matrix element of the term.
     * \param[in] orbital \f$\alpha\f$ - first orbital affected by this term.
     * \param[in] orbital \f$\alpha'\f$ - second orbital affected by this term.
     * \param[in] spin1 \f$\sigma\f$ - first affected spin component. By default set to \f$\uparrow\f$.
     * \param[in] spin2 \f$\sigma'\f$ - second affected spin component. By default set to \f$\downarrow\f$.
     */
    static Term* Spinflip ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1 = up, unsigned short spin2 = down);


    /** Generates a pair-hopping \f$ J c^\dagger_{i\alpha\sigma}c^\dagger_{i\alpha\sigma'}c_{i\alpha'\sigma}c_{i\alpha'\sigma'}, \alpha \neq \alpha', \sigma \neq \sigma' \f$ term.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$J\f$ - matrix element of the term.
     * \param[in] orbital \f$\alpha\f$ - first orbital affected by this term.
     * \param[in] orbital \f$\alpha'\f$ - second orbital affected by this term.
     * \param[in] spin1 \f$\sigma\f$ - first affected spin component. By default set to \f$\uparrow\f$.
     * \param[in] spin2 \f$\sigma'\f$ - second affected spin component. By default set to \f$\downarrow\f$.
     */
    static Term* PairHopping (const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1 = up, unsigned short spin2 = down);


    static Term* SplusSminus ( const std::string& label1, const std::string& label2, MelemType Value, unsigned short orbital); 
    static Term* SminusSplus ( const std::string& label1, const std::string& label2, MelemType Value, unsigned short orbital);

    /** Exception: wrong indices. */
    class exWrongIndices : public std::exception { 
    public:
        virtual const char* what() const throw();
    };
};

/** A set of presets to fill the Lattice::TermStorage and Lattice::Sites for some commonly used examples. */
class LatticePresets {
public:

    /** Adds an interaction with the hamiltonian \f[ \sum\limits_{\alpha, \sigma > \sigma'} Un_{i\alpha\sigma}Un_{i\alpha\sigma'} + \sum\limits_{\alpha,\sigma} \varepsilon n_{i\alpha\sigma} to a specified site. \f] 
     * \param[in] L A pointer to the Lattice to add the site.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] U \f$U\f$ - value of the onsite Coulomb interaction.
     * \param[in] Level \f$\varepsilon\f$ - the local energy level on the site.
     */
    static void addCoulombS(Lattice *L, const std::string& label, MelemType U, MelemType Level);

    /** Adds an interaction with the hamiltonian \f[ U \sum_{\alpha, \sigma > \sigma'} n_{i\alpha\sigma}n_{i\alpha\sigma'} + U' \sum_{\alpha\neq\alpha',\sigma > \sigma'} n_{i\alpha\sigma} n_{i\alpha'\sigma'} + \frac{U'-J}{2} \sum_{\alpha\neq\alpha',\sigma} n_{i\alpha\sigma} n_{i\alpha'\sigma} - J \sum_{\alpha\neq\alpha',\sigma > \sigma'} (c^\dagger_{i\alpha \sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'} + c^\dagger_{i\alpha'\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha\sigma}c_{i\alpha\sigma'}) to the specified site. \f] 
     * \param[in] L A pointer to the Lattice to add the site.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] U \f$U\f$ - Kanamori \f$U\f$,  value of the onsite Coulomb interaction.
     * \param[in] U_p \f$U'\f$ - Kanamori \f$U'\f$. 
     * \param[in] J \f$J\f$ - Kanamori J, value of the Hund's coupling.
     * \param[in] Level \f$\varepsilon\f$ - the local energy level on the site.
     */
    static void addCoulombP(Lattice *L, const std::string& label, MelemType U, MelemType U_p, MelemType J, MelemType Level);
    /** A shortcut to Lattice::Presets::addPSite with \f$U'=U-2J\f$, i.e. U_p = U - 2.0* J */
    static void addCoulombP(Lattice *L, const std::string& label, MelemType U, MelemType J, MelemType Level);

    /** Adds a magnetic \f$ \sum\limits_\alpha mH \frac{1}{2} (n_{i\alpha\uparrow} - n_{i\alpha\downarrow}) \f$ splitting to a given site. Valid only for 2 spins.
     * \param[in] L A pointer to the Lattice to add the terms. 
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] Magnetization \f$mH\f$ - magnetization to add.
     */
    static void addMagnetization( Lattice *L, const std::string& label, MelemType Magnetization);

    /** Adds a level \f$ \sum\limits_{\alpha, \sigma} \varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma} \f$.
     * \param[in] L A pointer to the Lattice to add the terms.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] Level \f$\varepsilon\f$ - energy level to add.
     */
    static void addLevel ( Lattice *L, const std::string& label, MelemType Level);

    /** Adds a SzSz \f[ \sum\limits_{\alpha} J \frac{1}{2}(n_{i\alpha\uparrow} - n_{i\alpha\downarrow})\frac{1}{2}(n_{j\alpha\uparrow} - n_{j\alpha\downarrow}) \f]
     * interaction terms. Valid only for 2 spins 
     * \param[in] L A pointer to the Lattice to add the terms. 
     * \param[in] Label1 \f$i\f$ - label of the first connected site.
     * \param[in] Label2 \f$j\f$ - label of the second connected site. Site can be choosen the same as the first site.
     * \param[in] ExchJ \f$J\f$ - magnetic exchange constant.
     * \param[in] Orbitals Total amount of orbitals on the site. By default equal to 1.
     * \param[in] Spins Total amount of spin components on the site. By default equal to 2. Works only for 2 spins.
     */
    static void addSzSz ( Lattice *L, const std::string& Label1, const std::string& Label2, MelemType ExchJ);

    /** Adds a spin-spin \f[ \sum\limits_{\alpha} J \hat S_{i\alpha} \hat S_{j\alpha} \f]
     * interaction terms. Valid only for 2 spins 
     * \param[in] L A pointer to the Lattice to add the terms. 
     * \param[in] Label1 \f$i\f$ - label of the first connected site.
     * \param[in] Label2 \f$j\f$ - label of the second connected site. Site can be choosen the same as the first site.
     * \param[in] ExchJ \f$J\f$ - magnetic exchange constant.
     */
    static void addSS ( Lattice *L, const std::string& Label1, const std::string& Label2, MelemType ExchJ);

    /** Adds a hopping \f$ t c^{\dagger}_{i\alpha\sigma}c_{j\alpha'\sigma'} \f$ term to the Lattice. This is a safe method : indices are checked to belong to the lattice.
     * \param[in] L A pointer to the Lattice to add the terms. 
     * \param[in] Label1 \f$i\f$ - label of the first connected site.
     * \param[in] Label2 \f$j\f$ - label of the second connected site. Site can be choosen the same as the first site.
     * \param[in] t \f$t\f$ - hopping matrix element.
     * \param[in] Orbital1 \f$\alpha\f$ - orbital of the first site.
     * \param[in] Orbital2 \f$\alpha\f$ - orbital of the second site.
     * \param[in] Spin1 \f$\sigma\f$ - spin component of the first site
     * \param[in] Spin2 \f$\sigma\f$ - spin component of the second site
     */
    static void addHopping ( Lattice *L, const std::string &label1, const std::string& label2, MelemType t, unsigned short Orbital1, unsigned short Orbital2, unsigned short spin1, unsigned short spin2 );
    /** A shortcut to addHopping \f$ t c^{\dagger}_{i\alpha\sigma}c_{j\alpha\sigma} \f$ */
    static void addHopping ( Lattice *L, const std::string &label1, const std::string& label2, MelemType t, unsigned short Orbital1, unsigned short Orbital2, unsigned short spin );
    /** A shortcut to addHopping \f$ \sum_{\sigma} t c^{\dagger}_{i\alpha\sigma}c_{j\alpha\sigma} \f$ */
    static void addHopping ( Lattice *L, const std::string &label1, const std::string& label2, MelemType t, unsigned short Orbital1, unsigned short Orbital2 );
    /** A shortcut to addHopping \f$ \sum_{\sigma\alpha} t c^{\dagger}_{i\alpha\sigma}c_{j\alpha\sigma} \f$ */
    static void addHopping ( Lattice *L, const std::string &label1, const std::string& label2, MelemType t );
};

}; // end of namespace Pomerol

#endif // endif : ifndef __INCLUDE_LATTICE_SITES_PRESETS_H
