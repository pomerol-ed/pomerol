/** \file include/pomerol/LatticePresets.h
** \brief A set of preset methods to simplify Pomerol::Lattice entering.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_LATTICEPRESETS_H
#define POMEROL_INCLUDE_POMEROL_LATTICEPRESETS_H

#include "Misc.hpp"
#include "Operators.hpp"

#include <ostream>
#include <string>

namespace Pomerol {

/** This is a set of presets of different Terms, most commonly used while writing a Hamiltonian. */
namespace LatticePresets {

    /** Possible spin projections are \b down and \b up */
    enum spin : short {undef = -1, down = 0, up = 1};

    std::ostream & operator<<(std::ostream & os, spin s);

    using RealExpr = Operators::expression<RealType, std::string, unsigned short, spin>;
    using ComplexExpr = Operators::expression<ComplexType, std::string, unsigned short, spin>;

    /** Generates a single energy level term \f$\varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma} \f$ on a local site for a given spin and orbital.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$\varepsilon\f$ - the energy level.
     * \param[in] orbital \f$\alpha\f$ - affected orbital of the site.
     * \param[in] spin \f$\sigma\f$ - affected spin component.
     */
    RealExpr Level(std::string const& Label, RealType Value, unsigned short Orbital, spin Spin);
    ComplexExpr Level(std::string const& Label, ComplexType Value, unsigned short Orbital, spin Spin);

    /** Adds a level \f$ \sum\limits_{\alpha, \sigma} \varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma} \f$.
     * \param[in] L A pointer to the Lattice to add the terms.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] Level \f$\varepsilon\f$ - energy level to add.
     */
    RealExpr Level(std::string const& Label, RealType Value, unsigned short NOrbitals = 1);
    ComplexExpr Level(std::string const& Label, ComplexType Value, unsigned short NOrbitals = 1);

    /** Generates a hopping term \f$ t c^{\dagger}_{i\alpha\sigma}c_{j\alpha'\sigma}, j \neq i \f$ between two sites.
     * \param[in] Label1 \f$i\f$ - the first site which is connected by this term.
     * \param[in] Label2 \f$j\f$ - the second site which is connected by this term.
     * \param[in] Value \f$t\f$ - matrix element of a term.
     * \param[in] orbital \f$\alpha\f$ - orbital of site \f$i\f$, which is connected by this term.
     * \param[in] orbital \f$\alpha'\f$ - orbital of site \f$j\f$, which is connected by this term.
     * \param[in] spin \f$\sigma\f$ - spins of sites, which are connected by this term.
     */
    RealExpr Hopping(std::string const& Label1, std::string const& Label2, RealType t, unsigned short Orbital1, unsigned short Orbital2, spin Spin1, spin Spin2);
    ComplexExpr Hopping(std::string const& Label1, std::string const& Label2, ComplexType t, unsigned short Orbital1, unsigned short Orbital2, spin Spin1, spin Spin2);
    /** A shortcut to hopping Lattice::Term \f$ t c^{\dagger}_{i\alpha\sigma}c_{j\alpha\sigma}, j \neq i \f$ */
    RealExpr Hopping(std::string const& Label1, std::string const& Label2, RealType t, unsigned short Orbital, spin Spin);
    ComplexExpr Hopping(std::string const& Label1, std::string const& Label2, ComplexType t, unsigned short Orbital, spin Spin);
    /** A shortcut to Hopping \f$ \sum_{\sigma} t c^{\dagger}_{i\alpha\sigma}c_{j\alpha\sigma} \f$ */
    RealExpr Hopping(std::string const& Label1, std::string const& Label2, RealType t, unsigned short Orbital1, unsigned short Orbital2);
    ComplexExpr Hopping(std::string const& Label1, std::string const& Label2, ComplexType t, unsigned short Orbital1, unsigned short Orbital2);
    /** A shortcut to Hopping \f$ \sum_{\sigma\alpha} t c^{\dagger}_{i\alpha\sigma}c_{j\alpha\sigma} \f$ */
    RealExpr Hopping(std::string const& Label1, std::string const& Label2, RealType t, unsigned short NOrbitals = 1);
    ComplexExpr Hopping(std::string const& Label1, std::string const& Label2, ComplexType t, unsigned short NOrbitals = 1);

    /** Generates a local density-density 4-point term \f$ U n_{i\alpha\sigma}n_{j\alpha'\sigma'} \f$.
     * \param[in] Label1 \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Label2 \f$j\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$U\f$ - matrix element of the term.
     * \param[in] orbital1 \f$\alpha\f$ - the orbital affected by the first density operator.
     * \param[in] orbital2 \f$\alpha'\f$ - the orbital affected by the second density operator.
     * \param[in] spin1 \f$\sigma\f$ - the spin component affected by the first density operator.
     * \param[in] spin2 \f$\sigma'\f$ - the spin component affected by the second density operator.
     */
    RealExpr NupNdown(std::string const& Label1, std::string const& Label2, RealType Value, unsigned short Orbital1, unsigned short Orbital2, spin Spin1, spin Spin2);
    ComplexExpr NupNdown(std::string const& Label1, std::string const& Label2, ComplexType Value, unsigned short Orbital1, unsigned short Orbital2, spin Spin1, spin Spin2);
    /** A shortcut to Pomerol::Lattice::Term::Presets::NupNdown \f$ U n_{i\alpha\uparrow}n_{j\alpha'\downarrow'} \f$ term for \f$i=j\f$. */
    RealExpr NupNdown(std::string const& Label, RealType Value, unsigned short Orbital1, unsigned short Orbital2, spin Spin1, spin Spin2);
    ComplexExpr NupNdown(std::string const& Label, ComplexType Value, unsigned short Orbital1, unsigned short Orbital2, spin Spin1, spin Spin2);
    /** A shortcut to Pomerol::Lattice::Term::Presets::NupNdown \f$ U n_{i\alpha\uparrow}n_{i\alpha'\downarrow'} \f$ term for spin1 = \f$\uparrow\f$, spin2 = \f$\downarrow\f$. */
    RealExpr NupNdown(std::string const& Label, RealType Value, unsigned short Orbital1, unsigned short Orbital2);
    ComplexExpr NupNdown(std::string const& Label, ComplexType Value, unsigned short Orbital1, unsigned short Orbital2);
    /** A shortcut to Pomerol::Lattice::Term::Presets::NupNdown \f$ U n_{i\alpha\uparrow}n_{i\alpha\downarrow'} \f$ term for the same orbital \f$m=m'\f$ and default parameters spin1 = \f$\uparrow\f$, spin2 = \f$\downarrow\f$. */
    RealExpr NupNdown(std::string const& Label, RealType Value, unsigned short Orbital, spin Spin1 = up, spin Spin2 = down);
    ComplexExpr NupNdown(std::string const& Label, ComplexType Value, unsigned short Orbital, spin Spin1 = up, spin Spin2 = down);

    /** Generates a spinflip \f$ J c^\dagger_{i\alpha\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'}, \alpha \neq \alpha', \sigma \neq \sigma' \f$ term.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$J\f$ - matrix element of the term.
     * \param[in] orbital \f$\alpha\f$ - first orbital affected by this term.
     * \param[in] orbital \f$\alpha'\f$ - second orbital affected by this term.
     * \param[in] spin1 \f$\sigma\f$ - first affected spin component. By default set to \f$\uparrow\f$.
     * \param[in] spin2 \f$\sigma'\f$ - second affected spin component. By default set to \f$\downarrow\f$.
     */
    RealExpr Spinflip(std::string const& Label, RealType Value, unsigned short Orbital1, unsigned short Orbital2, spin Spin1 = up, spin Spin2 = down);
    ComplexExpr Spinflip(std::string const& Label, ComplexType Value, unsigned short Orbital1, unsigned short Orbital2, spin Spin1 = up, spin Spin2 = down);

    /** Generates a pair-hopping \f$ J c^\dagger_{i\alpha\sigma}c^\dagger_{i\alpha\sigma'}c_{i\alpha'\sigma}c_{i\alpha'\sigma'}, \alpha \neq \alpha', \sigma \neq \sigma' \f$ term.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$J\f$ - matrix element of the term.
     * \param[in] orbital \f$\alpha\f$ - first orbital affected by this term.
     * \param[in] orbital \f$\alpha'\f$ - second orbital affected by this term.
     * \param[in] spin1 \f$\sigma\f$ - first affected spin component. By default set to \f$\uparrow\f$.
     * \param[in] spin2 \f$\sigma'\f$ - second affected spin component. By default set to \f$\downarrow\f$.
     */
    RealExpr PairHopping(std::string const& Label, RealType Value, unsigned short Orbital1, unsigned short Orbital2, spin Spin1 = up, spin Spin2 = down);
    ComplexExpr PairHopping(std::string const& Label, ComplexType Value, unsigned short Orbital1, unsigned short Orbital2, spin Spin1 = up, spin Spin2 = down);

    RealExpr SplusSminus(std::string const& Label1, std::string const& Label2, RealType Value, unsigned short Orbital);
    ComplexExpr SplusSminus(std::string const& Label1, std::string const& Label2, ComplexType Value, unsigned short Orbital);

    RealExpr SminusSplus(std::string const& label1, std::string const& Label2, RealType Value, unsigned short Orbital);
    ComplexExpr SminusSplus(std::string const& label1, std::string const& Label2, ComplexType Value, unsigned short Orbital);

    /** Adds an interaction with the hamiltonian \f[ \sum\limits_{\alpha, \sigma > \sigma'} Un_{i\alpha\sigma}Un_{i\alpha\sigma'} + \sum\limits_{\alpha,\sigma} \varepsilon n_{i\alpha\sigma} to a specified site. \f]
     * \param[in] L A pointer to the Lattice to add the site.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] U \f$U\f$ - value of the onsite Coulomb interaction.
     * \param[in] Level \f$\varepsilon\f$ - the local energy level on the site.
     */
    RealExpr CoulombS(std::string const& Label, RealType U, RealType Level, unsigned short NOrbitals = 1);
    ComplexExpr CoulombS(std::string const& Label, ComplexType U, ComplexType Level, unsigned short NOrbitals = 1);

    /** Adds an interaction with the hamiltonian \f[ U \sum_{\alpha, \sigma > \sigma'} n_{i\alpha\sigma}n_{i\alpha\sigma'} + U' \sum_{\alpha\neq\alpha',\sigma > \sigma'} n_{i\alpha\sigma} n_{i\alpha'\sigma'} + \frac{U'-J}{2} \sum_{\alpha\neq\alpha',\sigma} n_{i\alpha\sigma} n_{i\alpha'\sigma} - J \sum_{\alpha\neq\alpha',\sigma > \sigma'} (c^\dagger_{i\alpha \sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha'\sigma}c_{i\alpha\sigma'} + c^\dagger_{i\alpha'\sigma}c^\dagger_{i\alpha'\sigma'}c_{i\alpha\sigma}c_{i\alpha\sigma'}) to the specified site. \f]
     * \param[in] L A pointer to the Lattice to add the site.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] U \f$U\f$ - Kanamori \f$U\f$,  value of the onsite Coulomb interaction.
     * \param[in] U_p \f$U'\f$ - Kanamori \f$U'\f$.
     * \param[in] J \f$J\f$ - Kanamori J, value of the Hund's coupling.
     * \param[in] Level \f$\varepsilon\f$ - the local energy level on the site.
     */
    RealExpr CoulombP(std::string const& Label, RealType U, RealType U_p, RealType J, RealType Level, unsigned short NOrbitals = 3);
    ComplexExpr CoulombP(std::string const& Label, ComplexType U, ComplexType U_p, ComplexType J, ComplexType Level, unsigned short NOrbitals = 3);
    /** A shortcut to Lattice::Presets::addPSite with \f$U'=U-2J\f$, i.e. U_p = U - 2.0* J */
    RealExpr CoulombP(std::string const& Label, RealType U, RealType J, RealType Level, unsigned short NOrbitals = 3);
    ComplexExpr CoulombP(std::string const& Label, ComplexType U, ComplexType J, ComplexType Level, unsigned short NOrbitals = 3);

    /** Adds a magnetic \f$ \sum\limits_\alpha mH \frac{1}{2} (n_{i\alpha\uparrow} - n_{i\alpha\downarrow}) \f$ splitting to a given site. Valid only for 2 spins.
     * \param[in] L A pointer to the Lattice to add the terms.
     * \param[in] label \f$i\f$ - label of the site.
     * \param[in] Magnetization \f$mH\f$ - magnetization to add.
     */
    RealExpr Magnetization(std::string const& Label, RealType Magnetization, unsigned short NOrbitals = 1);
    ComplexExpr Magnetization(std::string const& Label, ComplexType Magnetization, unsigned short NOrbitals = 1);

    /** Adds a SzSz \f[ \sum\limits_{\alpha} J \frac{1}{2}(n_{i\alpha\uparrow} - n_{i\alpha\downarrow})\frac{1}{2}(n_{j\alpha\uparrow} - n_{j\alpha\downarrow}) \f]
     * interaction terms. Valid only for 2 spins
     * \param[in] L A pointer to the Lattice to add the terms.
     * \param[in] Label1 \f$i\f$ - label of the first connected site.
     * \param[in] Label2 \f$j\f$ - label of the second connected site. Site can be choosen the same as the first site.
     * \param[in] ExchJ \f$J\f$ - magnetic exchange constant.
     * \param[in] Orbitals Total amount of orbitals on the site. By default equal to 1.
     * \param[in] Spins Total amount of spin components on the site. By default equal to 2. Works only for 2 spins.
     */
    RealExpr SzSz(std::string const& Label1, std::string const& Label2, RealType ExchJ, unsigned short NOrbitals = 1);
    ComplexExpr SzSz(std::string const& Label1, std::string const& Label2, ComplexType ExchJ, unsigned short NOrbitals = 1);

    /** Adds a spin-spin \f[ \sum\limits_{\alpha} J \hat S_{i\alpha} \hat S_{j\alpha} \f]
     * interaction terms. Valid only for 2 spins
     * \param[in] L A pointer to the Lattice to add the terms.
     * \param[in] Label1 \f$i\f$ - label of the first connected site.
     * \param[in] Label2 \f$j\f$ - label of the second connected site. Site can be choosen the same as the first site.
     * \param[in] ExchJ \f$J\f$ - magnetic exchange constant.
     */
    RealExpr SS(std::string const& Label1, std::string const& Label2, RealType ExchJ, unsigned short NOrbitals = 1);
    ComplexExpr SS(std::string const& Label1, std::string const& Label2, ComplexType ExchJ, unsigned short NOrbitals = 1);

    //
    // Bosons
    //

    /** Generates a single energy level term \f$\varepsilon c^{\dagger}_{i\alpha\sigma}c_{i\alpha\sigma} \f$ on a local site for a given spin and orbital.
     * \param[in] Label \f$i\f$ - site affected by this Lattice::Term.
     * \param[in] Value \f$\varepsilon\f$ - the energy level.
     * \param[in] orbital \f$\alpha\f$ - affected orbital of the site.
     * \param[in] spin \f$\sigma\f$ - affected spin component.
     */
    RealExpr BosonLevel(std::string const& Label, RealType Value, unsigned short ExtraIndex);
    ComplexExpr BosonLevel(std::string const& Label, ComplexType Value, unsigned short ExtraIndex);

    /**
     * TODO: Bose-Hubbard interaction term
     */
    RealExpr BosonInteraction(std::string const& Label, RealType Value, unsigned short ExtraIndex);
    ComplexExpr BosonInteraction(std::string const& Label, ComplexType Value, unsigned short ExtraIndex);

    /**
     * TODO: Holstein coupling
     */
    RealExpr HolsteinInteraction(std::string const& Label, RealType Value, unsigned short Orbital, unsigned short BosonExtraIndex);
    ComplexExpr HolsteinInteraction(std::string const& Label, ComplexType Value, unsigned short Orbital, unsigned short BosonExtraIndex);

} // namespace Pomerol::LatticePresets
} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_LATTICEPRESETS_H
