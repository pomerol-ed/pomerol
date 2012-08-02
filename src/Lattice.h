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
** \author Andrey Antipov (antipov@ct-qmc.org)
** \author Igor Krivenko (igor@shg.ru)
*/

#ifndef __INCLUDE_LATTICE_H
#define __INCLUDE_LATTICE_H

#include "Misc.h"
#include "Logger.h"
#include <json/json.h>

namespace Pomerol{

/** This class stores the information about a lattice. 
 */
class Lattice
{
public:
    /** This holds the information about a given site of the lattice, namely it's label, number of spins and orbitals. */ 
    struct Site;
    /** This structure holds the information, about a written term in a formula - it's matrix element, corresponding site labels, spins and orbitals. */
    struct Term;
    /** A typedef for a list of pointers to the terms. */
    typedef std::list<Term*> TermList;
    /** A typedef for a map between the label and the corresponding site */
    typedef std::map<std::string, Site*> SiteMap;
    /** A storage for all the terms. Realized as a map between the order of the Lattice::Terms and the corresponding Lattice::TermList. */
    class TermStorage;
    /** A set of presets to fill the TermStorage and Sites for some commonly used examples. Look at the LatticePresets.h . */
    class Presets;
    
    /** Add a Site to list of Sites 
     * \param[in] S A site to add.
     */
    void addSite(Lattice::Site* S);
    /** Add a Site to list of Sites 
     * \param[in] Label Label of Site.
     * \param[in] Orbitals Amount of orbitals on Site.
     * \param[in] Spins Amount of spins on Site. By default 2
     */
    void addSite(const std::string &Label, unsigned short orbitals, unsigned short spins=2);

    /** Add a Term to the Lattice::Storage.
     * \param[in] T The Term to add.
     */
    void addTerm(const Term* T);
    /** Print all terms of the given order
     * \param[in] order The order of Terms to print 
     */
    void printTerms(unsigned int order);

    /** Print all Sites */
    void printSites();

protected:
    /** A map between the particular Lattice::Site and it's label. */
    SiteMap Sites; 
    //Lattice::TermStorage Terms1;
    TermStorage* Terms;
public:
    /** Empty constructor. */
    Lattice();
    /** Destructor. */
    ~Lattice();
    /** Returns a Lattice::Site for a given Label. */
    const Lattice::Site& getSite(const std::string& Label) const;
    /** Returns the map of Sites*/
    const Lattice::SiteMap& getSiteMap() const;
    /** Returns all stored terms */
    const Lattice::TermStorage& getTermStorage() const;
    /** An exception, which is thrown when a wrong Term added. */
    class exWrongLabel : public std::exception { 
        virtual const char* what() const throw();
    };
};


/** This structure holds the information about a given site of the lattice, namely it's label, number of spins and orbitals. */ 
struct Lattice::Site{
friend class Lattice;
public:
    /** Site label. */
    const std::string Label;
    /** Amount of orbitals on a site. */
    const unsigned short OrbitalSize;
    /** Amount of spins on a site. */
    const unsigned short SpinSize;
    /** Full constructor 
     * \param[in] Label Site label
     * \param[in] OrbitalSize Number of Orbitals for current site 
     * \param[in] SpinSize Number of spins for current site 
     * */
    Site(const std::string& Label, unsigned short OrbitalSize, unsigned short SpinSize );
/** Make the object printable. */
friend std::ostream& operator<<(std::ostream& output, const Site& out);
};


/** This structure holds the information, about a written term in a formula - it's matrix element, corresponding site labels, spins and orbitals. 
 */
struct Lattice::Term { 
friend class Lattice;
private:
    /** Total amount of operators in Lattice::Term. */
    unsigned int N;
public:
    /** A set of presets to simplify term generation */
    class Presets;
    /** The order of the creation/annihilation operator in the Lattice::Term. */
    std::vector<bool> OperatorSequence; 
    /** An array with labels of sites, connected by this Lattice::Term. */
    std::vector<std::string> SiteLabels;
    /** An array of spins on the sites, which are connected by this Lattice::Term. */
    std::vector<unsigned short> Spins;
    /** An array of orbitals on the sites, which are connected by this Lattice::Term. */
    std::vector<unsigned short> Orbitals;
    /** The matrix element of the Lattice::Term. */
    RealType Value;
    /** This returns the order of Term. Also the inheritance from TermPointer is provided by this method. */
    unsigned int getOrder() const;
    /** Constructor */
    Term(unsigned int N);

    /** Full constructor */
    Term(unsigned int N, bool OperatorSequence[ ], RealType Value, std::string SiteLabels[ ], unsigned short Orbitals[ ], unsigned short Spins[ ]);

    /** Copy-constuctor 
     * \param[in] in A Lattice::Term to copy.
     */
    Term(const Lattice::Term &in);
/** Make the Term printable */
friend std::ostream& operator<< (std::ostream& output, const Term& out);
};


/** A storage for all the terms. Realized as a map between the order of the Lattice::Terms and the corresponding Lattice::TermList */
class Lattice::TermStorage { 
friend class Lattice;
protected:
    /** A storage for the TermLists for the corresponding order */
    std::map<unsigned int, Lattice::TermList> Terms;
    /** Stores the maximum total number of operators in all Terms */
    unsigned int MaxTermOrder;
public:
    /** Add a Term to the storage.
     * \param[in] T The Term to add.
     */
    int addTerm(const Term* T);
    /** Get a List of Terms of a given order.
     * \param[in] N The required order of Terms. */
    const Lattice::TermList &getTerms (unsigned int N) const;
    
    /** Returns largest the number of operators in all stored Terms */
    const unsigned int getMaxTermOrder() const;
    /** Empty constructor */
    TermStorage();
};

/** This class stores the information about a lattice and reads it from a provided JSON file. */
class JSONLattice : public Lattice
{
    class JSONPresets;
    //typedef void (addSite)(Lattice *L ) JSONPreset;
    //std::map<std::string, JSONPreset> Meth; 
    /** Read and store the information about the sites of the lattice. This also add some local Terms<2> to the Terms map.
     * \param[in] JSONSites A "Sites" section of the dictionary from the JSON file.
     */
    void readSites(Json::Value &JSONSites);
    /** Read and store the information about the Terms between the sites of the lattice.
     * \param[in] JSONTerm A "Terms" section of the dictionary from the JSON file.
     */
    void readTerms(Json::Value &JSONTerms);
    public:
    /** Read the contents of a dictionary from an external JSON file. */
    int readin (const std::string &filename);
    /** Empty constructor. */
    JSONLattice();
};


} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_LATTICE_H
