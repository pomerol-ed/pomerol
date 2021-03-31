/** \file include/pomerol/Lattice.h
** \brief A lattice handler.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_LATTICE_H
#define __INCLUDE_LATTICE_H

#include "Misc.h"

namespace Pomerol{

/** This class stores the information about a lattice.
 */
class Lattice
{
friend class LatticePresets;
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
    void printSites() const;

protected:
    /** A map between the particular Lattice::Site and it's label. */
    SiteMap Sites;
    //Lattice::TermStorage Terms1;
    TermStorage* Terms;
public:
    /** Empty constructor. */
    Lattice();
    /** Copy-constructor. */
    Lattice(const Lattice &l);
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
    enum op_type  {annihilation, creation };
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
    ComplexType Value;
    /** This returns the order of Term. Also the inheritance from TermPointer is provided by this method. */
    unsigned int getOrder() const;
    /** Constructor */
    Term(unsigned int N);

    /** Full constructor */
    Term(unsigned int N, bool OperatorSequence[ ], ComplexType Value, std::string SiteLabels[ ], unsigned short Orbitals[ ], unsigned short Spins[ ]);

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

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_LATTICE_H
