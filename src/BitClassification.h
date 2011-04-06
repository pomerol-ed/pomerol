// pomerol/trunk/src/BitClassification.h
// This file is a part of pomerol diagonalization code

/** \file BitClassification.h
**  \brief Declaration of BitInfo, TermsList, BitClassification classes.
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __BIT_CLASSIFICATION__
#define __BIT_CLASSIFICATION__

#include "config.h"
#include "LatticeAnalysis.h"
#include "Term.h"
#include <json/json.h>
#include <map>
#include <string>
#include <list>


/**
 * BitInfo class is an abstract class to handle all info about current bit. 
 * The bit is a site+spin index. The class also handles the number of the bit and the type of the site which it belongs to
 */ 

class BitInfo
{
public:
    unsigned short site;         //!< The number of site
    unsigned short spin;         //!< Spin up(1) or Spin down(0)
    unsigned short type;        //!< The type of site which handles this spin. Can be s,p,d,f (0,1,2,3).
    unsigned short bitNumber;    //!< The number of bit
    RealType LocalMu;            //!< The energy in the -mu*N term for the site ( Useful for adjusting a filling in the orbital )
    void setBitNumber(const unsigned short &in); 
    virtual void print_to_screen()=0;
};

/**
 * sBitInfo is a class, inherited from BitInfo. It provides full info about the bit in s-orbital
 */
class sBitInfo : public BitInfo
{
public:
    RealType U;
    sBitInfo(unsigned short site_, unsigned short type_, unsigned short spin_, RealType U_, RealType LocalMu_):U(U_){site=site_;type=type_;spin=spin_;LocalMu=LocalMu_;};
    friend std::ostream& operator<<(std::ostream& output, const sBitInfo& out);
    void print_to_screen(){cout << *this << endl;};
};

/**
 * pBitInfo is a class, inherited from BitInfo. It provides full info about the bit in p-orbital
 */
class pBitInfo : public BitInfo
{
public:
    RealType U;
    RealType J;
    string basis;
    unsigned short orbital; 
    pBitInfo(unsigned short site_, unsigned short type_, unsigned short spin_,unsigned short orbital_, const string &basis_, RealType U_, RealType J_):U(U_),J(J_){site=site_;type=type_;spin=spin_;orbital=orbital_;basis=basis_;};
    friend std::ostream& operator<<(std::ostream& output, const pBitInfo& out);
    void print_to_screen(){cout << *this << endl;};
};

/** 
 * A TermsList is a class to handle all the terms of different orders.
 * It consists of several lists of pointers to Terms, with an unsigned int map between the order and the list.
 */
class TermsList 
{
   std::map <unsigned short, std::list<Term*> > TermsMap;            //!< The map between the order of the terms and corresponding list
   unsigned short maxOrder;

   public:
   void addTerm(Term* in);                            //!< Add a new term to the List
   TermsList(){maxOrder=0;};
   std::list<Term*> &getTerms(unsigned short order);                 //!< Get a list of terms of a given order
   friend std::ostream& operator<<(std::ostream& output,TermsList& out);
};

/**
 * The BitClassification class handles all Terms and Bits for a given Lattice.
 * It reads the structure of the lattice from file "Lattice.json".
 * The input file should be written in JSON format and contain all the info about the sites of the system.
 * BitClassification parses this file and creates a bit classification instead of site in order to associate a unique index "bit" to a site + spin configuration.
 * Then the HoppingMatrix to reproduce all the hoppings of the first order between different bits is created.
 * Finally, the list of formula for the model is written on a current lattice through a list of Terms to instruct the entering of the matrix elements of the hamiltonian.
 */
class BitClassification
{
    LatticeAnalysis &Lattice;
    int N_bit;                                          //!< The length of BitInfoList. Defines the number of states in system as 2^N_bit.
    RealMatrixType HoppingMatrix;                       //!< The matrix to show all hoppings between different bits.
    vector<BitInfo*> BitInfoList;                       //!< A list of all Bits.
    TermsList Terms;                                    //!< The list of all terms for the current Lattice
public:
    BitClassification(LatticeAnalysis &Lattice);
    int prepare();                                      //!< Reads all info from Lattice (LatticeAnalysis class)
    void printBitInfoList();                            //!< Print BitInfoList to screen
    void printHoppingMatrix();                          //!< Print HoppingMatrix to screen
    void printTerms();                                  //!< Print TermsList to screen
    RealMatrixType& getHoppingMatrix();                 //!< Returns a Hopping Matrix
    std::vector<BitInfo*> &getBitInfoList();            //!< Returns a BitInfoList
    const int& getBitSize() const;                      //!< Returns N_bit
    bool checkIndex(ParticleIndex in);                  //!< Returns true if current Index may exist, i.e. if it is between 0 and N_bit-1
    TermsList& getTermsList();                          //!< Returns Terms
private:
    void defineBits();                                  //!< Define the bit classification from the info from file
    void defineHopping();                               //!< Define Hopping from classification
    void defineTerms();                                 //!< Define all terms
    void definePorbitalSphericalTerms(pBitInfo **Bits); //!< The bunch of methods to add Terms for p-orbital written in Spherical basis
    void definePorbitalNativeTerms(pBitInfo **Bits);    //!< The bunch of methods to add Terms for p-orbital written in Native basis
    vector<unsigned short>& findBits(const unsigned short &site); //!< A method to find all bits for a corresponding site
    unsigned short findBit(const unsigned short &site, const unsigned short &orbital, const unsigned short &spin); //!< A method to find bit, which corresponds to given site, orbital and spin
};

#endif

