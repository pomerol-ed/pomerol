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


/** \file IndexClassification.h
**  \brief Declaration of SingleIndex, TermsList, IndexClassification classes.
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_INDEXCLASSIFICATION_H
#define __INCLUDE_INDEXCLASSIFICATION_H

#include"Misc.h"
#include"LatticeAnalysis.h"
#include"Term.h"
#include<json/json.h>

/**
 * SingleIndex class is an abstract class to handle all info about current ParticleIndex. 
 *  is a site+spin index. The class also handles the number of the bit and the type of the site which it belongs to
 */ 

class SingleIndex
{
public:
    unsigned short site;         //!< The number of site
    unsigned short spin;         //!< Spin up(1) or Spin down(0)
    unsigned short type;         //!< The type of site which handles this spin. Can be s,p,d,f (0,1,2,3).
    ParticleIndex bitNumber;    //!< The number of bit
    RealType LocalMu;            //!< The energy in the -mu*N term for the site ( Useful for adjusting a filling in the orbital )
    void setIndexNumber(const unsigned short &in); 
    virtual void print_to_screen()=0;
};

/**
 * sSingleIndex is a class, inherited from SingleIndex. It provides full info about an index in s-orbital
 */
class sSingleIndex : public SingleIndex
{
public:
    RealType U;
    sSingleIndex(unsigned short site_, unsigned short type_, unsigned short spin_, RealType U_, RealType LocalMu_):U(U_){site=site_;type=type_;spin=spin_;LocalMu=LocalMu_;};
    friend std::ostream& operator<<(std::ostream& output, const sSingleIndex& out);
    void print_to_screen(){std::cout << *this << std::endl;};
};

/**
 * pSingleIndex is a class, inherited from SingleIndex. It provides full info about an index in p-orbital
 */
class pSingleIndex : public SingleIndex
{
public:
    RealType U;
    RealType J;
    std::string basis;
    unsigned short orbital; 
    pSingleIndex(unsigned short site_, unsigned short type_, unsigned short spin_,unsigned short orbital_, const std::string &basis_, RealType U_, RealType J_):U(U_),J(J_){site=site_;type=type_;spin=spin_;orbital=orbital_;basis=basis_;};
    friend std::ostream& operator<<(std::ostream& output, const pSingleIndex& out);
    void print_to_screen(){std::cout << *this << std::endl;};
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
   void addTerm(Term* in);                                           //!< Add a new term to the List
   TermsList(){maxOrder=0;};
   std::list<Term*> &getTerms(unsigned short order);                 //!< Get a list of terms of a given order
   friend std::ostream& operator<<(std::ostream& output,TermsList& out);
};

/**
 * The IndexClassification class handles all Terms and Indexs for a given Lattice.
 * It reads the structure of the lattice from external lattice json file 
 * The lattice file should be written in JSON format and contain all the info about the sites of the system.
 * IndexClassification parses this file and creates an index classification instead of site in order to associate a unique index "bit" to a site + spin configuration.
 * Then the HoppingMatrix to reproduce all the hoppings of the first order between different bits is created.
 * Finally, the list of formula for the model is written on a current lattice through a list of Terms to instruct the entering of the matrix elements of the hamiltonian.
 */
class IndexClassification
{
public:
    struct IndexPermutation;
private:
    LatticeAnalysis &Lattice;
    int N_bit;                                          //!< The length of SingleIndexList. Defines the number of states in system as 2^N_bit.
    RealMatrixType HoppingMatrix;                       //!< The matrix to show all hoppings between different bits.
    std::vector<SingleIndex*> SingleIndexList;          //!< A list of all SingleIndex pointers -> contains all information on each index
    TermsList Terms;                                    //!< The list of all terms for the current Lattice
    std::list<IndexPermutation> EquivalentPermutations; //!< A list of all equivalent IndexPermutations
public:
    IndexClassification(LatticeAnalysis &Lattice);
    int prepare();                                      //!< Reads all info from Lattice (LatticeAnalysis class)
    void printIndexList();                              //!< Print SingleIndexList to screen
    void printHoppingMatrix();                          //!< Print HoppingMatrix to screen
    void printTerms();                                  //!< Print TermsList to screen
    void printEquivalentPermutations();                 //!< Print all equivalent ParticleIndex permutations
    RealMatrixType& getHoppingMatrix();                 //!< Returns a Hopping Matrix
    std::vector<SingleIndex*> &getSingleIndexList();    //!< Returns a SingleIndexList
    const int& getIndexSize() const;                    //!< Returns N_bit
    bool checkIndex(ParticleIndex in);                  //!< Returns true if current Index may exist, i.e. if it is between 0 and N_bit-1
    TermsList& getTermsList();                    //!< Returns Terms
    const std::list<IndexClassification::IndexPermutation> &getEquivalentIndexPermutations();  //!< Returns a vector of equivalent index permutation
private:
    void defineIndices();                               //!< Define the bit classification from the info from file
    void defineHopping();                               //!< Define Hopping from classification
    void defineTerms();                                 //!< Define all terms
    void defineEquivalentIndexPermutations();           //!< Define all equivalent index permutations
    void definePorbitalSphericalTerms(pSingleIndex **Indexs); //!< The bunch of methods to add Terms for p-orbital written in Spherical basis
    void definePorbitalNativeTerms(pSingleIndex **Indexs);    //!< The bunch of methods to add Terms for p-orbital written in Native basis
    std::vector<unsigned short>& findIndexs(const unsigned short &site); //!< A method to find all bits for a corresponding site
    unsigned short findIndex(const unsigned short &site, const unsigned short &orbital, const unsigned short &spin); //!< A method to find bit, which corresponds to given site, orbital and spin
};

/** IndexClassification::IndexPermutation - a container of permutation of indices of the system
 */
class IndexClassification::IndexPermutation
{
friend void IndexClassification::defineEquivalentIndexPermutations();
protected:
    IndexPermutation(std::vector<ParticleIndex>& in);
public:
    std::vector<ParticleIndex> Permutation;
    bool isCorrect();
    ~IndexPermutation();
    friend std::ostream& operator<<(std::ostream& output, const IndexPermutation& out);
};


#endif // endif :: #ifndef #__INCLUDE_INDEXCLASSIFICATION_H

