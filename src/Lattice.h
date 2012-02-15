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

#ifndef __INCLUDE_LATTICE_H
#define __INCLUDE_LATTICE_H

#include "Misc.h"
#include "Logger.h"
#include <json/json.h>

namespace Pomerol{

/** This class stores the information about a lattice. 
 *  It can be read from a JSON file.
 */
class Lattice
{
public:

    struct Site;
    template <unsigned short> struct Term;

private:
    /* A map between the particular Site and it's label */
    std::map<unsigned short, Site> Sites; 
    //std::list<Term<unsigned short> > Terms;

public:
    /** Empty constructor */
    Lattice ();
};

class JSONLattice : public Lattice
{
    Json::Value *root;
    void readSites(Json::Value &JSONSites);
    void readTerms(Json::Value &JSONTerms);
    public:
    /** Read the contents of a dictionary from an external JSON file */
    int readin (const std::string &filename);
    JSONLattice();
};

struct Lattice::Site{
friend class Lattice;
protected:
    std::string label;
    unsigned short OrbitalSize;
    unsigned short SpinSize;
friend std::ostream& operator<<(std::ostream& output, const Site& out);
};

template <unsigned short N> struct Lattice::Term{
friend class Lattice;
protected:
   unsigned short type;
   bool OperatorOrder[N]; 
   std::string ConnectedSites[N];
   unsigned short Spins[N];
   unsigned short Orbitals[N];
   RealType Value;
public:
    Term();
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_LATTICEREADER_H
