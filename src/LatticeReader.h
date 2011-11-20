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

/** \file src/LatticeReader.h
** \brief A lattice file dictionary unit.
** 
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_LATTICEREADER_H
#define __INCLUDE_LATTICEREADER_H

#include"Misc.h"
#include<json/json.h>

namespace Pomerol{

/** This class reads data from input JSON file, defined by a string and returns all of its entries by a request.
 *  The data is stored as a jsoncpp dictonary. */
class LatticeReader
{
    /** The jsoncpp dictionary variable */
    Json::Value *root;

public:
    /** Empty constructor */
    LatticeReader ();
    /** Read the contents of a dictionary from an external JSON file */
    int readinFromJSON(const std::string filename);
    /** Get the contents of the dictionary */
    const Json::Value& getDictionary();
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_LATTICEREADER_H
