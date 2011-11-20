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
//

/** \file src/ComputableObject.h
** \brief A parent abstract class for all complex computable objects, i.e. Greens Function, Density matrix, Two Particle GF etc.
** 
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#include "Misc.h"

#ifndef __INCLUDE_COMPUTABLEOBJECT_H
#define __INCLUDE_COMPUTABLEOBJECT_H

namespace Pomerol{

/** This abstract class is a prototype for any computable object in the current code, i.e. Greens Function, Two Particle GF, Density Matrix etc...
 *  It defines the current state of calculation and declares virtual methods "compute" and "dump". 
 *  The enumeration of statuses is listed in the Misc.h file
 */
class ComputableObject{
protected:
    /** Current status of an object */
    ObjectStatus  Status;
public:
    /** Returns the current status of an object */
    ObjectStatus getStatus(){return Status;};

    /** Constructor - set status to Constructed */
    ComputableObject():Status(Constructed){};            
    //virtual void prepare();                   //!< Prepare all the containers of the object, do preliminary fast routines
    //virtual void compute();                   //!< Do the most expensive calculation. Finishing it means the object has finished all calculation jobs
};

} // end of namespace Pomerol
#endif
