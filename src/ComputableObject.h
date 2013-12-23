//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2011 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "Misc.h"

#ifndef __INCLUDE_COMPUTABLEOBJECT_H
#define __INCLUDE_COMPUTABLEOBJECT_H

namespace Pomerol{

struct ComputableObject {
    /** Computation statuses of the object. */
    enum {Constructed, Prepared, Computed};
protected:
    /** Current status of an object */
    unsigned int Status = Constructed;
public:
    /** Returns the current status of an object */
    unsigned int getStatus(){return Status;};
    class exStatusMismatch : public std::exception { virtual const char* what() const throw() { return "Object status mismatch"; } };
};

} // end of namespace Pomerol
#endif
