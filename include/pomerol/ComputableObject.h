/** \file include/pomerol/ComputableObject.h
** \brief A parent abstract class for all complex computable objects, i.e. Greens Function, Density matrix, Two Particle GF etc.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_COMPUTABLEOBJECT_H
#define __INCLUDE_COMPUTABLEOBJECT_H

#include "Misc.h"

#include <exception>

namespace Pomerol {

struct ComputableObject {
    /** Computation statuses of the object. */
    enum {Constructed, Prepared, Computed};
protected:
    /** Current status of an object */
    unsigned int Status = Constructed;
public:
    ComputableObject() = default;

    /** Returns the current status of an object */
    unsigned int getStatus() const { return Status; }
    void setStatus(unsigned int Status_in){
        if(Status_in > Computed) throw exStatusMismatch();
        Status = Status_in;
    }

    class exStatusMismatch : public std::exception {
        virtual const char* what() const noexcept override {
            return "Object status mismatch";
        }
    };
};

} // namespace Pomerol
#endif
