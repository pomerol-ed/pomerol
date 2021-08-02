/** \file include/pomerol/ComputableObject.h
** \brief A parent abstract class for all complex computable objects, i.e. Greens Function, Density matrix, Two Particle GF etc.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_H
#define POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_H

#include <exception>

namespace Pomerol {

struct ComputableObject {
    /** Computation statuses of the object. */
    enum StatusEnum {Constructed, Prepared, Computed};
protected:
    /** Current status of an object */
    StatusEnum Status = Constructed;
public:
    ComputableObject() = default;

    /** Returns the current status of an object */
    StatusEnum getStatus() const { return Status; }
    void setStatus(StatusEnum Status_in){
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

#endif // #ifndef POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_H
