/** \file include/pomerol/ComputableObject.h
** \brief A parent abstract class for all complex computable objects, i.e. Greens Function, Density matrix, Two Particle GF etc.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_H
#define POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_H

#include <stdexcept>
#include <string>

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
        if(Status_in > Computed)
            throw StatusMismatch("Invalid computable object status");
        Status = Status_in;
    }

    class StatusMismatch : public std::runtime_error {
    public:
        explicit StatusMismatch(std::string const& str) : std::runtime_error(str) {}
    };
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_H
