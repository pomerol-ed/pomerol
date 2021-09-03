//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
    enum StatusEnum { Constructed, Prepared, Computed };

protected:
    /** Current status of an object */
    StatusEnum Status = Constructed;

public:
    ComputableObject() = default;

    /** Returns the current status of an object */
    StatusEnum getStatus() const { return Status; }
    void setStatus(StatusEnum Status_in) {
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
