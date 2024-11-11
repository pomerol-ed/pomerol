//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/ComputableObject.hpp
/// \brief A base class for computable objects, e.g. Green's function, two-particle GF, etc.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_HPP
#define POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_HPP

#include <stdexcept>
#include <string>

namespace Pomerol {

/// \addtogroup Misc
///@{

/// A base class for computable objects.
struct ComputableObject {

    /// Computation status of the object.
    enum StatusEnum {
        Constructed, ///< Object has been constructed.
        Prepared,    ///< Object has been prepared for computation (usually means memory allocation).
        Computed     ///< Object has been computed.
    };

protected:
    /// Current computation status
    StatusEnum Status = Constructed;

public:
    ComputableObject() = default;

    /// Return the current computation status.
    StatusEnum getStatus() const { return Status; }
    /// Set the computation status.
    /// \param[in] Status_in New computation status.
    void setStatus(StatusEnum Status_in) {
        if(Status_in > Computed)
            throw StatusMismatch("Invalid computable object status");
        Status = Status_in;
    }

    /// Exception: Unexpected computation status of a computable object
    class StatusMismatch : public std::runtime_error {
    public:
        explicit StatusMismatch(std::string const& str) : std::runtime_error(str) {}
    };
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_COMPUTABLEOBJECT_HPP
