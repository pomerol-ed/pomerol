//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/MatsubaraContainers.hpp
/// \brief Container class designed to store values of functions of
/// multiple fermionic Matsubara frequencies.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_POMEROL_MATSUBARACONTAINERS_HPP
#define POMEROL_INCLUDE_POMEROL_MATSUBARACONTAINERS_HPP

#include "Misc.hpp"

#include <cmath>
#include <vector>

namespace Pomerol {

/// \addtogroup Misc
///@{

/// \brief Container for functions of three Matsubara frequencies.
///
/// Container class that stores values of a function \f$f(i\omega_1,i\omega_2,i\omega_3)\f$,
/// where \f$\omega_1,\omega_2,\omega_3\f$ are three fermionic Matsubara frequencies.
/// \tparam SourceObject Type of the source function object \f$f\f$.
template <typename SourceObject> class MatsubaraContainer4 {

    /// Stored elements are created by calling
    /// Source.createElement(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3).
    SourceObject const& Source;

    /// Number of bosonic Matsubara frequencies \f$\Omega = \omega_1 + \omega_2\f$,
    /// for which values are precomputed and stored.
    long NumberOfMatsubaras = 0;
    /// Stored precomputed values. Each element of this vector corresponds to one bosonic
    /// Matsubara frequency \f$\Omega = \omega_1 + \omega_2\f$, its matrix elements
    /// correspond to fermionic frequencies \f$\nu = \omega_1\f$ and \f$\omega_3 = \nu'\f$.
    std::vector<ComplexMatrixType> Values;
    /// Index offsets between fermionic indices \f$\nu, \nu'\f$ and matrix elements of \ref Values,
    /// one offsets per bosonic frequency.
    std::vector<long> FermionicIndexOffset;

public:
    /// Construct from a source function object.
    /// The container is initially empty and shall be populated with values
    /// by a subsequent call to \ref fill().
    /// \param[in] Source Source function object used to compute stored values.
    explicit MatsubaraContainer4(SourceObject const& Source) : Source(Source) {}

    /// Return a value of the function for a combination of Matsubara frequencies \f$\omega_1,\omega_2,\omega_3\f$.
    /// If value has not been precomputed for the given combination of frequencies, then
    /// Source.value(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3) will be called to obtain it.
    /// \param[in] MatsubaraNumber1 Index \f$n_1\f$ of Matsubara frequency \f$\omega_1 = \pi(2n_1+1)/\beta.\f$.
    /// \param[in] MatsubaraNumber2 Index \f$n_2\f$ of Matsubara frequency \f$\omega_2 = \pi(2n_2+1)/\beta\f$.
    /// \param[in] MatsubaraNumber3 Index \f$n_3\f$ of Matsubara frequency \f$\omega_3 = \pi(2n_3+1)/\beta\f$.
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /// Fill the container with precomputed values from the source function object.
    /// Each value is created by calling Source.value(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3).
    /// \param[in] NumberOfMatsubaras Number of positive fermionic Matsubara frequencies
    ///            \f$\omega_1\f$ and \f$\omega_2\f$ for which values are precomputed and stored.
    void fill(long NumberOfMatsubaras);

    /// Get the number of positive fermionic Matsubara frequencies \f$\omega_1\f$ and \f$\omega_2\f$
    /// for which values are precomputed and stored.
    long getNumberOfMatsubaras() const { return NumberOfMatsubaras; }
};

///@}

template <typename SourceObject> inline void MatsubaraContainer4<SourceObject>::fill(long NumberOfMatsubaras) {
    this->NumberOfMatsubaras = NumberOfMatsubaras;

    if(NumberOfMatsubaras == 0) {
        Values.clear();
        FermionicIndexOffset.clear();
        return;
    } else {
        Values.resize(4 * NumberOfMatsubaras - 1);
        FermionicIndexOffset.resize(4 * NumberOfMatsubaras - 1);
    }

    // \omega_1 = \nu, \omega_3 = \nu', \omega_1+\omega_2 = \Omega
    for(long BosonicIndexV = 0; BosonicIndexV <= 4 * NumberOfMatsubaras - 2; ++BosonicIndexV) {
        long BosonicIndex = BosonicIndexV - 2 * NumberOfMatsubaras;
        long FermionicMatrixSize = 2 * NumberOfMatsubaras - std::abs(BosonicIndex + 1);
        Values[BosonicIndexV].resize(FermionicMatrixSize, FermionicMatrixSize);
        FermionicIndexOffset[BosonicIndexV] = (BosonicIndex < 0 ? 0 : BosonicIndex + 1) - NumberOfMatsubaras;

        for(long NuIndexM = 0; NuIndexM < FermionicMatrixSize; ++NuIndexM)
            for(long NupIndexM = 0; NupIndexM < FermionicMatrixSize; ++NupIndexM) {
                long MatsubaraNumber1 = NuIndexM + FermionicIndexOffset[BosonicIndexV];
                long MatsubaraNumber2 = BosonicIndex - MatsubaraNumber1;
                long MatsubaraNumber3 = NupIndexM + FermionicIndexOffset[BosonicIndexV];
                Values[BosonicIndexV](NuIndexM, NupIndexM) =
                    Source.value(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3);
            }
    }
}

template <typename SourceObject>
inline ComplexType MatsubaraContainer4<SourceObject>::operator()(long MatsubaraNumber1,
                                                                 long MatsubaraNumber2,
                                                                 long MatsubaraNumber3) const {
    long BosonicIndexV = MatsubaraNumber2 + MatsubaraNumber1 + 2 * NumberOfMatsubaras;
    if(BosonicIndexV >= 0 && BosonicIndexV <= 2 * (2 * NumberOfMatsubaras - 1)) {
        long NuIndexM = MatsubaraNumber1 - FermionicIndexOffset[BosonicIndexV];
        long NupIndexM = MatsubaraNumber3 - FermionicIndexOffset[BosonicIndexV];
        if(NuIndexM >= 0 && NuIndexM < Values[BosonicIndexV].rows() && NupIndexM >= 0 &&
           NupIndexM < Values[BosonicIndexV].cols())
            return Values[BosonicIndexV](NuIndexM, NupIndexM);
    }

    DEBUG("MatsubaraContainer4 at " << this << ": "
                                    << "cache miss for n1 = " << MatsubaraNumber1 << ", "
                                    << "n2 = " << MatsubaraNumber2 << ", "
                                    << "n3 = " << MatsubaraNumber3 << " (NumberOfMatsubaras = " << NumberOfMatsubaras
                                    << "), fetching a raw value from " << &Source);

    return Source.value(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3);
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_MATSUBARACONTAINERS_HPP
