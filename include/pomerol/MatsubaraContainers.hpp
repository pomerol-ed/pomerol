//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/MatsubaraContainers.h
** \brief Templated container classes designed to store values of functions of Matsubara frequencies.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_MATSUBARACONTAINERS_HPP
#define POMEROL_INCLUDE_POMEROL_MATSUBARACONTAINERS_HPP

#include "Misc.hpp"

#include <cmath>
#include <vector>

namespace Pomerol {

// class MatsubaraContainer4
template <typename SourceObject> class MatsubaraContainer4 {
    SourceObject const& Source;
    long NumberOfMatsubaras = 0;
    std::vector<ComplexMatrixType> Values;
    std::vector<long> FermionicIndexOffset;

public:
    explicit MatsubaraContainer4(SourceObject const& Source) : Source(Source) {}

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    // calls Source.value(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3)
    void fill(long NumberOfMatsubaras);

    long getNumberOfMatsubaras() const { return NumberOfMatsubaras; }
};

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
