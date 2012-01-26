//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
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

/** \file src/MatsubaraContainers.h
** \brief Templated container classes designed to store values of functions of Matsubara frequencies.
**
** \author Andrey Antipov (antipov@ct-qmc.org)
** \author Igor Krivenko (igor@shg.ru)
*/
#ifndef __INCLUDE_MATSUBARACONTAINERS_H 
#define __INCLUDE_MATSUBARACONTAINERS_H 

#include "Misc.h"

namespace Pomerol{

// class MatsubaraContainer1
template<typename ObjectToFillFrom>
class MatsubaraContainer1 {
    VectorType Values;

public:

    const long NumberOfMatsubaras;

    MatsubaraContainer1(long NumberOfMatsubaras);

    bool isInContainer(long MatsubaraNumber) const;
    ComplexType operator()(long MatsubaraNumber) const;

    void fill(const ObjectToFillFrom* pObject);	// calls pObject->rawValue(MatsubaraNumber)
};

template<typename ObjectToFillFrom> inline
MatsubaraContainer1<ObjectToFillFrom>::MatsubaraContainer1(long NumberOfMatsubaras) :
    NumberOfMatsubaras(NumberOfMatsubaras),
    Values(2*NumberOfMatsubaras)
{}

template<typename ObjectToFillFrom> inline
bool MatsubaraContainer1<ObjectToFillFrom>::isInContainer(long MatsubaraNumber) const
{
    return MatsubaraNumber < NumberOfMatsubaras && MatsubaraNumber >= -NumberOfMatsubaras;
}

template<typename ObjectToFillFrom> inline
ComplexType MatsubaraContainer1<ObjectToFillFrom>::operator()(long MatsubaraNumber) const
{
    return Values(NumberOfMatsubaras+MatsubaraNumber);
}

template<typename ObjectToFillFrom> inline
void MatsubaraContainer1<ObjectToFillFrom>::fill(const ObjectToFillFrom* pObject)
{
    for(long MatsubaraNum=-NumberOfMatsubaras; MatsubaraNum<NumberOfMatsubaras; MatsubaraNum++)
        Values[MatsubaraNum+NumberOfMatsubaras] = pObject->rawValue(MatsubaraNum);
}

// class MatsubaraContainer4
template<typename ObjectToFillFrom>
class MatsubaraContainer4 {

    std::vector<MatrixType> Values;
    std::vector<long> FermionicIndexOffset;

public:

    const long NumberOfMatsubaras;

    MatsubaraContainer4(long NumberOfMatsubaras);

    bool isInContainer(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    void fill(const ObjectToFillFrom* pObject);	// calls pObject->rawValue(MatsubaraNumber)
};

template<typename ObjectToFillFrom> inline
MatsubaraContainer4<ObjectToFillFrom>::MatsubaraContainer4(long NumberOfMatsubaras) :
    NumberOfMatsubaras(NumberOfMatsubaras),
    Values(4*NumberOfMatsubaras-1),
    FermionicIndexOffset(4*NumberOfMatsubaras-1)
{
    for(long BosonicIndex=-2*NumberOfMatsubaras; BosonicIndex<=2*NumberOfMatsubaras-2; ++BosonicIndex){
        long FermionicMatrixSize = 2*NumberOfMatsubaras - abs(BosonicIndex+1);
        Values[BosonicIndex+2*NumberOfMatsubaras].resize(FermionicMatrixSize,FermionicMatrixSize);
        FermionicIndexOffset[BosonicIndex+2*NumberOfMatsubaras] =
            (BosonicIndex < 0 ? 0 : BosonicIndex+1) -NumberOfMatsubaras;
    }
}

template<typename ObjectToFillFrom> inline
void MatsubaraContainer4<ObjectToFillFrom>::fill(const ObjectToFillFrom* pObject)
{
    // \omega_1 = \nu, \omega_3 = \nu', \omega_1+\omega_2 = \Omega
    for(long BosonicIndexV=0; BosonicIndexV<=4*NumberOfMatsubaras-2; ++BosonicIndexV){
        long FermionicMatrixSize = Values[BosonicIndexV].rows();
        for(long NuIndexM=0; NuIndexM<FermionicMatrixSize; ++NuIndexM)
        for(long NupIndexM=0; NupIndexM<FermionicMatrixSize; ++NupIndexM){
            long MatsubaraNumber1 = NuIndexM+FermionicIndexOffset[BosonicIndexV];
            long MatsubaraNumber2 = BosonicIndexV - MatsubaraNumber1 - 2*NumberOfMatsubaras;
            long MatsubaraNumber3 = NupIndexM+FermionicIndexOffset[BosonicIndexV];;
            Values[BosonicIndexV](NuIndexM,NupIndexM) =
                pObject->rawValue(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
        }
    }
}

template<typename ObjectToFillFrom> inline
bool MatsubaraContainer4<ObjectToFillFrom>::isInContainer(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    long BosonicIndexV = MatsubaraNumber1+MatsubaraNumber2+2*NumberOfMatsubaras;
    if(BosonicIndexV >=0 && BosonicIndexV <= 2*(2*NumberOfMatsubaras-1)){
        long NuIndexM = MatsubaraNumber1-FermionicIndexOffset[BosonicIndexV];
        long NupIndexM = MatsubaraNumber3-FermionicIndexOffset[BosonicIndexV];
        if(NuIndexM >= 0 && NuIndexM < Values[BosonicIndexV].size() &&
           NupIndexM >= 0 && NupIndexM <= Values[BosonicIndexV].size())
            return true;
    }

    return false;
}

template<typename ObjectToFillFrom> inline
ComplexType MatsubaraContainer4<ObjectToFillFrom>::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    long BosonicIndexV = MatsubaraNumber2 + MatsubaraNumber1 + 2*NumberOfMatsubaras;
    long NuIndexM = MatsubaraNumber1-FermionicIndexOffset[BosonicIndexV];
    long NupIndexM = MatsubaraNumber3-FermionicIndexOffset[BosonicIndexV];

    return Values[BosonicIndexV](NuIndexM,NupIndexM);
}

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_MATSUBARACONTAINERS_H 
