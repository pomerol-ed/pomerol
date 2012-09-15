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

#include"Misc.h"
#include"Logger.h"

namespace Pomerol{

// class MatsubaraContainer1
template<typename SourceObject>
class MatsubaraContainer1 {
    const SourceObject* pSource;
    ComplexVectorType Values;
    long NumberOfMatsubaras;

public:

    MatsubaraContainer1(const SourceObject* pSource);

    ComplexType operator()(long MatsubaraNumber) const;

    void fill(long NumberOfMatsubaras);	// calls pSource->value(MatsubaraNumber)
    long getNumberOfMatsubaras(void) const;
};

template<typename SourceObject> inline
MatsubaraContainer1<SourceObject>::MatsubaraContainer1(const SourceObject* pSource) :
    pSource(pSource),
    NumberOfMatsubaras(0)
{
}

template<typename SourceObject> inline
ComplexType MatsubaraContainer1<SourceObject>::operator()(long MatsubaraNumber) const
{
    if(MatsubaraNumber < NumberOfMatsubaras && MatsubaraNumber >= -NumberOfMatsubaras)
        return Values(NumberOfMatsubaras+MatsubaraNumber);
    else{
        DEBUG("MatsubaraContainer1 at " << this << ": " <<
              "cache miss for n = " << MatsubaraNumber <<
              " (NumberOfMatsubaras = " << NumberOfMatsubaras <<
              "), fetching a raw value from " << pSource
        )
        return pSource->value(MatsubaraNumber);
    }
}

template<typename SourceObject> inline
void MatsubaraContainer1<SourceObject>::fill(long NumberOfMatsubaras)
{
    this->NumberOfMatsubaras = NumberOfMatsubaras;
    Values.resize(2*NumberOfMatsubaras);
    for(long MatsubaraNum=-NumberOfMatsubaras; MatsubaraNum<NumberOfMatsubaras; MatsubaraNum++)
        Values[MatsubaraNum+NumberOfMatsubaras] = pSource->value(MatsubaraNum);
}

template<typename SourceObject> inline
long MatsubaraContainer1<SourceObject>::getNumberOfMatsubaras(void) const
{
    return NumberOfMatsubaras;
}

// class MatsubaraContainer4
template<typename SourceObject>
class MatsubaraContainer4 {
    const SourceObject* pSource;
    long NumberOfMatsubaras;
    std::vector<ComplexMatrixType> Values;
    std::vector<long> FermionicIndexOffset;

public:

    MatsubaraContainer4(void);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    // calls pSource->value(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3)
    void fill(const SourceObject* pSource, long NumberOfMatsubaras);
    long getNumberOfMatsubaras(void) const;
};

template<typename SourceObject> inline
MatsubaraContainer4<SourceObject>::MatsubaraContainer4(void) :
    pSource(NULL),
    NumberOfMatsubaras(0)
{}

template<typename SourceObject> inline
void MatsubaraContainer4<SourceObject>::fill(const SourceObject* pSource, long NumberOfMatsubaras)
{
    this->NumberOfMatsubaras = NumberOfMatsubaras;
    this->pSource = pSource;

    if(NumberOfMatsubaras==0){
        Values.resize(0);
        FermionicIndexOffset.resize(0);
        return;
    }else{
        Values.resize(4*NumberOfMatsubaras-1);
        FermionicIndexOffset.resize(4*NumberOfMatsubaras-1);
    }

    // \omega_1 = \nu, \omega_3 = \nu', \omega_1+\omega_2 = \Omega
    for(long BosonicIndexV=0; BosonicIndexV<=4*NumberOfMatsubaras-2; ++BosonicIndexV){
        long BosonicIndex = BosonicIndexV - 2*NumberOfMatsubaras;
        long FermionicMatrixSize = 2*NumberOfMatsubaras - abs(BosonicIndex+1);
        Values[BosonicIndexV].resize(FermionicMatrixSize,FermionicMatrixSize);
        FermionicIndexOffset[BosonicIndexV] =
            (BosonicIndex < 0 ? 0 : BosonicIndex+1) -NumberOfMatsubaras;

        for(long NuIndexM=0; NuIndexM<FermionicMatrixSize; ++NuIndexM)
        for(long NupIndexM=0; NupIndexM<FermionicMatrixSize; ++NupIndexM){
            long MatsubaraNumber1 = NuIndexM+FermionicIndexOffset[BosonicIndexV];
            long MatsubaraNumber2 = BosonicIndex - MatsubaraNumber1;
            long MatsubaraNumber3 = NupIndexM+FermionicIndexOffset[BosonicIndexV];
            Values[BosonicIndexV](NuIndexM,NupIndexM) =
                pSource->value(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
        }
    }
}

template<typename SourceObject> inline
ComplexType MatsubaraContainer4<SourceObject>::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    long BosonicIndexV = MatsubaraNumber2 + MatsubaraNumber1 + 2*NumberOfMatsubaras;
    if(BosonicIndexV >=0 && BosonicIndexV <= 2*(2*NumberOfMatsubaras-1)){
        long NuIndexM = MatsubaraNumber1-FermionicIndexOffset[BosonicIndexV];
        long NupIndexM = MatsubaraNumber3-FermionicIndexOffset[BosonicIndexV];
        if(NuIndexM >= 0 && NuIndexM < Values[BosonicIndexV].rows() &&
           NupIndexM >= 0 && NupIndexM < Values[BosonicIndexV].cols())
            return Values[BosonicIndexV](NuIndexM,NupIndexM);
    }

    DEBUG("MatsubaraContainer4 at " << this << ": " <<
          "cache miss for n1 = " << MatsubaraNumber1 << ", " <<
          "n2 = " << MatsubaraNumber2 << ", " <<
          "n3 = " << MatsubaraNumber3 <<
          " (NumberOfMatsubaras = " << NumberOfMatsubaras <<
          "), fetching a raw value from " << pSource
        )

    return pSource->value(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
}

template<typename SourceObject> inline
long MatsubaraContainer4<SourceObject>::getNumberOfMatsubaras(void) const
{
    return NumberOfMatsubaras;
}

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_MATSUBARACONTAINERS_H 
