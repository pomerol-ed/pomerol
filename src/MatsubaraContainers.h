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
    const long NumberOfMatsubaras;
    VectorType Values;

public:
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


} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_MATSUBARACONTAINERS_H 
