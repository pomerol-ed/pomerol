//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2011 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


/** \file src/FourIndexObject.h
** \brief A prototype class for all objects, depending on four indices
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "pomerol/MatsubaraContainer.h"

namespace Pomerol{

//
// Matsubara Container
//
MatsubaraContainer::MatsubaraContainer(RealType beta):MatsubaraSpacing(I*M_PI/beta),BosonicMin_(0),BosonicMax_(-1),FermionicMin_(-1),FermionicMax_(-1){};

void MatsubaraContainer::prepare(int BosonicMin, int BosonicMax, int FermionicMin, int FermionicMax)
{
    BosonicMin_ = BosonicMin;
    BosonicMax_ = BosonicMax;
    if (BosonicMax_ < BosonicMin_) throw std::logic_error("MatsubaraContainer : BosonicMin > BosonicMax"); 
    FermionicMin_ = FermionicMin;
    FermionicMax_ = FermionicMax;
    
    int NBosonic = this->NBosonic();
    int NFermionic = this->NFermionic();

    Data.resize(std::max(0, NBosonic));
    for (int BosonicIndex= 0; BosonicIndex < Data.size(); BosonicIndex++)
    {
        Data[BosonicIndex].resize(NFermionic,NFermionic);
        Data[BosonicIndex].setZero();
    };
};

MatsubaraContainer& MatsubaraContainer::operator+= (const MatsubaraContainer& rhs)
{
    for (long BosonicIndex=0;BosonicIndex < Data.size();BosonicIndex++){
        Data[BosonicIndex]+=rhs.Data[BosonicIndex];
    }
    return (*this);
};

void MatsubaraContainer::clear()
{
    for (long BosonicIndex=0;BosonicIndex < Data.size(); BosonicIndex++){
        Data[BosonicIndex].resize(0,0);
    }
}


} // end of namespace Pomerol

