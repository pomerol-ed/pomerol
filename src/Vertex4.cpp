//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


#include "Vertex4.h"
#include <Eigen/LU> 

namespace Pomerol{

Vertex4::Vertex4(TwoParticleGF& Chi4,
                 GreensFunction& G13, GreensFunction& G24,
                 GreensFunction& G14, GreensFunction& G23) :
    Thermal(Chi4.beta), ComputableObject(),
    Chi4(Chi4), G13(G13), G24(G24), G14(G14), G23(G23)
{}

void Vertex4::compute(long NumberOfMatsubaras)
{
    Storage.fill(this,NumberOfMatsubaras);
    Status = Computed;
}

ComplexType Vertex4::value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    ComplexType Value = Chi4(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);

    if(MatsubaraNumber1 == MatsubaraNumber3)
        Value += beta*  G13(MatsubaraNumber1)*G24(MatsubaraNumber2);
    if(MatsubaraNumber2 == MatsubaraNumber3)
        Value -= beta*  G14(MatsubaraNumber1)*G23(MatsubaraNumber2);

    return Value;
}

ComplexType Vertex4::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    //if(isVanishing())
    //    return 0.0;
    //else
        return Storage(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
}

bool Vertex4::isVanishing(void) const
{
    // TODO: We need a smarter mechanism to detect this.
    return false;
}

} // end of namespace Pomerol
