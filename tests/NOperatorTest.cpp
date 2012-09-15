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

/** \file tests/IndexTerm.cpp
** \brief Test of the Symmetrizer::IndexPermutation.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "Misc.h"
#include "Logger.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Operator.h"
#include "OperatorPresets.h"
#include "Logger.h"
#include <boost/shared_ptr.hpp>

using namespace Pomerol;

int main(int argc, char* argv[])
{

    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));
    L.addSite(new Lattice::Site("B",1,2));
    L.addSite(new Lattice::Site("C",1,2));
    Log.setDebugging(true);

    IndexClassification Indices(L.getSiteMap());
    Indices.prepare();

    ParticleIndex IndexSize = Indices.getIndexSize();
    OperatorPresets::N N(IndexSize);
    Operator *NN = &N;

    FockState ket(IndexSize,3);
    auto map1=NN->actRight(ket);
    for ( auto it1=map1.begin(); it1!=map1.end(); it1++) {
        FockState bra = it1->first;
        MelemType Value= it1->second; 
        INFO("<" << bra << "| N |" << ket << "> = " << Value ); 
        }
    if ( map1.size()!=1 && N.getMatrixElement(ket,ket)!=MelemType(2)) return EXIT_FAILURE;
    ket = FockState(IndexSize, 7);
    INFO ( "<" << ket << "| N |" << ket << "> = " << N.getMatrixElement(ket));
    if ( std::abs(N.getMatrixElement(ket) - MelemType(3)) >std::numeric_limits<double>::epsilon()) return EXIT_FAILURE;
    ket = FockState(IndexSize, 8);
    INFO ( "<" << ket << "| N |" << ket << "> = " << N.getMatrixElement(ket));
    if ( std::abs(N.getMatrixElement(ket) - MelemType(1)) >std::numeric_limits<double>::epsilon()) return EXIT_FAILURE;
    ket = FockState(IndexSize, 10);
    INFO ( "<" << ket << "| N |" << ket << "> = " << N.getMatrixElement(ket));
    INFO ( "<" << ket << "| N |" << ket << "> = " << NN->getMatrixElement(ket,ket));
    {
    };



    return EXIT_SUCCESS;
}

