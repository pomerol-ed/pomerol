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
#include "Index.h"
#include "IndexClassification.h"
#include "Operator.h"
#include "OperatorPresets.h"
#include <boost/shared_ptr.hpp>

using namespace Pomerol;

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;

    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));
    L.addSite(new Lattice::Site("B",1,2));
    L.addSite(new Lattice::Site("C",1,2));
    

    IndexClassification Indices(L.getSiteMap());
    Indices.prepare();

    ParticleIndex IndexSize = Indices.getIndexSize();
    OperatorPresets::Cdag Cdag_op(3);
    Operator *O = &Cdag_op;

    FockState ket(IndexSize,2);
    FockState ket2(IndexSize,2);
    std::map<FockState, MelemType> map1=Cdag_op.actRight(ket);
    for ( std::map<FockState, MelemType>::iterator it1=map1.begin(); it1!=map1.end(); it1++) {
        FockState bra = it1->first;
        MelemType Value= it1->second; 
        INFO("<" << bra << "| c^+_3 |" << ket << "> = " << Value ); 
        ket2=bra;
        }

    map1=O->actRight(ket);
    for ( std::map<FockState, MelemType>::iterator it1=map1.begin(); it1!=map1.end(); it1++) {
        FockState bra = it1->first;
        MelemType Value= it1->second; 
        INFO("<" << bra << "| c^+_3 |" << ket << "> = " << Value ); 
        ket2=bra;
        }

    ket=ket2;
    OperatorPresets::C C_op(1);
    map1=C_op.actRight(ket);
    for ( std::map<FockState, MelemType>::iterator it1=map1.begin(); it1!=map1.end(); it1++) {
        FockState bra = it1->first;
        MelemType Value= it1->second; 
        INFO("<" << bra << "| c_1 |" << ket << "> = " << Value ); 
        }

    return EXIT_SUCCESS;
}

