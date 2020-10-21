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
#include "Lattice.h"
#include "LatticePresets.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Operator.h"
#include "OperatorPresets.h"
#include "IndexHamiltonian.h"
#include "Symmetrizer.h"
#include "StatesClassification.h"
#include "Hamiltonian.h"
#include "FieldOperator.h"

using namespace Pomerol;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));
    L.addSite(new Lattice::Site("B",1,2));
    LatticePresets::addCoulombS(&L, "A", 1.0, -0.5);
    LatticePresets::addCoulombS(&L, "B", 2.0, -1.0);
    LatticePresets::addHopping(&L, "A", "B", -1.0);

    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare();

    IndexHamiltonian Storage(&L,IndexInfo);
    Storage.prepare();

    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute();

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    Hamiltonian H(IndexInfo, Storage, S);
    H.prepare();
    H.compute(MPI_COMM_WORLD);

    ParticleIndex op_index=3;
    CreationOperator Cdag(IndexInfo, S, H, op_index);
    Cdag.prepare();
    Cdag.compute();

    DEBUG(Cdag.getParts()[1]->getRowMajorValue());

    AnnihilationOperator C(IndexInfo, S, H, op_index);
    C.prepare();
    C.compute();

    DEBUG(C.getParts()[1]->getColMajorValue());

    DEBUG((C.getParts()[1]->getColMajorValue().transpose() - Cdag.getParts()[1]->getRowMajorValue()).norm());

    MPI_Finalize();
    return EXIT_SUCCESS;
}

