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
#include "HamiltonianPart.h"
#include <boost/shared_ptr.hpp>

using namespace Pomerol;

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;

    

    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));
    L.addSite(new Lattice::Site("B",1,2));
    LatticePresets::addCoulombS(&L, "A", 1.0, -0.5);
    LatticePresets::addCoulombS(&L, "B", 2.0, -1.0);
    LatticePresets::addHopping(&L, "A", "B", -1.0);

    IndexClassification Indices(L.getSiteMap());
    Indices.prepare();
    Indices.printIndices();

    ParticleIndex IndexSize = Indices.getIndexSize();
    INFO(IndexSize);

    IndexHamiltonian Storage(&L,Indices);
    Storage.prepare();

    Symmetrizer Symm(Indices, Storage);
    Symm.compute();

    StatesClassification S(Indices,Symm);
    S.compute();

    QuantumNumbers Q = Symm.getQuantumNumbers();
    Q.set(0,2); Q.set(1,0);
    HamiltonianPart Hpart(Indices, Storage, S, S.getBlockNumber(Q));
    RealMatrixType hmatrix(4,4);
    hmatrix << -0.402764, 0        , -0.707107, -0.581189,
               -0.581189, 0.707107 , 0        , 0.402764,
               -0.581189, -0.707107, 0        , 0.402764,
               -0.402764, 0        , 0.707107 , -0.581189;


    Hpart.prepare();
    Hpart.print_to_screen();
    Hpart.diagonalize();
    Hpart.print_to_screen();
    RealVectorType E_calc = Hpart.getEigenValues();
    RealVectorType E_real(4);
    E_real << -2.886, -1.5, 0, 1.386;
    INFO(Hpart.getEigenValues());

    if (std::abs((E_calc.sum()-E_real.sum())) > 1e-5) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}

