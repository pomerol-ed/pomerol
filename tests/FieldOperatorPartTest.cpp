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
#include "FieldOperatorPart.h"
#include "Logger.h"

using namespace Pomerol;

int main(int argc, char* argv[])
{

    Log.setDebugging(true);

    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));
    L.addSite(new Lattice::Site("B",1,2));
    LatticePresets::addCoulombS(&L, "A", 1.0, -0.5);
    LatticePresets::addCoulombS(&L, "B", 2.0, -1.0);
    LatticePresets::addHopping(&L, "A", "B", -1.0);

    IndexClassification Indices(L.getSiteMap());
    Indices.prepare();

    IndexHamiltonian Storage(&L,Indices);
    Storage.prepare();

    Symmetrizer Symm(Indices, Storage);
    Symm.compute();

    StatesClassification S(Indices,Symm);
    S.compute();

    BlockNumber test_block=4;
    ParticleIndex op_index=3;
    QuantumNumbers QN1( S.getQuantumNumbers(test_block));
    OperatorPresets::Cdag op(op_index);
    // This returns the quantum number of the block, which is a result of the cdag acting on the block test_block. 
    QuantumNumbers QN2( S.getQuantumNumbers(S.getBlockNumber(op.actRight(S.getFockState(test_block,0)).begin()->first)));
    BlockNumber result_block = S.getBlockNumber(QN2);
    INFO( "Acting with rotated cdag_" << op_index << " on block " << QN1 << " and receiving " << QN2);

    HamiltonianPart HpartRHS(Indices, Storage, S, test_block);
    HpartRHS.prepare();
    HpartRHS.diagonalize();
    HpartRHS.print_to_screen();

    HamiltonianPart HpartLHS(Indices, Storage, S, result_block);
    HpartLHS.prepare();
    HpartLHS.diagonalize();
    HpartLHS.print_to_screen();

    CreationOperatorPart Cdag1 (Indices, S, HpartRHS, HpartLHS, op_index); 
    Cdag1.compute();
    Cdag1.print_to_screen();

    ColMajorMatrixType cmatrix(2,4);
    cmatrix.coeffRef(0,0) = 0.67513198;
    cmatrix.coeffRef(1,0) = 0.21023036;
    cmatrix.coeffRef(0,1) = -0.43516215;
    cmatrix.coeffRef(1,1) = -0.55734541;
    cmatrix.coeffRef(0,2) = 0.55734541;
    cmatrix.coeffRef(1,2) = -0.43516215;
    cmatrix.coeffRef(0,3) = 0.21023036;
    cmatrix.coeffRef(1,3) = -0.67513198;

    ColMajorMatrixType cmatrix_result=Cdag1.getColMajorValue();
    # warning no good test condition for eigenfuctions defined with an arbitrary phase.
    //if ( std::abs(cmatrix_result.sum()) - std::abs(cmatrix.sum()) > 1e-6) return EXIT_FAILURE;

    // Check transposition

    AnnihilationOperatorPart C1 ( Cdag1.transpose() );
    AnnihilationOperatorPart C2 ( Indices, S, HpartLHS, HpartRHS, op_index );
    C1.compute(); // makes nothing
    C2.compute(); 
    
    if ( std::abs ((C1.getRowMajorValue() - C2.getRowMajorValue()).sum()) > 1e-6) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

