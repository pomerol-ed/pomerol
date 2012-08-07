//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
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


/** \file src/HDF5Storage.h
** \brief HDF5 Storage.
**
** \author Igor Krivenko (igor@shg.ru)
*/
#ifndef __INCLUDE_HDF5STORAGE_H
#define __INCLUDE_HDF5STORAGE_H

#include<H5Cpp.h>
#include"Misc.h"
#include"Logger.h"

namespace Pomerol{

class HDF5Storage;
class HDF5Storable {
public:
    virtual void save(H5::CommonFG* RootGroup) const = 0;
    virtual void load(const H5::CommonFG* RootGroup) = 0;
};

class HDF5Storage : public H5::H5File {
    static const unsigned int* HDF5Version;
    static const unsigned int* initHDF5(void);

    H5::CompType ComplexDataType;

    const H5::CompType initCompexDataType(void);

    static bool fileExists(const std::string& FileName);

public:
    HDF5Storage(const std::string& FileName);
    ~HDF5Storage(void);

    void save(const HDF5Storable& Object);
    void load(HDF5Storable& Object) const;

    static void saveInt(H5::CommonFG* FG, const std::string& Name, int x);
    static void saveReal(H5::CommonFG* FG, const std::string& Name, RealType x);
    static void saveComplex(H5::CommonFG* FG, const std::string& Name, ComplexType C);
    static void saveRealVector(H5::CommonFG* FG, const std::string& Name, const RealVectorType& V);
    static void saveRealMatrix(H5::CommonFG* FG, const std::string& Name, const RealMatrixType& M);
    static void saveMatrix(H5::CommonFG* FG, const std::string& Name, const MatrixType& M);
    static void saveColMajorMatrix(H5::CommonFG* FG, const std::string& Name, const ColMajorMatrixType& CMSM);
    static void saveRowMajorMatrix(H5::CommonFG* FG, const std::string& Name, const RowMajorMatrixType& RMSM);

    static int loadInt(const H5::CommonFG* FG, const std::string& Name);
    static RealType loadReal(const H5::CommonFG* FG, const std::string& Name);
    static ComplexType loadComplex(const H5::CommonFG* FG, const std::string& Name);
    static void loadRealVector(const H5::CommonFG* FG, const std::string& Name, RealVectorType& V);
    static void loadRealMatrix(const H5::CommonFG* FG, const std::string& Name, RealMatrixType& M);
    static void loadMatrix(const H5::CommonFG* FG, const std::string& Name, MatrixType& M);
    static void loadColMajorMatrix(H5::CommonFG* FG, const std::string& Name, ColMajorMatrixType& CMSM);
    static void loadRowMajorMatrix(H5::CommonFG* FG, const std::string& Name, RowMajorMatrixType& RMSM);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HDF5STORAGE_H

