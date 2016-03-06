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


/** \file src/HDF5Storage.cpp
** \brief HDF5 Storage.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "pomerol/HDF5Storage.h"

namespace Pomerol{

const unsigned int* HDF5Storage::HDF5Version = HDF5Storage::initHDF5();

const unsigned int* HDF5Storage::initHDF5()
{
    unsigned int* V = new unsigned int[3];
    H5::H5Library::getLibVersion(V[0],V[1],V[2]);

    //info() << "Initializing HDF5 Library (version " << V[0] << "." << V[1] << "." << V[2] << ")...";
    H5::H5Library::open();
    //H5::Exception::dontPrint();

    return V;
}

// Choose H5::PredType corresponding to RealType
#ifdef REALTYPE_FLOAT
    #define H5_REAL_TYPE	H5::PredType::NATIVE_FLOAT
#endif
#ifdef REALTYPE_DOUBLE
    #define H5_REAL_TYPE	H5::PredType::NATIVE_DOUBLE
#endif
#ifdef REALTYPE_LDOUBLE
    #define H5_REAL_TYPE	H5::PredType::NATIVE_LDOUBLE
#endif
#ifndef H5_REAL_TYPE
    #error We do not know how to choose an HDF5 DataType for this RealType.
#endif

const H5::CompType HDF5Storage::initCompexDataType()
{
    H5::CompType Type(sizeof(ComplexType));
    // FIXME: we should NOT rely on the internal structure of a complex type.
    // But currently there is no choice (perhaps HDF5 1.10 with native complex datatypes
    // will change things for the better).
    Type.insertMember("real",0,H5_REAL_TYPE);
    Type.insertMember("imag",sizeof(RealType),H5_REAL_TYPE);

    try {
      Type.commit(*this,"complex");
    } catch(H5::DataTypeIException){};

    return Type;
}

bool HDF5Storage::fileExists(const std::string& FileName)
{
    try {
      return isHdf5(FileName.c_str());
    } catch(H5::FileIException){
      return false;
    }
}

HDF5Storage::HDF5Storage(const std::string& FileName) :
H5File(FileName.c_str(), fileExists(FileName) ? H5F_ACC_RDWR : H5F_ACC_TRUNC),
ComplexDataType(initCompexDataType())
{
    INFO("Opened HDF5 file " << FileName)
}

HDF5Storage::~HDF5Storage(void)
{
    INFO("Closed HDF5 file " << getFileName())
    close();
}

void HDF5Storage::save(const HDF5Storable& Object)
{
    Object.save(this);
    flush(H5F_SCOPE_LOCAL);
}

void HDF5Storage::load(HDF5Storable& Object) const
{
    Object.load(this);
}

void HDF5Storage::saveInt(H5::CommonFG* FG, const std::string& Name, int x)
{
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),H5::PredType::NATIVE_INT,H5::DataSpace());
    DataSet.write(&x,H5::PredType::NATIVE_INT);
}

int HDF5Storage::loadInt(const H5::CommonFG* FG, const std::string& Name)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());
    int x;
    DataSet.read(&x,H5::PredType::NATIVE_INT);
    return x;
}



// Save/load RealType
void HDF5Storage::saveReal(H5::CommonFG* FG, const std::string& Name, RealType x)
{
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),H5_REAL_TYPE,H5::DataSpace());
    DataSet.write(&x,H5_REAL_TYPE);
}

RealType HDF5Storage::loadReal(const H5::CommonFG* FG, const std::string& Name)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());
    RealType x;
    DataSet.read(&x,H5_REAL_TYPE);
    return x;
}

// Save/load ComplexType
void HDF5Storage::saveComplex(H5::CommonFG* FG, const std::string& Name, ComplexType C)
{
    H5::CompType ComplexDataType = FG->openCompType("complex");
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),ComplexDataType,H5::DataSpace());
    RealType RealImag[2] = {real(C),imag(C)};
    DataSet.write(RealImag,ComplexDataType);
}

ComplexType HDF5Storage::loadComplex(const H5::CommonFG* FG, const std::string& Name)
{
    H5::CompType ComplexDataType = FG->openCompType("complex");
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());
    ComplexType C;
    RealType RealImag[2];
    DataSet.read(RealImag,ComplexDataType);
    return ComplexType(RealImag[0],RealImag[1]);
}

// Save/load RealVectorType
void HDF5Storage::saveRealVector(H5::CommonFG* FG, const std::string& Name, const RealVectorType& V)
{
    hsize_t Dim[1] = {hsize_t(V.size())};
    H5::DataSpace DataSpace(1,Dim);
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),H5_REAL_TYPE,DataSpace);
    DataSet.write(V.data(),H5_REAL_TYPE);
}

void HDF5Storage::loadRealVector(const H5::CommonFG* FG, const std::string& Name, RealVectorType& V)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());

    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1)
	throw(H5::DataSpaceIException("HDF5Storage::loadRealVector()","Unexpected multidimentional dataspace."));

    V.resize(DataSpace.getSimpleExtentNpoints());
    DataSet.read(V.data(),H5_REAL_TYPE);
}

// Save/load RealMatrixType
void HDF5Storage::saveRealMatrix(H5::CommonFG* FG, const std::string& Name, const RealMatrixType& M)
{
    hsize_t Dims[2] = {hsize_t(M.rows()),hsize_t(M.cols())};
    H5::DataSpace DataSpace(2,Dims);
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),H5_REAL_TYPE,DataSpace);
    DataSet.write(M.data(),H5_REAL_TYPE);
}

void HDF5Storage::loadRealMatrix(const H5::CommonFG* FG, const std::string& Name, RealMatrixType& M)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());

    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 2)
	throw(H5::DataSpaceIException("HDF5Storage::loadRealMatrix()",
				      "A dataspace must be precisely two-dimensional."));

    hsize_t Dims[2];
    DataSpace.getSimpleExtentDims(Dims);

    M.resize(Dims[0],Dims[1]);
    DataSet.read(M.data(),H5_REAL_TYPE);
}

// Save/load MatrixType
void HDF5Storage::saveMatrix(H5::CommonFG* FG, const std::string& Name, const MatrixType& M)
{
    H5::CompType ComplexDataType = FG->openCompType("complex");
    hsize_t Dims[2] = {hsize_t(M.rows()),hsize_t(M.cols())};
    H5::DataSpace DataSpace(2,Dims);
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),ComplexDataType,DataSpace);
    DataSet.write(M.data(),ComplexDataType);
}

void HDF5Storage::loadMatrix(const H5::CommonFG* FG, const std::string& Name, MatrixType& M)
{
    H5::CompType ComplexDataType = FG->openCompType("complex");  
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());

    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 2)
	throw(H5::DataSpaceIException("HDF5Storage::loadMatrix()","A dataspace must be precisely two-dimensional."));

    hsize_t Dims[2];
    DataSpace.getSimpleExtentDims(Dims);

    M.resize(Dims[0],Dims[1]);
    DataSet.read(M.data(),ComplexDataType);
}

// Save/load ColMajorMatrix
void HDF5Storage::saveColMajorMatrix(H5::CommonFG* FG, const std::string& Name, const ColMajorMatrixType& CMSM)
{
    // WARNING: This method uses undocumented methods of SparseMatrix
    H5::Group Group = FG->createGroup(Name.c_str());

    // Save outer indices array
    int outerSize = CMSM.outerSize();
    hsize_t outerDim[1] = {hsize_t(outerSize)};
    H5::DataSpace outerDataSpace(1,outerDim);
    Group.createDataSet("outerIndex",H5::PredType::NATIVE_INT,outerDataSpace)
	.write(CMSM.outerIndexPtr(),H5::PredType::NATIVE_INT);

    // Save inner size
    int innerSize = CMSM.innerSize();
    Group.createDataSet("innerSize",H5::PredType::NATIVE_INT,H5::DataSpace())
	.write(&innerSize,H5::PredType::NATIVE_INT);

    // Save data
    hsize_t innerDim[1] = {hsize_t(CMSM.nonZeros())};
    H5::DataSpace innerDataSpace(1,innerDim);
    Group.createDataSet("innerIndex",H5::PredType::NATIVE_INT,innerDataSpace)
	.write(CMSM.innerIndexPtr(),H5::PredType::NATIVE_INT);
    Group.createDataSet("values",H5_REAL_TYPE,innerDataSpace)
	.write(CMSM.valuePtr(),H5_REAL_TYPE);
}

void HDF5Storage::loadColMajorMatrix(H5::CommonFG* FG, const std::string& Name, ColMajorMatrixType& CMSM)
{
    // WARNING: This method uses undocumented methods of SparseMatrix
    H5::Group Group = FG->openGroup(Name.c_str());

    H5::DataSet outerIndexDataSet = Group.openDataSet("outerIndex");
    H5::DataSpace outerIndexDataSpace = outerIndexDataSet.getSpace();
    if(outerIndexDataSpace.getSimpleExtentNdims() != 1)
	throw(H5::DataSpaceIException("HDF5Storage::loadColMajorMatrix()",
				      "Unexpected multidimentional dataspace."));
    int outerSize = outerIndexDataSpace.getSimpleExtentNpoints();

    int innerSize;
    Group.openDataSet("innerSize").read(&innerSize,H5::PredType::NATIVE_INT);;

    CMSM.resize(innerSize,outerSize);

    outerIndexDataSet.read(CMSM.outerIndexPtr(),H5::PredType::NATIVE_INT);

    H5::DataSet innerIndexDataSet = Group.openDataSet("innerIndex");
    H5::DataSet valuesDataSet = Group.openDataSet("values");
    H5::DataSpace innerIndexDataSpace = innerIndexDataSet.getSpace();
    H5::DataSpace valuesDataSpace = valuesDataSet.getSpace();
    if(	innerIndexDataSpace.getSimpleExtentNdims() != 1 ||
	valuesDataSpace.getSimpleExtentNdims() != 1)
	throw(H5::DataSpaceIException("HDF5Storage::loadColMajorMatrix()",
				      "Unexpected multidimentional dataspace."));
    int innerIndexPoints = innerIndexDataSpace.getSimpleExtentNpoints();
    int valuesPoints = valuesDataSpace.getSimpleExtentNpoints();
    if(innerIndexPoints != valuesPoints)
	throw(H5::DataSpaceIException("HDF5Storage::loadColMajorMatrix()",
				      "innerIndex and values arrays must have the same number of elements."));
    if(innerIndexPoints > innerSize*outerSize)
	throw(H5::DataSpaceIException("HDF5Storage::loadColMajorMatrix()",
				      "Number of nonzero elements must not exceed "
				      "innerSize*outerSize."));

    CMSM.resizeNonZeros(valuesPoints);
    innerIndexDataSet.read(CMSM.innerIndexPtr(),H5::PredType::NATIVE_INT);
    valuesDataSet.read(CMSM.valuePtr(),H5_REAL_TYPE);
    CMSM.finalize();
}

// Save/load RowMajorMatrix
void HDF5Storage::saveRowMajorMatrix(H5::CommonFG* FG, const std::string& Name, const RowMajorMatrixType& RMSM)
{
    // WARNING: This method uses undocumented methods of SparseMatrix
    H5::Group Group = FG->createGroup(Name.c_str());

    // Save outer indices array
    int outerSize = RMSM.outerSize();
    hsize_t outerDim[1] = {hsize_t(outerSize)};
    H5::DataSpace outerDataSpace(1,outerDim);
    Group.createDataSet("outerIndex",H5::PredType::NATIVE_INT,outerDataSpace)
	.write(RMSM.outerIndexPtr(),H5::PredType::NATIVE_INT);

    // Save inner size
    int innerSize = RMSM.innerSize();
    Group.createDataSet("innerSize",H5::PredType::NATIVE_INT,H5::DataSpace())
	.write(&innerSize,H5::PredType::NATIVE_INT);

    // Save data
    hsize_t innerDim[1] = {hsize_t(RMSM.nonZeros())};
    H5::DataSpace innerDataSpace(1,innerDim);
    Group.createDataSet("innerIndex",H5::PredType::NATIVE_INT,innerDataSpace)
	.write(RMSM.innerIndexPtr(),H5::PredType::NATIVE_INT);
    Group.createDataSet("values",H5_REAL_TYPE,innerDataSpace)
	.write(RMSM.valuePtr(),H5_REAL_TYPE);
}

void HDF5Storage::loadRowMajorMatrix(H5::CommonFG* FG, const std::string& Name, RowMajorMatrixType& RMSM)
{
    // WARNING: This method uses undocumented methods of SparseMatrix
    H5::Group Group = FG->openGroup(Name.c_str());

    H5::DataSet outerIndexDataSet = Group.openDataSet("outerIndex");
    H5::DataSpace outerIndexDataSpace = outerIndexDataSet.getSpace();
    if(outerIndexDataSpace.getSimpleExtentNdims() != 1)
	throw(H5::DataSpaceIException("HDF5Storage::loadRowMajorMatrix()",
				      "Unexpected multidimentional dataspace."));
    int outerSize = outerIndexDataSpace.getSimpleExtentNpoints();

    int innerSize;
    Group.openDataSet("innerSize").read(&innerSize,H5::PredType::NATIVE_INT);;

    RMSM.resize(outerSize,innerSize);

    outerIndexDataSet.read(RMSM.outerIndexPtr(),H5::PredType::NATIVE_INT);

    H5::DataSet innerIndexDataSet = Group.openDataSet("innerIndex");
    H5::DataSet valuesDataSet = Group.openDataSet("values");
    H5::DataSpace innerIndexDataSpace = innerIndexDataSet.getSpace();
    H5::DataSpace valuesDataSpace = valuesDataSet.getSpace();
    if(	innerIndexDataSpace.getSimpleExtentNdims() != 1 ||
	valuesDataSpace.getSimpleExtentNdims() != 1)
	throw(H5::DataSpaceIException("HDF5Storage::loadRowMajorMatrix()",
				      "Unexpected multidimentional dataspace."));
    int innerIndexPoints = innerIndexDataSpace.getSimpleExtentNpoints();
    int valuesPoints = valuesDataSpace.getSimpleExtentNpoints();
    if(innerIndexPoints != valuesPoints)
	throw(H5::DataSpaceIException("HDF5Storage::loadRowMajorMatrix()",
				      "innerIndex and values arrays must have the same number of elements."));
    if(innerIndexPoints > innerSize*outerSize)
	throw(H5::DataSpaceIException("HDF5Storage::loadRowMajorMatrix()",
				      "Number of nonzero elements must not exceed "
				      "innerSize*outerSize."));

    RMSM.resizeNonZeros(valuesPoints);
    innerIndexDataSet.read(RMSM.innerIndexPtr(),H5::PredType::NATIVE_INT);
    valuesDataSet.read(RMSM.valuePtr(),H5_REAL_TYPE);
    RMSM.finalize();
}

} // end of namespace Pomerol
