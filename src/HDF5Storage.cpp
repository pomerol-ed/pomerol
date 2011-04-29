/** \file src/HDF5Storage.cpp
** \brief HDF5 Storage.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#include"HDF5Storage.h"

const unsigned int* HDF5Storage::HDF5Version = HDF5Storage::initHDF5();
const H5::CompType HDF5Storage::ComplexDataType = HDF5Storage::initCompexDataType();

const unsigned int* HDF5Storage::initHDF5()
{
    unsigned int* V = new unsigned int[3];
    H5::H5Library::getLibVersion(V[0],V[1],V[2]);

    INFO("Initializing HDF5 Library (version " << V[0] << "." << V[1] << "." << V[2] << ")...")
    H5::H5Library::open();

    return V;
}

const H5::CompType HDF5Storage::initCompexDataType()
{
    H5::CompType CType(sizeof(ComplexType));
    CType.insertMember("real",0,H5::PredType::NATIVE_DOUBLE);
    CType.insertMember("imag",sizeof(RealType),H5::PredType::NATIVE_DOUBLE);
    return CType;
}

HDF5Storage::HDF5Storage(const std::string& FileName) : H5File(FileName.c_str(),H5F_ACC_TRUNC)
{  
    INFO("Opened HDF5 file " << FileName)
}

HDF5Storage::~HDF5Storage(void )
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

void HDF5Storage::saveReal(H5::CommonFG* FG, const std::string& Name, RealType x)
{
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
    DataSet.write(&x,H5::PredType::NATIVE_DOUBLE);
}

void HDF5Storage::saveComplex(H5::CommonFG* FG, const std::string& Name, ComplexType C)
{

    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),ComplexDataType, H5::DataSpace());
    RealType RealImag[2] = {real(C),imag(C)};
    DataSet.write(RealImag,ComplexDataType);
}

void HDF5Storage::saveRealVector(H5::CommonFG* FG, const std::string& Name, const RealVectorType& V)
{
    hsize_t Dim[1] = {V.size()};
    H5::DataSpace DataSpace(1,Dim);
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE,DataSpace);
    DataSet.write(V.data(),H5::PredType::NATIVE_DOUBLE);
}

void HDF5Storage::saveRealMatrix(H5::CommonFG* FG, const std::string& Name, const RealMatrixType& M)
{
    hsize_t Dims[2] = {M.rows(),M.cols()};
    H5::DataSpace DataSpace(2,Dims);
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE,DataSpace);
    DataSet.write(M.data(),H5::PredType::NATIVE_DOUBLE);
}

void HDF5Storage::saveMatrix(H5::CommonFG* FG, const std::string& Name, const MatrixType& M)
{
    hsize_t Dims[2] = {M.rows(),M.cols()};
    H5::DataSpace DataSpace(2,Dims);
    H5::DataSet DataSet = FG->createDataSet(Name.c_str(),ComplexDataType,DataSpace);
    DataSet.write(M.data(),ComplexDataType);
}

RealType HDF5Storage::loadReal(const H5::CommonFG* FG, const std::string& Name)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());
    RealType x;
    DataSet.read(&x,H5::PredType::NATIVE_DOUBLE);
    return x;
}

ComplexType HDF5Storage::loadComplex(const H5::CommonFG* FG, const std::string& Name)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());
    ComplexType C;
    DataSet.read(&C,ComplexDataType);
    return C;
}

void HDF5Storage::loadRealVector(const H5::CommonFG* FG, const std::string& Name, RealVectorType& V)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());

    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1)
	throw(H5::DataSpaceIException("HDF5Storage::loadRealVector()","Unexpected multidimentional dataspace."));

    V.resize(DataSpace.getSimpleExtentNpoints());
    DataSet.read(V.data(),H5::PredType::NATIVE_DOUBLE);
}

void HDF5Storage::loadRealMatrix(const H5::CommonFG* FG, const std::string& Name, RealMatrixType& M)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());

    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 2)
	throw(H5::DataSpaceIException("HDF5Storage::loadRealMatrix()","A dataspace must be precisely two-dimensional."));

    hsize_t Dims[2];
    DataSpace.getSimpleExtentDims(Dims);

    M.resize(Dims[0],Dims[1]);
    DataSet.read(M.data(),H5::PredType::NATIVE_DOUBLE);
}

void HDF5Storage::loadMatrix(const H5::CommonFG* FG, const std::string& Name, MatrixType& M)
{
    H5::DataSet DataSet = FG->openDataSet(Name.c_str());

    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 2)
	throw(H5::DataSpaceIException("HDF5Storage::loadMatrix()","A dataspace must be precisely two-dimensional."));

    hsize_t Dims[2];
    DataSpace.getSimpleExtentDims(Dims);

    M.resize(Dims[0],Dims[1]);
    DataSet.read(M.data(),ComplexDataType);
}