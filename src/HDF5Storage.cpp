/** \file src/Dumper.cpp
** \brief HDF5 Dumper object..
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#ifdef pomerolHDF5
#include"HDF5Storage.h"

const unsigned int* HDF5Storage::HDF5Version = HDF5Storage::initHDF5();

const unsigned int* HDF5Storage::initHDF5()
{
    unsigned int* V = new unsigned int[3];
    H5::H5Library::getLibVersion(V[0],V[1],V[2]);

    INFO("Initializing HDF5 Library (version " << V[0] << "." << V[1] << "." << V[2] << ")...")
    H5::H5Library::open();

    return V;
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

void HDF5Storage::saveReal(H5::CommonFG& FG, const std::string& Name, RealType x)
{
    H5::DataSet DataSet = FG.createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
    DataSet.write(&x,H5::PredType::NATIVE_DOUBLE);
}

RealType HDF5Storage::loadReal(const H5::CommonFG& FG, const std::string& Name)
{
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    RealType x;
    DataSet.read(&x,H5::PredType::NATIVE_DOUBLE);
    return x;
}

inline
const H5::CompType& HDF5Storage::ComplexCType(void)
{
    static H5::CompType CType(sizeof(ComplexType));
    do_once
	CType.insertMember("real",0,H5::PredType::NATIVE_DOUBLE);
	CType.insertMember("imag",sizeof(RealType),H5::PredType::NATIVE_DOUBLE);
	CType.lock();
    end_do_once
    return CType;
}

void HDF5Storage::saveComplex(H5::CommonFG& FG, const std::string& Name, ComplexType C)
{

    H5::DataSet DataSet = FG.createDataSet(Name.c_str(),ComplexCType(), H5::DataSpace());
    RealType RealImag[2] = {real(C),imag(C)};
    DataSet.write(RealImag,ComplexCType());
}

void HDF5Storage::saveRealVector(H5::CommonFG& FG, const std::string& Name, const RealVectorType& V)
{
    hsize_t Dim[1] = {V.size()};
    H5::DataSpace DataSpace(1,Dim);
    H5::DataSet DataSet = FG.createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE,DataSpace);
    DataSet.write(V.data(),H5::PredType::NATIVE_DOUBLE);
}

void HDF5Storage::loadRealVector(const H5::CommonFG& FG, const std::string& Name, RealVectorType& V)
{
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());

    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1){
#warning TODO: throw an exception?
    }

    V.resize(DataSpace.getSimpleExtentNpoints());
    DataSet.read(V.data(),H5::PredType::NATIVE_DOUBLE);
}

#endif // endif :: #ifdef pomerolHDF5
