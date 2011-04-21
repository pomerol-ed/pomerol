/** \file src/Dumper.cpp
** \brief HDF5 Dumper object..
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#ifdef pomerolHDF5
#include "Dumper.h"
#include <complex>

const unsigned int* Dumper::HDF5Version = Dumper::initHDF5();

const unsigned int* Dumper::initHDF5()
{
    unsigned int* V = new unsigned int[3];
    H5::H5Library::getLibVersion(V[0],V[1],V[2]);

    INFO("Initializing HDF5 Library (version " << V[0] << "." << V[1] << "." << V[2] << ")...")
    H5::H5Library::open();
    
    return V;
}

Dumper::Dumper(std::string FileName) : H5File(FileName.c_str(),H5F_ACC_EXCL)
{  
    INFO("Opened HDF5 file " << FileName)
}

Dumper::~Dumper(void )
{
    INFO("Closed HDF5 file " << getFileName())
    close();
}

void Dumper::dump(const Dumpable& Object)
{
    Object.dumpIt(this);
    flush(H5F_SCOPE_LOCAL);
}

void Dumper::dumpReal(H5::CommonFG& FG, const std::string& Name, RealType x)
{
    H5::DataSet DataSet = FG.createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
    DataSet.write(&x,H5::PredType::NATIVE_DOUBLE);
}

inline
const H5::CompType& Dumper::ComplexCType(void)
{
    static H5::CompType CType(sizeof(ComplexType));
    do_once
#warning Dirty hack: we rely on the internal structure of std::complex which we do not know and must not know.\
Must look for a better solution...
	CType.insertMember("real",0,H5::PredType::NATIVE_DOUBLE);
	CType.insertMember("imag",sizeof(RealType),H5::PredType::NATIVE_DOUBLE);
    end_do_once
    return CType;
}

void Dumper::dumpComplex(H5::CommonFG& FG, const std::string& Name, ComplexType C)
{  
    
    H5::DataSet DataSet = FG.createDataSet(Name.c_str(),ComplexCType(), H5::DataSpace());
    DataSet.write(&C,ComplexCType());
}

void Dumper::dumpRealVector(H5::CommonFG& FG, const std::string& Name, RealVectorType V)
{
    hsize_t Dim[1] = {V.size()};
    H5::DataSpace DataSpace(1,Dim);
    H5::DataSet DataSet = FG.createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE,DataSpace);
    DataSet.write(V.data(),H5::PredType::NATIVE_DOUBLE);
}

#endif // endif :: #ifdef pomerolHDF5
