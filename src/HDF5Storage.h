/** \file src/HDF5Storage.h
** \brief HDF5 Storage.
**
** \author Igor Krivenko (igor@shg.ru)
*/
#ifndef __INCLUDE_HDF5STORAGE_H
#define __INCLUDE_HDF5STORAGE_H

#include<H5Cpp.h>
#include"Misc.h"

class HDF5Storage;
class HDF5Storable {
public:
    virtual void save(H5::CommonFG* FG) const = 0;
    virtual void load(const H5::CommonFG* FG) = 0;
};

class HDF5Storage : public H5::H5File {
    static const unsigned int* HDF5Version;
    static const unsigned int* initHDF5(void);

    static const H5::CompType& ComplexCType(void);
public:
    HDF5Storage(const std::string& FileName);
    ~HDF5Storage(void);

    void save(const HDF5Storable& Object);
    void load(HDF5Storable& Object) const;

    static void saveReal(H5::CommonFG& FG, const std::string& Name, RealType x);
    static void saveComplex(H5::CommonFG& FG, const std::string& Name, ComplexType C);
    static void saveRealVector(H5::CommonFG& FG, const std::string& Name, const RealVectorType& V);

    static RealType loadReal(const H5::CommonFG& FG, const std::string& Name);
    static void loadRealVector(const H5::CommonFG& FG, const std::string& Name, RealVectorType& V);
};

#endif // endif :: #ifndef __INCLUDE_HDF5STORAGE_H