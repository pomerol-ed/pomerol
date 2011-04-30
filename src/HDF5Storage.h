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
    virtual void save(H5::CommonFG* RootGroup, HDF5Storage const* const Storage) const = 0;
    virtual void load(const H5::CommonFG* RootGroup, HDF5Storage const* const Storage) = 0;
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

    void saveReal(H5::CommonFG* FG, const std::string& Name, RealType x) const;
    void saveComplex(H5::CommonFG* FG, const std::string& Name, ComplexType C) const;
    void saveRealVector(H5::CommonFG* FG, const std::string& Name, const RealVectorType& V) const;
    void saveRealMatrix(H5::CommonFG* FG, const std::string& Name, const RealMatrixType& M) const;
    void saveMatrix(H5::CommonFG* FG, const std::string& Name, const MatrixType& M) const;

    RealType loadReal(const H5::CommonFG* FG, const std::string& Name) const;
    ComplexType loadComplex(const H5::CommonFG* FG, const std::string& Name) const;
    void loadRealVector(const H5::CommonFG* FG, const std::string& Name, RealVectorType& V) const;
    void loadRealMatrix(const H5::CommonFG* FG, const std::string& Name, RealMatrixType& M) const;
    void loadMatrix(const H5::CommonFG* FG, const std::string& Name, MatrixType& M) const;
};

#endif // endif :: #ifndef __INCLUDE_HDF5STORAGE_H

