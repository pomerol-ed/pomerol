/** \file src/Dumper.h
** \brief HDF5 Dumper object.
**
** \author Igor Krivenko (igor@shg.ru)
*/
#ifdef pomerolHDF5
#ifndef __INCLUDE_DUMPER_H
#define __INCLUDE_DUMPER_H

#include<H5Cpp.h>
#include"Misc.h"

class Dumper;
class Dumpable {
public:
    virtual void dumpIt(H5::CommonFG* FG) const = 0;
};

class Dumper : public H5::H5File {
    static const unsigned int* HDF5Version;
    static const unsigned int* initHDF5(void);
    
    static const H5::CompType& ComplexCType(void);
public:
    Dumper(std::string FileName);
    ~Dumper(void);
    
    void dump(const Dumpable& Object);
    
    static void dumpReal(H5::CommonFG& FG, const std::string& Name, RealType x);
    static void dumpComplex(H5::CommonFG& FG, const std::string& Name, ComplexType C);
    static void dumpRealVector(H5::CommonFG& FG, const std::string& Name, RealVectorType V);
};

#endif // endif :: #ifndef __INCLUDE_DUMPER_H
#endif // endif :: #ifndef pomerolHDF5
