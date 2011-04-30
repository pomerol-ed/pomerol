/** \file src/Thermal.h
** \brief Thermal object (an object which has sense only for a finite temperature).
**
** \author Igor Krivenko (igor@shg.ru)
*/
#ifndef __INCLUDE_THERMAL_H
#define __INCLUDE_THERMAL_H

#include "Misc.h"
#include "HDF5Storage.h"

struct Thermal {
    const RealType beta;

    Thermal(RealType beta);
};

#endif // endif :: #ifndef __INCLUDE_THERMAL_H

