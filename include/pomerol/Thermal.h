/** \file src/Thermal.h
** \brief Thermal object (an object which has sense only for a finite temperature).
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef __INCLUDE_THERMAL_H
#define __INCLUDE_THERMAL_H

#include "Misc.h"

namespace Pomerol{

struct Thermal {
    const RealType beta;
    const ComplexType MatsubaraSpacing;

    Thermal(RealType beta);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_THERMAL_H

