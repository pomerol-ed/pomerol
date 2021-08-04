/** \file include/pomerol/Thermal.h
** \brief Thermal object (an object which has sense only for a finite temperature).
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_THERMAL_H
#define POMEROL_INCLUDE_THERMAL_H

#include "Misc.hpp"

#include <cmath>

namespace Pomerol {

struct Thermal {
    const RealType beta;
    const ComplexType MatsubaraSpacing;

    explicit Thermal(RealType beta) : beta(beta), MatsubaraSpacing(I * M_PI / beta) {}
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_THERMAL_H
