#include"pomerol/Thermal.h"

namespace Pomerol{

Thermal::Thermal(RealType beta) : beta(beta), MatsubaraSpacing(I*M_PI/beta) {}

} // end of namespace Pomerol

