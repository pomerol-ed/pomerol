#include "pomerol/Vertex4.h"
#include <Eigen/LU>

namespace Pomerol{

template<bool Complex>
Vertex4<Complex>::Vertex4(TwoParticleGF<Complex>& Chi4,
                 GreensFunction<Complex>& G13, GreensFunction<Complex>& G24,
                 GreensFunction<Complex>& G14, GreensFunction<Complex>& G23) :
    Thermal(Chi4.beta), ComputableObject(),
    Chi4(Chi4), G13(G13), G24(G24), G14(G14), G23(G23)
{}

template<bool Complex>
void Vertex4<Complex>::compute(long NumberOfMatsubaras)
{
    Storage.fill(this,NumberOfMatsubaras);
    Status = Computed;
}

template<bool Complex>
ComplexType Vertex4<Complex>::value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    ComplexType Value = Chi4(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);

    if(MatsubaraNumber1 == MatsubaraNumber3)
        Value += beta*  G13(MatsubaraNumber1)*G24(MatsubaraNumber2);
    if(MatsubaraNumber2 == MatsubaraNumber3)
        Value -= beta*  G14(MatsubaraNumber1)*G23(MatsubaraNumber2);

    return Value;
}

template<bool Complex>
ComplexType Vertex4<Complex>::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    //if(isVanishing())
    //    return 0.0;
    //else
        return Storage(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
}

template<bool Complex>
bool Vertex4<Complex>::isVanishing(void) const
{
    // TODO: We need a smarter mechanism to detect this.
    return false;
}

template class Vertex4<false>;
template class Vertex4<true>;

} // end of namespace Pomerol
