#include "pomerol/Vertex4.hpp"

namespace Pomerol {

Vertex4::Vertex4(TwoParticleGF const& Chi4,
                 GreensFunction const& G13, GreensFunction const& G24,
                 GreensFunction const& G14, GreensFunction const& G23) :
    Thermal(Chi4.beta), ComputableObject(),
    Chi4(Chi4), G13(G13), G24(G24), G14(G14), G23(G23), Storage(*this)
{}

void Vertex4::compute(long NumberOfMatsubaras)
{
    if(getStatus() >= Computed) return;
    Storage.fill(NumberOfMatsubaras);
    setStatus(Computed);
}

ComplexType Vertex4::value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    ComplexType Value = Chi4(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);

    if(MatsubaraNumber1 == MatsubaraNumber3)
        Value += beta*  G13(MatsubaraNumber1)*G24(MatsubaraNumber2);
    if(MatsubaraNumber2 == MatsubaraNumber3)
        Value -= beta*  G14(MatsubaraNumber1)*G23(MatsubaraNumber2);

    return Value;
}

ComplexType Vertex4::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    return Storage(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
}

} // namespace Pomerol
