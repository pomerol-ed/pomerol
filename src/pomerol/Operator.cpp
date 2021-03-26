#include "pomerol/Operator.h"
#include <algorithm>
#include <iterator>

namespace Pomerol{

template<bool Complex>
const char* Operator<Complex>::exWrongLabel::what() const throw(){
    return "Wrong labels";
};

template<bool Complex>
const char* Operator<Complex>::exMelemVanishes::what() const throw(){
    return "Matrix element vanishes";
};

template<bool Complex>
bool Operator<Complex>::commutes(const Operator &rhs) const
{
    return ( Operator((*this)*rhs) == Operator(rhs*(*this)));
}

template<bool Complex>
Operator<Complex> Operator<Complex>::getCommutator(const Operator &rhs) const
{
    return (*this)*rhs - rhs*(*this);
}

template<bool Complex>
Operator<Complex> Operator<Complex>::getAntiCommutator(const Operator &rhs) const
{
    return (*this)*rhs + rhs*(*this);
}

template<bool Complex>
auto Operator<Complex>::actRight(const monomial_t &in, const FockState &ket) -> std::tuple<FockState, MelemT>
{
    if (in.size()==0) return std::make_tuple(ket, MelemT(1));
    //DEBUG(in << "|" << ket << ">");
    //ParticleIndex prev_pos_ = ket.size(); // Here we'll store the index of the last operator to speed up sign counting
    ParticleIndex prev_pos_ = 0;
    //DEBUG(prev_pos_ <<"|" << ket);
    int sign=1;
    FockState bra = ket;
    unsigned int N=in.size();
    for (int i=N-1; i>=0; i--) // Is the number of operator in OperatorSequence. Now we need to count them from back.
        {
            bool op; ParticleIndex ind;
            std::tie(op,ind)=in[i];
            //DEBUG(op << "_" << ind);
            if ((op == creation && bra[ind]) || (op == annihilation && !bra[ind]) ) return std::make_tuple(ERROR_FOCK_STATE, 0); // This is Pauli principle.
            if (ind > prev_pos_)
                for (ParticleIndex j=prev_pos_; j<ind; ++j) { if (bra[j]) sign*=-1; }
            else
                for (ParticleIndex j=prev_pos_; j>ind; j--) { if (bra[j]) sign*=-1; }
            bra[ind] = (op == creation); // This is c or c^+ acting
            //prev_pos_ = 0;

        }
    return std::make_tuple(bra, MelemT(sign));
}


template <class R>
inline bool __is_zero(std::pair<FockState,R> in){return (std::abs(in.second)<std::numeric_limits<RealType>::epsilon());};

template<bool Complex>
auto Operator<Complex>::actRight(const FockState &ket) const -> std::map<FockState, MelemT>
{
    std::map<FockState, MelemT> result1;
    for (auto it = monomials.begin(); it!=monomials.end(); it++)
        {
            FockState bra;
            MelemT melem;
            std::tie(bra,melem) = actRight(it->first,ket);
            //if (std::abs(melem)>1e-8) DEBUG(bra << "|*" << melem);
            if (bra!=ERROR_FOCK_STATE && std::abs(melem)>std::numeric_limits<RealType>::epsilon())
                result1[bra]+=melem*(it->second);
        };
    // C++11 remove_if has a different behaviour, so this is a hck around. */
    std::map<FockState, MelemT> result2;
    std::remove_copy_if(result1.begin(), result1.end(), std::inserter(result2, result2.end()), __is_zero<MelemT> );
    return result2;
}

template<bool Complex>
auto Operator<Complex>::getMatrixElement( const FockState & bra, const FockState &ket) const -> MelemT
{
    std::map<FockState, MelemT> output = this->actRight(ket);
    if (output.find(bra)==output.end())
        return 0;
    else {
        return output[bra];
        }
}

template<bool Complex>
auto Operator<Complex>::getMatrixElement(const VectorType<Complex> &bra,
                                         const VectorType<Complex> &ket,
                                         const std::vector<FockState> &states) const -> MelemT
{
    if (bra.size()!=ket.size() || bra.size()!=states.size()) throw (exMelemVanishes());
    MelemT melem = 0.0;
    for (int i=0; i<ket.size(); ++i) {
        FockState current_state = states[i];
        MelemT overlap = ket[i];
        if (std::abs(overlap)>std::numeric_limits<RealType>::epsilon()) {
            //DEBUG(overlap << "," << current_state);
            std::map<FockState, MelemT> map1 = this->actRight(current_state);
            for (typename std::map<FockState, MelemT>::const_iterator it = map1.begin(); it!= map1.end(); it++) {
                FockState result_state = it->first;
                MelemT melem2 = it->second;
                //DEBUG("\t<" << result_state << "|" << melem2);
                std::vector<FockState>::const_iterator it1 = std::find(states.begin(), states.end(), result_state);
                MelemT overlap2;
                if (it1 != states.end() ) {
                    size_t j = std::distance(states.begin(), it1);
                    overlap2 = conj(bra(j));
                }
                else overlap2 = 0.0;
                //DEBUG(overlap2);
                //DEBUG("<" << result_state << "|" << overlap2 << "*" << melem << "*" << overlap << "|" << current_state << ">");
                melem += overlap2 * melem2 * overlap;
                }
            };
        };
    return melem;
}

template<bool Complex>
bool operator==(const typename Operator<Complex>::monomials_map_t::value_type& lhs,
                const typename Operator<Complex>::monomials_map_t::value_type& rhs)
{
    return (std::equal(lhs.first.begin(), lhs.first.end(), rhs.first.begin()) &&
            std::abs(rhs.second - lhs.second)<100*std::numeric_limits<RealType>::epsilon());
}

template<bool Complex>
bool operator==(const Operator<Complex> &lhs, const Operator<Complex> &rhs)
{
    return (lhs.monomials.size() == rhs.monomials.size() &&
           std::equal(lhs.begin(), lhs.end(), rhs.begin()));
}

// Explicit instantiations

template class Operator<false>;
template class Operator<true>;

namespace OperatorPresets {
template Operator<false> c<false>(ParticleIndex);
template Operator<false> c_dag<false>(ParticleIndex);
template Operator<false> n<false>(ParticleIndex);
template Operator<false> n_offdiag<false>(ParticleIndex, ParticleIndex);

template Operator<true> c<true>(ParticleIndex);
template Operator<true> c_dag<true>(ParticleIndex);
template Operator<true> n<true>(ParticleIndex);
template Operator<true> n_offdiag<true>(ParticleIndex, ParticleIndex);
}

} // end of namespace Pomerol
