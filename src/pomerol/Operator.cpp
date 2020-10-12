#include "pomerol/Operator.h"
#include <algorithm>
#include <iterator>
//#include <boost/lambda/lambda.hpp>
//#include <boost/lambda/bind.hpp>
#include "boost/utility/swap.hpp"
#include <boost/functional/hash.hpp>

namespace Pomerol{

const char* Operator::exWrongLabel::what() const throw(){
    return "Wrong labels";
};

const char* Operator::exMelemVanishes::what() const throw(){
    return "Matrix element vanishes";
};

bool Operator::commutes(const Operator &rhs) const
{
    return ( Operator((*this)*rhs) == Operator(rhs*(*this)));
}

Operator Operator::getCommutator(const Operator &rhs) const
{
    return (*this)*rhs - rhs*(*this);
}

Operator Operator::getAntiCommutator(const Operator &rhs) const
{
    return (*this)*rhs + rhs*(*this);
}

std::tuple<FockState,MelemType> Operator::actRight(const monomial_t &in, const FockState &ket)
{
    if (in.size()==0) return std::make_tuple(ket, MelemType(1));
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
    return std::make_tuple(bra, MelemType(sign));
}


template <class R>
inline bool __is_zero(std::pair<FockState,R> in){return (std::abs(in.second)<std::numeric_limits<RealType>::epsilon());};

std::map<FockState, MelemType> Operator::actRight(const FockState &ket) const
{
    std::map<FockState, MelemType> result1;
    for (std::map<monomial_t,MelemType>::const_iterator it = monomials.begin(); it!=monomials.end(); it++)
        {
            FockState bra;
            MelemType melem;
            std::tie(bra,melem) = actRight(it->first,ket);
            //if (std::abs(melem)>1e-8) DEBUG(bra << "|*" << melem);
            if (bra!=ERROR_FOCK_STATE && std::abs(melem)>std::numeric_limits<RealType>::epsilon())
                result1[bra]+=melem*(it->second);
        };
    // C++11 remove_if has a different behaviour, so this is a hck around. */
    std::map<FockState, MelemType> result2;
    std::remove_copy_if(result1.begin(), result1.end(), std::inserter(result2, result2.end()), __is_zero<MelemType> );//__is_zero<MelemType>);
    return result2;
}


MelemType Operator::getMatrixElement( const FockState & bra, const FockState &ket) const
{
    std::map<FockState, MelemType> output = this->actRight(ket);
    if (output.find(bra)==output.end())
        return 0;
    else {
        return output[bra];
        }
}

MelemType Operator::getMatrixElement( const VectorType & bra, const VectorType &ket, const std::vector<FockState> &states) const
{
    if (bra.size()!=ket.size() || bra.size()!=states.size()) throw (exMelemVanishes());
    MelemType melem = 0.0;
    for (int i=0; i<ket.size(); ++i) {
        FockState current_state = states[i];
        MelemType overlap = ket[i];
        if (std::abs(overlap)>std::numeric_limits<RealType>::epsilon()) {
            //DEBUG(overlap << "," << current_state);
            std::map<FockState, MelemType> map1 = this->actRight(current_state);
            for (std::map<FockState, MelemType>::const_iterator it = map1.begin(); it!= map1.end(); it++) {
                FockState result_state = it->first;
                MelemType melem2 = it->second;
                //DEBUG("\t<" << result_state << "|" << melem2);
                std::vector<FockState>::const_iterator it1 = std::find(states.begin(), states.end(), result_state);
                MelemType overlap2;
                if (it1 != states.end() ) {
                    size_t j = std::distance(states.begin(), it1);
                #ifdef POMEROL_COMPLEX_MATRIX_ELEMENTS
                    overlap2 = std::conj(bra(j));
                #else
                    overlap2 = bra(j);
                #endif
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

bool operator==(const Operator::monomials_map_t::value_type& lhs, const Operator::monomials_map_t::value_type& rhs)
{
    return (std::equal(lhs.first.begin(), lhs.first.end(), rhs.first.begin()) && std::abs(rhs.second - lhs.second)<100*std::numeric_limits<RealType>::epsilon());
}

bool operator==(const Operator &lhs, const Operator &rhs)
{
    return (lhs.monomials.size() == rhs.monomials.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin()));
}

} // end of namespace Pomerol
