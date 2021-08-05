/** \file IndexClassification.h
**  \brief Declaration of IndexClassification class.
**
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_INDEXCLASSIFICATION_H
#define POMEROL_INCLUDE_POMEROL_INDEXCLASSIFICATION_H

#include "Misc.hpp"
#include "Operators.hpp"

#include <libcommute/utility.hpp>

#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace Pomerol {

/** This class handles all the indices classification, it allocates the indices to particular Site+Spin+Orbital configuration.
 *  It also returns the information about current ParticleIndex on request. */
template<typename... IndexTypes> class IndexClassification {
public:
    /** A structure, which holds the site label, orbital and spin of a ParticleIndex. */
    using IndexInfo = std::tuple<IndexTypes...>;

private:
    /** A map of each ParticleIndex to the information about it. */
    std::map<IndexInfo, ParticleIndex> InfoToIndices;
    /** A vector of IndexInfo - each element corresponds to its number. */
    std::vector<IndexInfo> IndicesToInfo;

public:

    /** Constructor
     * \param[in] L A pointer to a Lattice Object.
     */
    template<typename ScalarType>
    explicit IndexClassification(Operators::expression<ScalarType, IndexTypes...> const& H) {
        // Collect indices of fermionic operators in the Hamiltonian
        for(auto const& mon : H) {
            for(auto const& g : mon.monomial) {
                if(libcommute::is_fermion(g))
                    InfoToIndices.emplace(g.indices(), 0);
            }
        }
        UpdateMaps();
    }
    IndexClassification() = default;

    void addInfo(IndexInfo const& info) {
        InfoToIndices.emplace(info, 0);
        UpdateMaps();
    }

    void addInfo(IndexTypes... indices) {
        InfoToIndices.emplace(IndexInfo{indices...}, 0);
        UpdateMaps();
    }

    /** Checks if the index belongs to the space of indices
     * \param[in] in Index to check. */
    bool checkIndex(ParticleIndex in) const { return in < InfoToIndices.size(); }

    /** Returns a ParticleIndex, which corresponds to a given site, orbital and spin. */
    ParticleIndex getIndex(IndexInfo const& info) const {
        auto it = InfoToIndices.find(info);
        if(it != InfoToIndices.end())
            return it->second;
        else {
            std::stringstream ss;
            libcommute::print_tuple(ss, info);
            throw std::runtime_error("Wrong indices " + ss.str());
        }
    }

    ParticleIndex getIndex(IndexTypes... info) const {
        return getIndex(std::make_tuple(info...));
    }

    /** Return all information about the given index. */
    IndexInfo const& getInfo(ParticleIndex in) const {
        if(in >= InfoToIndices.size())
            throw std::runtime_error("Wrong particle index " + std::to_string(in));
        return IndicesToInfo[in];
    }

    /** Returns total number of ParticleIndices. */
    ParticleIndex getIndexSize() const { return InfoToIndices.size(); }

    /** Print all Indices to the information stream */
    friend std::ostream & operator<<(std::ostream & os, IndexClassification const& ic) {
        for(ParticleIndex i=0; i<ic.InfoToIndices.size(); ++i) {
            os << "Index " << i << " = (";
            libcommute::print_tuple(os, ic.IndicesToInfo[i]);
            os << ")" << std::endl;
        }
        return os;
    }

private:

    void UpdateMaps() {
        IndicesToInfo.clear();
        IndicesToInfo.reserve(InfoToIndices.size());
        for(auto & ii : InfoToIndices) {
            ii.second = IndicesToInfo.size();
            IndicesToInfo.emplace_back(ii.first);
        }
    }
};

template<typename ScalarType, typename... IndexTypes>
IndexClassification<IndexTypes...>
MakeIndexClassification(Operators::expression<ScalarType, IndexTypes...> const& H) {
  return IndexClassification<IndexTypes...>(H);
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_INDEXCLASSIFICATION_H
