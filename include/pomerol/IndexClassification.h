/** \file IndexClassification.h
**  \brief Declaration of IndexClassification class.
**
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_INDEXCLASSIFICATION_H
#define __INCLUDE_INDEXCLASSIFICATION_H

#include "Misc.h"
#include "Operators.h"

#include <libcommute/utility.hpp>

#include <map>
#include <iostream>
#include <string>
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
    IndexClassification(Operators::expression<ScalarType, IndexTypes...> const& H) {
        for(auto const& mon : H) {
            for(auto const& g : mon.monomial)
                InfoToIndices.emplace(g.indices(), 0);
        }

        IndicesToInfo.reserve(InfoToIndices.size());
        for(auto it = InfoToIndices.begin(); it != InfoToIndices.end(); ++it) {
            it->second = IndicesToInfo.size();
            IndicesToInfo.push_back(it->first);
        }
    }
    IndexClassification() = default;

    /** Checks if the index belongs to the space of indices
     * \param[in] in Index to check. */
    bool checkIndex(ParticleIndex in) const { return in < InfoToIndices.size(); }

    /** Returns a ParticleIndex, which corresponds to a given site, orbital and spin. */
    ParticleIndex getIndex(const IndexClassification::IndexInfo& info) const {
        auto it = InfoToIndices.find(info);
        if(it != InfoToIndices.end())
            return it->second;
        else
            throw exWrongIndex();
    }

    /** Return all information about the given index. */
    IndexInfo const& getInfo(ParticleIndex in) const {
        if(in >= InfoToIndices.size()) throw exWrongIndex(in);
        return IndicesToInfo[in];
    }

    /** Returns total number of ParticleIndices. */
    const ParticleIndex getIndexSize() const { return InfoToIndices.size(); }

    /** Print all Indices to the information stream */
    void printIndices() {
        for(ParticleIndex i=0; i<InfoToIndices.size(); ++i) {
            std::cout << "Index " << i << " = (";
            libcommute::print_tuple(std::cout, IndicesToInfo[i]);
            std::cout << ")" << std::endl;
        }
    }

    /** Exception - wrong index. */
    class exWrongIndex : public std::exception {
        std::string msg;
    public:
        exWrongIndex(ParticleIndex index) :
          msg("Wrong index " + std::to_string(index)) {}
        virtual const char* what() const noexcept override {
          return msg.c_str();
        }
    };
};

template<typename ScalarType, typename... IndexTypes>
IndexClassification<IndexTypes...>
MakeIndexClassification(Operators::expression<ScalarType, IndexTypes...> const& H) {
  return IndexClassification<IndexTypes...>(H);
}

} // end of namespace Pomerol
#endif // endif :: #ifndef #__INCLUDE_INDEXCLASSIFICATION_H
