//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/IndexClassification.hpp
/// \brief Classification of indices of fermionic creation/annihilation operators.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_INDEXCLASSIFICATION_HPP
#define POMEROL_INCLUDE_POMEROL_INDEXCLASSIFICATION_HPP

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

/// \defgroup ED Exact diagonalization: Hilbert space, Hamiltonian, monomial operators and density matrix
///@{

/// \brief Contiguous list of operator index tuples.
///
/// This class establishes correspondence between index tuples of fermionic creation/annihilation
/// operators and values of a contiguous integer index (\ref ParticleIndex).
/// \tparam IndexTypes Types of indices carried by a single creation/annihilation operator.
template <typename... IndexTypes> class IndexClassification {
public:
    /// Tuple of indices carried by a single creation/annihilation operator.
    using IndexInfo = std::tuple<IndexTypes...>;

private:
    /// The map from operator index tuples to \ref ParticleIndex.
    std::map<IndexInfo, ParticleIndex> InfoToIndices;
    /// A reverse map from \ref ParticleIndex to the operator index tuples.
    std::vector<IndexInfo> IndicesToInfo;

public:
    /// Populate the index map by extracting all index tuples from a given
    /// \ref Operators "polynomial expression". Mapped \ref ParticleIndex values
    /// are assigned according to the order that keys (index tuples) are stored in.
    /// \tparam ScalarType Coefficient type of expression H.
    /// \param[in] H Expression to be analyzed.
    template <typename ScalarType>
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
    /// Construct an empty map.
    IndexClassification() = default;

    /// Add an operator index tuple to the map.
    /// \param[in] info Index tuple to be added.
    void addInfo(IndexInfo const& info) {
        InfoToIndices.emplace(info, 0);
        UpdateMaps();
    }

    /// Create an operator index tuple and add it to the map.
    /// \param[in] indices A pack of operator indices.
    void addInfo(IndexTypes... indices) {
        InfoToIndices.emplace(IndexInfo{indices...}, 0);
        UpdateMaps();
    }

    /// Check if a given \ref ParticleIndex has a corresponding index tuple in the map.
    /// \param[in] in Particle index to check.
    bool checkIndex(ParticleIndex in) const { return in < InfoToIndices.size(); }

    /// Return the \ref ParticleIndex corresponding to a given operator index tuple.
    /// \param[in] info Index tuple.
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

    /// Return the \ref ParticleIndex corresponding to a given operator index tuple.
    /// Elements of the tuple are provided as a sequence of argument.
    /// \param[in] indices A pack of operator indices.
    ParticleIndex getIndex(IndexTypes... indices) const { return getIndex(std::make_tuple(indices...)); }

    /// Return the operator index tuple corresponding to a given \ref ParticleIndex.
    /// \param[in] in Particle index to retrieve information for.
    IndexInfo const& getInfo(ParticleIndex in) const {
        if(in >= InfoToIndices.size())
            throw std::runtime_error("Wrong particle index " + std::to_string(in));
        return IndicesToInfo[in];
    }

    /// Return the total number of elements in the map.
    ParticleIndex getIndexSize() const { return InfoToIndices.size(); }

    /// Print an entire \ref IndexClassification map into an output stream.
    /// \param[out] os Output stream.
    /// \param[in] ic The \ref IndexClassification object to be printed.
    /// \return Reference to the output stream.
    friend std::ostream& operator<<(std::ostream& os, IndexClassification const& ic) {
        for(ParticleIndex i = 0; i < ic.InfoToIndices.size(); ++i) {
            os << "Index " << i << " = (";
            libcommute::print_tuple(os, ic.IndicesToInfo[i]);
            os << ")" << std::endl;
        }
        return os;
    }

private:
    /// Re-enumerate all index tuples stored in the map.
    void UpdateMaps() {
        IndicesToInfo.clear();
        IndicesToInfo.reserve(InfoToIndices.size());
        for(auto& ii : InfoToIndices) {
            ii.second = IndicesToInfo.size();
            IndicesToInfo.emplace_back(ii.first);
        }
    }
};

/// A factory function for \ref IndexClassification, which populates the index map by
/// extracting all index tuples from a given \ref Operators "polynomial expression".
/// \tparam ScalarType Coefficient type of expression H.
/// \tparam IndexTypes Types of indices carried by a single creation/annihilation operator.
/// \param[in] H Expression to be analyzed.
template <typename ScalarType, typename... IndexTypes>
IndexClassification<IndexTypes...> MakeIndexClassification(Operators::expression<ScalarType, IndexTypes...> const& H) {
    return IndexClassification<IndexTypes...>(H);
}

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_INDEXCLASSIFICATION_HPP
