//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file Symmetrizer.h
**  \brief Declaration of the Symmetrizer class - a class to store and get the information about the symmetries of the system.
**
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_HILBERTSPACE_HPP
#define POMEROL_INCLUDE_POMEROL_HILBERTSPACE_HPP

#include "ComputableObject.hpp"
#include "IndexClassification.hpp"
#include "Misc.hpp"
#include "Operators.hpp"

#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/loperator.hpp>
#include <libcommute/loperator/space_partition.hpp>

#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace Pomerol {

template <typename... IndexTypes> class HilbertSpace : public ComputableObject {

public:
    using FullHilbertSpaceType = libcommute::hilbert_space<IndexTypes...>;
    using SpacePartitionType = libcommute::space_partition;

private:
    IndexClassification<IndexTypes...> const& IndexInfo;

    FullHilbertSpaceType FullHilbertSpace;

    bool HamiltonianComplex;

    std::shared_ptr<void> HOp;

    std::unique_ptr<SpacePartitionType> Partition = nullptr;

    struct boson_es_constructor_from_map {
        std::map<std::tuple<IndexTypes...>, unsigned int> const& bits_per_boson;

        explicit boson_es_constructor_from_map(std::map<std::tuple<IndexTypes...>, unsigned int> const& bits_per_boson)
            : bits_per_boson(bits_per_boson) {}

        inline std::unique_ptr<libcommute::elementary_space<IndexTypes...>>
        operator()(libcommute::generator<IndexTypes...> const& g) const {
            using namespace libcommute;
            if(is_fermion(g)) {
                return make_unique<elementary_space_fermion<IndexTypes...>>(g.indices());
            } else if(is_boson(g)) {
                auto it = bits_per_boson.find(g.indices());
                unsigned int bits = (it == bits_per_boson.end()) ? 1 : it->second;
                return make_unique<elementary_space_boson<IndexTypes...>>(bits, g.indices());
            } else {
                std::stringstream ss;
                ss << "Unexpected algebra generator " << g;
                throw std::runtime_error(ss.str());
            }
        }
    };

public:
    template <typename ScalarType>
    HilbertSpace(IndexClassification<IndexTypes...> const& IndexInfo,
                 Operators::expression<ScalarType, IndexTypes...> const& Hamiltonian,
                 unsigned int bits_per_boson = 1)
        : IndexInfo(IndexInfo),
          FullHilbertSpace(Hamiltonian, libcommute::boson_es_constructor(bits_per_boson)),
          HamiltonianComplex(std::is_same<ScalarType, ComplexType>::value),
          HOp(std::make_shared<LOperatorType<ScalarType>>(Hamiltonian, FullHilbertSpace)) {}

    template <typename ScalarType>
    HilbertSpace(IndexClassification<IndexTypes...> const& IndexInfo,
                 Operators::expression<ScalarType, IndexTypes...> const& Hamiltonian,
                 std::map<std::tuple<IndexTypes...>, unsigned int> const& bits_per_boson_map)
        : IndexInfo(IndexInfo),
          FullHilbertSpace(Hamiltonian, boson_es_constructor_from_map(bits_per_boson_map)),
          HamiltonianComplex(std::is_same<ScalarType, ComplexType>::value),
          HOp(std::make_shared<LOperatorType<ScalarType>>(Hamiltonian, FullHilbertSpace)) {}

    void compute() {
        if(Status >= Computed)
            return;

        // Phase I of auto-partition algorithm
        if(HamiltonianComplex) {
            auto const& op = *std::static_pointer_cast<LOperatorTypeRC<true>>(HOp);
            Partition.reset(new libcommute::space_partition(op, FullHilbertSpace));
        } else {
            auto const& op = *std::static_pointer_cast<LOperatorTypeRC<false>>(HOp);
            Partition.reset(new libcommute::space_partition(op, FullHilbertSpace));
        }
        // Phase II of auto-partition algorithm
        for(ParticleIndex in = 0; in < IndexInfo.getIndexSize(); ++in) {
            auto const& info = IndexInfo.getInfo(in);
            using Operators::c_dag;
            using Operators::c;
            auto Cd =
                LOperatorType<RealType>(Operators::Detail::apply(c_dag<double, IndexTypes...>, info), FullHilbertSpace);
            auto C =
                LOperatorType<RealType>(Operators::Detail::apply(c<double, IndexTypes...>, info), FullHilbertSpace);
            Partition->merge_subspaces(Cd, C, FullHilbertSpace, false);
        }

        Status = Computed;
    }

    IndexClassification<IndexTypes...> const& getIndexInfo() const { return IndexInfo; }
    FullHilbertSpaceType const& getFullHilbertSpace() const { return FullHilbertSpace; }
    SpacePartitionType const& getSpacePartition() const {
        if(Status < Computed)
            throw std::runtime_error("Hilbert space partition has not been computed");
        return *Partition;
    }
};

template <typename ScalarType, typename... IndexTypes>
HilbertSpace<IndexTypes...> MakeHilbertSpace(IndexClassification<IndexTypes...> const& IndexInfo,
                                             Operators::expression<ScalarType, IndexTypes...> const& Hamiltonian,
                                             unsigned int bits_per_boson = 1) {
    return HilbertSpace<IndexTypes...>(IndexInfo, Hamiltonian, bits_per_boson);
}

template <typename ScalarType, typename... IndexTypes>
HilbertSpace<IndexTypes...>
MakeHilbertSpace(IndexClassification<IndexTypes...> const& IndexInfo,
                 Operators::expression<ScalarType, IndexTypes...> const& Hamiltonian,
                 std::map<std::tuple<IndexTypes...>, unsigned int> const& bits_per_boson_map) {
    return HilbertSpace<IndexTypes...>(IndexInfo, Hamiltonian, bits_per_boson_map);
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_HILBERTSPACE_HPP
