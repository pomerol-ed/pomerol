//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/HilbertSpace.hpp
/// \brief Hilbert space of a system and invariant subspaces of its Hamiltonian.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

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

/// \addtogroup ED
///@{

/// \brief Hilbert space of a quantum system.
///
/// A thin wrapper around libcommute's
/// <a href="https://krivenko.github.io/libcommute/loperator/hilbert_space.html">hilbert_space</a>
/// (information about a finite-dimensional state space)
/// and <a href="https://krivenko.github.io/libcommute/loperator/space_partition.html">space_partition</a>
/// (partition of the full state space into invariant subspaces of a Hamiltonian).
/// \tparam IndexTypes Types of indices carried by operators acting in this Hilbert space.
template <typename... IndexTypes> class HilbertSpace : public ComputableObject {

public:
    /// Type of the full Hilbert space.
    using FullHilbertSpaceType = libcommute::hilbert_space<IndexTypes...>;
    /// Type of the partition into invariant subspaces.
    using SpacePartitionType = libcommute::space_partition<FullHilbertSpaceType>;

private:
    /// Parent operator index tuple map.
    IndexClassification<IndexTypes...> const& IndexInfo;

    /// Full Hilbert space of the problem.
    FullHilbertSpaceType FullHilbertSpace;

    /// Has a complex-valued Hamiltonian been used to construct this object?
    bool HamiltonianComplex;

    /// A type-erased real- or complex-valued linear operator corresponding to system's Hamiltonian.
    std::shared_ptr<void> HOp;

    /// A Hilbert space partition object.
    std::unique_ptr<SpacePartitionType> Partition = nullptr;

    /// A libcommute-compatible bosonic elementary space constructor that allows automatic
    /// creation of spaces with different sizes for different indices of respective
    /// operators \f$a^+\f$/\f$a\f$.
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
                return make_unique<elementary_space_boson<IndexTypes...>>(1 << bits, g.indices());
            } else {
                std::stringstream ss;
                ss << "Unexpected algebra generator " << g;
                throw std::runtime_error(ss.str());
            }
        }
    };

public:
    /// Construct a full Hilbert space from an \ref IndexClassification object and
    /// a \ref Operators "polynomial expression" of system's Hamiltonian. The Hilbert space is constructed
    /// as a direct product of elementary spaces, each associated with a single fermionic or bosonic
    /// degree of freedom (an index tuple carried by a boson creation/annihilation operator).
    /// \tparam ScalarType Coefficient type of expression H.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] H Hamiltonian of the system.
    /// \param[in] bits_per_boson Each bosonic degree of freedom  will result in a truncated elementary bosonic space
    ///            of dimension \f$2^\verb|bits_per_boson|\f$.
    template <typename ScalarType>
    HilbertSpace(IndexClassification<IndexTypes...> const& IndexInfo,
                 Operators::expression<ScalarType, IndexTypes...> const& H,
                 unsigned int bits_per_boson = 1)
        : IndexInfo(IndexInfo),
          FullHilbertSpace(H, libcommute::boson_es_constructor(1 << bits_per_boson)),
          HamiltonianComplex(std::is_same<ScalarType, ComplexType>::value),
          HOp(std::make_shared<LOperatorType<ScalarType>>(H, FullHilbertSpace)) {
        if(FullHilbertSpace.is_sparse())
            throw std::runtime_error("libcommute's sparse Hilbert spaces are not supported by Pomerol");
    }

    /// Construct a full Hilbert space from an \ref IndexClassification object and
    /// a \ref Operators "polynomial expression" of system's Hamiltonian. The Hilbert space is constructed
    /// as a direct product of elementary spaces, each associated with a single fermionic or bosonic
    /// degree of freedom (an index tuple carried by a boson creation/annihilation operator).
    /// \tparam ScalarType Coefficient type of expression H.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] H Hamiltonian of the system.
    /// \param[in] bits_per_boson_map A bosonic degree of freedom with a certain operator index tuple
    ///            will result in a truncated elementary bosonic space
    ///            of dimension \f$2^b\f$, where \f$b\f$ is the value in this map corresponding to
    ///            the index tuple. If the tuple is missing from the map, \f$b\f$ is taken to be 1.
    template <typename ScalarType>
    HilbertSpace(IndexClassification<IndexTypes...> const& IndexInfo,
                 Operators::expression<ScalarType, IndexTypes...> const& H,
                 std::map<std::tuple<IndexTypes...>, unsigned int> const& bits_per_boson_map)
        : IndexInfo(IndexInfo),
          FullHilbertSpace(H, boson_es_constructor_from_map(bits_per_boson_map)),
          HamiltonianComplex(std::is_same<ScalarType, ComplexType>::value),
          HOp(std::make_shared<LOperatorType<ScalarType>>(H, FullHilbertSpace)) {}

    /// Find a partition of the full Hilbert space into invariant subspaces of the Hamiltonian.
    /// The partition fulfills an additional requirement that all fermionic creation/annihilation
    /// operators connect one invariant subspace to at most one subspace.
    void compute() {
        if(getStatus() >= Computed)
            return;

        // Phase I of auto-partition algorithm
        if(HamiltonianComplex) {
            auto const& op = *std::static_pointer_cast<LOperatorTypeRC<true>>(HOp);
            Partition.reset(new SpacePartitionType(op, FullHilbertSpace));
        } else {
            auto const& op = *std::static_pointer_cast<LOperatorTypeRC<false>>(HOp);
            Partition.reset(new SpacePartitionType(op, FullHilbertSpace));
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
            Partition->merge_subspaces(Cd, C, false);
        }

        setStatus(Computed);
    }

    /// Access the \ref IndexClassification object used to construct this Hilbert space.
    IndexClassification<IndexTypes...> const& getIndexInfo() const { return IndexInfo; }

    /// Access the full Hilbert space object.
    FullHilbertSpaceType const& getFullHilbertSpace() const { return FullHilbertSpace; }

    /// Access the space partition object.
    /// \pre \ref compute() has been called.
    SpacePartitionType const& getSpacePartition() const {
        if(getStatus() < Computed)
            throw std::runtime_error("Hilbert space partition has not been computed");
        return *Partition;
    }
};

/// A factory function for \ref HilbertSpace that constructs it from an \ref IndexClassification object and
/// a \ref Operators "polynomial expression" of system's Hamiltonian. The Hilbert space is constructed
/// as a direct product of elementary spaces, each associated with a single fermionic or bosonic
/// degree of freedom (an index tuple carried by a boson creation/annihilation operator).
/// \tparam ScalarType Coefficient type of expression H.
/// \param[in] IndexInfo Map for fermionic operator index tuples.
/// \param[in] H Hamiltonian of the system.
/// \param[in] bits_per_boson Each bosonic degree of freedom  will result in a truncated elementary bosonic space
///            of dimension \f$2^\verb|bits_per_boson|\f$.
template <typename ScalarType, typename... IndexTypes>
HilbertSpace<IndexTypes...> MakeHilbertSpace(IndexClassification<IndexTypes...> const& IndexInfo,
                                             Operators::expression<ScalarType, IndexTypes...> const& H,
                                             unsigned int bits_per_boson = 1) {
    return HilbertSpace<IndexTypes...>(IndexInfo, H, bits_per_boson);
}

/// A factory function for \ref HilbertSpace that constructs it from an \ref IndexClassification object and
/// a \ref Operators "polynomial expression" of system's Hamiltonian. The Hilbert space is constructed
/// as a direct product of elementary spaces, each associated with a single fermionic or bosonic
/// degree of freedom (an index tuple carried by a boson creation/annihilation operator).
/// \tparam ScalarType Coefficient type of expression H.
/// \param[in] IndexInfo Map for fermionic operator index tuples.
/// \param[in] H Hamiltonian of the system.
/// \param[in] bits_per_boson_map A bosonic degree of freedom with a certain operator index tuple
///            will result in a truncated elementary bosonic space of dimension \f$2^b\f$,
///            where \f$b\f$ is the value in this map corresponding to
///            the index tuple. If the tuple is missing from the map, \f$b\f$ is taken to be 1.
template <typename ScalarType, typename... IndexTypes>
HilbertSpace<IndexTypes...>
MakeHilbertSpace(IndexClassification<IndexTypes...> const& IndexInfo,
                 Operators::expression<ScalarType, IndexTypes...> const& H,
                 std::map<std::tuple<IndexTypes...>, unsigned int> const& bits_per_boson_map) {
    return HilbertSpace<IndexTypes...>(IndexInfo, H, bits_per_boson_map);
}

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_HILBERTSPACE_HPP
