//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/MonomialOperator.hpp
/// \brief Storage for an operator that is a product of creation/annihilation operators.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_MONOMIALOPERATOR_HPP
#define POMEROL_INCLUDE_MONOMIALOPERATOR_HPP

#include "ComputableObject.hpp"
#include "Hamiltonian.hpp"
#include "HilbertSpace.hpp"
#include "Misc.hpp"
#include "MonomialOperatorPart.hpp"
#include "Operators.hpp"
#include "StatesClassification.hpp"

#include "mpi_dispatcher/misc.hpp"
#include "mpi_dispatcher/mpi_skel.hpp"

#include <boost/bimap.hpp>

#include <cassert>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Pomerol {

/// \addtogroup ED
///@{

/// A pair of invariant subspace indices.
using BlockMapping = std::pair<BlockNumber, BlockNumber>;

/// \brief Monomial quantum operator.
///
/// This class stores an operator \f$\hat M\f$, which is a monomial, i.e. a product
/// of fermionic and/or bosonic creation/annihilation operators. The operator is stored as
/// a list of matrix blocks (\ref MonomialOperatorPart), each connecting a pair of
/// invariant subspaces of the Hamiltonian. For a given right invariant subspace,
/// there exists at most one part connecting it to a left subspace (and the other way around).
class MonomialOperator : public ComputableObject {
public:
    /// A bi-map container for connections between invariant subspaces established by a monomial operator.
    using BlocksBimap = boost::bimaps::bimap<boost::bimaps::set_of<BlockNumber>, boost::bimaps::set_of<BlockNumber>>;
    /// A single subspace-to-subspace connection established by a monomial operator.
    using BlockMapping = BlocksBimap::value_type;

protected:
    friend class FieldOperatorContainer;

    /// Whether the \ref MOp object is complex-valued.
    bool MOpComplex;
    /// A type-erased real/complex-valued \p libcommute::loperator object.
    std::shared_ptr<void> MOp;

    /// Return a constant reference to the stored \p libcommute::loperator object.
    /// \tparam Complex Request a reference to the complex-valued linear operator.
    /// \pre The compile-time value of \ref Complex must agree with the result of \ref isComplex().
    template <bool Complex> LOperatorTypeRC<Complex> const& getMOp() const {
        assert(MOpComplex == Complex);
        return *std::static_pointer_cast<LOperatorTypeRC<Complex>>(MOp);
    }

    /// Whether the stored parts are complex-valued.
    bool Complex;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;
    /// The Hamiltonian.
    Hamiltonian const& H;

    /// A map between positions of parts in the \ref parts list and the respective right subspace indices.
    std::unordered_map<std::size_t, BlockNumber> mapPartsFromRight;
    /// A map between positions of parts in the \ref parts list and the respective left subspace indices.
    std::unordered_map<std::size_t, BlockNumber> mapPartsFromLeft;

    /// Left-to-right connections between invariant subspaces established by this monomial operator.
    BlocksBimap LeftRightBlocks;

    /// List of parts (matrix blocks).
    std::vector<MonomialOperatorPart> parts;

public:
    /// Constructor.
    /// \tparam ScalarType Scalar type (either double or std::complex<double>) of the expression \p MO.
    /// \tparam IndexTypes Types of indices carried by operators in the expression \p MO.
    /// \param[in] MO Expression of the monomial operator \f$\hat M\f$.
    /// \param[in] HS Hilbert space.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \pre \p MO is a monomial operator.
    template <typename ScalarType, typename... IndexTypes>
    MonomialOperator(libcommute::expression<ScalarType, IndexTypes...> const& MO,
                     HilbertSpace<IndexTypes...> const& HS,
                     StatesClassification const& S,
                     Hamiltonian const& H)
        : MOpComplex(std::is_same<ScalarType, ComplexType>::value),
          MOp(std::make_shared<LOperatorType<ScalarType>>(MO, HS.getFullHilbertSpace())),
          Complex(MOpComplex || H.isComplex()),
          S(S),
          H(H) {
        if(MO.size() > 1)
            throw std::runtime_error("Only monomial expressions are supported");
    }

    /// Is the monomial operator a complex-valued matrix?
    bool isComplex() const { return Complex; }

    /// Return a reference to the part by a given left invariant subspace.
    /// \param[in] LeftIndex Index of the left invariant subspace
    /// \pre \ref prepare() has been called.
    MonomialOperatorPart& getPartFromLeftIndex(BlockNumber LeftIndex);
    /// Return a constant reference to the part by a given left invariant subspace.
    /// \param[in] LeftIndex Index of the left invariant subspace
    /// \pre \ref prepare() has been called.
    MonomialOperatorPart const& getPartFromLeftIndex(BlockNumber LeftIndex) const;
    /// Return a reference to the part by a given right invariant subspace.
    /// \param[in] RightIndex Index of the right invariant subspace
    /// \pre \ref prepare() has been called.
    MonomialOperatorPart& getPartFromRightIndex(BlockNumber RightIndex);
    /// Return a constant reference to the part by a given right invariant subspace.
    /// \param[in] RightIndex Index of the right invariant subspace
    /// \pre \ref prepare() has been called.
    MonomialOperatorPart const& getPartFromRightIndex(BlockNumber RightIndex) const;

    /// For a given right invariant subspace, return the corresponding left invariant subspace.
    /// \param[in] RightIndex Index of the right subspace.
    /// \return Index of the left subspace.
    /// \pre \ref prepare() has been called.
    BlockNumber getLeftIndex(BlockNumber RightIndex) const;
    /// For a given left invariant subspace, return the corresponding right invariant subspace.
    /// \param[in] LeftIndex Index of the left subspace.
    /// \return Index of the right subspace.
    /// \pre \ref prepare() has been called.
    BlockNumber getRightIndex(BlockNumber LeftIndex) const;

    /// Return a constant reference to the left-to-right connection map.
    /// \pre \ref prepare() has been called.
    BlocksBimap const& getBlockMapping() const;

    /// Allocate memory for all parts.
    /// \tparam IndexTypes Types of indices carried by operators acting in the Hilbert space \p HS.
    /// \param[in] HS The Hilbert space.
    template <typename... IndexTypes> void prepare(HilbertSpace<IndexTypes...> const& HS) {
        if(getStatus() >= Prepared)
            return;

        if(HS.getStatus() != Computed) { // Hilbert space has not been partitioned
            mapPartsFromRight.emplace(0, 0);
            mapPartsFromLeft.emplace(0, 0);
            LeftRightBlocks.insert(BlockMapping(0, 0));

            if(MOpComplex)
                parts.emplace_back(getMOp<true>(), S, H.getPart(0), H.getPart(0));
            else
                parts.emplace_back(getMOp<false>(), S, H.getPart(0), H.getPart(0));

            Status = Prepared;
            return;
        }

        auto const& Partition = HS.getSpacePartition();
        auto Connections =
            MOpComplex ? Partition.find_connections(getMOp<true>()) : Partition.find_connections(getMOp<false>());

        parts.reserve(Connections.size());
        for(auto const& Conn : Connections) {
            mapPartsFromRight.emplace(Conn.first, parts.size());
            mapPartsFromLeft.emplace(Conn.second, parts.size());
            LeftRightBlocks.insert(BlockMapping(Conn.second, Conn.first));

            if(MOpComplex)
                parts.emplace_back(getMOp<true>(), S, H.getPart(Conn.first), H.getPart(Conn.second));
            else
                parts.emplace_back(getMOp<false>(), S, H.getPart(Conn.first), H.getPart(Conn.second));
        }

        Status = Prepared;
    }

    /// Compute matrix elements of all parts in parallel.
    /// \param[in] Tolerance Matrix elements with the absolute value equal or below this threshold
    ///                      are considered negligible.
    /// \param[in] comm MPI communicator used to parallelize the computation.
    /// \pre \ref prepare() has been called.
    void compute(RealType Tolerance = 1e-8, MPI_Comm const& comm = MPI_COMM_WORLD);

private:
    // Implementation details
    void checkPrepared() const;
};

/// A special case of a monomial operator: A single fermion creation or annihilation operator \f$\hat F_i\f$.
class FieldOperator : public MonomialOperator {
    /// The single-particle index corresponding to the field operator.
    ParticleIndex Index;

public:
    /// Constructor.
    /// \tparam ScalarType Scalar type (either double or std::complex<double>) of the expression \p FO.
    /// \tparam IndexTypes Types of indices carried by operators in the expression \p FO.
    /// \param[in] FO Expression of the field operator \f$\hat F_i\f$.
    /// \param[in] HS Hilbert space.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] Index The single-particle index \f$i\f$.
    /// \pre \p FO is a field operator.
    template <typename ScalarType, typename... IndexTypes>
    FieldOperator(libcommute::expression<ScalarType, IndexTypes...> const& FO,
                  HilbertSpace<IndexTypes...> const& HS,
                  StatesClassification const& S,
                  Hamiltonian const& H,
                  ParticleIndex Index)
        : MonomialOperator(FO, HS, S, H), Index(Index) {
        auto const& mon = FO.get_monomials().cbegin()->first;
        if(mon.size() != 1 || !is_fermion(mon[0]))
            throw std::runtime_error("Expected a single-fermion monomial expression");
    }

    /// Return the single-particle index \f$i\f$.
    ParticleIndex getIndex() const { return Index; }
};

/// A special case of a monomial operator: A single fermion creation operator \f$c^\dagger_i\f$.
class CreationOperator : public FieldOperator {
public:
    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by operators acting in the Hilbert space \p HS.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] HS Hilbert space.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] Index The single-particle index \f$i\f$.
    template <typename... IndexTypes>
    CreationOperator(IndexClassification<IndexTypes...> const& IndexInfo,
                     HilbertSpace<IndexTypes...> const& HS,
                     StatesClassification const& S,
                     Hamiltonian const& H,
                     ParticleIndex Index)
        : FieldOperator(Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index)),
                        HS,
                        S,
                        H,
                        Index) {}
};

/// A special case of a monomial operator: A single fermion annihilation operator \f$c_i\f$.
class AnnihilationOperator : public FieldOperator {
public:
    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by operators acting in the Hilbert space \p HS.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] HS Hilbert space.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] Index The single-particle index \f$i\f$.
    template <typename... IndexTypes>
    AnnihilationOperator(IndexClassification<IndexTypes...> const& IndexInfo,
                         HilbertSpace<IndexTypes...> const& HS,
                         StatesClassification const& S,
                         Hamiltonian const& H,
                         ParticleIndex Index)
        : FieldOperator(Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index)),
                        HS,
                        S,
                        H,
                        Index) {}
};

/// A special case of a monomial operator: A product of two fermionic operators \f$O_i O_j\f$.
/// Each of \f$O_i\f$ and \f$O_j\f$ can be either a creation or annihilation operator.
class QuadraticOperator : public MonomialOperator {
    /// The single-particle index corresponding to the first operator.
    ParticleIndex Index1;
    /// The single-particle index corresponding to the second operator.
    ParticleIndex Index2;
    /// Indicates whether each of the two operators is a creator.
    std::tuple<bool, bool> Dagger;

public:
    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by operators acting in the Hilbert space \p HS.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] HS Hilbert space.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] Index1 The single-particle index \f$i\f$ of the first operator.
    /// \param[in] Index2 The single-particle index \f$j\f$ of the second operator.
    /// \param[in] Dagger Indicates whether each of the two operators is a creator.
    template <typename... IndexTypes>
    QuadraticOperator(IndexClassification<IndexTypes...> const& IndexInfo,
                      HilbertSpace<IndexTypes...> const& HS,
                      StatesClassification const& S,
                      Hamiltonian const& H,
                      ParticleIndex Index1,
                      ParticleIndex Index2,
                      std::tuple<bool, bool> const& Dagger = std::tuple<bool, bool>(true, false))
        : MonomialOperator(
              (std::get<0>(Dagger) ?
                   Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index1)) :
                   Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index1))) *
                  (std::get<1>(Dagger) ?
                       Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index2)) :
                       Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index2))),
              HS,
              S,
              H),
          Index1(Index1),
          Index2(Index2),
          Dagger(Dagger) {}

    /// Return the single-particle index \f$i\f$
    ParticleIndex getIndex1() const { return Index1; }
    /// Return the single-particle index \f$j\f$
    ParticleIndex getIndex2() const { return Index2; }

    /// Return the single-particle index \f$i\f$ under the assumption that O_i is a creation operator.
    ParticleIndex getCXIndex() const {
        assert(std::get<0>(Dagger));
        return Index1;
    }
    /// Return the single-particle index \f$j\f$ under the assumption that O_j is an annihilation operator.
    ParticleIndex getCIndex() const {
        assert(!std::get<1>(Dagger));
        return Index2;
    }

    /// Return the creation/annihilation type of each of the two operators.
    std::tuple<bool, bool> const& getDagger() const { return Dagger; }
};

/// A special case of a monomial operator: A product of four fermionic operators \f$O_i O_j O_k O_l\f$.
/// Each of \f$O_i, O_j, O_k\f$ and \f$O_l\f$ can be either a creation or annihilation operator.
class QuarticOperator : public MonomialOperator {
    /// The single-particle index corresponding to the first operator.
    ParticleIndex Index1;
    /// The single-particle index corresponding to the second operator.
    ParticleIndex Index2;
    /// The single-particle index corresponding to the third operator.
    ParticleIndex Index3;
    /// The single-particle index corresponding to the fourth operator.
    ParticleIndex Index4;
    /// Indicates whether each of the four operators is a creator.
    std::tuple<bool, bool, bool, bool> Dagger;

public:
    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by operators acting in the Hilbert space \p HS.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] HS Hilbert space.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] Index1 The single-particle index \f$i\f$ of the first operator.
    /// \param[in] Index2 The single-particle index \f$j\f$ of the second operator.
    /// \param[in] Index3 The single-particle index \f$k\f$ of the third operator.
    /// \param[in] Index4 The single-particle index \f$l\f$ of the fourth operator.
    /// \param[in] Dagger Indicates whether each of the four operators is a creator.
    template <typename... IndexTypes>
    QuarticOperator(
        IndexClassification<IndexTypes...> const& IndexInfo,
        HilbertSpace<IndexTypes...> const& HS,
        StatesClassification const& S,
        Hamiltonian const& H,
        ParticleIndex Index1,
        ParticleIndex Index2,
        ParticleIndex Index3,
        ParticleIndex Index4,
        std::tuple<bool, bool, bool, bool> const& Dagger = std::tuple<bool, bool, bool, bool>(true, true, false, false))
        : MonomialOperator(
              (std::get<0>(Dagger) ?
                   Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index1)) :
                   Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index1))) *
                  (std::get<1>(Dagger) ?
                       Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index2)) :
                       Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index2))) *
                  (std::get<2>(Dagger) ?
                       Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index3)) :
                       Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index3))) *
                  (std::get<3>(Dagger) ?
                       Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index4)) :
                       Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index4))),
              HS,
              S,
              H),
          Index1(Index1),
          Index2(Index2),
          Index3(Index3),
          Index4(Index4),
          Dagger(Dagger) {}

    /// Return the single-particle index \f$i\f$
    ParticleIndex getIndex1() const { return Index1; }
    /// Return the single-particle index \f$j\f$
    ParticleIndex getIndex2() const { return Index2; }
    /// Return the single-particle index \f$k\f$
    ParticleIndex getIndex3() const { return Index3; }
    /// Return the single-particle index \f$l\f$
    ParticleIndex getIndex4() const { return Index4; }

    /// Return the single-particle index \f$i\f$ under the assumption that O_i is a creation operator.
    ParticleIndex getCX1Index() const {
        assert(std::get<0>(Dagger));
        return Index1;
    }
    /// Return the single-particle index \f$j\f$ under the assumption that O_j is a creation operator.
    ParticleIndex getCX2Index() const {
        assert(std::get<1>(Dagger));
        return Index2;
    }
    /// Return the single-particle index \f$k\f$ under the assumption that O_k is an annihilation operator.
    ParticleIndex getC1Index() const {
        assert(!std::get<2>(Dagger));
        return Index3;
    }
    /// Return the single-particle index \f$l\f$ under the assumption that O_l is an annihilation operator.
    ParticleIndex getC2ndex() const {
        assert(!std::get<3>(Dagger));
        return Index4;
    }

    /// Return the creation/annihilation type of each of the four operators.
    std::tuple<bool, bool, bool, bool> const& getDagger() const { return Dagger; }
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_MONOMIALOPERATOR_HPP
