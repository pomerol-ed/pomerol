//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/MonomialOperator.hpp
/// \brief Storage for an operator that is a product of creation/annihilation operators.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)
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

#include <cstddef>
#include <memory>
#include <stdexcept>
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

        auto const& FullHS = HS.getFullHilbertSpace();
        auto const& Partition = HS.getSpacePartition();
        auto Connections = MOpComplex ? Partition.find_connections(getMOp<true>(), FullHS) :
                                        Partition.find_connections(getMOp<false>(), FullHS);

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
    /// \param[in] comm MPI communicator used to parallelize the computation.
    /// \pre \ref prepare() has been called.
    void compute(MPI_Comm const& comm = MPI_COMM_WORLD);

private:
    // Implementation details
    void checkPrepared() const;
};

/// A special case of a monomial operator: A single fermion creation operator \f$c^\dagger_i\f$.
class CreationOperator : public MonomialOperator {
    /// The single-particle index corresponding to the creation operator.
    ParticleIndex Index;

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
        : MonomialOperator(Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index)),
                           HS,
                           S,
                           H),
          Index(Index) {}

    /// Return the single-particle index \f$i\f$.
    ParticleIndex getIndex() const { return Index; }
};

/// A special case of a monomial operator: A single fermion annihilation operator \f$c_i\f$.
class AnnihilationOperator : public MonomialOperator {
    /// The single-particle index corresponding to the annihilation operator.
    ParticleIndex Index;

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
        : MonomialOperator(Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index)),
                           HS,
                           S,
                           H),
          Index(Index) {}

    /// Return the single-particle index \f$i\f$.
    ParticleIndex getIndex() const { return Index; }
};

/// A special case of a monomial operator: A single quadratic fermionic operator \f$c^\dagger_i c_j\f$.
class QuadraticOperator : public MonomialOperator {
    /// The single-particle index corresponding to the creation operator.
    ParticleIndex Index1;
    /// The single-particle index corresponding to the annihilation operator.
    ParticleIndex Index2;

public:
    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by operators acting in the Hilbert space \p HS.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] HS Hilbert space.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] Index1 The single-particle index \f$i\f$ of the creation operator.
    /// \param[in] Index2 The single-particle index \f$j\f$ of the annihilation operator.
    template <typename... IndexTypes>
    QuadraticOperator(IndexClassification<IndexTypes...> const& IndexInfo,
                      HilbertSpace<IndexTypes...> const& HS,
                      StatesClassification const& S,
                      Hamiltonian const& H,
                      ParticleIndex Index1,
                      ParticleIndex Index2)
        : MonomialOperator(
              Operators::Detail::apply(Operators::c_dag<double, IndexTypes...>, IndexInfo.getInfo(Index1)) *
                  Operators::Detail::apply(Operators::c<double, IndexTypes...>, IndexInfo.getInfo(Index2)),
              HS,
              S,
              H),
          Index1(Index1),
          Index2(Index2) {}

    /// Return the single-particle index \f$i\f$.
    ParticleIndex getCXXIndex() const { return Index1; }
    /// Return the single-particle index \f$j\f$.
    ParticleIndex getCIndex() const { return Index2; }
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_MONOMIALOPERATOR_HPP
