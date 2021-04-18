/** \file include/pomerol/MonomialOperator.h
** \brief Declaration of field operators : creation and annihilation operators.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_FIELDOPERATOR_H
#define __INCLUDE_FIELDOPERATOR_H

#include"Misc.h"
#include"StatesClassification.h"
#include"HilbertSpace.h"
#include"Hamiltonian.h"
#include"MonomialOperatorPart.h"
#include"Operators.h"

#include "mpi_dispatcher/misc.hpp"
#include "mpi_dispatcher/mpi_skel.hpp"

#include<boost/bimap.hpp>

#include <memory>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Pomerol {

/** \typedef
 * A pair of left and right indices of a part in a Field Operator. Each part is a non-vanishing worldline in an operator
 */
typedef std::pair<BlockNumber, BlockNumber> BlockMapping;

/** This class is a parent class for creation/annihilation operators which act
 * on all blocks of quantum states */
class MonomialOperator : public ComputableObject
{
public:

    typedef boost::bimaps::bimap<
        boost::bimaps::set_of<BlockNumber>,
        boost::bimaps::set_of<BlockNumber>
    > BlocksBimap;
    typedef BlocksBimap::value_type BlockMapping;

protected:

    friend class FieldOperatorContainer;

    bool Complex;

    /** A reference an Operator object (OperatorPresets::C or Cdag). */
    template<bool Complex>
    using LOperatorType = libcommute::loperator<MelemType<Complex>, libcommute::fermion>;
    std::shared_ptr<void> MOp;

    template<bool Complex>
    const LOperatorType<Complex>& getMOp() const {
        return *std::static_pointer_cast<LOperatorType<Complex>>(MOp);
    }

    /** A reference to a StatesClassification object */
    const StatesClassification &S;
    /** A reference to a Hamiltonian object */
    const Hamiltonian &H;

    const libcommute::space_partition &Partition;
    /** A map between non-vanishing parts (internal numbering) and their R.H.S. BlockNumbers  */
    std::unordered_map<size_t, BlockNumber> mapPartsFromRight;
    /** A map between non-vanishing parts (internal numbering) and their L.H.S. BlockNumbers  */
    std::unordered_map<size_t, BlockNumber> mapPartsFromLeft;

    BlocksBimap LeftRightBlocks;

    /** A vector of parts */
    std::vector<MonomialOperatorPart> parts;

public:

    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    template<typename ScalarType, typename... IndexTypes>
    MonomialOperator(libcommute::expression<ScalarType, IndexTypes...> const& MO,
                     HilbertSpace<IndexTypes...> const& HS,
                     const StatesClassification &S,
                     const Hamiltonian &H) :
        Complex(std::is_same<ScalarType, ComplexType>::value || H.isComplex()),
        MOp(std::make_shared<libcommute::loperator<ScalarType, libcommute::fermion>>(MO, HS.getFullHilbertSpace())),
        S(S), H(H), Partition(HS.getSpacePartition()) {
            if(MO.size() > 1)
                throw std::runtime_error("Only monomial expressions are supported");
        }

    bool isComplex() const { return Complex; }

    /** Returns a FieldOperatorPart based on its left BlockNumber */
    MonomialOperatorPart& getPartFromLeftIndex(BlockNumber in);
    const MonomialOperatorPart& getPartFromLeftIndex(BlockNumber in) const;
    /** Returns a FieldOperatorPart based on its right BlockNumber */
    MonomialOperatorPart& getPartFromRightIndex(BlockNumber out);
    const MonomialOperatorPart& getPartFromRightIndex(BlockNumber out) const;

    /** Returns a left BlockNumber for a given right BlockNumber */
    BlockNumber getLeftIndex(BlockNumber RightIndex) const;
    /** Returns a right BlockNumber for a given left BlockNumber */
    BlockNumber getRightIndex(BlockNumber LeftIndex) const;
    /** Returns a reference to BlockMapping */
    BlocksBimap const& getBlockMapping() const;

    /** Virtual method for assigning world-lines */
    template<typename... IndexTypes>
    void prepare(HilbertSpace<IndexTypes...> const& HS) {
        if (Status >= Prepared) return;

        auto const& FullHS = HS.getFullHilbertSpace();
        auto Connections = Complex ?
            Partition.find_connections(getMOp<true>(), FullHS) :
            Partition.find_connections(getMOp<false>(), FullHS);

        parts.reserve(Connections.size());
        for(auto const& Conn : Connections) {
            mapPartsFromRight.emplace(Conn.first, parts.size());
            mapPartsFromLeft.emplace(Conn.second, parts.size());
            LeftRightBlocks.insert(BlockMapping(Conn.second, Conn.first));

            if(Complex)
                parts.emplace_back(getMOp<true>(), S, H.getPart(Conn.first), H.getPart(Conn.second));
            else
                parts.emplace_back(getMOp<false>(), S, H.getPart(Conn.first), H.getPart(Conn.second));
        }

        Status = Prepared;
    }
    /** Computes all world-lines */
    void compute(const MPI_Comm& comm = MPI_COMM_WORLD);
};

class CreationOperator : public MonomialOperator
{
    ParticleIndex Index;

public:
    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    template<typename... IndexTypes>
    CreationOperator(const IndexClassification<IndexTypes...> &IndexInfo,
                     const HilbertSpace<IndexTypes...> &HS,
                     const StatesClassification &S,
                     const Hamiltonian &H,
                     ParticleIndex Index) :
    MonomialOperator(
        Operators::detail::apply(Operators::c_dag<double, IndexTypes...>,
                                 HS.getIndexInfo().getInfo(Index)
                                 ),
        HS, S, H
    ), Index(Index)
    {}

    ParticleIndex getIndex() const { return Index; }
};

class AnnihilationOperator : public MonomialOperator
{
    ParticleIndex Index;

public:

    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    template<typename... IndexTypes>
    AnnihilationOperator(const IndexClassification<IndexTypes...> &IndexInfo,
                         const HilbertSpace<IndexTypes...> &HS,
                         const StatesClassification &S,
                         const Hamiltonian &H,
                         ParticleIndex Index) :
    MonomialOperator(
        Operators::detail::apply(Operators::c<double, IndexTypes...>,
                                 HS.getIndexInfo().getInfo(Index)
                                 ),
        HS, S, H
    ), Index(Index)
    {}

    ParticleIndex getIndex() const { return Index; }
};

/** A quadratic operator, c_1^+ c_2, in the eigenbasis of a Hamiltonian */
class QuadraticOperator : public MonomialOperator
{
    ParticleIndex Index1, Index2;

public:
    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index1 An index of a creation operator
     * \param[in] Index2 An index of an annihilation operator
     */
    template<typename... IndexTypes>
    QuadraticOperator(const IndexClassification<IndexTypes...> &IndexInfo,
                      const HilbertSpace<IndexTypes...> &HS,
                      const StatesClassification &S,
                      const Hamiltonian &H,
                      ParticleIndex Index1, ParticleIndex Index2) :
    MonomialOperator(
        Operators::detail::apply(Operators::c_dag<double, IndexTypes...>,
                                 HS.getIndexInfo().getInfo(Index1)
                                 ) *
        Operators::detail::apply(Operators::c<double, IndexTypes...>,
                                 HS.getIndexInfo().getInfo(Index2)
                                 ),
        HS, S, H
    ), Index1(Index1), Index2(Index2)
    {}
};

} // end of namespace Pomerol
#endif // endif :: #ifdef __INCLUDE_FIELDOPERATOR_H