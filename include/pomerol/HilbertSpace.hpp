/** \file Symmetrizer.h
**  \brief Declaration of the Symmetrizer class - a class to store and get the information about the symmetries of the system.
**
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_HILBERTSPACE_H
#define POMEROL_INCLUDE_POMEROL_HILBERTSPACE_H

#include "Misc.hpp"
#include "IndexClassification.hpp"
#include "ComputableObject.hpp"
#include "Operators.hpp"

#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/loperator.hpp>
#include <libcommute/loperator/space_partition.hpp>

#include <memory>
#include <stdexcept>
#include <type_traits>

namespace Pomerol {

template<typename... IndexTypes>
class HilbertSpace : public ComputableObject {

public:

    using FullHilbertSpaceType = libcommute::hilbert_space<IndexTypes...>;
    using SpacePartitionType = libcommute::space_partition;

private:

    IndexClassification<IndexTypes...> const& IndexInfo;

    FullHilbertSpaceType FullHilbertSpace;

    bool HamiltonianComplex;

    std::shared_ptr<void> HOp;

    std::unique_ptr<SpacePartitionType> Partition = nullptr;

public:

    template<typename ScalarType>
    HilbertSpace(const IndexClassification<IndexTypes...> &IndexInfo,
                 const Operators::expression<ScalarType, IndexTypes...>& Hamiltonian) :
      IndexInfo(IndexInfo),
      FullHilbertSpace(InitFullHilbertSpace(IndexInfo)),
      HamiltonianComplex(std::is_same<ScalarType, ComplexType>::value),
      HOp(std::make_shared<LOperatorType<ScalarType>>(Hamiltonian, FullHilbertSpace))
    {}

    void compute() {
        if(Status >= Computed) return;

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
            using Operators::Detail::apply;
            auto Cd = LOperatorType<RealType>(
                apply(c_dag<double, IndexTypes...>, info),
                FullHilbertSpace
            );
            auto C = LOperatorType<RealType>(
                apply(c<double, IndexTypes...>, info),
                FullHilbertSpace
            );
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

private:

  FullHilbertSpaceType InitFullHilbertSpace(const IndexClassification<IndexTypes...> &IndexInfo) {
      FullHilbertSpaceType FullHS;
      for(ParticleIndex p = 0; p < IndexInfo.getIndexSize(); ++p) {
          FullHS.add(libcommute::elementary_space_fermion<IndexTypes...>(IndexInfo.getInfo(p)));
      }
      return FullHS;
  }
};

template<typename ScalarType, typename... IndexTypes>
HilbertSpace<IndexTypes...>
MakeHilbertSpace(const IndexClassification<IndexTypes...> &IndexInfo,
                 const Operators::expression<ScalarType, IndexTypes...>& Hamiltonian) {
  return HilbertSpace<IndexTypes...>(IndexInfo, Hamiltonian);
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_HILBERTSPACE_H
