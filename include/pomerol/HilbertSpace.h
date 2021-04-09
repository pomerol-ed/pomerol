/** \file Symmetrizer.h
**  \brief Declaration of the Symmetrizer class - a class to store and get the information about the symmetries of the system.
**
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_SYMMETRIZER_H
#define __INCLUDE_SYMMETRIZER_H

#include "Misc.h"
#include "IndexClassification.h"
#include "ComputableObject.h"
#include "Operators.h"

#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/loperator.hpp>
#include <libcommute/loperator/space_partition.hpp>

#include <memory>
#include <stdexcept>

namespace Pomerol {

template<typename ScalarType, typename... IndexTypes>
class HilbertSpace : public ComputableObject {

public:

    using OperatorType = Operators::expression<ScalarType, IndexTypes...>;
    using FullHilbertSpaceType = libcommute::hilbert_space<IndexTypes...>;
    using SpacePartitionType = libcommute::space_partition;

private:

    IndexClassification<IndexTypes...> const& IndexInfo;

    FullHilbertSpaceType FullHilbertSpace;

    libcommute::loperator<ScalarType, libcommute::fermion> HamiltonianLOp;

    std::unique_ptr<SpacePartitionType> partition = nullptr;

public:

    HilbertSpace(const IndexClassification<IndexTypes...> &IndexInfo,
                 const OperatorType& Hamiltonian) :
      IndexInfo(IndexInfo),
      FullHilbertSpace(InitFullHilbertSpace(IndexInfo)),
      HamiltonianLOp(Hamiltonian, FullHilbertSpace)
    {}

    void compute() {
        if(Status >= Computed) return;
        partition.reset(new libcommute::space_partition(HamiltonianLOp, FullHilbertSpace));
        Status = Computed;
    }

    FullHilbertSpaceType const& getFullHilbertSpace() const { return FullHilbertSpace; }
    SpacePartitionType const& getSpacePartition() const {
        if(Status < Computed)
            throw std::runtime_error("Hilbert space partition has not been computed");
        return *partition;
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
HilbertSpace<ScalarType, IndexTypes...>
MakeHilbertSpace(const IndexClassification<IndexTypes...> &IndexInfo,
                 const Operators::expression<ScalarType, IndexTypes...>& Hamiltonian) {
  return HilbertSpace<ScalarType, IndexTypes...>(IndexInfo, Hamiltonian);
}

}; // end of namespace Pomerol

#endif //  endif :: #ifndef __INCLUDE_SYMMETRIZER_H
