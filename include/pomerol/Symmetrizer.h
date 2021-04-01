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
class Symmetrizer : public ComputableObject {

public:

    using OperatorType = Operators::expression<ScalarType, IndexTypes...>;
    using HilbertSpaceType = libcommute::hilbert_space<IndexTypes...>;
    using SpacePartitionType = libcommute::space_partition;

private:

    IndexClassification<IndexTypes...> const& IndexInfo;

    HilbertSpaceType HilbertSpace;

    libcommute::loperator<ScalarType, libcommute::fermion> HamiltonianLOp;

    std::unique_ptr<SpacePartitionType> partition = nullptr;

public:

    Symmetrizer(const IndexClassification<IndexTypes...> &IndexInfo,
                const OperatorType& Hamiltonian) :
      IndexInfo(IndexInfo),
      HilbertSpace(InitHilbertSpace(IndexInfo)),
      HamiltonianLOp(Hamiltonian, HilbertSpace)
    {}

    void compute() {
        if(Status >= Computed) return;
        partition.reset(new libcommute::space_partition(HamiltonianLOp, HilbertSpace));
        Status = Computed;
    }

    HilbertSpaceType const& getHilbertSpace() const { return HilbertSpace; }
    SpacePartitionType const& getSpacePartition() const {
        if(Status < Computed)
            throw std::runtime_error("Hilbert space partition has not been computed");
        return *partition;
    }

private:

  HilbertSpaceType InitHilbertSpace(const IndexClassification<IndexTypes...> &IndexInfo) {
      HilbertSpaceType hs;
      for(ParticleIndex p = 0; p < IndexInfo.getIndexSize(); ++p) {
          hs.add(libcommute::elementary_space_fermion<IndexTypes...>(IndexInfo.getInfo(p)));
      }
      return hs;
  }
};

template<typename ScalarType, typename... IndexTypes>
Symmetrizer<ScalarType, IndexTypes...>
MakeSymmetrizer(const IndexClassification<IndexTypes...> &IndexInfo,
                const Operators::expression<ScalarType, IndexTypes...>& Hamiltonian) {
  return Symmetrizer<ScalarType, IndexTypes...>(IndexInfo, Hamiltonian);
}

}; // end of namespace Pomerol

#endif //  endif :: #ifndef __INCLUDE_SYMMETRIZER_H
