//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/GFContainer.h
** \brief Storage of GF for multiple indices (obsolete, remove)
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_GFCONTAINER_HPP
#define POMEROL_INCLUDE_POMEROL_GFCONTAINER_HPP

#include "DensityMatrix.hpp"
#include "FieldOperatorContainer.hpp"
#include "GreensFunction.hpp"
#include "Hamiltonian.hpp"
#include "Index.hpp"
#include "IndexClassification.hpp"
#include "IndexContainer2.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

#include <memory>
#include <set>

namespace Pomerol {

class GFContainer : public IndexContainer2<GreensFunction, GFContainer>, public Thermal {

public:
    template <typename... IndexTypes>
    GFContainer(IndexClassification<IndexTypes...> const& IndexInfo,
                StatesClassification const& S,
                Hamiltonian const& H,
                DensityMatrix const& DM,
                FieldOperatorContainer const& Operators)
        : IndexContainer2<GreensFunction, GFContainer>(*this, IndexInfo),
          Thermal(DM),
          S(S),
          H(H),
          DM(DM),
          Operators(Operators) {}

    void prepareAll(std::set<IndexCombination2> const& InitialIndices = {});
    void computeAll();

protected:
    friend class IndexContainer2<GreensFunction, GFContainer>;
    std::shared_ptr<GreensFunction> createElement(IndexCombination2 const& Indices) const;

    StatesClassification const& S;

    Hamiltonian const& H;
    DensityMatrix const& DM;
    FieldOperatorContainer const& Operators;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_GFCONTAINER_HPP
