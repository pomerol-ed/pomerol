//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/EnsembleAverage.h
** \brief Ensemble average.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
*/
#ifndef POMEROL_INCLUDE_POMEROL_ENSEMBLEAVERAGE_HPP
#define POMEROL_INCLUDE_POMEROL_ENSEMBLEAVERAGE_HPP

#include "ComputableObject.hpp"
#include "DensityMatrix.hpp"
#include "Hamiltonian.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

namespace Pomerol {

/** This class represents the ensemble average of a quadratic operator.
 *
 * Exact definition:
 *
 * \f[
 *      \langle A \rangle = \langle c_i^{\dag} c_j \rangle
 * \f]
 *
 * How to use:
 *   EnsembleAverage EA(A, DM);
 *   EA.prepare()
 *   EA.getResult()
 */
class EnsembleAverage : public Thermal, public ComputableObject {

    /** A reference to a bosonic operator. */
    MonomialOperator const& A;
    /** A reference to a density matrix. */
    DensityMatrix const& DM;

    ComplexType Result = 0;

    /** Returns the contribution to the ensemble average from a part. Called in prepare() */
    template <bool Complex> ComplexType computeImpl(MonomialOperatorPart const& Apart, DensityMatrixPart const& DMpart);

public:
    /** Constructor.
     * \param[in] A A reference to a quadratic operator.
     * \param[in] DM A reference to a density matrix.
     */
    EnsembleAverage(MonomialOperator const& A, DensityMatrix const& DM);
    /** Copy-constructor.
     * \param[in] EA EnsembleAverage object to be copied.
     */
    EnsembleAverage(EnsembleAverage const& EA);

    /** Compute the ensemble average of A by choosing relevant parts of A and sum up each contribution. */
    void compute();

    /** Returns the ensemble average */
    ComplexType operator()() const { return Result; };
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_ENSEMBLEAVERAGE_HPP
