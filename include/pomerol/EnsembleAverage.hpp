//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/EnsembleAverage.hpp
/// \brief Ensemble average of a monomial operator representing a physical observable.
/// \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

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

/// \addtogroup Susc
///@{

/// \brief Canonical ensemble average of a monomial operator.
///
/// This class represents the ensemble average of a monomial operator \f$\hat A\f$,
/// \f[
///   \langle A \rangle = Tr[\hat\rho \hat A].
/// \f]
///
/// Usage example:
/// \code{.cpp}
///   EnsembleAverage EA(A /* Monomial operator */, DM /* Density matrix */);
///
///   EA.prepare();
///   auto average = EA();
/// \endcode
class EnsembleAverage : public Thermal, public ComputableObject {

    /// The monomial operator \f$\hat A\f$.
    MonomialOperator const& A;
    /// Many-body density matrix \f$\hat\rho\f$.
    DensityMatrix const& DM;

    /// Computed result
    ComplexType Result = 0;

    /// Implementation detail of prepare().
    template <bool Complex> ComplexType computeImpl(MonomialOperatorPart const& Apart, DensityMatrixPart const& DMpart);

public:
    /// Constructor.
    /// \param[in] A Monomial operator \f$\hat A\f$.
    /// \param[in] DM Many-body density matrix \f$\hat\rho\f$.
    EnsembleAverage(MonomialOperator const& A, DensityMatrix const& DM);
    /// Copy-constructor.
    /// \param[in] EA EnsembleAverage object to be copied.
    EnsembleAverage(EnsembleAverage const& EA);

    /// Compute the ensemble average of \f$\hat A\f$.
    void compute();

    /// Return the ensemble average.
    ComplexType operator()() const { return Result; };
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_ENSEMBLEAVERAGE_HPP
