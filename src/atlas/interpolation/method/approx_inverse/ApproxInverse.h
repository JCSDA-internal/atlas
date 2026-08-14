/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/functionspace/FunctionSpace.h"

#include "atlas/interpolation/method/Method.h"
#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"
#include "atlas/linalg/sparse/SparseMatrixView.h"

#include "eckit/config/Configuration.h"

namespace atlas {

class Grid;

namespace interpolation {

class Cache;

namespace method {

/// @brief Abstract base class for approximate inverse interpolation methods.
///
/// Given a forward interpolation from a @e source grid to a @e target grid,
/// this class constructs an approximate inverse interpolation operator that
/// maps from the target grid back to the source grid.
///
/// Setup builds an ancillary forward interpolation in the reverse direction
/// (target → source), redistributes its matrix across MPI ranks via a halo
/// exchange, delegates to the subclass-specific @c approx_inverse_transform,
/// and finally normalises each row of the resulting matrix so that its
/// weights sum to one.
///
/// Concrete subclasses must implement @c approx_inverse_transform to provide
/// the specific inversion strategy (e.g. binning or local pseudoinverse).
class ApproxInverse : public Method {
public:
    using Method::do_setup;
    using Value               = double;
    using Index               = int;
    using Triplet             = linalg::Triplet<Value, Index>;
    using SparseMatrixStorage = linalg::SparseMatrixStorage;
    using SparseMatrixView    = linalg::SparseMatrixView<Value, Index>;

    /// @brief Construct from an eckit configuration.
    /// @param config  Must be (or derive from) @c eckit::LocalConfiguration;
    ///                must contain a @c "scheme" sub-configuration that
    ///                describes the ancillary forward interpolation.
    ApproxInverse(const Config& config);
    ~ApproxInverse() override {}

    void print(std::ostream&) const override;
    /// @brief The inverse-interpolation source function space (original target).
    const FunctionSpace& source() const override { return inverse_interp_source_; }
    /// @brief The inverse-interpolation target function space (original source).
    const FunctionSpace& target() const override { return inverse_interp_target_; }

protected:
    

private:
    /// @brief Compute the approximate inverse of the ancillary interpolation matrix.
    /// @param interp_matrix  The halo-exchanged ancillary forward interpolation
    ///                       matrix (target → source direction).
    /// @return               An approximate inverse interpolation matrix
    ///                       (source → target direction), before row normalisation.
    virtual SparseMatrixStorage approx_inverse_transform(const SparseMatrixStorage& interp_matrix) const = 0;

    /// @brief Redistribute the interpolation matrix rows across MPI ranks.
    ///
    /// Each rank of the forward interpolation may reference columns (source
    /// points) owned by remote ranks.  This function uses MPI all-to-all
    /// communication to ensure that every rank receives all matrix rows whose
    /// columns it owns, making the subsequent transpose-based inversion
    /// local.
    /// @param interp_matrix  The local portion of the ancillary interpolation matrix.
    /// @return               A redistributed matrix where each rank holds the
    ///                       rows relevant to its owned source points.
    SparseMatrixStorage halo_exchange(const SparseMatrixStorage& interp_matrix) const;

    /// @brief Normalise each row of the inverse interpolation matrix to sum to one.
    /// @param inverse_interp_matrix  Un-normalised inverse interpolation matrix.
    /// @return                       Row-normalised matrix.
    SparseMatrixStorage normalise_rows(const SparseMatrixStorage& inverse_interp_matrix) const;

    void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
    
    eckit::LocalConfiguration interp_ancillary_scheme_{};

    FunctionSpace inverse_interp_source_{};
    FunctionSpace inverse_interp_target_{};

    FunctionSpace interp_source_{};
    FunctionSpace interp_target_{};
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
