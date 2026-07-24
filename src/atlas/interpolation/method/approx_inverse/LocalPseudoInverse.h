/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/interpolation/method/approx_inverse/ApproxInverse.h"


namespace atlas {
namespace interpolation {
namespace method {


/// @brief Approximate inverse interpolation via local pseudoinverse.
///
/// Groups the rows of the ancillary forward interpolation matrix by their
/// set of non-zero column indices (i.e. rows that share the same set of
/// source stencil points are grouped together).  For each group a dense
/// sub-matrix is formed and its Moore–Penrose pseudoinverse is computed
/// using Eigen's complete orthogonal decomposition.  The resulting
/// sub-matrix inverses are assembled into the global inverse interpolation
/// matrix.
///
/// A warning is emitted (on MPI rank 0) when any sub-matrix has fewer rows
/// than columns (underdetermined system), as the pseudoinverse may be
/// inaccurate in that case.
///
/// Registered under the factory key @c "local-pseudoinverse".
class LocalPseudoInverse : public ApproxInverse {
 public:
  using ApproxInverse::ApproxInverse;
  ~LocalPseudoInverse() override {}

 private:
    /// @brief Build the approximate inverse matrix using per-stencil pseudoinverses.
    /// @param interp_matrix  Halo-exchanged ancillary forward interpolation matrix.
    /// @return               Assembled pseudoinverse interpolation matrix.
    SparseMatrixStorage approx_inverse_transform(
        const SparseMatrixStorage& interp_matrix) const override;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
