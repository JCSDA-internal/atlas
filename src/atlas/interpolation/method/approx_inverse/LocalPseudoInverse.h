/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/interpolation/method/approx_inverse/ApproxInverse.h"

#if ATLAS_HAVE_EIGEN

namespace atlas {
namespace interpolation {
namespace method {


/// @brief Compact-support approximate inverse interpolation via local pseudoinverses.
///
/// This method builds an approximate left inverse of the forward interpolation
/// operator by working stencil-by-stencil rather than globally. Rows of the
/// ancillary forward interpolation matrix are first grouped by their identical
/// sparsity pattern (that is, by the set of non-zero source columns they touch).
/// For each group a dense sub-matrix is formed, its Moore-Penrose pseudoinverse
/// is computed with Eigen's complete orthogonal decomposition, and the result
/// is truncated back into the original compact-support pattern.
///
/// The resulting operator is sparse and local by construction. On each block,
/// and provided that the block is not underdetermined, the construction behaves
/// like an exact left inverse for the corresponding forward stencil. In the
/// underdetermined case the pseudoinverse is still formed, but it is only an
/// approximate inverse and a warning is emitted during assembly.
///
/// After each block pseudoinverse is assembled, the block is scaled so that its
/// Frobenius norm is proportional to @f$\sqrt{m/n}@f$, where @f$m@f$ is the
/// number of rows in the block and @f$n@f$ is the number of columns. This
/// provides a simple blockwise energy normalisation and helps keep the assembled
/// rows visually coherent when a block spans mixed element types such as points,
/// lines, triangles, and quadrilaterals.
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

#endif
