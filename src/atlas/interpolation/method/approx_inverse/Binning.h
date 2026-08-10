/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "atlas/interpolation/method/approx_inverse/ApproxInverse.h"


namespace atlas {
namespace interpolation {
namespace method {


/// @brief Approximate inverse interpolation via area-weighted binning.
///
/// Computes the approximate inverse of a forward interpolation matrix by
/// transposing it and weighting each entry by the area of the corresponding
/// source node.  For cubed-sphere grids with a @c NodeColumns function space
/// the nodal area is obtained from the mesh; for all other grid types equal
/// weights are used.
///
/// @deprecated Prefer @c LocalPseudoInverse for better variance preservation.
///
/// Registered under the factory key @c "binning".
class Binning : public ApproxInverse {
 public:
  using ApproxInverse::ApproxInverse;
  ~Binning() override {}

 private:
    /// @brief Build the approximate inverse matrix by transposing the
    ///        interpolation matrix and multiplying each weight by the
    ///        source-node area.
    /// @param interp_matrix  Halo-exchanged ancillary forward interpolation matrix.
    /// @return               Area-weighted transposed (binned) matrix.
    SparseMatrixStorage approx_inverse_transform(
        const SparseMatrixStorage& interp_matrix) const override;

    /// @brief Retrieve per-node area weights for the source function space.
    ///
    /// Returns nodal areas from the cubed-sphere mesh when the source is a
    /// @c NodeColumns function space on a @c CubedSphereGrid, otherwise
    /// returns a vector of ones (equal weights).
    /// @return  Vector of length @c source().size() containing the area weight
    ///          for each source node.
    std::vector<double> get_area_weights() const;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
