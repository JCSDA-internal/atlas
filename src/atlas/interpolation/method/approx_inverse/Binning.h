/*
 * (C) Crown Copyright 2024 Met Office
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


class Binning : public ApproxInverse {
 public:
  using ApproxInverse::ApproxInverse;
  ~Binning() override {}

 private:
    SparseMatrixStorage approx_inverse_transform(
        const SparseMatrixStorage& interp_matrix) const override;
    std::vector<double> get_area_weights() const;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
