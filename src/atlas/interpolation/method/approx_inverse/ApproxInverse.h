/*
 * (C) Crown Copyright 2024 Met Office
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

class ApproxInverse : public Method {
public:
    using Value               = double;
    using Index               = int;
    using Triplet             = linalg::Triplet<Value, Index>;
    using SparseMatrixStorage = linalg::SparseMatrixStorage;
    using SparseMatrixView    = linalg::SparseMatrixView<Value, Index>;

    ApproxInverse(const Config& config);
    ~ApproxInverse() override {}

    void print(std::ostream&) const override;
    const FunctionSpace& source() const override { return inverse_interp_source_; }
    const FunctionSpace& target() const override { return inverse_interp_target_; }

protected:
    

private:
    virtual SparseMatrixStorage approx_inverse_transform(const SparseMatrixStorage& interp_matrix) const = 0;
    SparseMatrixStorage halo_exchange(const SparseMatrixStorage& interp_matrix) const;
    SparseMatrixStorage normalise_rows(const SparseMatrixStorage& inverse_interp_matrix) const;

    void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
    void do_setup(const Grid& source, const Grid& target, const Cache&) override;
    void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) override;

    
    
    eckit::LocalConfiguration interp_ancillary_scheme_{};

    FunctionSpace inverse_interp_source_{};
    FunctionSpace inverse_interp_target_{};

    FunctionSpace interp_source_{};
    FunctionSpace interp_target_{};
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
