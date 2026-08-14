/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/approx_inverse/LocalPseudoInverse.h"

#include "atlas/library/defines.h"

#if ATLAS_HAVE_EIGEN

#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "eckit/maths/Eigen.h"


namespace atlas {
namespace interpolation {
namespace method {

namespace {

MethodBuilder<LocalPseudoInverse> builder("local-pseudoinverse");

constexpr auto tiny_weight = 1e-12;

// Class to handle a sub-matrix of the interpolation matrix, which is used to compute part of the pseudoinverse.
// Sub matrix contains all rows of interpolation matrix which share a common set of column indices.
class SubMatrix {
public:
    using Index   = LocalPseudoInverse::Index;
    using Value   = LocalPseudoInverse::Value;
    using Triplet = LocalPseudoInverse::Triplet;

    SubMatrix(const std::vector<Index>& col_indices): col_indices_(col_indices) {}

    void add_row(Index row_idx, const std::vector<Value>& row_values) {
        ATLAS_ASSERT(row_values.size() == col_indices_.size());
        row_indices_.emplace_back(row_idx);
        values_.emplace_back(row_values);
    }

    void compute_pseudoinverse_triplets(std::vector<Triplet>& tripletBuffer) const {
        // Create a matrix from the row vectors.
        const auto mat_rows = static_cast<Eigen::Index>(num_rows());
        const auto mat_cols = static_cast<Eigen::Index>(num_cols());

        using MatrixType = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        MatrixType sub_matrix(mat_rows, mat_cols);

        for (Eigen::Index i = 0; i < mat_rows; ++i) {
            for (Eigen::Index j = 0; j < mat_cols; ++j) {
                sub_matrix(i, j) = values_.at(static_cast<std::size_t>(i)).at(static_cast<std::size_t>(j));
            }
        }

        // Compute the pseudoinverse using complete orthogonal decomposition.
        const MatrixType sub_matrix_pseudoinverse = sub_matrix.completeOrthogonalDecomposition().pseudoInverse();

        const auto inv_mat_rows = sub_matrix_pseudoinverse.rows();
        const auto inv_mat_cols = sub_matrix_pseudoinverse.cols();
        const auto squared_norm = sub_matrix_pseudoinverse.squaredNorm();
        const auto norm = squared_norm > 0.
                      ? std::sqrt(inv_mat_cols / (inv_mat_rows * squared_norm))
                      : 0.;

        // Write the pseudoinverse values to triplets.
        for (Eigen::Index i = 0; i < inv_mat_rows; ++i) {
            const auto triplet_row = col_indices_.at(static_cast<std::size_t>(i));
            for (Eigen::Index j = 0; j < inv_mat_cols; ++j) {
                const auto value       = sub_matrix_pseudoinverse(i, j);
                const auto triplet_col = row_indices_.at(static_cast<std::size_t>(j));

                tripletBuffer.emplace_back(triplet_row, triplet_col, value * norm);
            }
        }
    }

    size_t num_rows() const { return row_indices_.size(); }
    size_t num_cols() const { return col_indices_.size(); }


private:
    std::vector<Index> row_indices_{};
    std::vector<Index> col_indices_{};
    std::vector<std::vector<Value>> values_{};
};

// Class to manage a collection of sub-matrices, which are used to compute a pseudoinverse of the interpolation matrix.
// Sub-matrices are stored in a map, keyed to their column indices.
class SubMatrices {
public:
    using Value               = LocalPseudoInverse::Value;
    using Index               = LocalPseudoInverse::Index;
    using SparseMatrixStorage = LocalPseudoInverse::SparseMatrixStorage;
    using Triplet             = LocalPseudoInverse::Triplet;

    explicit SubMatrices(const SparseMatrixStorage& interp_matrix) {
        const auto interp_matrix_view = linalg::make_host_view<Value, Index>(interp_matrix);

        for (std::size_t row_idx = 0; row_idx < interp_matrix_view.rows(); ++row_idx) {
            // Use a map to order columns by index and accumulate values for duplicate columns.
            // Should already be true for interpolation matrices, but not guaranteed.
            auto col_value_map = std::map<Index, Value>{};

            linalg::sparse_matrix_row_for_each(row_idx, interp_matrix_view, [&](Index col, Value value) {
                // Skip tiny weights to avoid numerical issues in the pseudoinverse computation.
                if (std::abs(value) < tiny_weight) {
                    return;
                }
                col_value_map[col] += value;
            });

            if (col_value_map.empty()) {
                continue;
            }

            auto col_indices = std::vector<Index>{};
            col_indices.reserve(col_value_map.size());

            auto values = std::vector<Value>{};
            values.reserve(col_value_map.size());

            for (const auto& [col_idx, value] : col_value_map) {
                col_indices.emplace_back(col_idx);
                values.emplace_back(value);
            }

            auto [it, inserted] = sub_matrices_.try_emplace(col_indices, col_indices);
            it->second.add_row(static_cast<Index>(row_idx), values);
            num_elems_ += values.size();
        }
    }

    void print_underdetermined_warnings() const {
        auto underdetermined_matrices = 0;
        for (const auto& element : sub_matrices_) {
            if (element.second.num_rows() < element.second.num_cols()) {
                ++underdetermined_matrices;
            }
        }

        auto total_elements = sub_matrices_.size();
        mpi::comm().reduceInPlace(underdetermined_matrices, eckit::mpi::Operation::SUM, 0);
        mpi::comm().reduceInPlace(total_elements, eckit::mpi::Operation::SUM, 0);
        if (mpi::rank() == 0 && underdetermined_matrices > 0) {
            Log::warning() << "LocalPseudoInverse: " << underdetermined_matrices << " of " << total_elements
                           << " mesh elements have fewer rows than columns. Consider checking accuracy of regridding."
                           << std::endl;
        }
    }

    std::vector<Triplet> compute_pseudoinverse_triplets() const {
        auto triplets = std::vector<Triplet>{};
        triplets.reserve(num_elems_);

        for (const auto& element : sub_matrices_) {
            element.second.compute_pseudoinverse_triplets(triplets);
        }

        return triplets;
    }

private:
    size_t num_elems_{};
    std::map<std::vector<Index>, SubMatrix> sub_matrices_{};
};

}  // namespace

LocalPseudoInverse::SparseMatrixStorage LocalPseudoInverse::approx_inverse_transform(
    const SparseMatrixStorage& interp_matrix) const {
    const auto sub_matrices = SubMatrices{interp_matrix};
    sub_matrices.print_underdetermined_warnings();
    auto triplets = sub_matrices.compute_pseudoinverse_triplets();
    return linalg::make_sparse_matrix_storage_from_triplets(target().size(), source().size(), triplets);
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas

#endif
