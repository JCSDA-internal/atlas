/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/approx_inverse/ApproxInverse.h"

#include <cmath>
#include <cstddef>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

#include "atlas/array/native/NativeArrayView.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "eckit/config/LocalConfiguration.h"

namespace atlas {
namespace interpolation {
namespace method {

ApproxInverse::ApproxInverse(const Config& config): Method(config) {
    const auto* conf = dynamic_cast<const eckit::LocalConfiguration*>(&config);
    ATLAS_ASSERT(conf, "config must be derived from eckit::LocalConfiguration");
    interp_ancillary_scheme_     = conf->getSubConfiguration("scheme");
}

void ApproxInverse::print(std::ostream&) const {
    ATLAS_NOTIMPLEMENTED;
}

void ApproxInverse::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    ATLAS_TRACE("atlas::interpolation::method::ApproxInverse::do_setup()");

    inverse_interp_source_ = source;
    inverse_interp_target_ = target;

    interp_source_ = target;
    interp_target_ = source;

    if (inverse_interp_target_.size() == 0) {
        return;
    }

    const auto interp = Interpolation(interp_ancillary_scheme_, interp_source_, interp_target_);

    auto matrix = interpolation::MatrixCache(interp).matrix();

    matrix = halo_exchange(matrix);
    matrix = approx_inverse_transform(matrix);
    matrix = normalise_rows(matrix);

    setMatrix(std::move(matrix));
}

void ApproxInverse::do_setup(const Grid& source, const Grid& target, const Cache&) {
    ATLAS_NOTIMPLEMENTED;
}

void ApproxInverse::do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) {
    ATLAS_NOTIMPLEMENTED;
}

ApproxInverse::SparseMatrixStorage ApproxInverse::halo_exchange(const SparseMatrixStorage& interp_matrix) const {
    const auto& comm              = mpi::comm();
    const auto interp_matrix_view = linalg::make_host_view<Value, Index>(interp_matrix);

    struct Views {
        explicit Views(const FunctionSpace& functionSpace):
            remote_index(array::make_indexview<idx_t, 1>(functionSpace.remote_index())),
            partition(array::make_view<int, 1>(functionSpace.partition())),
            ghost(array::make_view<int, 1>(functionSpace.ghost())) {}

        array::IndexView<idx_t, 1> remote_index;
        array::ArrayView<int, 1> partition;
        array::ArrayView<int, 1> ghost;
    };

    const auto source_views = Views{interp_source_};
    const auto target_views = Views{interp_target_};

    using PartRidx       = std::pair<idx_t, idx_t>;
    using PartRidxBuffer = std::vector<std::vector<PartRidx>>;
    using ValueBuffer    = std::vector<std::vector<Value>>;

    // Note: in practice (partition, remote_index) provides a more reliable unique identifier than global_index.
    auto global_row_send_buffer = PartRidxBuffer{comm.size()};
    auto global_col_send_buffer = PartRidxBuffer{comm.size()};
    auto value_send_buffer      = ValueBuffer{comm.size()};

    for (std::size_t row = 0; row < interp_matrix_view.rows(); ++row) {
        // Ghost rows should not be present in interpolation matrix.
        if (target_views.ghost(row)) {
            continue;
        }

        // Find all column ranks in row.
        auto col_rank_set = std::set<int>{};
        linalg::sparse_matrix_row_for_each(row, interp_matrix_view,
                                           [&](Index col, Value) { col_rank_set.insert(source_views.partition(col)); });

        // Send full row to all ranks in col_rank_set.
        linalg::sparse_matrix_row_for_each(row, interp_matrix_view, [&](Index col, Value weight) {
            for (auto rank : col_rank_set) {
                global_row_send_buffer[rank].emplace_back(target_views.partition(row), target_views.remote_index(row));
                global_col_send_buffer[rank].emplace_back(source_views.partition(col), source_views.remote_index(col));
                value_send_buffer[rank].emplace_back(weight);
            }
        });
    }

    // Perform MPI communication.
    const auto do_all_to_all = [&](const auto& send_buffer) {
        auto recv_buffer = decltype(send_buffer){comm.size()};
        comm.allToAll(send_buffer, recv_buffer);
        return recv_buffer;
    };

    const auto global_row_recv_buffer = do_all_to_all(global_row_send_buffer);
    const auto global_col_recv_buffer = do_all_to_all(global_col_send_buffer);
    const auto value_recv_buffer      = do_all_to_all(value_send_buffer);

    // Create a map of local indices, keyed to (partition, remote_index) pair.
    const auto make_part_ridx_to_lidx_map = [&](const Views& views) {
        auto part_ridx_to_lidx_map = std::map<PartRidx, idx_t>{};

        for (idx_t lidx = 0; lidx < views.partition.size(); ++lidx) {
            const auto part_ridx = PartRidx{views.partition(lidx), views.remote_index(lidx)};
            part_ridx_to_lidx_map.emplace(part_ridx, lidx);
        }
        return part_ridx_to_lidx_map;
    };

    const auto source_part_ridx_to_lidx_map = make_part_ridx_to_lidx_map(source_views);
    const auto target_part_ridx_to_lidx_map = make_part_ridx_to_lidx_map(target_views);

    auto triplets = std::vector<Triplet>{};
    triplets.reserve(std::accumulate(global_row_recv_buffer.begin(), global_row_recv_buffer.end(), size_t{0},
                                     [](std::size_t sum, const auto& buffer) { return sum + buffer.size(); }));

    std::size_t missing_rows = 0;
    std::size_t missing_cols = 0;

    for (std::size_t rank = 0; rank < comm.size(); ++rank) {
        const auto& row_sub_buffer   = global_row_recv_buffer[rank];
        const auto& col_sub_buffer   = global_col_recv_buffer[rank];
        const auto& value_sub_buffer = value_recv_buffer[rank];
        ATLAS_ASSERT(row_sub_buffer.size() == col_sub_buffer.size());
        ATLAS_ASSERT(row_sub_buffer.size() == value_sub_buffer.size());

        // Look up local indices for each (partition, remote_index) pair and accumulate triplets.
        for (std::size_t i = 0; i < row_sub_buffer.size(); ++i) {
            const auto row_part_ridx = row_sub_buffer[i];
            const auto col_part_ridx = col_sub_buffer[i];
            const auto value         = value_sub_buffer[i];

            const auto row_iter = target_part_ridx_to_lidx_map.find(row_part_ridx);
            const auto col_iter = source_part_ridx_to_lidx_map.find(col_part_ridx);

            auto triplet_found = true;
            if (row_iter == target_part_ridx_to_lidx_map.end()) {
                ++missing_rows;
                triplet_found = false;
            }
            if (col_iter == source_part_ridx_to_lidx_map.end()) {
                ++missing_cols;
                triplet_found = false;
            }
            if (!triplet_found) {
                continue;
            }

            triplets.emplace_back(row_iter->second, col_iter->second, value);
        }
    }

    comm.allReduceInPlace(missing_rows, eckit::mpi::Operation::SUM);
    comm.allReduceInPlace(missing_cols, eckit::mpi::Operation::SUM);

    if (missing_rows > 0 || missing_cols > 0) {
        throw_Exception("ApproxInverse::halo_exchange: Missing halo points. Missed rows=" +
                            std::to_string(missing_rows) + ", missed columns=" + std::to_string(missing_cols) +
                            ". Try increase source functionspace halo size for rows and target functionspace halo "
                            "size for columns.",
                        Here());
    }

    return linalg::make_sparse_matrix_storage_from_triplets(interp_target_.size(), interp_source_.size(), triplets);
}

ApproxInverse::SparseMatrixStorage ApproxInverse::normalise_rows(
    const ApproxInverse::SparseMatrixStorage& inverse_interp_matrix) const {
    const auto interp_matrix_view = linalg::make_host_view<Value, Index>(inverse_interp_matrix);

    auto triplets = std::vector<Triplet>{};
    triplets.reserve(interp_matrix_view.nnz());

    const auto target_ghost_view = array::make_view<int, 1>(inverse_interp_target_.ghost());

    for (std::size_t row_idx = 0; row_idx < interp_matrix_view.rows(); ++row_idx) {
        if (target_ghost_view(row_idx)) {
            continue;
        }

        Value row_sum = 0.;
        linalg::sparse_matrix_row_for_each(row_idx, interp_matrix_view, [&](Value value) { row_sum += value; });

        ATLAS_ASSERT(std::isfinite(row_sum) && row_sum > 0., "Row sum for row " + std::to_string(row_idx) +
                                                                 " should be positive and finite, but is " +
                                                                 std::to_string(row_sum));

        const auto norm = 1. / row_sum;
        linalg::sparse_matrix_row_for_each(row_idx, interp_matrix_view, [&](Index row, Index col, Value value) {
            triplets.emplace_back(row, col, value * norm);
        });
    }

    return linalg::make_sparse_matrix_storage_from_triplets(inverse_interp_target_.size(),
                                                            inverse_interp_source_.size(), triplets);
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
