/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/binning/Binning.h"

#include <climits>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>

#include "atlas/array/native/NativeArrayView.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/library/config.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"
#include "atlas/mesh/actions/GetCubedSphereNodalArea.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/maths/Eigen.h"


namespace atlas {
namespace interpolation {
namespace method {

namespace {

MethodBuilder<Binning> __builder("binning");

using TripletType = linalg::Triplet<Binning::ValueType, Binning::IndexType>;
static_assert(std::is_trivially_copyable_v<TripletType>);

class MpiBuffer {
public:
    MpiBuffer(): buffer_{mpi::comm().size()} {}

    std::size_t size(std::size_t rank) const {
        const auto n = buffer_.at(rank).size();
        ATLAS_ASSERT(n % sizeof(TripletType) == 0, "Buffer size not a multiple of Triplet size");
        return n / sizeof(TripletType);
    }

    void pushBack(std::size_t rank, const TripletType& triplet) {
        auto& subBuffer                        = buffer_.at(rank);
        alignas(TripletType) auto appendBuffer = std::array<char, sizeof(TripletType)>{};
        std::memcpy(appendBuffer.data(), &triplet, sizeof(TripletType));
        subBuffer.insert(subBuffer.end(), appendBuffer.begin(), appendBuffer.end());
    }

    auto get(std::size_t rank) const {
        const auto& subBuffer = buffer_.at(rank);
        return [&subBuffer](size_t index) -> TripletType {
            TripletType triplet{};
            std::memcpy(&triplet, subBuffer.data() + index * sizeof(TripletType), sizeof(TripletType));
            return triplet;
        };
    }

    MpiBuffer allToAll() const {
        auto recvBuffer = MpiBuffer{};
        mpi::comm().allToAll(buffer_, recvBuffer.buffer_);
        return recvBuffer;
    }

private:
    std::vector<std::vector<char>> buffer_{};
};

}  // namespace

Binning::Binning(const Config& config): Method(config) {
    const auto* conf = dynamic_cast<const eckit::LocalConfiguration*>(&config);
    ATLAS_ASSERT(conf, "config must be derived from eckit::LocalConfiguration");
    interpAncillaryScheme_ = conf->getSubConfiguration("scheme");
}


void Binning::do_setup(const Grid& source, const Grid& target, const Cache&) {
    ATLAS_NOTIMPLEMENTED;
}

void Binning::do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) {
    ATLAS_NOTIMPLEMENTED;
}


void Binning::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    ATLAS_TRACE("atlas::interpolation::method::Binning::do_setup()");

    source_ = source;
    target_ = target;

    if (target_.size() == 0) {
        return;
    }

    const auto interp = Interpolation(interpAncillaryScheme_, target_, source_);

    auto matrix = interpolation::MatrixCache(interp).matrix();

    // matrix = cleanInterpMatrix(matrix);
    matrix = haloExchange(matrix);
    matrix = approxInverseTransform(matrix);
    matrix = normaliseRows(matrix);

    setMatrix(std::move(matrix));
}


void Binning::print(std::ostream&) const {
    ATLAS_NOTIMPLEMENTED;
}


Binning::SparseMatrixStorage Binning::cleanInterpMatrix(const Binning::SparseMatrixStorage& interpMatrix) const {
    const auto interpMatrixView = linalg::make_host_view<ValueType, IndexType>(interpMatrix);

    auto triplets = std::vector<TripletType>{};
    triplets.reserve(interpMatrixView.nnz());

    const auto targetGhostView     = array::make_view<int, 1>(target_.ghost());
    const auto sourcePartitionView = array::make_view<int, 1>(source_.partition());
    const auto sourceRidxView      = array::make_indexview<idx_t, 1>(source_.remote_index());

    // Resolve (part, ridx) -> local index dengeneracy in source function space.
    auto sourceRemoteToLocalMap = std::map<std::pair<int, idx_t>, idx_t>{};
    for (idx_t localIdx = 0; localIdx < source_.size(); ++localIdx) {
        const int part   = sourcePartitionView(localIdx);
        const idx_t ridx = sourceRidxView(localIdx);
        sourceRemoteToLocalMap.insert({{part, ridx}, localIdx});
    }


    // Remove rows that map to ghost elements in the target function space.
    for (std::size_t rowIdx = 0; rowIdx < interpMatrixView.rows(); ++rowIdx) {
        // Skip rows that map to ghost elements. They shouldn't be here.
        if (targetGhostView(rowIdx)) {
            continue;
        }

        linalg::sparse_matrix_for_each_row(rowIdx, interpMatrixView,
                                           [&](IndexType row, IndexType col, ValueType value) {
                                               const int part   = sourcePartitionView(col);
                                               const idx_t ridx = sourceRidxView(col);

                                               const auto lidx = sourceRemoteToLocalMap.at({part, ridx});
                                               triplets.emplace_back(row, lidx, value);
                                           });
    }

    return linalg::make_sparse_matrix_storage_from_triplets(static_cast<IndexType>(source_.size()),
                                                            static_cast<IndexType>(target_.size()), triplets);
}

Binning::SparseMatrixStorage Binning::haloExchange(const Binning::SparseMatrixStorage& interpMatrix) const {
    const auto interpMatrixView = linalg::make_host_view<ValueType, IndexType>(interpMatrix);

    struct Views {
        Views(const FunctionSpace& functionSpace):
            remote_index{array::make_indexview<idx_t, 1>(functionSpace.remote_index())},
            partition{array::make_view<int, 1>(functionSpace.partition())},
            ghost{array::make_view<int, 1>(functionSpace.ghost())} {}
        array::IndexView<idx_t, 1> remote_index;
        array::ArrayView<int, 1> partition;
        array::ArrayView<int, 1> ghost;
    };

    // Note: source_ and target_ for interp matrix are swapped.
    const auto interpSource = target_;
    const auto interpTarget = source_;


    const Views sourceViews{interpSource};
    const Views targetViews{interpTarget};

    using PartRidxPair = std::pair<int, idx_t>;

    auto global_row_send_buffer   = std::vector<std::vector<PartRidxPair>>{mpi::comm().size()};
    auto global_col_send_buffer   = std::vector<std::vector<PartRidxPair>>{mpi::comm().size()};
    auto global_value_send_buffer = std::vector<std::vector<ValueType>>{mpi::comm().size()};

    linalg::sparse_matrix_for_each(interpMatrixView, [&](IndexType row, IndexType col, ValueType weight) {
        // Ghost rows should not be present in interpolation matrix.
        if (targetViews.ghost(row)) {
            return;
        }

        // Find all column ranks in row
        auto col_ranks = std::set<int>{};
        linalg::sparse_matrix_for_each_row(row, interpMatrixView, [&](IndexType row, IndexType col, ValueType weight) {
            col_ranks.insert(sourceViews.partition(col));
        });

        // Send entire row to all columns.
        linalg::sparse_matrix_for_each_row(row, interpMatrixView, [&](IndexType row, IndexType col, ValueType weight) {
            for (auto rank : col_ranks) {
                global_row_send_buffer[rank].push_back({targetViews.partition(row), targetViews.remote_index(row)});
                global_col_send_buffer[rank].push_back({sourceViews.partition(col), sourceViews.remote_index(col)});
                global_value_send_buffer[rank].push_back(weight);
            }
        });
    });

    const auto do_all_to_all = [&](const auto& send_buffer) {
        auto recv_buffer = decltype(send_buffer){mpi::comm().size()};
        mpi::comm().allToAll(send_buffer, recv_buffer);
        return recv_buffer;
    };

    const auto global_row_recv_buffer   = do_all_to_all(global_row_send_buffer);
    const auto global_col_recv_buffer   = do_all_to_all(global_col_send_buffer);
    const auto global_value_recv_buffer = do_all_to_all(global_value_send_buffer);


    const auto make_global_to_local_map = [&](const Views& views) {
        auto global_to_local_map = std::map<PartRidxPair, idx_t>{};
        for (idx_t localIdx = 0; localIdx < views.partition.size(); ++localIdx) {
            const auto part_ridx = PartRidxPair{views.partition(localIdx), views.remote_index(localIdx)};
            global_to_local_map.insert({part_ridx, localIdx});
        }
        return global_to_local_map;
    };

    const auto source_global_to_local_map = make_global_to_local_map(sourceViews);
    const auto target_global_to_local_map = make_global_to_local_map(targetViews);


    ATLAS_ASSERT(global_row_recv_buffer.size() == mpi::comm().size());
    ATLAS_ASSERT(global_col_recv_buffer.size() == mpi::comm().size());
    ATLAS_ASSERT(global_value_recv_buffer.size() == mpi::comm().size());

    auto triplets = std::vector<TripletType>{};

    for (std::size_t rank = 0; rank < mpi::comm().size(); ++rank) {
        const auto& row_buffer   = global_row_recv_buffer[rank];
        const auto& col_buffer   = global_col_recv_buffer[rank];
        const auto& value_buffer = global_value_recv_buffer[rank];

        ATLAS_ASSERT(row_buffer.size() == col_buffer.size());
        ATLAS_ASSERT(row_buffer.size() == value_buffer.size());

        for (std::size_t i = 0; i < row_buffer.size(); ++i) {
            const auto part_ridx_row = row_buffer[i];
            const auto part_ridx_col = col_buffer[i];
            const auto weight        = value_buffer[i];

            const auto row_iter = target_global_to_local_map.find(part_ridx_row);
            ATLAS_ASSERT_MSG(
                row_iter != target_global_to_local_map.end(),
                "Source local index not found for global index. Try increasing source functionspace halo size.");
            const auto row = row_iter->second;

            const auto col_iter = source_global_to_local_map.find(part_ridx_col);
            ATLAS_ASSERT_MSG(
                col_iter != source_global_to_local_map.end(),
                "Target local index not found for global index. Try increasing target functionspace halo size.");
            const auto col = col_iter->second;

            triplets.emplace_back(row, col, weight);
        }
    }

    // Note: source_ and target_ for interp matrix are swapped.
    return linalg::make_sparse_matrix_storage_from_triplets(static_cast<IndexType>(interpTarget.size()),
                                                            static_cast<IndexType>(interpSource.size()), triplets);
}

// Binning::SparseMatrixStorage Binning::approxInverseTransform(const Binning::SparseMatrixStorage& interpMatrix) const {
//     const auto interpMatrixView = linalg::make_host_view<ValueType, IndexType>(interpMatrix);

//     auto triplets = std::vector<TripletType>{};
//     triplets.reserve(interpMatrixView.nnz());

//     const auto areaWeights = getAreaWeights();

//     // Approximate inverse transform by transposing the interpolation matrix.
//     linalg::sparse_matrix_for_each(interpMatrixView, [&](IndexType row, IndexType col, ValueType weight) {
//         triplets.emplace_back(col, row, weight * areaWeights.at(col));
//     });

//     return linalg::make_sparse_matrix_storage_from_triplets(static_cast<IndexType>(target_.size()),
//                                                             static_cast<IndexType>(source_.size()), triplets);

// }

Binning::SparseMatrixStorage Binning::approxInverseTransform(const Binning::SparseMatrixStorage& interpMatrix) const {
    const auto interpMatrixView = linalg::make_host_view<ValueType, IndexType>(interpMatrix);

    auto triplets = std::vector<TripletType>{};
    triplets.reserve(interpMatrixView.nnz());

    using Value         = Binning::ValueType;
    using ColIndex      = Binning::IndexType;
    using ColIndices    = std::vector<ColIndex>;
    using RowVectorData = std::array<Binning::ValueType, 4>;
    using RowIndex      = Binning::IndexType;
    using RowVector     = std::pair<RowIndex, RowVectorData>;
    using RowVectors    = std::vector<RowVector>;
    using RowVectorMap  = std::map<ColIndices, RowVectors>;

    // Store a map of accumulated row-vectors, keyed by vector of sorted column indices.
    auto row_vector_map = RowVectorMap{};

    for (std::size_t row_idx = 0; row_idx < interpMatrixView.rows(); ++row_idx) {
        // use a std::map to order columns by index and accumulate values for duplicate columns
        using RowMap = std::map<ColIndex, Value>;
        auto row_map = RowMap{};

        linalg::sparse_matrix_for_each_row(row_idx, interpMatrixView,
                                           [&](RowIndex, ColIndex col, Value value) { row_map[col] += value; });

        if (row_map.empty()) {
            continue;
        }

        if (row_map.size() > 4 || row_map.size() < 3) {
            // print row, col, value
            for (const auto& [col, value] : row_map) {
                std::cout << "Row: " << row_idx << ", Col: " << col << ", Value: " << value << std::endl;
            }

            ATLAS_ASSERT(false, "Row " + std::to_string(row_idx) + " has " + std::to_string(row_map.size()) +
                                    " non-zero entries. Expected between 3 and 4.");
        }


        // Make sure no more that 4 non-zero entries in a row.
        ATLAS_ASSERT(row_map.size() <= 4, "Row " + std::to_string(row_idx) + " has more than 4 non-zero entries");

        // Make sure there are at least 3 non-zero entries in a row.
        ATLAS_ASSERT(row_map.size() >= 3, "Row " + std::to_string(row_idx) + " has less than 3 non-zero entries");


        // Extract column indices and values from the row_map
        auto col_indices = ColIndices{};
        col_indices.reserve(row_map.size());
        std::transform(row_map.begin(), row_map.end(), std::back_inserter(col_indices),
                       [](const auto& pair) { return pair.first; });

        auto row_vector_data = RowVectorData{};
        std::transform(row_map.begin(), row_map.end(), row_vector_data.data(),
                       [](const auto& pair) { return pair.second; });

        row_vector_map[col_indices].emplace_back(static_cast<RowIndex>(row_idx), row_vector_data);
    }

    auto pseudoinverse_triplets = std::vector<linalg::Triplet<Value, Binning::IndexType>>{};
    pseudoinverse_triplets.reserve(interpMatrixView.nnz());

    // Now calculate the pseudo-inverse for each unique set of column indices and accumulate the triplets.
    for (const auto& [col_indices, row_vectors] : row_vector_map) {
        // Get row indices:
        auto row_indices = std::vector<RowIndex>{};
        row_indices.reserve(row_vectors.size());
        std::transform(row_vectors.begin(), row_vectors.end(), std::back_inserter(row_indices),
                       [](const auto& pair) { return pair.first; });


        // Create a matrix from the row vectors
        const auto num_rows = static_cast<Eigen::Index>(row_vectors.size());
        const auto num_cols = static_cast<Eigen::Index>(col_indices.size());
        Eigen::Matrix<Binning::ValueType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> sub_matrix(num_rows,
                                                                                                      num_cols);

        for (Eigen::Index i = 0; i < num_rows; ++i) {
            for (Eigen::Index j = 0; j < num_cols; ++j) {
                sub_matrix(i, j) = row_vectors.at(i).second.at(j);
            }
        }

        // Compute the pseudoinverse of interpolation polygon using Eigen's complete orthogonal decomposition
        const auto poly_pseudoinverse = sub_matrix.completeOrthogonalDecomposition().pseudoInverse();

        ATLAS_ASSERT(poly_pseudoinverse.rows() == num_cols, "Pseudoinverse rows do not match number of columns");
        ATLAS_ASSERT(poly_pseudoinverse.cols() == num_rows, "Pseudoinverse cols do not match number of rows");

        const auto target_ghost_view = array::make_view<int, 1>(target_.ghost());

        for (Eigen::Index i = 0; i < poly_pseudoinverse.rows(); ++i) {
            const auto triplet_row = col_indices.at(i);
            if (target_ghost_view(triplet_row)) {
                continue;
            }

            // Normalise by row RMS, to minimise total energy of pseudoinverse matrix.
            const auto norm = std::sqrt(poly_pseudoinverse.rows()) / poly_pseudoinverse.row(i).norm();

            for (Eigen::Index j = 0; j < poly_pseudoinverse.cols(); ++j) {
                const auto value       = poly_pseudoinverse(i, j);
                const auto triplet_col = row_indices.at(j);

                pseudoinverse_triplets.emplace_back(triplet_row, triplet_col, value * norm);
            }
        }
    }
    return linalg::make_sparse_matrix_storage_from_triplets(target_.size(), source_.size(), pseudoinverse_triplets);
}

Binning::SparseMatrixStorage Binning::normaliseRows(const Binning::SparseMatrixStorage& invInterpMatrix) const {
    const auto interpMatrixView = linalg::make_host_view<ValueType, IndexType>(invInterpMatrix);

    auto triplets = std::vector<TripletType>{};
    triplets.reserve(interpMatrixView.nnz());

    // Normalise rows.
    for (std::size_t rowIdx = 0; rowIdx < interpMatrixView.rows(); ++rowIdx) {
        ValueType rowSum = 0.;
        linalg::sparse_matrix_for_each_row(rowIdx, interpMatrixView,
                                           [&](IndexType row, IndexType col, ValueType value) { rowSum += value; });

        linalg::sparse_matrix_for_each_row(
            rowIdx, interpMatrixView,
            [&](IndexType row, IndexType col, ValueType value) { triplets.emplace_back(row, col, value / rowSum); });
    }

    return linalg::make_sparse_matrix_storage_from_triplets(static_cast<IndexType>(target_.size()),
                                                            static_cast<IndexType>(source_.size()), triplets);
};


std::vector<double> Binning::getAreaWeights() const {
    const auto csGrid    = CubedSphereGrid(source_.grid());
    auto ncFunctionSpace = functionspace::NodeColumns(source_);

    if (!(csGrid && ncFunctionSpace)) {
        // If not a cubed sphere grid and a node columns function space, return equal weights.
        return std::vector<double>(source_.size(), 1.);
    }

    auto mesh                  = ncFunctionSpace.mesh();
    const auto nodeWeights     = mesh::actions::GetCubedSphereNodalArea()(mesh);
    const auto nodeWeightsView = array::make_view<double, 1>(nodeWeights);

    auto weightVector = std::vector<double>(source_.size());
    for (idx_t i = 0; i < source_.size(); ++i) {
        weightVector[i] = nodeWeightsView(i);
    }
    return weightVector;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
