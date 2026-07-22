/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/approx_inverse/Binning.h"

#include <vector>

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"
#include "atlas/mesh/actions/GetCubedSphereNodalArea.h"


namespace atlas {
namespace interpolation {
namespace method {

namespace {

MethodBuilder<Binning> builder("binning");

using TripletType = linalg::Triplet<Binning::Value, Binning::Index>;

}  // namespace

Binning::SparseMatrixStorage Binning::approx_inverse_transform(
    const Binning::SparseMatrixStorage& interp_matrix) const {
    const auto interp_matrix_view = linalg::make_host_view<Value, Index>(interp_matrix);

    auto triplets = std::vector<TripletType>{};
    triplets.reserve(interp_matrix_view.nnz());

    const auto area_weights = get_area_weights();

    // Approximate inverse transform by transposing the interpolation matrix.
    linalg::sparse_matrix_for_each(interp_matrix_view, [&](Index row, Index col, Value weight) {
        triplets.emplace_back(col, row, weight * area_weights.at(col));
    });

    return linalg::make_sparse_matrix_storage_from_triplets(static_cast<Index>(target().size()),
                                                             static_cast<Index>(source().size()), triplets);
}

std::vector<double> Binning::get_area_weights() const {
    const auto cs_grid    = CubedSphereGrid(source().grid());
    auto nc_function_space = functionspace::NodeColumns(source());

    if (!(cs_grid && nc_function_space)) {
        // If not a cubed sphere grid and a node columns function space, return equal weights.
        return std::vector<double>(source().size(), 1.);
    }

    auto mesh                  = nc_function_space.mesh();
    const auto node_weights     = mesh::actions::GetCubedSphereNodalArea()(mesh);
    const auto node_weights_view = array::make_view<double, 1>(node_weights);

    auto weight_vector = std::vector<double>(source().size());
    for (idx_t i = 0; i < source().size(); ++i) {
        weight_vector[i] = node_weights_view(i);
    }
    return weight_vector;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
