/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/ElementKDTree.h"

#include "eckit/config/Resource.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"
#include "atlas/util/Topology.h"

namespace atlas {
namespace interpolation {
namespace method {

util::IndexKDTree3D build_element_kdtree(const Mesh& mesh, const Field& field_centres) {
    ATLAS_TRACE();
    const auto centres = array::make_view<double, 2>(field_centres);
    const auto flags   = array::make_view<int, 1>(mesh.cells().flags());
    auto include_element = [&](idx_t e) {
        using util::Topology;
        return not Topology::view(flags(e)).check(Topology::INVALID);
    };

    // When false, values are inserted one by one without reserve(), which is slower for the eckit backend
    static bool fastBuildKDTrees = eckit::Resource<bool>("$ATLAS_FAST_BUILD_KDTREES", true);

    util::IndexKDTree3D tree;
    const idx_t nb_cells = centres.shape(0);
    if (fastBuildKDTrees) {
        tree.reserve(nb_cells);
    }
    for (idx_t j = 0; j < nb_cells; ++j) {
        if (include_element(j)) {
            tree.insert(PointXYZ{centres(j, XX), centres(j, YY), centres(j, ZZ)}, j);
        }
    }
    if (fastBuildKDTrees) {
        tree.build();
    }
    return tree;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
