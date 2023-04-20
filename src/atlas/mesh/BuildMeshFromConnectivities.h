/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/mesh/Mesh.h"

namespace atlas {
class Grid;
}  // namespace atlas

namespace atlas {

//-----------------------------------------------------------------------------

struct TriangConnectivityData {
    const idx_t local_cell_index;
    const gidx_t global_cell_index;
    const std::array<idx_t, 3> boundary_nodes_of_cell;
};

struct QuadConnectivityData {
    const idx_t local_cell_index;
    const gidx_t global_cell_index;
    const std::array<idx_t, 4> boundary_nodes_of_cell;
};

Mesh build_mesh_from_coordinates_and_connectivities(
    const std::vector<double>& lons,
    const std::vector<double>& lats,
    const std::vector<int>& ghosts,
    const std::vector<int>& global_indices,
    const std::vector<int>& partitions,
    const std::vector<TriangConnectivityData>& triang_connectivities,
    const std::vector<QuadConnectivityData>& quad_connectivities
    );

//-----------------------------------------------------------------------------

}  // namespace atlas
