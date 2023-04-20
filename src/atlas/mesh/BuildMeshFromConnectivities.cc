/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/BuildMeshFromConnectivities.h"

#include "eckit/utils/Hash.h"

//#include "atlas/array/ArrayView.h"
//#include "atlas/array/MakeView.h"
//#include "atlas/grid/Distribution.h"
//#include "atlas/grid/Grid.h"
//#include "atlas/grid/Iterator.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/runtime/Log.h"

using atlas::Mesh;

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

Mesh build_mesh_from_coordinates_and_connectivities(
    const std::vector<double>& lons,
    const std::vector<double>& lats,
    const std::vector<int>& ghosts,
    const std::vector<int>& global_indices,
    const std::vector<int>& partitions,
    const std::vector<TriangConnectivityData>& triang_connectivities,
    const std::vector<QuadConnectivityData>& quad_connectivities
    )
{
  const size_t nb_points = lons.size();
  ATLAS_ASSERT(nb_points == lats.size());
  ATLAS_ASSERT(nb_points == ghosts.size());
  ATLAS_ASSERT(nb_points == global_indices.size());
  ATLAS_ASSERT(nb_points == partitions.size());

  Mesh mesh;

  // nodes

  mesh.nodes().resize(nb_points);
  auto xy        = array::make_view<double, 2>(mesh.nodes().xy());
  auto lonlat    = array::make_view<double, 2>(mesh.nodes().lonlat());
  auto ghost     = array::make_view<int, 1>(mesh.nodes().ghost());
  auto gidx      = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
  auto partition = array::make_view<idx_t, 1>(mesh.nodes().partition());
  auto halo      = array::make_view<int, 1>(mesh.nodes().halo());

  for (size_t i = 0; i < nb_points; ++i) {
    xy(i, size_t(XX)) = lons[i];
    xy(i, size_t(YY)) = lats[i];
    // identity projection, therefore (lon,lat) = (x,y)
    lonlat(i, size_t(LON)) = lons[i];
    lonlat(i, size_t(LAT)) = lats[i];
    ghost(i) = ghosts[i];
    gidx(i) = global_indices[i];
    partition(i) = partitions[i];
  }
  halo.assign(0);

  // cells

  // first, count how many cells of each type are on this processor
  const size_t ntriangs = triang_connectivities.size();
  const size_t nquads = quad_connectivities.size();
  // TODO: accumulate over MPI to optimize away if globally ntriangs or nquads is zero

  mesh.cells().add(new mesh::temporary::Triangle(), ntriangs);
  mesh.cells().add(new mesh::temporary::Quadrilateral(), nquads);

  atlas::mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();
  auto cells_part = array::make_view<int, 1>(mesh.cells().partition());
  auto cells_gidx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
  auto cells_ridx = array::make_view<idx_t, 1>(mesh.cells().remote_index());

  size_t index = 0;
  for (const auto& conn : triang_connectivities) {
    ATLAS_ASSERT(conn.local_cell_index < ntriangs);  // check triangle cells come first
    node_connectivity.set(index, conn.boundary_nodes_of_cell.data());  // triangs and quads
    cells_gidx(index) = conn.global_cell_index;
    cells_ridx(index) = conn.local_cell_index;
    index++;
  }
  for (const auto& conn : quad_connectivities) {
    ATLAS_ASSERT(conn.local_cell_index >= ntriangs);  // check quad cells come last
    node_connectivity.set(index, conn.boundary_nodes_of_cell.data());  // triangs and quads
    cells_gidx(index) = conn.global_cell_index;
    cells_ridx(index) = conn.local_cell_index;
    index++;
  }
  ATLAS_ASSERT(index == ntriangs + nquads);

  cells_part.assign(atlas::mpi::comm().rank());

  return mesh;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
