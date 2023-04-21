/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iterator>

#include "atlas/mesh/BuildMeshFromConnectivities.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/output/Gmsh.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::output;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_cs_c2_mesh_parallel") {

    ATLAS_ASSERT(mpi::comm().size() == 6);
    const int rank = mpi::comm().rank();

    // coordinates of C2 lfric cubed-sphere grid: grid("CS-LFR-2");
    const std::vector<double> global_lons = {{
      337.5, 22.5, 337.5, 22.5,    // +x
      67.5, 112.5, 67.5, 112.5,    // +y
      202.5, 202.5, 157.5, 157.5,  // -x
      292.5, 292.5, 247.5, 247.5,  // -y
      315, 45, 225, 135,           // +z
      315, 225, 45, 135}};         // -z

    const std::vector<double> global_lats = {{
      -20.941, -20.941, 20.941, 20.941,
      -20.941, -20.941, 20.941, 20.941,
      -20.941, 20.941, -20.941, 20.941,
      -20.941, 20.941, -20.941, 20.941,
      59.6388, 59.6388, 59.6388, 59.6388,
      -59.6388, -59.6388, -59.6388, -59.6388}};

    const std::vector<int> global_partitions = {{
      0, 0, 0, 0,
      1, 1, 1, 1,
      2, 2, 2, 2,
      3, 3, 3, 3,
      4, 4, 4, 4,
      5, 5, 5, 5}};

    const auto map_tri_to_local = [&rank](const std::vector<int>& global_indices,
                                     const std::array<int, 3>& global_to_map) {
      std::array<int, 3> result;
      for (size_t i = 0; i < 3; ++i) {
        const auto& it = find(global_indices.begin(), global_indices.end(), global_to_map[i]);
        ATLAS_ASSERT(it != global_indices.end());
        result[i] = std::distance(global_indices.begin(), it);
      }
      return result;
    };
    const auto map_quad_to_local = [&rank](const std::vector<int>& global_indices,
                                      const std::array<int, 4>& global_to_map) {
      std::array<int, 4> result;
      for (size_t i = 0; i < 4; ++i) {
        const auto& it = find(global_indices.begin(), global_indices.end(), global_to_map[i]);
        ATLAS_ASSERT(it != global_indices.end());
        result[i] = std::distance(global_indices.begin(), it);
      }
      return result;
    };

    std::vector<int> global_indices;
    std::vector<atlas::TriangConnectivityData> tris;
    std::vector<atlas::QuadConnectivityData> quads;
    if (rank == 0) {
      // global indices for points and cells are 1-based
      global_indices = {{1, 2, 3, 4, 5, 7, 17, 18}};
      tris.push_back({0, 2, map_tri_to_local(global_indices, {{18, 4, 7}})});
      quads.push_back({1, 9, map_quad_to_local(global_indices, {{1, 2, 4, 3}})});
      quads.push_back({2, 15, map_quad_to_local(global_indices, {{2, 5, 7, 4}})});
      quads.push_back({3, 22, map_quad_to_local(global_indices, {{3, 4, 18, 17}})});
    } else if (rank == 1) {
      global_indices = {{5, 6, 7, 8, 11, 12, 18, 20}};
      tris.push_back({0, 3, map_tri_to_local(global_indices, {{20, 8, 12}})});
      quads.push_back({1, 10, map_quad_to_local(global_indices, {{5, 6, 8, 7}})});
      quads.push_back({2, 16, map_quad_to_local(global_indices, {{6, 11, 12, 8}})});
      quads.push_back({3, 19, map_quad_to_local(global_indices, {{7, 8, 20, 18}})});
    } else if (rank == 2) {
      global_indices = {{9, 10, 11, 12, 15, 16, 22, 24}};
      tris.push_back({0, 8, map_tri_to_local(global_indices, {{22, 15, 9}})});
      quads.push_back({1, 11, map_quad_to_local(global_indices, {{11, 9, 10, 12}})});
      quads.push_back({2, 17, map_quad_to_local(global_indices, {{9, 15, 16, 10}})});
      quads.push_back({3, 24, map_quad_to_local(global_indices, {{24, 22, 9, 11}})});
    } else if (rank == 3) {
      global_indices = {{1, 3, 13, 14, 15, 16, 21, 22}};
      tris.push_back({0, 5, map_tri_to_local(global_indices, {{21, 1, 13}})});
      quads.push_back({1, 12, map_quad_to_local(global_indices, {{15, 13, 14, 16}})});
      quads.push_back({2, 18, map_quad_to_local(global_indices, {{13, 1, 3, 14}})});
      quads.push_back({3, 25, map_quad_to_local(global_indices, {{22, 21, 13, 15}})});
    } else if (rank == 4) {
      global_indices = {{3, 10, 12, 14, 16, 17, 18, 19, 20}};
      tris.push_back({0, 1, map_tri_to_local(global_indices, {{17, 14, 3}})});
      tris.push_back({1, 4, map_tri_to_local(global_indices, {{19, 10, 16}})});
      quads.push_back({2, 13, map_quad_to_local(global_indices, {{17, 18, 20, 19}})});
      quads.push_back({3, 20, map_quad_to_local(global_indices, {{12, 10, 19, 20}})});
      quads.push_back({4, 21, map_quad_to_local(global_indices, {{16, 14, 17, 19}})});
    } else {  // rank == 5
      global_indices = {{1, 2, 5, 6, 11, 21, 22, 23, 24}};
      tris.push_back({0, 6, map_tri_to_local(global_indices, {{23, 5, 2}})});
      tris.push_back({1, 7, map_tri_to_local(global_indices, {{24, 11, 6}})});
      quads.push_back({2, 14, map_quad_to_local(global_indices, {{21, 22, 24, 23}})});
      quads.push_back({3, 23, map_quad_to_local(global_indices, {{23, 24, 6, 5}})});
      quads.push_back({4, 26, map_quad_to_local(global_indices, {{21, 23, 2, 1}})});
    }

    // compute lon, lat, ghosts, and partitions from global indices
    std::vector<double> lons;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(lons),
        [&global_lons](const int index) { return global_lons[index - 1]; });

    std::vector<double> lats;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(lats),
        [&global_lats](const int index) { return global_lats[index - 1]; });

    std::cout << "FH DEBUG -- max lons = " << *std::max_element(lons.begin(), lons.end()) << std::endl;
    std::cout << "FH DEBUG -- max lats = " << *std::max_element(lats.begin(), lats.end()) << std::endl;
    std::cout << "FH DEBUG -- min lons = " << *std::min_element(lons.begin(), lons.end()) << std::endl;
    std::cout << "FH DEBUG -- min lats = " << *std::min_element(lats.begin(), lats.end()) << std::endl;

    std::vector<int> ghosts;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(ghosts),
        [&global_partitions, &rank](const int index) {
          return static_cast<int>(global_partitions[index - 1] != rank);
        });

    std::vector<int> partitions;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(partitions),
        [&global_partitions](const int index) { return global_partitions[index - 1]; });

    Mesh mesh = build_mesh_from_coordinates_and_connectivities(
        lons,
        lats,
        ghosts,
        global_indices,
        partitions,
        tris,
        quads
        );

    Gmsh gmsh("fh_debug_out.msh", util::Config("coordinates", "xyz"));
    gmsh.write(mesh);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
  return atlas::test::run(argc, argv);
}
