/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <numeric>

#include "atlas/mesh/BuildMeshFromConnectivities.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/output/Gmsh.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::output;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_tiny_mesh") {

    // small regional grid whose cell-centers are connected as:
    //
    //   0 - 4 ----- 5
    //   |     \   / |  <-- cells 3,1,2 respectively
    //   1 ----- 2 - 3
    //
    std::vector<double> lons{{0.0, 0.0, 10.0, 15.0, 5.0, 15.0}};
    std::vector<double> lats{{5.0, 0.0, 0.0, 0.0, 5.0, 5.0}};

    std::vector<int> ghosts(6, 0);  // all points owned
    std::vector<int> global_indices(6);
    std::iota(global_indices.begin(), global_indices.end(), 1);  // 1-based numbering
    std::vector<int> partitions(6, 0);  // all points on proc 0

    // triangles
    std::vector<atlas::TriangConnectivityData> tris = {{
      // local cell index is 0-based; global cell index is 1-based
      {0, 1, {{2, 5, 4}}},
      {1, 2, {{2, 3, 5}}}
    }};
    // quads
    std::vector<atlas::QuadConnectivityData> quads = {{
      {2, 3, {{0, 1, 2, 4}}}
    }};

    Mesh mesh = build_mesh_from_coordinates_and_connectivities(
        lons,
        lats,
        ghosts,
        global_indices,
        partitions,
        tris,
        quads
        );

    //Grid grid = mesh.grid();
    //std::cout << grid.spec() << std::endl;  // doesn't work yet because mesh.grid is not set

    //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

CASE("test_cs_c2_mesh_serial") {

    // coordinates of C2 lfric cubed-sphere grid: grid("CS-LFR-2");
    std::vector<double> lons = {{
      337.5, 22.5, 337.5, 22.5,    // +x
      67.5, 112.5, 67.5, 112.5,    // +y
      202.5, 202.5, 157.5, 157.5,  // -x
      292.5, 292.5, 247.5, 247.5,  // -y
      315, 45, 225, 135,           // +z
      315, 225, 45, 135}};         // -z

    std::vector<double> lats = {{
      -20.941, -20.941, 20.941, 20.941,
      -20.941, -20.941, 20.941, 20.941,
      -20.941, 20.941, -20.941, 20.941,
      -20.941, 20.941, -20.941, 20.941,
      59.6388, 59.6388, 59.6388, 59.6388,
      -59.6388, -59.6388, -59.6388, -59.6388}};

    std::vector<int> ghosts(24, 0);
    std::vector<int> global_indices(24);
    std::iota(global_indices.begin(), global_indices.end(), 1);
    std::vector<int> partitions(24, 0);

    // triangles
    std::vector<atlas::TriangConnectivityData> tris = {{
      // corners
      {0, 1, {{16, 13, 2}}},
      {1, 2, {{17, 3, 6}}},
      {2, 3, {{19, 7, 11}}},
      {3, 4, {{18, 9, 15}}},
      {4, 5, {{20, 0, 12}}},
      {5, 6, {{22, 4, 1}}},
      {6, 7, {{23, 10, 5}}},
      {7, 8, {{21, 14, 8}}}
    }};
    // quads
    std::vector<atlas::QuadConnectivityData> quads = {{
      // faces
      {8, 9, {{0, 1, 3, 2}}},
      {9, 10, {{4, 5, 7, 6}}},
      {10, 11, {{10, 8, 9, 11}}},
      {11, 12, {{14, 12, 13, 15}}},
      {12, 13, {{16, 17, 19, 18}}},
      {13, 14, {{20, 21, 23, 22}}},
      // edges between faces
      {14, 15, {{1, 4, 6, 3}}},
      {15, 16, {{5, 10, 11, 7}}},
      {16, 17, {{8, 14, 15, 9}}},
      {17, 18, {{12, 0, 2, 13}}},
      {18, 19, {{6, 7, 19, 17}}},
      {19, 20, {{11, 9, 18, 19}}},
      {20, 21, {{15, 13, 16, 18}}},
      {21, 22, {{2, 3, 17, 16}}},
      {22, 23, {{22, 23, 5, 4}}},
      {23, 24, {{23, 21, 8, 10}}},
      {24, 25, {{21, 20, 12, 14}}},
      {25, 26, {{20, 22, 1, 0}}}
    }};

    Mesh mesh = build_mesh_from_coordinates_and_connectivities(
        lons,
        lats,
        ghosts,
        global_indices,
        partitions,
        tris,
        quads
        );

    //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
