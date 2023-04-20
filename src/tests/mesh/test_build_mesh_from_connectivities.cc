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
    //   |     \   / |  <-- cells 2,0,1 respectively
    //   1 ----- 2 - 3
    //
    std::vector<double> lons{{0.0, 0.0, 10.0, 15.0, 5.0, 15.0}};
    std::vector<double> lats{{5.0, 0.0, 0.0, 0.0, 5.0, 5.0}};

    std::vector<int> ghosts(6, 0);  // all points owned
    std::vector<int> global_indices({0, 1, 2, 3, 4, 5});
    std::vector<int> partitions(6, 0);  // all points on proc 0

    // triangles
    std::vector<atlas::TriangConnectivityData> tris = {{
      {0, 0, {{2, 5, 4}}},  // cell 0
      {1, 1, {{2, 3, 5}}}   // cell 1
    }};
    // quads
    std::vector<atlas::QuadConnectivityData> quads = {{
      {2, 2, {{0, 1, 2, 4}}}  // cell 2
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

    //Gmsh gmsh("fh_debug_out.msh", util::Config("coordinates", "xyz"));
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

    std::vector<int> ghosts(24, 0);  // all points owned
    std::vector<int> global_indices(24);
    std::iota(global_indices.begin(), global_indices.end(), 0);
    std::vector<int> partitions(24, 0);  // all points on proc 0

    // triangles
    std::vector<atlas::TriangConnectivityData> tris = {{
      // corners
      {0, 0, {{16, 13, 2}}},
      {1, 1, {{17, 3, 6}}},
      {2, 2, {{19, 7, 11}}},
      {3, 3, {{18, 9, 15}}},
      {4, 4, {{20, 0, 12}}},
      {5, 5, {{22, 4, 1}}},
      {6, 6, {{23, 10, 5}}},
      {7, 7, {{21, 14, 8}}}
    }};
    // quads
    std::vector<atlas::QuadConnectivityData> quads = {{
      // faces
      {8, 8, {{0, 1, 3, 2}}},
      {9, 9, {{4, 5, 7, 6}}},
      {10, 10, {{10, 8, 9, 11}}},
      {11, 11, {{14, 12, 13, 15}}},
      {12, 12, {{16, 17, 19, 18}}},
      {13, 13, {{20, 21, 23, 22}}},
      // edges between faces
      {14, 14, {{1, 4, 6, 3}}},
      {15, 15, {{5, 10, 11, 7}}},
      {16, 16, {{8, 14, 15, 9}}},
      {17, 17, {{12, 0, 2, 13}}},
      {18, 18, {{6, 7, 19, 17}}},
      {19, 19, {{11, 9, 18, 19}}},
      {20, 20, {{15, 13, 16, 18}}},
      {21, 21, {{2, 3, 17, 16}}},
      {22, 22, {{22, 23, 5, 4}}},
      {23, 23, {{23, 21, 8, 10}}},
      {24, 24, {{21, 20, 12, 14}}},
      {25, 25, {{20, 22, 1, 0}}}
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

    //Gmsh gmsh("fh_debug_out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
