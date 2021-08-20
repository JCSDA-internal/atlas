/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Tiles.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/output/Gmsh.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/CubedSpherePartitioner.h"
#include "atlas/option.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

CASE("cubedsphere_mesh_jacobian_test") {

  using namespace meshgenerator::detail::cubedsphere;

  // Set N.
  const idx_t N = 2;

  // Set grid.
  const auto grid = atlas::Grid("CS-LFR-C-" + std::to_string(N));

  // Set Jacobian
  const auto jacobian = NeighbourJacobian(CubedSphereGrid(grid));


  // Play around with some grids.
  for (idx_t t = 0; t < 6; ++t) {

    std::cout << t << std::endl;
    for (idx_t j = - 1; j < N + 2; ++j) {
      for (idx_t i = - 1; i < N + 2; ++i) {

        const auto ij = PointIJ(i, j);

        if (!jacobian.ijCross(ij)) continue;
        if (jacobian.ijInterior(ij)) continue;

        const auto xyLocal = jacobian.xy(ij, t);
        const auto xytGlobal = jacobian.xyLocalToGlobal(xyLocal, t);
        const auto ijGlobal = jacobian.ijLocalToGlobal(ij, t);

        std::cout << ij << " ";
        std::cout << ijGlobal.first << " " << ijGlobal.second << "   ";


      }
    std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
  }
}

CASE("cubedsphere_mesh_test") {


  // Set grid.
  const auto grid = atlas::Grid("CS-LFR-C-24");

  // Set mesh config.
  auto meshConfig = util::Config("partitioner", "equal_regions");

  // Set mesh generators.
  const auto csMeshgen = atlas::MeshGenerator("cubedsphere"); // defaults to cubed sphere partitioner.
  const auto erMeshgen = atlas::MeshGenerator("cubedsphere", meshConfig); // Equal regions partitioner.

  const auto csMesh = csMeshgen.generate(grid);
  const auto erMesh = erMeshgen.generate(grid);


  const auto erFuncSpace = atlas::functionspace::NodeColumns(erMesh);
  auto erFieldSet = atlas::FieldSet{};
  erFieldSet.add(erMesh.nodes().xy());
  erFieldSet.add(erMesh.nodes().lonlat());
  erFieldSet.add(erMesh.nodes().ghost());
  erFieldSet.add(erMesh.nodes().halo());
  erFieldSet.add(erMesh.nodes().remote_index());
  erFieldSet.add(erMesh.nodes().partition());
  erFieldSet.add(erMesh.nodes().global_index());
  erFieldSet.add(erMesh.nodes().field("ijt"));


  const auto csFuncSpace = atlas::functionspace::NodeColumns(csMesh);
  auto csFieldSet = atlas::FieldSet{};
  csFieldSet.add(csMesh.nodes().xy());
  csFieldSet.add(csMesh.nodes().lonlat());
  csFieldSet.add(csMesh.nodes().ghost());
  csFieldSet.add(csMesh.nodes().halo());
  csFieldSet.add(csMesh.nodes().remote_index());
  csFieldSet.add(csMesh.nodes().partition());
  csFieldSet.add(csMesh.nodes().global_index());
  csFieldSet.add(csMesh.nodes().field("ijt"));


  // Set gmsh config.
  auto gmshConfigXy = atlas::util::Config("coordinates", "xy");
  auto gmshConfigXyz = atlas::util::Config("coordinates", "xyz");
  auto gmshConfigLonLat = atlas::util::Config("coordinates", "lonlat");

  gmshConfigXy.set("ghost", true);

  gmshConfigXyz.set("ghost", true);

  gmshConfigLonLat.set("ghost", true);

  // Set gmsh objects.
  auto gmshXy = atlas::output::Gmsh("cs_xy_mesh.msh", gmshConfigXy);
  auto gmshXyz = atlas::output::Gmsh("cs_xyz_mesh.msh", gmshConfigXyz);

  // Write gmsh.
  gmshXy.write(csMesh);
  gmshXy.write(csFieldSet, csFuncSpace);
  gmshXyz.write(csMesh);
  gmshXyz.write(csFieldSet, csFuncSpace);


  // Set gmsh objects.
  gmshXy = atlas::output::Gmsh("er_xy_mesh.msh", gmshConfigXy);
  gmshXyz = atlas::output::Gmsh("er_xyz_mesh.msh", gmshConfigXyz);

  // Write gmsh.
  gmshXy.write(erMesh);
  gmshXy.write(erFieldSet, erFuncSpace);
  gmshXyz.write(erMesh);
  gmshXyz.write(erFieldSet, erFuncSpace);


}



}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
  return atlas::test::run( argc, argv );
}
