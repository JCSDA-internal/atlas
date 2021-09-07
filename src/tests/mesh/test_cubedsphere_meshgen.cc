/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/CubedSphereCellColumns.h"
#include "atlas/functionspace/CubedSphereNodeColumns.h"
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

  // Set grid an N2 grid with a halo size of 1.
  const auto grid = atlas::Grid("CS-LFR-C-2");

  // Set Jacobian
  const auto jacobian = NeighbourJacobian(CubedSphereGrid(grid));

  // Set vectors of known good outputs.
  const auto xyLocalKgoVec = std::vector<PointXY>{
    {0,-90}, {45,-90}, {90,-90}, {-45,-45}, {135,-45}, {-45,0}, {135,0},
    {-45,45}, {135,45}, {0,90}, {45,90}, {90,90}, {90,-90}, {135,-90},
    {180,-90}, {45,-45}, {225,-45}, {45,0}, {225,0}, {45,45}, {225,45}, {90,90},
    {135,90}, {180,90}, {315,-45}, {315,0}, {315,45}, {270,-90}, {270,90},
    {225,-90}, {225,90}, {180,-90}, {180,90}, {135,-45}, {135,0}, {135,45},
    {405,-45}, {405,0}, {405,45}, {360,-90}, {360,90}, {315,-90}, {315,90},
    {270,-90}, {270,90}, {225,-45}, {225,0}, {225,45}, {0,0}, {45,0}, {90,0},
    {-45,45}, {135,45}, {-45,90}, {135,90}, {-45,135}, {135,135}, {0,180},
    {45,180}, {90,180}, {-45,-45}, {-45,-90}, {-45,-135}, {0,0}, {0,-180},
    {45,0}, {45,-180}, {90,0}, {90,-180}, {135,-45}, {135,-90}, {135,-135}};

  const auto xyGlobalKgoVec = std::vector<PointXY>{
    {315,-45}, {45,-90}, {135,-45}, {315,-45}, {135,-45}, {315,0}, {135,0},
    {0,90}, {90,90}, {0,90}, {45,90}, {90,90}, {45,-45}, {45,-90}, {225,-45},
    {45,-45}, {225,-45}, {45,0}, {225,0}, {45,45}, {45,135}, {45,45}, {45,90},
    {45,135}, {315,-45}, {315,0}, {0,90}, {315,-45}, {0,90}, {45,-90}, {45,90},
    {135,-45}, {90,90}, {135,-45}, {135,0}, {90,90}, {45,-45}, {45,0}, {45,45},
    {45,-45}, {45,45}, {45,-90}, {45,90}, {225,-45}, {45,135}, {225,-45},
    {225,0}, {45,135}, {0,0}, {45,0}, {90,0}, {0,0}, {90,0}, {315,0}, {135,0},
    {270,0}, {180,0}, {270,0}, {225,0}, {180,0}, {0,0}, {315,0}, {270,0}, {0,0},
    {270,0}, {45,0}, {225,0}, {90,0}, {180,0}, {90,0}, {135,0}, {180,0}};

  const auto ijGlobalKgoVec = std::vector<PointIJ>{
    {0,1}, {1,1}, {1,0}, {0,1}, {1,0}, {1,1}, {1,1}, {0,1}, {2,1}, {0,1}, {1,1},
    {2,1}, {1,0}, {1,1}, {0,1}, {1,0}, {0,1}, {1,1}, {1,1}, {1,0}, {1,2}, {1,0},
    {1,1}, {1,2}, {0,1}, {1,1}, {0,1}, {0,1}, {0,1}, {1,1}, {1,1}, {1,0}, {2,1},
    {1,0}, {1,1}, {2,1}, {1,0}, {1,1}, {1,0}, {1,0}, {1,0}, {1,1}, {1,1}, {0,1},
    {1,2}, {0,1}, {1,1}, {1,2}, {0,1}, {1,1}, {0,1}, {0,1}, {0,1}, {1,1}, {1,1},
    {1,2}, {1,2}, {1,2}, {1,1}, {1,2}, {0,1}, {1,1}, {1,2}, {0,1}, {1,2}, {1,1},
    {1,1}, {0,1}, {1,2}, {0,1}, {1,1}, {1,2}};

  const auto tKgoVec = std::vector<idx_t>{
    3, 5, 1, 3, 1, 3, 1, 4, 4, 4, 4, 4, 0, 5, 2, 0, 2, 0, 2, 4, 4, 4, 4, 4, 3,
    3, 4, 3, 4, 5, 4, 1, 4, 1, 1, 4, 0, 0, 4, 0, 4, 5, 4, 2, 4, 2, 2, 4, 0, 0,
    1, 0, 1, 3, 1, 3, 2, 3, 2, 2, 0, 3, 3, 0, 3, 0, 2, 1, 2, 1, 1, 2};

  // Set kgo iterators.
  auto xyLocalKgoIt = xyLocalKgoVec.cbegin();
  auto xyGlobalKgoIt = xyGlobalKgoVec.cbegin();
  auto ijGlobalKgoIt = ijGlobalKgoVec.cbegin();
  auto tKgoIt = tKgoVec.cbegin();


  // Play around with some grids.
  for (idx_t t = 0; t < 6; ++t) {
    for (idx_t j = - 1; j < 4; ++j) {
      for (idx_t i = - 1; i < 4; ++i) {

        // Set ij object.
        const auto ij = PointIJ(i, j);

        // Only look at halo values.
        if (!jacobian.ijCross(ij)) continue;
        if (jacobian.ijInterior(ij)) continue;

        // Get known good outputs.
        const auto xyLocalKgo = *xyLocalKgoIt++;
        const auto xyGlobalKgo = *xyGlobalKgoIt++;
        const auto ijGlobalKgo = *ijGlobalKgoIt++;
        const auto tKgo = *tKgoIt++;

        // Test known good output.
        const auto xyLocal = jacobian.xy(ij, t);
        const auto xytGlobal = jacobian.xyLocalToGlobal(xyLocal, t);
        const auto ijtGlobal = jacobian.ijLocalToGlobal(ij, t);

        ATLAS_ASSERT(xyLocal == xyLocalKgo);
        ATLAS_ASSERT(xytGlobal.first == xyGlobalKgo);
        ATLAS_ASSERT(xytGlobal.second == tKgo);
        ATLAS_ASSERT(ijtGlobal.first == ijGlobalKgo);
        ATLAS_ASSERT(ijtGlobal.second == tKgo);

        // Check xy and ij transforms are consistent.
        ATLAS_ASSERT(jacobian.ij((jacobian.xyLocalToGlobal(xyLocal, t))) ==
          jacobian.ijLocalToGlobal(ij, t).first);

        ATLAS_ASSERT(jacobian.xy((jacobian.ijLocalToGlobal(ij, t))) ==
          jacobian.xyLocalToGlobal(xyLocal, t).first);

      }
    }
  }
}


double testFunction(double lon, double lat) {
  return std::sin(3 * lon * M_PI / 180) * std::sin(2 * lat * M_PI / 180);
}


void testHaloExchange(const std::string& gridStr, const std::string& partitionerStr,
                     idx_t halo) {

  // Set grid.
  const auto grid = Grid(gridStr);

  // Set mesh config.
  const auto meshConfig =
    util::Config("partitioner", partitionerStr) |
    util::Config("halo", halo);

  // Set mesh generator.
  const auto meshGen = MeshGenerator("cubedsphere", meshConfig);

  // Set mesh
  const auto mesh = meshGen.generate(grid);

  // Set function space.
  const auto nodeColumns = functionspace::NodeColumns(mesh);
  const auto cellColumns = functionspace::CellColumns(mesh);

  // Make a field set of useful diagnostic quantities.

  auto fieldSet = atlas::FieldSet{};
  fieldSet.add(mesh.nodes().xy());
  fieldSet.add(mesh.nodes().lonlat());
  fieldSet.add(mesh.nodes().ghost());
  fieldSet.add(mesh.nodes().halo());
  fieldSet.add(mesh.nodes().remote_index());
  fieldSet.add(mesh.nodes().partition());
  fieldSet.add(mesh.nodes().global_index());
  fieldSet.add(mesh.nodes().field("ijt"));

  // Set gmsh config.
  const auto gmshConfigXy =
    util::Config("coordinates", "xy") |
    util::Config("ghost", true);

  const auto gmshConfigXyz =
    util::Config("coordinates", "xyz") |
    util::Config("ghost", true);

  // Set gmsh objects.
  const auto fileStr =
   gridStr + "_" + partitionerStr + "_halo" + std::to_string(halo);
  const auto gmshXy =
   output::Gmsh(fileStr + "_xy.msh", gmshConfigXy);
  const auto gmshXyz =
   output::Gmsh(fileStr + "_xyz.msh", gmshConfigXyz);

  // Write outputs.
  gmshXy.write(mesh);
  gmshXy.write(fieldSet, nodeColumns);

  gmshXyz.write(mesh);
  gmshXyz.write(fieldSet, nodeColumns);

  // ---------------------------------------------------------------------------
  // Test node columns halo exchange.
  // ---------------------------------------------------------------------------

  // make a test field.
  auto testField1 = nodeColumns.createField<double>(util::Config("name", "test field (node columns)"));

  // Make some field views.
  auto testView1 = array::make_view<double, 1>(testField1);
  auto lonLatView = array::make_view<double, 2>(nodeColumns.lonlat());
  auto ghostView = array::make_view<idx_t, 1>(nodeColumns.ghost());

  // Set non-ghost values.
  idx_t testFuncCallCount = 0;
  for (idx_t i = 0; i < nodeColumns.size(); ++i) {

    // Stop once we've reached ghost points.
    if (ghostView(i)) break;

    testView1(i) = testFunction(lonLatView(i, LON), lonLatView(i, LAT));
    ++testFuncCallCount;

  }

  // Make sure that some of the field values were ghosts.
  if (halo > 0) {
    ATLAS_ASSERT(testFuncCallCount < nodeColumns.size());
  }

  nodeColumns.haloExchange(testField1);

  // Check all values after halo exchange.
  for (idx_t i = 0; i < nodeColumns.size(); ++i) {

    // Test field and test function should be the same.
    ATLAS_ASSERT(is_approximately_equal(
      testView1(i), testFunction(lonLatView(i, LON), lonLatView(i, LAT))));

  }

  // Write fields.
  gmshXy.write(testField1, nodeColumns);
  gmshXyz.write(testField1, nodeColumns);

  // ---------------------------------------------------------------------------
  // Test cell columns halo exchange.
  // ---------------------------------------------------------------------------

  // make a test field.
  auto testField2 = cellColumns.createField<double>(util::Config("name", "test field (cell columns)"));

  // Make some field views.
  auto testView2 = array::make_view<double, 1>(testField2);
  auto haloView = array::make_view<idx_t, 1>(cellColumns.mesh().cells().halo());
  lonLatView = array::make_view<double, 2>(cellColumns.mesh().cells().field("lonlat"));

  // Set non-halo values.
  testFuncCallCount = 0;
  for (idx_t i = 0; i < cellColumns.size(); ++i) {


    if (haloView(i)) break;

    testView2(i) = testFunction(lonLatView(i, LON), lonLatView(i, LAT));

  }

  // Make sure that some of the field values were ghosts.
  if (halo > 0) {
    ATLAS_ASSERT(testFuncCallCount < cellColumns.size());
  }

  cellColumns.haloExchange(testField2);

  // Check all values after halo exchange.
  for (idx_t i = 0; i < cellColumns.size(); ++i) {

    // Test field and test function should be the same.
    ATLAS_ASSERT(is_approximately_equal(
      testView2(i), testFunction(lonLatView(i, LON), lonLatView(i, LAT))));

  }


}

CASE("cubedsphere_mesh_test") {

  SECTION("halo = 0") {
    testHaloExchange("CS-LFR-C-12", "equal_regions", 0);
    testHaloExchange("CS-LFR-C-12", "cubed_sphere", 0);
  }
  SECTION("halo = 1") {
    testHaloExchange("CS-LFR-C-12", "equal_regions", 1);
    testHaloExchange("CS-LFR-C-12", "cubed_sphere", 1);
  }
  SECTION("halo = 2") {
    testHaloExchange("CS-LFR-C-12", "equal_regions", 2);
    testHaloExchange("CS-LFR-C-12", "cubed_sphere", 2);
  }
  SECTION("halo = 3") {
    testHaloExchange("CS-LFR-C-12", "equal_regions", 3);
    testHaloExchange("CS-LFR-C-12", "cubed_sphere", 3);
  }
}

template<typename FSpace>
void testFunctionSpace(const FSpace& functionspace) {

  // Make field.
  auto field = functionspace.template createField<double>(
    util::Config("name", "test field"));
  auto fieldView = array::make_view<double, 1>(field);

  // Get view of lonlat.
  const auto lonLatView = array::make_view<double, 2>(functionspace.lonlat());

  // Get view of halo/ghosts.
  const auto ghostView = array::make_view<idx_t, 1>(functionspace.ghost());

  // Loop over all non halo elements of test field.
  idx_t testFuncCallCount = 0;
  functionspace.for_each(
    [&](idx_t index, idx_t i, idx_t j, idx_t t) {

      // Make sure index matches ijt.
      ATLAS_ASSERT(index == functionspace.index(i, j, t));

      // Check that indices of "+" stencil are valid.
      const auto badIdx = functionspace.invalid_index();
      ATLAS_ASSERT(functionspace.index(i - 1, j    , t) != badIdx);
      ATLAS_ASSERT(functionspace.index(i + 1, j    , t) != badIdx);
      ATLAS_ASSERT(functionspace.index(i    , j - 1, t) != badIdx);
      ATLAS_ASSERT(functionspace.index(i    , j + 1, t) != badIdx);

      // Make sure we're avoiding halos.
      ATLAS_ASSERT(!ghostView(index));

      // Set field values.
      fieldView(index) = testFunction(lonLatView(index, LON), lonLatView(index, LAT));
      ++testFuncCallCount;

  });

  // Make sure call count is less than functionspace.size() as we skipped halos.
  ATLAS_ASSERT(testFuncCallCount < functionspace.size());

  // Perform halo exchange.
  functionspace.haloExchange(field);

  // Loop over elements including halo
  testFuncCallCount = 0;
  functionspace.for_each(
    [&](idx_t index, idx_t i, idx_t j, idx_t t) {

      // Make sure index matches ijt.
      ATLAS_ASSERT(index == functionspace.index(i, j, t));

      // Set field values.
      ATLAS_ASSERT(is_approximately_equal(
        fieldView(index), testFunction(lonLatView(index, LON), lonLatView(index, LAT))));
      ++testFuncCallCount;

  }, true);

  // Make sure call count is equal to functionspace.size().
  ATLAS_ASSERT(testFuncCallCount == functionspace.size());


}

CASE("cubedsphere_mesh_functionspace") {

  // Set grid.
  const auto grid = Grid("CS-LFR-C-12");

  // Set mesh config.
  const auto meshConfigEqualRegions =
    util::Config("partitioner", "equal_regions") |
    util::Config("halo", 1);
  const auto meshConfigCubedSphere =
    util::Config("partitioner", "cubed_sphere") |
    util::Config("halo", 1);

  // Set mesh generator.
  const auto meshGenEqualRegions = MeshGenerator("cubedsphere", meshConfigEqualRegions);
  const auto meshGenCubedSphere = MeshGenerator("cubedsphere", meshConfigCubedSphere);

  // Set mesh
  const auto meshEqualRegions = meshGenEqualRegions.generate(grid);
  const auto meshCubedSphere = meshGenCubedSphere.generate(grid);

  // Set functionspace.
  const auto equalRegionsCellColumns =
    functionspace::CubedSphereCellColumns(meshEqualRegions);
  const auto cubedSphereCellColumns =
    functionspace::CubedSphereCellColumns(meshCubedSphere);
  const auto equalRegionsNodeColumns =
    functionspace::CubedSphereNodeColumns(meshEqualRegions);
  const auto cubedSphereNodeColumns =
    functionspace::CubedSphereNodeColumns(meshCubedSphere);

  // test functionspaces.
  SECTION("CellColumns: equal_regions") {
    testFunctionSpace(equalRegionsCellColumns);
  }
  SECTION("CellColumns: cubed_sphere") {
    testFunctionSpace(cubedSphereCellColumns);
  }
  SECTION("NodeColumns: equal_regions") {
    testFunctionSpace(equalRegionsNodeColumns);
  }
  SECTION("NodeColumns: cubed_sphere") {
    testFunctionSpace(cubedSphereNodeColumns);
  }


}



}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
  return atlas::test::run( argc, argv );
}

