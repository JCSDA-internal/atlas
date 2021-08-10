/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "eckit/utils/Hash.h"

#include "atlas/field/Field.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Tiles.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/detail/CubedSphereMeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

namespace atlas {
namespace meshgenerator {

// -----------------------------------------------------------------------------

CubedSphereMeshGenerator::CubedSphereMeshGenerator(const eckit::Parametrisation& p) {

  configure_defaults();

  // Get number of partitions.
  size_t nb_parts;
  if (p.get("nb_parts", nb_parts)) options.set("nb_parts", nb_parts);

  // Get this partition.
  size_t part;
  if (p.get("part", part)) options.set("part", part);

  // Get partitioner.
  std::string partitioner;
  if (p.get("partitioner", partitioner)) options.set("partitioner", partitioner);

}

// -----------------------------------------------------------------------------


void CubedSphereMeshGenerator::configure_defaults() {
  // This option sets number of parts the mesh will be split in.
  options.set( "nb_parts", mpi::size() );

  // This option sets the part that will be generated.
  options.set("part", mpi::rank());

  // This options sets the default partitioner.
  options.set<std::string>("partitioner", "cubed_sphere");
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate(const Grid& grid, Mesh& mesh) const {

  // Get partitioner type and number of partitions from config.
  const auto nParts = static_cast<idx_t>(options.get<size_t>("nb_parts"));
  const auto partType = options.get<std::string>("partitioner");

  auto partConfig = util::Config{};
  partConfig.set("type", partType);
  partConfig.set("partitions", nParts);

  // Use lonlat instead of xy for equal_regions partitioner.
  if (partType == "equal_regions") partConfig.set("coordinates", "lonlat");

  // Set distribution.
  const auto partitioner = grid::Partitioner(partConfig);
  const auto distribution = grid::Distribution(grid, partitioner);

  generate(grid, distribution, mesh);

}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate(const Grid& grid, const grid::Distribution& distribution, Mesh& mesh) const {

  // Check for correct grid and need for mesh
  ATLAS_ASSERT(!mesh.generated());
  if (!CubedSphereGrid(grid)) {
    throw_Exception("CubedSphereMeshGenerator can only work "
    "with a cubedsphere grid.", Here());
  }

  // Check for correct stagger
  const auto gridName = grid.name();
  const auto gridStagger = gridName.substr(gridName.rfind("-") - 1, 1);

  if (gridStagger != "C") {
    throw_Exception("CubedSphereMeshGenerator can only work with a"
    "cell-centroid grid. Try NodalCubedSphereMeshGenerator instead.");
  }

  // Cast grid to cubed sphere grid.
  const auto csGrid = CubedSphereGrid(grid);

  // Clone some grid properties.
  setGrid(mesh, csGrid, distribution);

  generate_mesh(csGrid, distribution, mesh);
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate_mesh(const CubedSphereGrid& csGrid,
  const grid::Distribution& distribution, Mesh& mesh) const {

  ATLAS_TRACE("CubedSphereMeshGenerator::generate");

  // ---------------------------------------------------------------------------
  // CUBED SPHERE MESH GENERATOR
  // ---------------------------------------------------------------------------
  //
  // Mesh generator creates a mesh by placing cells with centroid positions at
  // grid xy coordinates. Nodes are then placed around cells by the generator.
  //
  // The node-ownership of cells is determined by the following rules:
  //    * node (i > 0, j > 0) is owned by cell (i - 1, j - 1)
  //    * node (i = 0, j > 0) is owned by cell (i    , j - 1)
  //    * node (i > 0, j = 0) is owned by cell (i - 1, j    )
  //    * node (i = 0, j = 0) is owned by cell (i    , j    )
  //
  // In addtion, there are two types of ghost node that need to be added to the
  // mesh:
  //    * Type-A ghost nodes. These are added to an edge between two tiles. One
  //      tile owns the nodes, the other tile has ghost nodes in the same place.
  //      The ownership rules are determined the CubedSphereTiles class. These
  //      nodes are defined globally and have their own unique global index.
  //      Their placement and local indexing is unaffected by the partitioner.
  //    * Type-B ghost nodes. These are added when a tile is partitioned. One
  //      side of the partation boundary has cells which own the nodes
  //      (see ownership rules above). The other side has ghost nodes. The
  //      placement and indexing of these ghost nodes may vary with different
  //      partitioners. The global index of these nodes are taken from their
  //      owned counterparts on the other side of the partition boundary.
  //
  // Global indices of cells are the same as the grid point indices. Global
  // indices of nodes first count the owned nodes, then the type-A ghost points.
  //
  // Local indices of cells follows the order of the local grid points. Local
  // indces of the nodes first count the owned nodes, followed by type-A ghost
  // nodes, followed by type-B ghost nodes.
  //
  // There are several stages to the mesh generator:
  //    1. Preamble.
  //    2. Define global cell distribution.
  //    3. Define global owned node and type-A ghost node distribution.
  //    4. Locate and count local type-B ghost nodes.
  //    5. Assign nodes to mesh.
  //    6. Assign cells to mesh.
  //    7. Finalise.

  // ---------------------------------------------------------------------------
  // 1. PREAMBLE
  // ---------------------------------------------------------------------------

  using Topology = atlas::mesh::Nodes::Topology;

  using namespace detail::cubedsphere;

  // Get dimensions of grid
  const auto N      = csGrid.N();
  const auto nTiles = csGrid.GetNTiles();

  const auto nNodesUnique = nTiles * N * N + 2;
  const auto nNodesAll    = nTiles * (N + 1) * (N + 1);
  const auto nCells       = nTiles * N * N;

  Log::debug() << "Number of cells per tile edge = "
    << std::to_string(N) << std::endl;

  // Define bad index values.
  constexpr auto undefinedIdx = -1;
  constexpr auto undefinedGlobalIdx = -1;

  // Projection and tiles
  const auto* const csProjection =
    dynamic_cast<const projection::detail::CubedSphereProjectionBase*>(
    csGrid.projection().get());
  ATLAS_ASSERT(csProjection);
  const auto csTiles = csProjection->getCubedSphereTiles();

  // Get partition information.
  const auto nParts =   mpi::comm().size();
  const auto thisPart = mpi::comm().rank();

  // Define an index counter.
  const auto idxSum = [](const std::vector<idx_t>& idxCounts) -> idx_t {
    return std::accumulate(idxCounts.begin(), idxCounts.end(), 0);};

  // Helper functions to get node and cell idx from (t, j, i).
  const auto getNodeIdx = [&](idx_t t, idx_t j, idx_t i){
    // Adjust bounds.
    t = std::max(std::min(t, nTiles - 1), 0);
    j = std::max(std::min(j, N), 0);
    i = std::max(std::min(i, N), 0);
    return idx2st(t * (N + 1) * (N + 1) + j * (N + 1) + i);
  };

  const auto getCellIdx = [&](idx_t t, idx_t j, idx_t i){
    // Adjust bounds.
    t = std::max(std::min(t, nTiles - 1), 0);
    j = std::max(std::min(j, N - 1), 0);
    i = std::max(std::min(i, N - 1), 0);
    return idx2st(t * N * N + j * N + i);
  };

  // ---------------------------------------------------------------------------
  // 2. GLOBAL CELL DISTRIBUTION
  // ---------------------------------------------------------------------------

  // Define cell record.
  // This currently keeps track of more information than we probably need.
  // This may be reduced once we've implemented halos.
  struct CellRecord {
    gidx_t   globalIdx{undefinedGlobalIdx}; // Global ID.
    idx_t   localIdx{undefinedIdx};         // Local ID.
    idx_t   part{undefinedIdx};             // Partition.
  };

  // Define ij bounding box for each face (this partition only).
  struct BoundingBox {
    idx_t iBegin{std::numeric_limits<idx_t>::max()};
    idx_t iEnd{std::numeric_limits<idx_t>::min()};
    idx_t jBegin{std::numeric_limits<idx_t>::max()};
    idx_t jEnd{std::numeric_limits<idx_t>::min()};
  };

  // Make list of all cells.
  auto globalCells = std::vector<CellRecord>(idx2st(nCells));

  // Initialise bounding box.
  auto cellBounds = std::vector<BoundingBox>(idx2st(nTiles));

  // Set xy and tji grid iterators.
  auto tjiIt = csGrid.tij().begin();

  // Set counters for cell local indices.
  auto cellIdxCount = std::vector<idx_t>(nParts, 0);

  for (gidx_t gridIdx = 0; gridIdx < csGrid.size(); ++gridIdx) {

    // Grid (t, j, i) order does not have to match mesh (t, j, i), although
    // in practice it probably does.

    // Get cell index.
    const auto t = (*tjiIt).t();
    const auto j = (*tjiIt).j();
    const auto i = (*tjiIt).i();
    const auto cellIdx = getCellIdx(t, j, i);
    auto& cell = globalCells[cellIdx];

    // Set global index.
    cell.globalIdx = gridIdx + 1;

    // Set partition and remote index.
    cell.part = distribution.partition(gridIdx);

    // Set local index.
    cell.localIdx = cellIdxCount[idx2st(cell.part)]++;

    if (cell.part == st2idx(thisPart)) {

      // Keep track of local (t, j, i) bounds.
      auto& bounds = cellBounds[idx2st(t)];
      bounds.iBegin = std::min(bounds.iBegin, i    );
      bounds.iEnd   = std::max(bounds.iEnd  , i + 1);
      bounds.jBegin = std::min(bounds.jBegin, j    );
      bounds.jEnd   = std::max(bounds.jEnd  , j + 1);

    }
    // Increment iterators.
    ++tjiIt;
  }

  ATLAS_ASSERT(idxSum(cellIdxCount) == nCells);

  // ---------------------------------------------------------------------------
  // 3. GLOBAL OWNED NODE AND TYPE-A GHOST NODE DISTRIBUTION
  // ---------------------------------------------------------------------------

  const auto jacs = CubedSphereJacobian(csGrid);

  // Define global node record.
  // This currently keeps track of more information than we probably need.
  // This may be reduced once we've implemented halos.
  struct GlobalNode {
    idx_t             globalIdx{undefinedIdx};  // Global ID,
    idx_t             part{undefinedIdx};       // Partition node is on.
    idx_t             localIdx{undefinedIdx};   // Local ID of owned node.
    const GlobalNode* remoteNode{};             // Pointer to remote node.
  };

  // Make list of all nodes.
  auto globalNodes = std::vector<GlobalNode>(idx2st(nNodesAll));

  // Set counters for local node indices.
  auto nodeLocalIdxCount = std::vector<idx_t>(nParts, 0);

  // Set counter for global indices.
  idx_t nodeGlobalOwnedIdxCount = 1;
  idx_t nodeGlobalGhostIdxCount = nNodesUnique + 1;

  // Make a list of type-A ghost nodes and process later.
  auto typeAGhostNodes = std::vector<GlobalNode*>{};

  for (idx_t t = 0; t < nTiles; ++t) {
    for (idx_t j = 0; j < N + 1; ++j) {
      for (idx_t i = 0; i < N + 1; ++i) {

        // Get this node.
        auto& node = globalNodes[getNodeIdx(t, j, i)];

        // Get owning node
        const auto ownerTij = jacs.tijLocalToGlobal(PointTIJ(t, i, j));
        const auto ownerIdx =
          getNodeIdx(ownerTij.t(), ownerTij.jRound(), ownerTij.iRound());
        const auto& ownerNode = globalNodes[ownerIdx];

        // Get owning cell.
        const auto& cell = globalCells[getCellIdx(t, j - 1, i - 1)];

        // Determine node ownership.
        if (&node == &ownerNode) {

          // Set partition.
          node.part = cell.part;

          // Owner node. Set local index.
          node.localIdx = nodeLocalIdxCount[idx2st(node.part)]++;

          // Set global index.
          node.globalIdx = nodeGlobalOwnedIdxCount++;

        } else {

          // Ghost node on tile edge. Point to owner.
          node.remoteNode = &ownerNode;

          // Set global index.
          node.globalIdx = nodeGlobalGhostIdxCount++;

        }

      }
    }
  }

  ATLAS_ASSERT(nodeGlobalOwnedIdxCount == nNodesUnique + 1);
  ATLAS_ASSERT(nodeGlobalGhostIdxCount == nNodesAll + 1);
  ATLAS_ASSERT(idxSum(nodeLocalIdxCount) == nNodesUnique);

  enum class NodeType {
    OWNED,
    GHOST_A,
    GHOST_B,
  };

  // Define local node record.
  struct LocalNode {

    NodeType      type{};
    idx_t         globalIdx{};
    idx_t         remoteIdx{undefinedGlobalIdx};
    idx_t         part{undefinedGlobalIdx};
    idx_t*        localIdx;
    PointXY       xyLocal{};
    PointXY       xyGlobal{};
  };

  auto localNodes = std::vector<LocalNode>{};

  // Loop over all possible local nodes.
  for (idx_t t = 0; t < nTiles; ++t) {

    // Limit range to bounds recorded earlier.
    const auto bounds = cellBounds[idx2st(t)];

    for (idx_t j = bounds.jBegin; j < bounds.jEnd + 1; ++j) {
      for (idx_t i = bounds.iBegin; i < bounds.iEnd + 1; ++i) {

        // Get node.
        auto& globalNode = globalNodes[getNodeIdx(t, j, i)];
        auto localNode = LocalNode{};

        // Get four neighbouring cells of node. (cell0 is owner of node).
        const auto& cell0 = globalCells[getCellIdx(t, j - 1, i - 1)];
        const auto& cell1 = globalCells[getCellIdx(t, j - 1, i    )];
        const auto& cell2 = globalCells[getCellIdx(t, j    , i    )];
        const auto& cell3 = globalCells[getCellIdx(t, j    , i - 1)];

        if (cell0.part == st2idx(thisPart)) {

          // Node is either owned, or a ghost on a tile edge.
          if (!globalNode.remoteNode) {

            // Node is owned.
            localNode.type = NodeType::OWNED;
            localNode.globalIdx = globalNode.globalIdx;
            localNode.remoteIdx = globalNode.localIdx;
            localNode.part      = globalNode.part;
            localNode.localIdx  = &globalNode.localIdx;

            const auto txyLocal = jacs.tijToTxy(PointTIJ(t, i, j));

            localNode.xyLocal   = txyLocal.xy();
            localNode.xyGlobal  = localNode.xyLocal;

            localNodes.push_back(localNode);

          } else {

            // Node is a ghost on a tile edge.
            localNode.type = NodeType::GHOST_A;
            localNode.globalIdx = globalNode.globalIdx;
            localNode.remoteIdx = globalNode.remoteNode->localIdx;
            localNode.part      = globalNode.remoteNode->part;
            localNode.localIdx  = &globalNode.localIdx;

            const auto txyLocal = jacs.tijToTxy(PointTIJ(t, i, j));

            localNode.xyLocal   = txyLocal.xy();
            localNode.xyGlobal  = jacs.txyLocalToGlobal(txyLocal).xy();

            localNodes.push_back(localNode);

          }

        } else if (cell1.part == st2idx(thisPart)
                || cell2.part == st2idx(thisPart)
                || cell3.part == st2idx(thisPart)) {

          // Node is a ghost on a partition edge.
          localNode.type = NodeType::GHOST_B;
          localNode.globalIdx = globalNode.globalIdx;
          localNode.remoteIdx = globalNode.localIdx;
          localNode.part      = globalNode.part;
          localNode.localIdx  = &globalNode.localIdx;

          const auto txyLocal = jacs.tijToTxy(PointTIJ(t, i, j));

          localNode.xyLocal   = txyLocal.xy();
          localNode.xyGlobal  = jacs.txyLocalToGlobal(txyLocal).xy();

          localNodes.push_back(localNode);

        }

      }
    }
  }

  // Ghost nodes are mixed in with owned nodes. Sort by node type to correct.
  // Note: value of localIdx for owned nodes should match element number of
  // localNodes vector. localIdx value will be overwritten for non-owned nodes.
  std::stable_sort(localNodes.begin(), localNodes.end(),
    [](const LocalNode& nodeA, const LocalNode& nodeB){
      return nodeA.type < nodeB.type;
    });

  // Check/Set local indices.
  idx_t nodeLocalIdx = 0;
  for (auto& node : localNodes) {

    if (node.type == NodeType::OWNED) {
      ATLAS_ASSERT(*node.localIdx == nodeLocalIdx);
    } else {
      *node.localIdx = nodeLocalIdx;
    }
    ++nodeLocalIdx;
  }

  // ---------------------------------------------------------------------------
  // 5. ASSIGN NODES TO MESH
  // ---------------------------------------------------------------------------

  // Resize nodes.
  mesh.nodes().resize(st2idx(localNodes.size()));

  // Get field views
  auto nodesGlobalIdx = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
  auto nodesRemoteIdx = array::make_indexview<idx_t, 1>(mesh.nodes().remote_index());
  auto nodesXy        = array::make_view<double, 2>(mesh.nodes().xy());
  auto nodesLonLat    = array::make_view<double, 2>(mesh.nodes().lonlat());
  auto nodesPart      = array::make_view<int, 1>(mesh.nodes().partition());
  auto nodesGhost     = array::make_view<int, 1>(mesh.nodes().ghost());
  auto nodesFlags     = array::make_view<int, 1>(mesh.nodes().flags());

  // Set fields.
  nodeLocalIdx = 0;
  for (const auto node : localNodes) {

    // Set global index.
    nodesGlobalIdx(nodeLocalIdx) = node.globalIdx;

    // Set node remote index.
    nodesRemoteIdx(nodeLocalIdx) = node.remoteIdx;

    // Set node partition.
    nodesPart(nodeLocalIdx) = node.part;

    // Set xy.
    nodesXy(nodeLocalIdx, XX) = node.xyLocal.x();
    nodesXy(nodeLocalIdx, YY) = node.xyLocal.y();

    // Set lon-lat.
    const auto lonLat = csProjection->lonlat(node.xyGlobal);
    nodesLonLat(nodeLocalIdx, LON) = lonLat.lon();
    nodesLonLat(nodeLocalIdx, LAT) = lonLat.lat();

    // Set flags.
    Topology::reset(nodesFlags(nodeLocalIdx));
    switch(node.type) {

      case NodeType::OWNED : {
        nodesGhost(nodeLocalIdx) = 0;
        break;
      }
      case NodeType::GHOST_A : {
        nodesGhost(nodeLocalIdx) = 1;
        Topology::set(nodesFlags(nodeLocalIdx), Topology::GHOST);
        break;
      }
      case NodeType::GHOST_B : {
        nodesGhost(nodeLocalIdx) = 1;
        Topology::set(nodesFlags(nodeLocalIdx), Topology::GHOST);
        break;
      }
    }

    ++nodeLocalIdx;
  }

  // ---------------------------------------------------------------------------
  // 6. ASSIGN CELLS TO MESH
  // ---------------------------------------------------------------------------

  // Resize cells.
  mesh.cells().add(new mesh::temporary::Quadrilateral(), cellIdxCount[thisPart]);

  // Set field views.
  auto cellsGlobalIdx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
  auto cellsPart      = array::make_view<int, 1>(mesh.cells().partition());

  // Set local cells.
  auto& nodeConnectivity = mesh.cells().node_connectivity();
  const auto idx0 = mesh.cells().elements(0).begin();

  for (idx_t t = 0; t < nTiles; ++t) {

    // Use bounds from before.
    const auto bounds = cellBounds[idx2st(t)];

    for (idx_t j = bounds.jBegin; j < bounds.jEnd; ++j) {
      for (idx_t i = bounds.iBegin; i < bounds.iEnd; ++i) {

        // Get cell.
        const auto& cell = globalCells[getCellIdx(t, j, i)];

        // Only add cells on this partition.
        if (cell.part == st2idx(thisPart)) {

          // Get local index.
          const auto localIdx = cell.localIdx + idx0;

          // Set quadrilateral.
          const auto& node0 = globalNodes[getNodeIdx(t, j    , i    )];
          const auto& node1 = globalNodes[getNodeIdx(t, j    , i + 1)];
          const auto& node2 = globalNodes[getNodeIdx(t, j + 1, i + 1)];
          const auto& node3 = globalNodes[getNodeIdx(t, j + 1, i    )];

          const auto quadNodeIdx = std::array<idx_t, 4> {
            node0.localIdx, node1.localIdx, node2.localIdx, node3.localIdx};

          // Set connectivity.
          nodeConnectivity.set(localIdx, quadNodeIdx.data());

          // Set global index.
          cellsGlobalIdx(localIdx) = cell.globalIdx;

          // Set partition.
          cellsPart(localIdx) = cell.part;

        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // 7. FINALISE
  // ---------------------------------------------------------------------------

  mesh.nodes().global_index().metadata().set("human_readable", true);
  mesh.nodes().global_index().metadata().set("min", 1);
  mesh.nodes().global_index().metadata().set("max", nNodesAll);
  mesh.nodes().metadata().set("parallel", true);

  mesh.cells().global_index().metadata().set("human_readable", true);
  mesh.cells().global_index().metadata().set("min", 1);
  mesh.cells().global_index().metadata().set("max", nCells);
  mesh.cells().metadata().set("parallel", true);

  return;
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::hash(eckit::Hash& h) const {
h.add("CubedSphereMeshGenerator");
options.hash(h);
}

// -----------------------------------------------------------------------------

namespace {
static MeshGeneratorBuilder<CubedSphereMeshGenerator> CubedSphereMeshGenerator(
CubedSphereMeshGenerator::static_type());
}

// -----------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
