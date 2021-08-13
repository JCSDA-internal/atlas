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

  // Get halo size.
  idx_t halo;
  if (p.get("halo", halo)) options.set("halo", halo);

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

  // This options sets the number of halo elements around each part.
  options.set("halo", idx_t{1});

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

  // Get size of halo.
  idx_t halo = 0;
  options.get("halo", halo);

  // Set bounds for (i, j) based loops.
  idx_t ijBegin = -halo;
  idx_t ijEnd = N + halo;

  const auto nNodesUnique = nTiles * N * N + 2;
  const auto nNodesGrid   = nTiles * (N + 2 * halo + 1) * (N + 2 * halo + 1);
  const auto nNodesTotal  = nNodesGrid - nTiles * 4 * halo * halo;
  const auto nCellsUnique = nTiles * N * N;
  const auto nCellsGrid   = nTiles * (N + 2 * halo) * (N + 2 * halo);
  const auto nCellsTotal  = nCellsGrid - nTiles * 4 * halo * halo;

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

  // Helper functions to get node idx from (t, j, i).
  const auto getNodeIdx = [&](idx_t t, idx_t j, idx_t i, bool includeHalo){

    // Adjust bounds.
    idx_t ijMin = includeHalo ? -halo : 0;
    idx_t ijMax = includeHalo ? N + halo : N;

    t = std::max(std::min(t, nTiles - 1), 0);
    j = std::max(std::min(j, ijMax), ijMin);
    i = std::max(std::min(i, ijMax), ijMin);

    const auto rowSize = (N + 2 * halo + 1);
    const auto tileSize = rowSize * rowSize;

    return idx2st(t * tileSize + (j + halo) * rowSize + i + halo);
  };

  // Helper functions to get cell idx from (t, j, i).
  const auto getCellIdx = [&](idx_t t, idx_t j, idx_t i, bool includeHalo){

    // Adjust bounds.
    idx_t ijMin = includeHalo ? -halo : 0;
    idx_t ijMax = includeHalo ? N + halo - 1: N - 1;

    t = std::max(std::min(t, nTiles - 1), 0);
    j = std::max(std::min(j, ijMax), ijMin);
    i = std::max(std::min(i, ijMax), ijMin);

    const auto rowSize = (N + 2 * halo);
    const auto tileSize = rowSize * rowSize;

    return idx2st(t * tileSize + (j + halo) * rowSize + i + halo);

  };

  // Return true for nodes interior to or on edge of tile.
  const auto interiorNode = [&](idx_t j, idx_t i){
    return i >= 0 and i < N + 1 and j >= 0 and j < N + 1;
  };

  // Return true for nodes on edge of tile.
  const auto edgeNode = [&](idx_t j, idx_t i){
    return interiorNode(j , i) and (i == 0 or i == N or j == 0 or j == N);
  };

  // Return true for cells interior to tile.
  const auto interiorCell = [&](idx_t j, idx_t i){
    return i >= 0 and i < N and j >= 0 and j < N;
  };

  // Return true for impossible combinations of (j, i) for nodes.
  const auto invalidNode = [&](idx_t j, idx_t i){

    return
      (i < 0 and j < 0) or // Bottom-left corner.
      (i > N and j < 0) or // Bottom-right corner.
      (i > N and j > N) or // Top-right corner.
      (i < 0 and j > N);   // Top-left corner.
  };

  // Return true for impossible combinations of (j, i) for nodes.
  const auto invalidCell = [&](idx_t j, idx_t i){

    return
      (i < 0     and j < 0    ) or // Bottom-left corner.
      (i > N - 1 and j < 0    ) or // Bottom-right corner.
      (i > N - 1 and j > N - 1) or // Top-right corner.
      (i < 0     and j > N - 1);   // Top-left corner.
  };

  const auto hackIndex = [&](idx_t i) -> double {
    if (i<=0) return i + 0.01;
    if (i >= N) return i -0.01;

    return static_cast<double>(i);
  };

  // ---------------------------------------------------------------------------
  // 2. GLOBAL CELL DISTRIBUTION
  // ---------------------------------------------------------------------------

  // Define cell record.
  // This currently keeps track of more information than we probably need.
  // This may be reduced once we've implemented halos.
  enum class CellType {
    UNDEFINED,
    OWNER,
    HALO
  };

  struct GlobalCell {
    CellType type{CellType::UNDEFINED};
    gidx_t  globalIdx{undefinedGlobalIdx}; // Global ID.
    idx_t   localIdx{undefinedIdx};         // Local ID.
    idx_t   part{undefinedIdx};             // Partition.
    const GlobalCell* remoteCell{};
  };

  // Define ij bounding box for each face (this partition only).
  struct BoundingBox {
    idx_t iBegin{std::numeric_limits<idx_t>::max()};
    idx_t iEnd{std::numeric_limits<idx_t>::min()};
    idx_t jBegin{std::numeric_limits<idx_t>::max()};
    idx_t jEnd{std::numeric_limits<idx_t>::min()};
  };

  // Make list of all cells.
  auto globalCells = std::vector<GlobalCell>(idx2st(nCellsGrid));

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
    const auto cellIdx = getCellIdx(t, j, i, false);
    auto& cell = globalCells[cellIdx];

    // cell is an owner.
    cell.type = CellType::OWNER;

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

  ATLAS_ASSERT(idxSum(cellIdxCount) == nCellsUnique);


  // Give possible edge-halo cells a unique global ID.
  const auto jacs = CubedSphereJacobian(csGrid);

  gidx_t cellGlobalIdxCount = nCellsUnique + 1;

  for (idx_t t = 0; t < nTiles; ++t) {
    for (idx_t j = ijBegin; j < ijEnd; ++j) {
      for (idx_t i = ijBegin; i < ijEnd; ++i) {

        // Skip if interior or invalid cell.
        if (interiorCell(j, i) or invalidCell(j, i)) continue;

        // Get cell.
        auto& cell = globalCells[getCellIdx(t, j, i, true)];

        // cell is and always will be a halo cell.
        cell.type = CellType::HALO;

        // Set global index.
        cell.globalIdx = cellGlobalIdxCount++;

        // Get owning cell
        const auto ownerTij = jacs.tijLocalToGlobal(PointTIJ(t, i, j));
        const auto ownerIdx =
          getCellIdx(ownerTij.t(), ownerTij.jRound(), ownerTij.iRound(), true);
        const auto& ownerCell = globalCells[ownerIdx];
        cell.remoteCell = &ownerCell;

      }
    }
  }

  ATLAS_ASSERT(cellGlobalIdxCount == nCellsTotal + 1);

  std::cout << "generated global cells" << std::endl;

  // ---------------------------------------------------------------------------
  // 3. GLOBAL OWNED NODE AND TYPE-A GHOST NODE DISTRIBUTION
  // ---------------------------------------------------------------------------

  // Define global node record.
  // This currently keeps track of more information than we probably need.
  // This may be reduced once we've implemented halos.

  enum class NodeType {
    UNDEFINED,
    OWNER,
    GHOST,
    HALO
  };

  struct GlobalNode {
    NodeType          type{NodeType::UNDEFINED};
    gidx_t            globalIdx{undefinedGlobalIdx};  // Global ID,
    idx_t             part{undefinedIdx};       // Partition node is on.
    idx_t             localIdx{undefinedIdx};   // Local ID of owned node.
    const GlobalNode* remoteNode{};             // Pointer to remote node.
    const GlobalCell* ownerCell{};              // Pointer to cell that owns node.
  };

  // Make list of all nodes.
  auto globalNodes = std::vector<GlobalNode>(idx2st(nNodesGrid));

  // Set counters for local node indices.
  auto nodeLocalIdxCount = std::vector<idx_t>(nParts, 0);

  // Set counter for global indices.
  gidx_t nodeGlobalOwnedIdxCount = 1;
  gidx_t nodeGlobalGhostIdxCount = nNodesUnique + 1;

  for (idx_t t = 0; t < nTiles; ++t) {
    for (idx_t j = ijBegin; j < ijEnd + 1; ++j) {
      for (idx_t i = ijBegin; i < ijEnd + 1; ++i) {

        // Skip if not a valid node.
       if (invalidNode(j, i)) continue;

        // Get this node.
        auto& node = globalNodes[getNodeIdx(t, j, i, true)];

        // Check for edge boundary ghosts.
        if (edgeNode(j, i) or !interiorNode(j, i)) {

          // Get owning node
          const auto ownerTij = jacs.tijLocalToGlobal(PointTIJ(t, i, j));
          const auto ownerIdx =
            getNodeIdx(ownerTij.t(), ownerTij.jRound(), ownerTij.iRound(), true);
          const auto& ownerNode = globalNodes[ownerIdx];


          const auto txy = jacs.tijToTxy(PointTIJ(t, i, j));
          const auto gtxy = jacs.txyLocalToGlobal(txy);

          std::cout << t << " " << j << " " << i << std::endl;
          std::cout << ownerTij.t() << " " << ownerTij.j() << " " << ownerTij.i() << std::endl;
          std::cout << ownerTij.t() << " " << ownerTij.jRound() << " " << ownerTij.iRound() << std::endl << std::endl;

          std::cout << txy.t() << " " << txy.y() << " " << txy.x() << std::endl;
          std::cout << gtxy.t() << " " << gtxy.y() << " " << gtxy.x() << std::endl <<std::endl;


          ATLAS_ASSERT(interiorNode(ownerTij.jRound(), ownerTij.iRound()));

          if (&node != &ownerNode) {

            // Node is and always will be a ghost.
            node.type = NodeType::GHOST;

            // Get owning cell.
            node.ownerCell = &globalCells[getCellIdx(t, j - 1, i - 1, true)];

            // Ghost node on tile edge. Point to owner.
            node.remoteNode = &ownerNode;

            // Set global index.
            node.globalIdx = nodeGlobalGhostIdxCount++;

            // Done with this loop iteration.
            continue;
          }

        }

        // Node is an owner.
        node.type = NodeType::OWNER;

        // Get owning cell.
        node.ownerCell = &globalCells[getCellIdx(t, j - 1, i - 1, false)];

        // Set partition.
        node.part = node.ownerCell->part;

        // Owner node. Set local index.
        node.localIdx = nodeLocalIdxCount[idx2st(node.part)]++;

        // Set global index.
        node.globalIdx = nodeGlobalOwnedIdxCount++;

      }
    }
  }

  ATLAS_ASSERT(nodeGlobalOwnedIdxCount == nNodesUnique + 1);
  ATLAS_ASSERT(nodeGlobalGhostIdxCount == nNodesTotal + 1);
  ATLAS_ASSERT(idxSum(nodeLocalIdxCount) == nNodesUnique);
  std::cout << "generated global nodes" << std::endl;

  // Find local cells.

  // Define local cell record.
  struct LocalCell {
    CellType      type{CellType::UNDEFINED};
    gidx_t        globalIdx{undefinedGlobalIdx};
    idx_t         remoteIdx{undefinedIdx};
    idx_t         part{undefinedIdx};
    idx_t         t{}, j{}, i{};
  };

  auto localCells = std::vector<LocalCell>{};

  const auto isOwner = [&](const GlobalCell& cell){
    return cell.type == CellType::OWNER and cell.part == st2idx(thisPart);
  };

  // Loop over all possible local cells.
  for (idx_t t = 0; t < nTiles; ++t) {

    // Limit range to bounds recorded earlier.
    const auto bounds = cellBounds[idx2st(t)];
    const auto iBegin = bounds.iBegin - halo;
    const auto iEnd   = bounds.iEnd   + halo;
    const auto jBegin = bounds.jBegin - halo;
    const auto jEnd   = bounds.jEnd   + halo;

    for (idx_t j = jBegin; j < jEnd; ++j) {
      for (idx_t i = iBegin; i < iEnd; ++i) {

        if (invalidCell(j, i)) continue;

        // Get cell
        auto& globalCell = globalCells[getCellIdx(t, j, i, true)];

        // Define data copying lambda.
        const auto copyCellData = [&](const CellType type) {

          auto localCell = LocalCell{};

          localCell.type = type;
          localCell.globalIdx = globalCell.globalIdx;
          if (globalCell.remoteCell) {
            // This is an exterior halo cell.
            localCell.remoteIdx = globalCell.remoteCell->localIdx;
            localCell.part = globalCell.remoteCell->part;
          } else {
            // This is an interior halo cell or owner.
            localCell.remoteIdx = globalCell.localIdx;
            localCell.part = globalCell.part;
          }
          localCell.t = t;
          localCell.i = i;
          localCell.j = j;
          localCells.push_back(localCell);

          return;

        };

        // Check if cell is an owner.
        if (isOwner(globalCell)) {

          copyCellData(CellType::OWNER);

        } else {

          // Check if cell is a halo.
          auto isHalo = false;
          for (idx_t jHalo = j - halo; jHalo < j + halo + 1; ++jHalo) {
            for (idx_t iHalo = i - halo; iHalo < i + halo + 1; ++iHalo) {

              if (invalidCell(jHalo, iHalo)) continue;

              // Is there a nearby owner cell?
              const auto ownerCell = globalCells[getCellIdx(t, jHalo, iHalo, true)];
              isHalo = isHalo or isOwner(ownerCell);
            }
          }

          if (isHalo) {

            copyCellData(CellType::HALO);

            // Overwrite globalCell properties to help find nodes later.
            globalCell.type = CellType::HALO;
            globalCell.part = st2idx(thisPart);

          }
        }

      }
    }
  }


  // Halo cells are mixed in with owner cells. Sort by cell type to correct.
  std::stable_sort(localCells.begin(), localCells.end(),
    [](const LocalCell& cellA, const LocalCell& cellB){
      return cellA.type < cellB.type;
    });

  std::cout << "generated local cells" << std::endl;

  // Define local node record.
  struct LocalNode {
    NodeType      type{NodeType::UNDEFINED};
    gidx_t        globalIdx{undefinedGlobalIdx};
    idx_t         remoteIdx{undefinedIdx};
    idx_t         part{undefinedIdx};
    idx_t         t{}, j{}, i{};
    idx_t*        localIdx{};
  };

  auto localNodes = std::vector<LocalNode>{};

  // Loop over all possible local nodes.
  for (idx_t t = 0; t < nTiles; ++t) {

      // Limit range to bounds recorded earlier.
      const auto bounds = cellBounds[idx2st(t)];
      const auto iBegin = bounds.iBegin - halo;
      const auto iEnd   = bounds.iEnd   + halo + 1;
      const auto jBegin = bounds.jBegin - halo;
      const auto jEnd   = bounds.jEnd   + halo + 1;

    for (idx_t j = jBegin; j < jEnd; ++j) {
      for (idx_t i = iBegin; i < iEnd; ++i) {

        if (invalidNode(j, i)) continue;

        // Get node.
        auto& globalNode = globalNodes[getNodeIdx(t, j, i, true)];

        // Define data copying lambda.
        const auto copyNodeData = [&](const NodeType type) {

          auto localNode = LocalNode{};

          localNode.type = type;
          localNode.globalIdx = globalNode.globalIdx;
          if (globalNode.remoteNode) {
            // This is an exterior ghost node.
            ATLAS_ASSERT(globalNode.remoteNode->type == NodeType::OWNER);
            localNode.remoteIdx = globalNode.remoteNode->localIdx;
            localNode.part = globalNode.remoteNode->part;
          } else {
            // This is an interior ghost node or owner.
            localNode.remoteIdx = globalNode.localIdx;
            localNode.part = globalNode.part;
          }
          localNode.t = t;
          localNode.i = i;
          localNode.j = j;

          // Need to overwrite this later for cell connectivity.
          localNode.localIdx = &globalNode.localIdx;

          localNodes.push_back(localNode);

          return;

        };

        // Check if node is an owner.
        ATLAS_ASSERT(globalNode.ownerCell);
        if (globalNode.ownerCell->part == st2idx(thisPart) and
            globalNode.ownerCell->type == CellType::OWNER and
            !globalNode.remoteNode) {

          // Node is owner.
          copyNodeData(NodeType::OWNER);
        } else {

          // Node is possibly a ghost or halo point.

          // Get four neighbouring cells of node.
          const auto& cell0 = globalCells[getCellIdx(t, j - 1, i - 1, true)];
          const auto& cell1 = globalCells[getCellIdx(t, j - 1, i    , true)];
          const auto& cell2 = globalCells[getCellIdx(t, j    , i    , true)];
          const auto& cell3 = globalCells[getCellIdx(t, j    , i - 1, true)];

          // Check if neighbour cells are on this partition.
          if (cell0.part == st2idx(thisPart) or
              cell1.part == st2idx(thisPart) or
              cell2.part == st2idx(thisPart) or
              cell3.part == st2idx(thisPart)) {

            // Node is flagged as a halo if none of the neighbouring cells are
            // owners. Otherwise it's a ghost node.
            const auto isGhost =
              (cell0.part == st2idx(thisPart) and cell0.type == CellType::OWNER) or
              (cell1.part == st2idx(thisPart) and cell1.type == CellType::OWNER) or
              (cell2.part == st2idx(thisPart) and cell2.type == CellType::OWNER) or
              (cell3.part == st2idx(thisPart) and cell3.type == CellType::OWNER);

            if (isGhost) copyNodeData(NodeType::GHOST);
            else copyNodeData(NodeType::HALO);


          }

        }

      }
    }
  }

  std::cout << "generated global nodes" << std::endl;

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

    if (node.type == NodeType::OWNER) {
      ATLAS_ASSERT(*node.localIdx == nodeLocalIdx);
    } else {
      *node.localIdx = nodeLocalIdx;
    }
    ++nodeLocalIdx;
  }

  std::cout << "generated local nodes" << std::endl;

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
  auto nodesHalo      = array::make_view<int, 1>(mesh.nodes().halo());
  auto nodesFlags     = array::make_view<int, 1>(mesh.nodes().flags());

  // Set fields.
  nodeLocalIdx = 0;
  for (const auto& node : localNodes) {

    // Set global index.
    nodesGlobalIdx(nodeLocalIdx) = node.globalIdx;

    // Set node remote index.
    nodesRemoteIdx(nodeLocalIdx) = node.remoteIdx;

    // Set node partition.
    nodesPart(nodeLocalIdx) = node.part;

    // Set xy.
    const auto txyLocal = jacs.tijToTxy(PointTIJ(node.t, node.i, node.j));
    const auto txyGlobal = jacs.txyLocalToGlobal(txyLocal);

    nodesXy(nodeLocalIdx, XX) = txyLocal.x();
    nodesXy(nodeLocalIdx, YY) = txyLocal.y();

    // Set lon-lat.
    const auto lonLat = csProjection->lonlat(txyGlobal.xy());
    nodesLonLat(nodeLocalIdx, LON) = lonLat.lon();
    nodesLonLat(nodeLocalIdx, LAT) = lonLat.lat();

    // Set flags.
    Topology::reset(nodesFlags(nodeLocalIdx));
    switch(node.type) {

      case NodeType::UNDEFINED : {
        ATLAS_ASSERT(NodeType::UNDEFINED);
        break;
      }
      case NodeType::OWNER : {
        nodesGhost(nodeLocalIdx) = 0;
        nodesHalo(nodeLocalIdx) = 0;
        break;
      }
      case NodeType::GHOST : {
        nodesGhost(nodeLocalIdx) = 1;
        nodesHalo(nodeLocalIdx) = 0;
        Topology::set(nodesFlags(nodeLocalIdx), Topology::GHOST);
        break;
      }
      case NodeType::HALO : {
        nodesGhost(nodeLocalIdx) = 1;
        nodesHalo(nodeLocalIdx) = 1;
        Topology::set(nodesFlags(nodeLocalIdx), Topology::GHOST);
        break;
      }
    }

    ++nodeLocalIdx;
  }

  std::cout << "nodes to fields" << std::endl;

  // ---------------------------------------------------------------------------
  // 6. ASSIGN CELLS TO MESH
  // ---------------------------------------------------------------------------

  // Resize cells.
  mesh.cells().add(new mesh::temporary::Quadrilateral(), st2idx(localCells.size()));

  // Set field views.
  auto cellsGlobalIdx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
  auto cellsPart      = array::make_view<int, 1>(mesh.cells().partition());

  // Set local cells.
  auto& nodeConnectivity = mesh.cells().node_connectivity();

  idx_t cellLocalIdx = mesh.cells().elements(0).begin();
  for (const auto& cell : localCells) {

    // Get four surroundings nodes.
    const auto& node0 = globalNodes[getNodeIdx(cell.t, cell.j    , cell.i    , true)];
    const auto& node1 = globalNodes[getNodeIdx(cell.t, cell.j    , cell.i + 1, true)];
    const auto& node2 = globalNodes[getNodeIdx(cell.t, cell.j + 1, cell.i + 1, true)];
    const auto& node3 = globalNodes[getNodeIdx(cell.t, cell.j + 1, cell.i    , true)];

    const auto quadNodeIdx = std::array<idx_t, 4> {
      node0.localIdx, node1.localIdx, node2.localIdx, node3.localIdx};

    // Set connectivity.
    nodeConnectivity.set(cellLocalIdx, quadNodeIdx.data());

    // Set global index.
    cellsGlobalIdx(cellLocalIdx) = cell.globalIdx;

    // Set partition.
    cellsPart(cellLocalIdx) = cell.part;

    ++cellLocalIdx;
  }

  std::cout << "cells to fields" << std::endl;

  // ---------------------------------------------------------------------------
  // 7. FINALISE
  // ---------------------------------------------------------------------------

  mesh.nodes().global_index().metadata().set("human_readable", true);
  mesh.nodes().global_index().metadata().set("min", 1);
  mesh.nodes().global_index().metadata().set("max", nNodesTotal);
  mesh.nodes().metadata().set("parallel", true);

  mesh.cells().global_index().metadata().set("human_readable", true);
  mesh.cells().global_index().metadata().set("min", 1);
  mesh.cells().global_index().metadata().set("max", nCellsTotal);
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

