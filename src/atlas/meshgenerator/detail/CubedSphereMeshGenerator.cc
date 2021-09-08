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

  // This option sets number of partitions.
  options.set( "nb_parts", mpi::size() );

  // This option sets the part that will be generated.
  options.set("part", mpi::rank());

  // This options sets the number of halo elements around each partition.
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

void CubedSphereMeshGenerator::generate(const Grid& grid,
  const grid::Distribution& distribution, Mesh& mesh) const {

  // Check for correct grid and need for mesh
  ATLAS_ASSERT(!mesh.generated());
  if (!CubedSphereGrid(grid)) {
    throw_Exception("CubedSphereMeshGenerator can only work "
    "with a cubedsphere grid.", Here());
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
  // Mesh generator creates a cubed sphere mesh by generating individual meshes
  // for each cubed sphere tile and then working out the correspondence between
  // overlapping nodes and cells.
  //
  // Meshgenerator places cell at each grid point on this partition. Halo
  // cells are added to the mesh if options.get("halo") > 0. The halo cells may
  // either be interior to the tile, || exterior. Interior halo cells share
  // their xy coordinate and global ID with the corresponding cells on other
  // partitions. Exterior halo cells have a unique global ID and their xy
  // coordinates are extrapolated from the interior tile. The global IDs of non-
  // halo cells match the global IDs of the cell-centre grid points.
  //
  // Nodes are added around cells. The nodes around interior cells are assigned
  // an owner by the following rules:
  //    * node (i > 0, j > 0) is owned by cell (i - 1, j - 1)
  //    * node (i = 0, j > 0) is owned by cell (i    , j - 1)
  //    * node (i > 0, j = 0) is owned by cell (i - 1, j    )
  //    * node (i = 0, j = 0) is owned by cell (i    , j    )
  //
  // The partition of the owning cell determines the partition of the node.
  // Ghost nodes are added to the mesh to complete cells at partition
  // boundaries, cells exterior to the tiles, oron tile edges which do not
  // own nodes.
  //
  // There are several stages to the mesh generator:
  //    1. Preamble.
  //    2. Define global cell distribution.
  //    3. Define global node distribution.
  //    4. Locate local cells.
  //    5. Locate local nodes
  //    6. Assign nodes to mesh.
  //    7. Assign cells to mesh.
  //    8. Finalise.

  // ---------------------------------------------------------------------------
  // 1. PREAMBLE
  //    Setup some general parameters of the mesh, as well as define useful
  //    lambda functions.
  // ---------------------------------------------------------------------------

  using Topology = atlas::mesh::Nodes::Topology;
  using atlas::array::make_datatype;
  using atlas::array::make_shape;

  using namespace detail::cubedsphere;


  // Check for correct grid stagger.
  const Stagger::s stagger = getStagger(csGrid.name());
  if (stagger != Stagger::CELL) {
    throw_Exception("CubedSphereMeshGenerator will only work with a"
    "cell-centroid grid. Try NodalCubedSphereMeshGenerator instead.", Here());
  }

  // Get dimensions of grid
  const idx_t N = csGrid.N();

  // Get size of halo.
  idx_t nHalo = 0;
  options.get("halo", nHalo);

  // Set default bounds for (i, j) based loops.
  idx_t ijBegin = -nHalo;
  idx_t ijEnd = N + nHalo;

  // Unique non-halo nodes and cells.
  const idx_t nNodesUnique = 6 * N * N + 2;
  const idx_t nCellsUnique = 6 * N * N;

  // Total array sizes (including invalid corner ijs).
  const idx_t nNodesArray   = 6 * (N + 2 * nHalo + 1) * (N + 2 * nHalo + 1);
  const idx_t nCellsArray   = 6 * (N + 2 * nHalo) * (N + 2 * nHalo);

  // Total number of possible cells and nodes (including halos).
  const idx_t nNodesTotal  = nNodesArray - 6 * 4 * nHalo * nHalo;
  const idx_t nCellsTotal  = nCellsArray - 6 * 4 * nHalo * nHalo;

  Log::debug() << "Number of cells per tile edge = "
    << std::to_string(N) << std::endl;

  // Define bad index values.
  constexpr idx_t undefinedIdx = -1;
  constexpr gidx_t undefinedGlobalIdx = -1;

  // Projection and tiles
  const auto* const csProjection = castProjection(csGrid.projection().get());

  // Get partition information.
  const auto nParts = options.get<size_t>("nb_parts");
  const auto thisPart = options.get<size_t>("part");

  // Define an index counter.
  const auto idxSum = [](const std::vector<idx_t>& idxCounts) -> idx_t {
    return std::accumulate(idxCounts.begin(), idxCounts.end(), 0);};

  // Helper functions to get node idx from (t, j, i).
  const auto getNodeIdx = [&](idx_t t, idx_t j, idx_t i, bool includeHalo) -> size_t {

    // Adjust bounds.
    idx_t ijMin = includeHalo ? -nHalo : 0;
    idx_t ijMax = includeHalo ? N + nHalo : N;

    t = std::max(std::min(t, 5), 0);
    j = std::max(std::min(j, ijMax), ijMin);
    i = std::max(std::min(i, ijMax), ijMin);

    const idx_t rowSize = (N + 2 * nHalo + 1);
    const idx_t tileSize = rowSize * rowSize;

    return static_cast<size_t>(t * tileSize + (j + nHalo) * rowSize + i + nHalo);
  };

  // Helper functions to get cell idx from (t, j, i).
  const auto getCellIdx = [&](idx_t t, idx_t j, idx_t i, bool includeHalo) -> size_t {

    // Adjust bounds.
    idx_t ijMin = includeHalo ? -nHalo : 0;
    idx_t ijMax = includeHalo ? N + nHalo - 1: N - 1;

    t = std::max(std::min(t, 5), 0);
    j = std::max(std::min(j, ijMax), ijMin);
    i = std::max(std::min(i, ijMax), ijMin);

    const idx_t rowSize = (N + 2 * nHalo);
    const idx_t tileSize = rowSize * rowSize;

    return static_cast<size_t>(t * tileSize + (j + nHalo) * rowSize + i + nHalo);

  };

  // Return true for nodes interior to or on edge of tile.
  const auto interiorNode = [&](idx_t j, idx_t i) -> bool {
    return i >= 0 && i < N + 1 && j >= 0 && j < N + 1;
  };

  // Return true for nodes on edge of tile.
  const auto edgeNode = [&](idx_t j, idx_t i) -> bool {
    return interiorNode(j , i) && (i == 0 ||  i == N ||  j == 0 ||  j == N);
  };

  // Return true for cells interior to tile.
  const auto interiorCell = [&](idx_t j, idx_t i) -> bool {
    return i >= 0 && i < N && j >= 0 && j < N;
  };

  // Return true for impossible combinations of (j, i) for nodes.
  const auto invalidNode = [&](idx_t j, idx_t i) -> bool {

    return
      (i < 0 && j < 0) ||  // Bottom-left corner.
      (i > N && j < 0) ||  // Bottom-right corner.
      (i > N && j > N) ||  // Top-right corner.
      (i < 0 && j > N);   // Top-left corner.
  };

  // Return true for impossible combinations of (j, i) for nodes.
  const auto invalidCell = [&](idx_t j, idx_t i) -> bool {

    return
      (i < 0     && j < 0    ) ||  // Bottom-left corner.
      (i > N - 1 && j < 0    ) ||  // Bottom-right corner.
      (i > N - 1 && j > N - 1) ||  // Top-right corner.
      (i < 0     && j > N - 1);   // Top-left corner.
  };

  // ---------------------------------------------------------------------------
  // 2. GLOBAL CELL DISTRIBUTION
  //    Need to construct a lightweight global mesh. This will help us pair up
  //    local and remote indices on different partitions.
  // ---------------------------------------------------------------------------

  // Define cell types.
  enum class CellType {
    UNDEFINED,
    OWNER,
    HALO
  };

  // Define cell record.
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

  // Make list of all possible cells.
  auto globalCells = std::vector<GlobalCell>(static_cast<size_t>(nCellsArray));

  // Initialise bounding box.
  auto cellBounds = std::vector<BoundingBox>(static_cast<size_t>(6));

  // Set xy and tji grid iterators.
  auto tjiIt = csGrid.tij().begin();

  // Set counters for cell local indices.
  auto cellLocalIdxCount = std::vector<idx_t>(nParts, 0);

  for (gidx_t gridIdx = 0; gridIdx < csGrid.size(); ++gridIdx) {

    // Grid (t, j, i) order does not have to match mesh (t, j, i), although
    // in practice it probably does.

    // Get cell index.
    const idx_t t = (*tjiIt).t();
    const idx_t j = (*tjiIt).j();
    const idx_t i = (*tjiIt).i();
    const size_t cellIdx = getCellIdx(t, j, i, false);
    GlobalCell& cell = globalCells[cellIdx];

    // cell is an owner.
    cell.type = CellType::OWNER;

    // Set global index.
    cell.globalIdx = gridIdx + 1;

    // Set partition and remote index.
    cell.part = distribution.partition(gridIdx);

    // Set local index.
    cell.localIdx = cellLocalIdxCount[static_cast<size_t>(cell.part)]++;

    if (cell.part == static_cast<idx_t>(thisPart)) {

      // Keep track of local (t, j, i) bounds.
      BoundingBox& bounds = cellBounds[static_cast<size_t>(t)];
      bounds.iBegin = std::min(bounds.iBegin, i    );
      bounds.iEnd   = std::max(bounds.iEnd  , i + 1);
      bounds.jBegin = std::min(bounds.jBegin, j    );
      bounds.jEnd   = std::max(bounds.jEnd  , j + 1);

    }
    // Increment iterators.
    ++tjiIt;
  }

  ATLAS_ASSERT(idxSum(cellLocalIdxCount) == nCellsUnique);


  // Give possible edge-halo cells a unique global ID.
  const auto jacobian = NeighbourJacobian(csGrid);

  gidx_t cellGlobalIdxCount = nCellsUnique + 1;

  for (idx_t t = 0; t < 6; ++t) {
    for (idx_t j = ijBegin; j < ijEnd; ++j) {
      for (idx_t i = ijBegin; i < ijEnd; ++i) {

        // Skip if interior or invalid cell.
        if (interiorCell(j, i) ||  invalidCell(j, i)) continue;

        // Get cell.
        GlobalCell& cell = globalCells[getCellIdx(t, j, i, true)];

        // cell is and always will be a halo cell.
        cell.type = CellType::HALO;

        // Set global index.
        cell.globalIdx = cellGlobalIdxCount++;

        // Get i, j and t of owning cell
        PointIJ ijGlobal;
        idx_t tGlobal;

        // Cells reside at half integer ij values.
        std::tie(ijGlobal, tGlobal) =
          jacobian.ijLocalToGlobal(PointIJ(i + 0.5, j + 0.5), t);

        const size_t ownerIdx =
          getCellIdx(tGlobal, ijGlobal.jCell(), ijGlobal.iCell(), true);
        const GlobalCell& ownerCell = globalCells[ownerIdx];
        cell.remoteCell = &ownerCell;

      }
    }
  }

  ATLAS_ASSERT(cellGlobalIdxCount == nCellsTotal + 1);

  // ---------------------------------------------------------------------------
  // 3. GLOBAL NODE DISTRIBUTION
  //    Construct a lightweight global distribution of nodes. Again, this will
  //    help us match up local and remote indices accross partitions.
  // ---------------------------------------------------------------------------


  // Define node types.
  enum class NodeType {
    UNDEFINED,
    OWNER,
    GHOST,
    HALO
  };

  // Define global node record.
  struct GlobalNode {
    NodeType          type{NodeType::UNDEFINED};
    gidx_t            globalIdx{undefinedGlobalIdx};  // Global ID,
    idx_t             part{undefinedIdx};       // Partition node is on.
    idx_t             localIdx{undefinedIdx};   // Local ID of owned node.
    const GlobalNode* remoteNode{};             // Pointer to remote node.
    const GlobalCell* ownerCell{};              // Pointer to cell that owns node.
  };

  // Make list of all nodes.
  auto globalNodes = std::vector<GlobalNode>(static_cast<size_t>(nNodesArray));

  // Set counters for local node indices.
  auto nodeLocalIdxCount = std::vector<idx_t>(nParts, 0);

  // Set counter for global indices.
  gidx_t nodeGlobalOwnedIdxCount = 1;
  gidx_t nodeGlobalGhostIdxCount = nNodesUnique + 1;

  for (idx_t t = 0; t < 6; ++t) {
    for (idx_t j = ijBegin; j < ijEnd + 1; ++j) {
      for (idx_t i = ijBegin; i < ijEnd + 1; ++i) {

        // Skip if not a valid node.
       if (invalidNode(j, i)) continue;

        // Get this node.
        GlobalNode& node = globalNodes[getNodeIdx(t, j, i, true)];

        // Check for edge boundary ghosts.
        if (edgeNode(j, i) ||  !interiorNode(j, i)) {

          // Get owning node
          PointIJ ijGlobal;
          idx_t tGlobal;

          // Nodes reside on interger ij values.
          std::tie(ijGlobal, tGlobal) = jacobian.ijLocalToGlobal(PointIJ(i, j), t);

          const size_t ownerIdx =
            getNodeIdx(tGlobal, ijGlobal.jNode(), ijGlobal.iNode(), true);
          const GlobalNode& ownerNode = globalNodes[ownerIdx];

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

        // Get partition from cell.
        node.part = node.ownerCell->part;

        // Owner node. Set local index.
        node.localIdx = nodeLocalIdxCount[static_cast<size_t>(node.part)]++;

        // Set global index.
        node.globalIdx = nodeGlobalOwnedIdxCount++;

      }
    }
  }

  ATLAS_ASSERT(nodeGlobalOwnedIdxCount == nNodesUnique + 1);
  ATLAS_ASSERT(nodeGlobalGhostIdxCount == nNodesTotal + 1);
  ATLAS_ASSERT(idxSum(nodeLocalIdxCount) == nNodesUnique);

  // ---------------------------------------------------------------------------
  // 4. LOCATE LOCAL CELLS.
  //    Now that we know where all the nodes and cells are, we can make a list
  //    of local cells to add the the mesh. "Owner" cells correspond to grid
  //    points on this parition. "Halo" cells are at most nHalo grid points
  //    away from an owner cell.
  // ---------------------------------------------------------------------------


  // Define local cell record.
  struct LocalCell {
    GlobalCell*   globalCellPtr{};                // Pointer to global cell.
    gidx_t        globalIdx{undefinedGlobalIdx};  // Global index.
    CellType      type{CellType::UNDEFINED};      // Cell Type.
    idx_t         remoteIdx{undefinedIdx};        // Remote index.
    idx_t         part{undefinedIdx};             // Partition.
    idx_t         t{}, j{}, i{};                  // t, i, and j.
  };

  auto localCells = std::vector<LocalCell>{};

  const auto isOwner = [&](const GlobalCell& cell) -> bool {
    return cell.type == CellType::OWNER && cell.part == static_cast<idx_t>(thisPart);
  };

  // Loop over all possible local cells.
  for (idx_t t = 0; t < 6; ++t) {

    // Limit range to bounds recorded earlier.
    const BoundingBox bounds = cellBounds[static_cast<size_t>(t)];
    const idx_t iBegin = bounds.iBegin - nHalo;
    const idx_t iEnd   = bounds.iEnd   + nHalo;
    const idx_t jBegin = bounds.jBegin - nHalo;
    const idx_t jEnd   = bounds.jEnd   + nHalo;

    for (idx_t j = jBegin; j < jEnd; ++j) {
      for (idx_t i = iBegin; i < iEnd; ++i) {

        if (invalidCell(j, i)) continue;

        // Get cell
        GlobalCell& globalCell = globalCells[getCellIdx(t, j, i, true)];

        // Define data copying lambda. This manipulates data via the capture list.
        const auto copyCellData = [&](const CellType type) -> void {

          // Add cell to list.
          localCells.emplace_back();
          LocalCell& localCell = localCells.back();

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


          if (type == CellType::HALO) {
            // Point to global data so we can overwrite members later.
            // This indirection is needed to find local nodes.
            localCell.globalCellPtr = &globalCell;
          }

        };

        // Check if cell is an owner.
        if (isOwner(globalCell)) {

          copyCellData(CellType::OWNER);

        } else {

          // Check if cell is a halo.
          bool isHalo = false;
          for (idx_t jHalo = j - nHalo; jHalo < j + nHalo + 1; ++jHalo) {
            for (idx_t iHalo = i - nHalo; iHalo < i + nHalo + 1; ++iHalo) {

              if (invalidCell(jHalo, iHalo)) continue;

              // Is there a nearby owner cell?
              const GlobalCell& ownerCell =
                globalCells[getCellIdx(t, jHalo, iHalo, true)];
              isHalo = isHalo ||  isOwner(ownerCell);
            }
          }

          if (isHalo) {

            copyCellData(CellType::HALO);

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

  // Overwrite partition and type for halos.
  for (LocalCell& cell : localCells) {

    if (cell.type == CellType::HALO) {
      cell.globalCellPtr->part = static_cast<idx_t>(thisPart);
      cell.globalCellPtr->type = CellType::HALO;
    }
  }

  // ---------------------------------------------------------------------------
  // 5. LOCATE LOCAL NODES.
  //    We can now locate the nodes surrounding the owner and halo cells that
  //    we are adding to the mesh. Nodes which are owned by a cell are marked
  //    as owners. Nodes which neighbour at least one owned cell are marked as
  //    ghosts. Nodes which only neighbour halo cells are marked as halo nodes.
  // ---------------------------------------------------------------------------

  // Define local node record.
  struct LocalNode {
    GlobalNode*   globalNodePtr{};                // Pointer to global node.
    gidx_t        globalIdx{undefinedGlobalIdx};  // Global index.
    NodeType      type{NodeType::UNDEFINED};      // Node type.
    idx_t         remoteIdx{undefinedIdx};        // Remote index.
    idx_t         part{undefinedIdx};             // Partition.
    idx_t         t{}, j{}, i{};                  // t, j and i.
  };

  auto localNodes = std::vector<LocalNode>{};

  // Loop over all possible local nodes.
  for (idx_t t = 0; t < 6; ++t) {

      // Limit range to bounds recorded earlier.
      const BoundingBox bounds = cellBounds[static_cast<size_t>(t)];
      const idx_t iBegin = bounds.iBegin - nHalo;
      const idx_t iEnd   = bounds.iEnd   + nHalo + 1;
      const idx_t jBegin = bounds.jBegin - nHalo;
      const idx_t jEnd   = bounds.jEnd   + nHalo + 1;

    for (idx_t j = jBegin; j < jEnd; ++j) {
      for (idx_t i = iBegin; i < iEnd; ++i) {

        if (invalidNode(j, i)) continue;

        // Get node.
        GlobalNode& globalNode = globalNodes[getNodeIdx(t, j, i, true)];

        // Define data copying lambda. This manipulates data via the capture list.
        const auto copyNodeData = [&](const NodeType type) -> void {

          // Add node to list.
          localNodes.emplace_back();
          LocalNode& localNode = localNodes.back();

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

          // Need to overwrite localIdx later for cell connectivity.
          localNode.globalNodePtr = &globalNode;

        };

        // Check if node is an owner.
        ATLAS_ASSERT(globalNode.ownerCell);
        if (globalNode.ownerCell->part == static_cast<idx_t>(thisPart) and
            globalNode.ownerCell->type == CellType::OWNER and
            !globalNode.remoteNode) {

          // Node is owner.
          copyNodeData(NodeType::OWNER);
        } else {

          // Node is possibly a ghost or halo point.

          // Get four neighbouring cells of node.
          const GlobalCell& cell0 = globalCells[getCellIdx(t, j - 1, i - 1, true)];
          const GlobalCell& cell1 = globalCells[getCellIdx(t, j - 1, i    , true)];
          const GlobalCell& cell2 = globalCells[getCellIdx(t, j    , i    , true)];
          const GlobalCell& cell3 = globalCells[getCellIdx(t, j    , i - 1, true)];

          // Check if neighbour cells are on this partition.
          if (cell0.part == static_cast<idx_t>(thisPart) or
              cell1.part == static_cast<idx_t>(thisPart) or
              cell2.part == static_cast<idx_t>(thisPart) or
              cell3.part == static_cast<idx_t>(thisPart)) {

            // Node is flagged as a halo if none of these neighbour cells is an
            // owner. Otherwise it's a ghost node.
            const bool isGhost =
              (cell0.part == static_cast<idx_t>(thisPart) && cell0.type == CellType::OWNER) ||
              (cell1.part == static_cast<idx_t>(thisPart) && cell1.type == CellType::OWNER) ||
              (cell2.part == static_cast<idx_t>(thisPart) && cell2.type == CellType::OWNER) ||
              (cell3.part == static_cast<idx_t>(thisPart) && cell3.type == CellType::OWNER);

            if (isGhost) copyNodeData(NodeType::GHOST);
            else copyNodeData(NodeType::HALO);

          }
        }
      }
    }
  }

  // Ghost nodes are mixed in with owner nodes. Sort by node type to correct.
  std::stable_sort(localNodes.begin(), localNodes.end(),
    [](const LocalNode& nodeA, const LocalNode& nodeB){
      return nodeA.type < nodeB.type;
    });

  // Check/Set local indices.
  // Note: value of localIdx for owner nodes should now match element number of
  // localNodes vector. localIdx value will be overwritten for non-owner nodes.
  idx_t nodeLocalIdx = 0;
  for (LocalNode& node : localNodes) {

    if (node.type == NodeType::OWNER) {
      ATLAS_ASSERT(node.globalNodePtr->localIdx == nodeLocalIdx);
    } else {
      node.globalNodePtr->localIdx = nodeLocalIdx;
    }
    ++nodeLocalIdx;
  }

  // ---------------------------------------------------------------------------
  // 6. ASSIGN NODES TO MESH
  //    In addition to the usual fields, we will add t, i and j to mesh.nodes().
  // ---------------------------------------------------------------------------

  // Resize nodes.
  auto& nodes = mesh.nodes();
  nodes.resize(static_cast<idx_t>(localNodes.size()));

  // Add extra field.
  Field ijtField =
    nodes.add(Field("ijt", make_datatype<idx_t>(), make_shape(nodes.size(), 3)));
  ijtField.set_variables(3);

  // Get field views
  auto nodesGlobalIdx = array::make_view<gidx_t, 1>(nodes.global_index());
  auto nodesRemoteIdx = array::make_indexview<idx_t, 1>(nodes.remote_index());
  auto nodesXy        = array::make_view<double, 2>(nodes.xy());
  auto nodesLonLat    = array::make_view<double, 2>(nodes.lonlat());
  auto nodesPart      = array::make_view<int, 1>(nodes.partition());
  auto nodesGhost     = array::make_view<int, 1>(nodes.ghost());
  auto nodesHalo      = array::make_view<int, 1>(nodes.halo());
  auto nodesFlags     = array::make_view<int, 1>(nodes.flags());
  auto nodesIjt       = array::make_view<idx_t, 2>(ijtField);


  // Set fields.
  nodeLocalIdx = 0;
  for (const LocalNode& node : localNodes) {

    // Set global index.
    nodesGlobalIdx(nodeLocalIdx) = node.globalIdx;

    // Set node remote index.
    nodesRemoteIdx(nodeLocalIdx) = node.remoteIdx;

    // Set node partition.
    nodesPart(nodeLocalIdx) = node.part;

    // Set xy.
    const PointXY xyLocal = jacobian.xy(PointIJ(node.i, node.j), node.t);
    const PointXY xyGlobal = jacobian.xyLocalToGlobal(xyLocal, node.t).first;

    nodesXy(nodeLocalIdx, XX) = xyLocal.x();
    nodesXy(nodeLocalIdx, YY) = xyLocal.y();

    // Set lon-lat.
    const PointLonLat lonLat = csProjection->lonlat(xyGlobal);
    nodesLonLat(nodeLocalIdx, LON) = lonLat.lon();
    nodesLonLat(nodeLocalIdx, LAT) = lonLat.lat();

    // Set tij.
    nodesIjt(nodeLocalIdx, Coordinates::I) = node.i;
    nodesIjt(nodeLocalIdx, Coordinates::J) = node.j;
    nodesIjt(nodeLocalIdx, Coordinates::T) = node.t;

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

  // ---------------------------------------------------------------------------
  // 7. ASSIGN CELLS TO MESH
  //    Again, we'll add t, i and j to mesh.cells().
  // ---------------------------------------------------------------------------

  auto& cells = mesh.cells();

  // Resize cells.
  cells.add(new mesh::temporary::Quadrilateral(), static_cast<idx_t>(localCells.size()));

  // Add extra fields.
  ijtField = cells.add(
      Field("ijt", make_datatype<idx_t>(), make_shape(cells.size(), 3)));
  ijtField.set_variables(3);

  Field xyField =
   cells.add(Field("xy", make_datatype<double>(), make_shape(cells.size(), 2)));
  xyField.set_variables(2);

  Field lonLatField =
    cells.add(Field("lonlat", make_datatype<double>(), make_shape(cells.size(), 2)));
  lonLatField.set_variables(2);

  // Set field views.
  auto cellsGlobalIdx = array::make_view<gidx_t, 1>(cells.global_index());
  auto cellsRemoteIdx = array::make_indexview<idx_t, 1>(cells.remote_index());
  auto cellsPart      = array::make_view<int, 1>(cells.partition());
  auto cellsHalo      = array::make_view<int, 1>(cells.halo());
  auto cellsFlags     = array::make_view<int, 1>(cells.flags());
  auto cellsIjt       = array::make_view<idx_t, 2>(ijtField);
  auto cellsXy        = array::make_view<double, 2>(xyField);
  auto cellsLonLat    = array::make_view<double, 2>(lonLatField);

  // Set local cells.
  auto& nodeConnectivity = cells.node_connectivity();
  const idx_t cellElemIdx0 = cells.elements(0).begin();

  idx_t cellLocalIdx = 0;
  for (const LocalCell& cell : localCells) {


    // Set xy.
    const PointXY xyLocal = jacobian.xy(PointIJ(cell.i + 0.5, cell.j + 0.5), cell.t);
    const PointXY xyGlobal = jacobian.xyLocalToGlobal(xyLocal, cell.t).first;

    cellsXy(cellLocalIdx, XX) = xyLocal.x();
    cellsXy(cellLocalIdx, YY) = xyLocal.y();

    // Set lon-lat.
    const PointLonLat lonLat = csProjection->lonlat(xyGlobal);
    cellsLonLat(cellLocalIdx, LON) = lonLat.lon();
    cellsLonLat(cellLocalIdx, LAT) = lonLat.lat();

    // Get four surroundings nodes.
    const size_t nodeIdx0 = getNodeIdx(cell.t, cell.j    , cell.i    , true);
    const size_t nodeIdx1 = getNodeIdx(cell.t, cell.j    , cell.i + 1, true);
    const size_t nodeIdx2 = getNodeIdx(cell.t, cell.j + 1, cell.i + 1, true);
    const size_t nodeIdx3 = getNodeIdx(cell.t, cell.j + 1, cell.i    , true);
    const GlobalNode& node0 = globalNodes[nodeIdx0];
    const GlobalNode& node1 = globalNodes[nodeIdx1];
    const GlobalNode& node2 = globalNodes[nodeIdx2];
    const GlobalNode& node3 = globalNodes[nodeIdx3];

    const auto quadNodeIdx = std::array<idx_t, 4> {
      node0.localIdx, node1.localIdx, node2.localIdx, node3.localIdx};

    // Set connectivity.
    nodeConnectivity.set(cellLocalIdx + cellElemIdx0, quadNodeIdx.data());

    // Set global index.
    cellsGlobalIdx(cellLocalIdx) = cell.globalIdx;

    // Set cell remote index.
    cellsRemoteIdx(cellLocalIdx) = cell.remoteIdx;

    // Set partition.
    cellsPart(cellLocalIdx) = cell.part;

    // Set ijt.
    cellsIjt(cellLocalIdx, Coordinates::I) = cell.i;
    cellsIjt(cellLocalIdx, Coordinates::J) = cell.j;
    cellsIjt(cellLocalIdx, Coordinates::T) = cell.t;

    // Set flags.
    Topology::reset(cellsFlags(cellLocalIdx));
    switch(cell.type) {

      case CellType::UNDEFINED : {
        ATLAS_ASSERT(CellType::UNDEFINED);
        break;
      }
      case CellType::OWNER : {
        cellsHalo(cellLocalIdx) = 0;
        break;
      }
      case CellType::HALO : {
        cellsHalo(cellLocalIdx) = 1;
        Topology::set(cellsFlags(cellLocalIdx), Topology::GHOST);
        break;
      }
    }

    ++cellLocalIdx;
  }


  // ---------------------------------------------------------------------------
  // 8. FINALISE
  //    Done. That was rather a lot of bookkeeping!
  // ---------------------------------------------------------------------------

  mesh.metadata().set("halo", nHalo);
  mesh.metadata().set("mesh_type", "cubed_sphere");

  mesh.nodes().metadata().set("parallel", true);
  mesh.nodes().metadata().set("nb_owned", nodeLocalIdxCount[thisPart]);

  mesh.cells().metadata().set("parallel", true);
  mesh.cells().metadata().set("nb_owned", cellLocalIdxCount[thisPart]);

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

