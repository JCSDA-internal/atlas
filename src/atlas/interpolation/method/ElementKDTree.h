/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/util/KDTree.h"

namespace atlas {
class Field;
class Mesh;
}  // namespace atlas

namespace atlas {
namespace interpolation {
namespace method {

//----------------------------------------------------------------------------------------------------------------------

/// @brief Build a k-d tree of mesh cell centres, with the cell index as payload.
///
/// Cells flagged as INVALID (util::Topology) are skipped.
/// The tree uses the default util::KDTree backend.
///
/// @param mesh           Mesh whose cells are indexed
/// @param field_centres  Cell-centre field of shape (cells, 3) with 3D cartesian coordinates (x,y,z)
util::IndexKDTree3D build_element_kdtree(const Mesh& mesh, const Field& field_centres);

//----------------------------------------------------------------------------------------------------------------------

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
