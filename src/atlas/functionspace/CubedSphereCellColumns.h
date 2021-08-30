/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/detail/CubedSphereUtility.h"
#include "atlas/mesh.h"
#include "atlas/mesh/HybridElements.h"

namespace atlas {
namespace functionspace {

class CubedSphereCellColumns : public functionspace::CellColumns,
  public detail::CubedSphereUtility {

  using CellColumns::CellColumns;

public:
  inline CubedSphereCellColumns(const FunctionSpace& functionSpace) :
    CellColumns(functionSpace), CubedSphereUtility(
      this->mesh().metadata().get<idx_t>("halo"),
      this->mesh().cells().field("ijt"),
      this->mesh().cells().halo()) {}

  inline CubedSphereCellColumns(const Mesh& mesh,
    const eckit::Configuration& configuration) :
    CellColumns(mesh, configuration), CubedSphereUtility(
      this->mesh().metadata().get<idx_t>("halo"),
      this->mesh().cells().field("ijt"),
      this->mesh().cells().halo()) {}

  inline CubedSphereCellColumns(const Mesh& mesh) :
    CellColumns(mesh), CubedSphereUtility(
      this->mesh().metadata().get<idx_t>("halo"),
      this->mesh().cells().field("ijt"),
      this->mesh().cells().halo()) {}

};

} // namespace functionspace
} // namespace atlas
