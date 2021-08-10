/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"

namespace atlas {
namespace detail {
namespace cubedsphere {

CubedSphereJacobian::CubedSphereJacobian(const CubedSphereGrid& csGrid) {

  // Get projection.
  csProjection_ = dynamic_cast<const CubedSphereProjectionBase*>
    (csGrid.projection().get());

  // Get number of tiles.
  const auto nTiles = idx2st(csGrid.GetNTiles());

  // Get xy of points (i = 0, j = 0), (i = 1, j = 0) and (i = 0, j = 0) on tiles.
  auto xy00 = std::vector<PointXY>(nTiles);   // (i = 0, j = 0)
  auto xy10 = std::vector<PointXY>(nTiles);   // (i = 1, j = 0)
  auto xy01 = std::vector<PointXY>(nTiles);   // (i = 0, j = 1)

  // Loop over grid points.
  auto tijIt = csGrid.tij().begin();
  for (const auto& xy : csGrid.xy()) {

    const auto t = idx2st((*tijIt).t());
    const auto i = (*tijIt).i();
    const auto j = (*tijIt).j();

    if (i == 0 && j == 0) xy00[t] = xy;
    if (i == 1 && j == 0) xy10[t] = xy;
    if (i == 0 && j == 1) xy01[t] = xy;

    ++tijIt;
  }

  // Calculate Jacobians.
  size_t t = 0;
  std::generate_n(std::back_inserter(jacobians_), nTiles, [&](){

      // Initialise element.
      auto jacElem = Jacobian{};

      // Calculate Jacobian.
      jacElem.dxyByDi = xy10[t] - xy00[t];
      jacElem.dxyByDj = xy01[t] - xy00[t];

      // Calculate inverse Jacobian.
      const auto invDet = 1./(jacElem.dxyByDi.x() * jacElem.dxyByDj.y()
        - jacElem.dxyByDj.x() * jacElem.dxyByDi.y());

      jacElem.dijByDx.i() =  jacElem.dxyByDj.y() * invDet;
      jacElem.dijByDy.i() = -jacElem.dxyByDj.x() * invDet;
      jacElem.dijByDx.j() = -jacElem.dxyByDi.y() * invDet;
      jacElem.dijByDy.j() =  jacElem.dxyByDi.x() * invDet;

      // Extrapolate cell(t, 0, 0) xy to get node(t, 0, 0) xy.
      jacElem.xy00 = PointXY{
        xy00[t] - jacElem.dxyByDi * 0.5  - jacElem.dxyByDj * 0.5};

      ++t;
      return jacElem;

    });
}

PointTXY CubedSphereJacobian::tijToTxy(const PointTIJ& tijLocal) const {

  // Get jacobian.
  const auto& jac = jacobians_[idx2st(tijLocal.t())];

  return PointTXY(tijLocal.t(),
    jac.xy00 + jac.dxyByDi * tijLocal.i() + jac.dxyByDj * tijLocal.j());
}

PointTIJ CubedSphereJacobian::txyToTij(const PointTXY& txyLocal) const {

  // Get jacobian.
  const auto& jac = jacobians_[idx2st(txyLocal.t())];

  // Set xy displacement
  PointXY dxy = txyLocal.xy() - jac.xy00;

  return PointTIJ(txyLocal.t(), jac.dijByDx * dxy.x() + jac.dijByDy * dxy.y());
}

PointTXY CubedSphereJacobian::txyLocalToGlobal(const PointTXY &txyLocal) const {

  // Get tiles object.
  const auto& csTiles = csProjection_->getCubedSphereTiles();

  // Get global xy
  const auto xyGlobal =
    csTiles.tileCubePeriodicity(txyLocal.xy(), txyLocal.t());

  // Get global t
  const auto tGlobal = csTiles.tileFromXY(xyGlobal.data());

  return PointTXY(tGlobal, xyGlobal);
}

PointTIJ CubedSphereJacobian::tijLocalToGlobal(const PointTIJ &tijLocal) const {

  // Convert ij to xy.
  const auto txyLocal = this->tijToTxy(tijLocal);

  // Convert local xy to global xy.
  const auto txyGlobal = this->txyLocalToGlobal(txyLocal);

  // Return global ij.
  return this->txyToTij(txyGlobal);
}

}
}
}
