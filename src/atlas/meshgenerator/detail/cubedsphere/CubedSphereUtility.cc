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
namespace meshgenerator {
namespace detail {
namespace cubedsphere {

// -----------------------------------------------------------------------------
// Jacobian2 class
// -----------------------------------------------------------------------------

Jacobian2::Jacobian2(
  double df0_by_dx0, double df0_by_dx1, double df1_by_dx0, double df1_by_dx1) {
  data_ << df0_by_dx0, df0_by_dx1, df1_by_dx0, df1_by_dx1;
}

Jacobian2::Jacobian2(const Point2& f00, const Point2& f10, const Point2& f01,
  double dx0, double dx1) : Jacobian2(
    (f10[0] - f00[0]) / dx0, (f01[0] - f00[0]) / dx1,
    (f10[1] - f00[1]) / dx0, (f01[1] - f00[1]) / dx1) {}

Jacobian2::Jacobian2(const Point2& f00, const Point2& f10, const Point2& f01) :
  Jacobian2(f00, f10, f01, 1., 1.) {}

Jacobian2::Jacobian2(const Eigen::Matrix2d& data) : data_(data) {}

Point2 Jacobian2::operator*(const Point2& dx) const {

  // Declare result.
  Point2 df;

  // Compute product using Eigen maps.
  Eigen::Map<Eigen::Vector2d>(df.data()) =
    data_ * Eigen::Map<const Eigen::Vector2d>(dx.data());

  return df;
}

Jacobian2 Jacobian2::operator*(const Jacobian2& Jb) const {
  return Jacobian2(data_ * Jb.data_);
}

Jacobian2 Jacobian2::inverse() const {
  return Jacobian2(data_.inverse());
}


// -----------------------------------------------------------------------------
// NeighbourJacobian class
// -----------------------------------------------------------------------------

NeighbourJacobian::NeighbourJacobian(const CubedSphereGrid& csGrid) {

  // Get projection.
  csProjection_ =
    dynamic_cast<const CubedSphereProjectionBase*>(csGrid.projection().get());
  ATLAS_ASSERT(csProjection_);

  // Get tiles.
  const auto& csTiles = csProjection_->getCubedSphereTiles();

  // Get grid size.
  N_ = csGrid.N();

  // Get xy of points (i = 0, j = 0), (i = 1, j = 0) and (i = 0, j = 0) on tiles.
  std::array<PointXY, 6> xy00;
  std::array<PointXY, 6> xy10;
  std::array<PointXY, 6> xy01;

  // Loop over grid points.
  auto tijIt = csGrid.tij().begin();
  for (const auto& xy : csGrid.xy()) {

    const auto t = idx2st((*tijIt).t());
    const auto i = (*tijIt).i();
    const auto j = (*tijIt).j();

    if      (i == 0 and j == 0) xy00[t] = xy;
    else if (i == 1 and j == 0) xy10[t] = xy;
    else if (i == 0 and j == 1) xy01[t] = xy;

    ++tijIt;
  }

  for (size_t t = 0; t < 6; ++t) {

    // Calculate tile Jacobians.
    dxy_by_dij_[t] = Jacobian2(xy00[t], xy10[t], xy01[t]);
    dij_by_dxy_[t]  = dxy_by_dij_[t].inverse();

    // Set xy00. Grid point needs moving to (i = 0, j = 0).
    xy00_[t] = xy00[t] + dxy_by_dij_[t] * PointIJ(-0.5, -0.5);

    // Neighbour assignment lambda.
    const auto neighbourAssignment = [&](Locations::k k){

      // Shift points in to neighbouring tiles.
      PointIJ ijDisplacement;
      switch (k) {
        case Locations::LEFT : {
          ijDisplacement = PointIJ(-2, 0);
          break;
          }
        case Locations::TOP : {
          ijDisplacement = PointIJ(0, N_);
          break;
        }
        case Locations::RIGHT : {
          ijDisplacement = PointIJ(N_, 0);
          break;
        }
        case Locations::BOTTOM : {
          ijDisplacement = PointIJ(0, -2);
          break;
        }
      }

      // Convert displacement from ij to xy.
      const PointXY xyDisplacement = dxy_by_dij_[t] * ijDisplacement;

      // Get neighbour xy points in xy space local to tile.
      const auto xy00Local = xy00[t] + xyDisplacement;
      const auto xy10Local = xy10[t] + xyDisplacement;
      const auto xy01Local = xy01[t] + xyDisplacement;

      // Convert from local xy to global xy.
      const auto xy00Global = csTiles.tileCubePeriodicity(xy00Local, st2idx(t));
      const auto xy10Global = csTiles.tileCubePeriodicity(xy10Local, st2idx(t));
      const auto xy01Global = csTiles.tileCubePeriodicity(xy01Local, st2idx(t));

      // Get neighbour tile ID.
      neighbours_[t].t_[k] = csTiles.tileFromXY(xy00Global.data());

      // Set Jacobian of global xy with respect to local ij.
      const auto dxyGlobal_by_dij = Jacobian2(
        xy00Global, xy10Global, xy01Global);

      // Chain rule to get Jacobian with respect to local xy.
      neighbours_[t].dxyGlobal_by_dxyLocal_[k] =
        dxyGlobal_by_dij * dij_by_dxy_[t];

      // Set local xy00
      neighbours_[t].xy00Local_[k] = xy00Local;

      // Set global xy00
      neighbours_[t].xy00Global_[k] = xy00Global;

    };

    // Assign neighbours (good job we put it all in a lambda!).
    neighbourAssignment(Locations::LEFT);
    neighbourAssignment(Locations::TOP);
    neighbourAssignment(Locations::RIGHT);
    neighbourAssignment(Locations::BOTTOM);

  }

}

PointXY NeighbourJacobian::xy(const PointIJ& ij, idx_t t) const {

  // Get jacobian.
  const auto& jac = dxy_by_dij_[idx2st(t)];
  const auto& xy00 = xy00_[idx2st(t)];

  // Return ij
  return xy00 + jac * ij;
}

PointXY NeighbourJacobian::xy(const std::pair<PointIJ, idx_t>& ijt) const {
  return xy(ijt.first, ijt.second);
}

PointIJ NeighbourJacobian::ij(const PointXY& xy, idx_t t) const {

  // Get jacobian.
  const auto& jac = dij_by_dxy_[idx2st(t)];
  const auto& xy00 = xy00_[idx2st(t)];

  // Return ij
  return jac * (xy - xy00);
}

PointIJ NeighbourJacobian::ij(const std::pair<PointXY, idx_t>& xyt) const {
  return ij(xyt.first, xyt.second);
}

std::pair<PointXY, idx_t> NeighbourJacobian::xyLocalToGlobal(
  const PointXY& xyLocal, idx_t tLocal) const {

  // The tileCubePeriodicity method fails when extrapolating along an unowned
  // tile edge. This method explicitly places an xy point on to a neighbouring
  // tile to avoid this. Once the correct xy position has been found,
  // tileCubePeriodicity will correcty find the "owned" xy position of a point
  // on an unowned tile edge.

  // Declare result.
  PointXY xyGlobal;
  idx_t tGlobal;

  // Get ij.
  const auto ijLocal = ij(xyLocal, tLocal);

  // Exclude invalid halo corners.
  ATLAS_ASSERT(ijCross(ijLocal));

  // Get tiles.
  const auto& csTiles = csProjection_->getCubedSphereTiles();

  if (ijInterior(ijLocal)) {
    // That was easy.
    xyGlobal = xyLocal;
    tGlobal = tLocal;
  }
  else {
    // Figure out which tile xy is on.
    Locations::k k;
    if      (ijLocal.i() < 0 ) k = Locations::LEFT;
    else if (ijLocal.j() > N_) k = Locations::TOP;
    else if (ijLocal.i() > N_) k = Locations::RIGHT;
    else if (ijLocal.j() < 0 ) k = Locations::BOTTOM;

    // Get reference points and jacobian.
    const auto& xy00Local_ = neighbours_[idx2st(tLocal)].xy00Local_[k];
    const auto& xy00Global_ = neighbours_[idx2st(tLocal)].xy00Global_[k];
    const auto& jac = neighbours_[idx2st(tLocal)].dxyGlobal_by_dxyLocal_[k];

    // Get t.
    tGlobal =  neighbours_[idx2st(tLocal)].t_[k];

    // Calculate global xy.
    xyGlobal = xy00Global_ + jac * (xyLocal - xy00Local_);
  }

  // Correct for edge-ownership rules.
  xyGlobal = csTiles.tileCubePeriodicity(xyGlobal, tGlobal);
  tGlobal = csTiles.tileFromXY(xyGlobal.data());

  return std::make_pair(xyGlobal, tGlobal);
}

std::pair<PointXY, idx_t> NeighbourJacobian::xyLocalToGlobal(
  const std::pair<PointXY, idx_t>& xytLocal) const {
  return xyLocalToGlobal(xytLocal.first, xytLocal.second);
}

std::pair<PointIJ, idx_t> NeighbourJacobian::ijLocalToGlobal(
  const PointIJ &ijLocal, idx_t tLocal) const {

  // Use xyLocalToGlobal method to take care of this.

  // Get global xyt.
  auto xytGlobal = xyLocalToGlobal(xy(ijLocal, tLocal), tLocal);

  // convert to ijt
  return std::make_pair(ij(xytGlobal), xytGlobal.second);
}

std::pair<PointIJ, idx_t> NeighbourJacobian::ijLocalToGlobal(
  const std::pair<PointIJ, idx_t>& ijtLocal) const {
  return ijLocalToGlobal(ijtLocal.first, ijtLocal.second);
}

bool NeighbourJacobian::ijInterior(const PointIJ& ij) const {
  return ij.i() >= 0 and ij.i() <= N_ and
         ij.j() >= 0 and ij.j() <= N_;
}

bool NeighbourJacobian::ijCross(const PointIJ& ij) const {

  const bool inCorner = (ij.i() < 0  and ij.j() < 0 ) or // bottom-left corner.
                        (ij.i() > N_ and ij.j() < 0 ) or // bottom-right corner.
                        (ij.i() > N_ and ij.j() > N_) or // top-right corner.
                        (ij.i() < 0  and ij.j() > N_);   // top-left corner.
  return !inCorner;
}

}
}
}
}
