/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/library/config.h"
#include "atlas/util/Object.h"
#include "atlas/util/Point.h"

namespace atlas {
class CubedSphereGrid;
}

namespace atlas {
namespace projection {
namespace detail {
class CubedSphereProjectionBase;
}
}
}

namespace atlas {
namespace detail {
namespace cubedsphere {

using namespace projection::detail;

// Shorthand static casts.
inline idx_t st2idx(size_t i) {return static_cast<idx_t>(i);}
inline size_t idx2st(idx_t i) {return static_cast<size_t>(i);}

class PointIJ : public Point2 {
public:

  using Point2::Point2;

  inline PointIJ(idx_t i, idx_t j) :
    Point2(static_cast<double>(i), static_cast<double>(j)) {}

  inline double i() const {return x_[0];}
  inline double j() const {return x_[1];}
  inline double& i() {return x_[0];}
  inline double& j() {return x_[1];}

  inline idx_t iRound() const {return static_cast<idx_t>(std::round(i()));}
  inline idx_t jRound() const {return static_cast<idx_t>(std::round(j()));}
  inline idx_t iFloor() const {return static_cast<idx_t>(std::floor(i()));}
  inline idx_t jFloor() const {return static_cast<idx_t>(std::floor(j()));}

};

class PointTIJ : public PointIJ {
public:

  inline PointTIJ(idx_t t, double i, double j) : PointIJ(i, j) {t_ = t;}
  inline PointTIJ(idx_t t, const PointIJ& ij) : PointIJ(ij) {t_ = t;}
  inline PointTIJ(idx_t t, idx_t i, idx_t j) : PointIJ(i, j) {t_ = t;}

  inline idx_t t() const {return t_;}
  inline idx_t& t() {return t_;}

  inline PointIJ ij() const {return PointIJ(*this);}

private:
  idx_t t_{};

};

class PointTXY : public PointXY {
public:

  inline PointTXY(idx_t t, double x, double y) : PointXY(x, y) {t_ = t;}
  inline PointTXY(idx_t t, const PointXY& xy) : PointXY(xy) {t_ = t;}

  inline idx_t t() const {return t_;}
  inline idx_t& t() {return t_;}

  inline PointXY xy() const {return PointXY(*this);}

private:
  idx_t t_{};

};

class CubedSphereJacobian {
public:

  // Constructor.
  CubedSphereJacobian(const CubedSphereGrid& csGrid);

  // linear coordinate transform from ij to xy (given t).
  PointTXY tijToTxy(const PointTIJ& tijLocal) const;

  // linear coordinate transform from xy to ij (given t).
  PointTIJ txyToTij(const PointTXY& txyLocal) const;

  // Convert local xy on tile t coordinate to global coordinate.
  PointTXY txyLocalToGlobal(const PointTXY& txyLocal) const;

  // Convert local ij to tile t to global coordinate (i.e., unique indices).
  PointTIJ tijLocalToGlobal(const PointTIJ& tijLocal) const;

private:

  // Pointer to grid projection.
  const CubedSphereProjectionBase* csProjection_;

  // Coordinate partial derivatives
  struct Jacobian {
    PointXY dxyByDi;  // (dx/di, dy/di)
    PointXY dxyByDj;  // (dx/dj, dy/dj)
    PointIJ dijByDx;  // (di/dx, dj/dx)
    PointIJ dijByDy;  // (di/dy, dj/dy)
    PointXY xy00{};   // xy coordinate of node(i = 0, j = 0).
  };
  std::vector<Jacobian> jacobians_{};


};

}
}
}
