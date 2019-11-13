/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Polygon.h"

namespace atlas {
namespace util {

namespace {

//------------------------------------------------------------------------------------------------------

double cross_product_analog( const PointLonLat& A, const PointLonLat& B, const PointLonLat& C ) {
    return ( A.lon() - C.lon() ) * ( B.lat() - C.lat() ) - ( A.lat() - C.lat() ) * ( B.lon() - C.lon() );
}

}  // namespace

//------------------------------------------------------------------------------------------------------

Polygon::Polygon() = default;

Polygon::Polygon( const Polygon::edge_set_t& edges ) {
    ATLAS_TRACE();
    // get external edges by attempting to remove reversed edges, if any
    edge_set_t extEdges;
    for ( const edge_t& e : edges ) {
        if ( !extEdges.erase( e.reverse() ) ) {
            extEdges.insert( e );
        }
    }
    ATLAS_ASSERT( extEdges.size() >= 2 );

    // set one polygon cycle, by picking next edge with first node same as second
    // node of last edge
    clear();
    reserve( extEdges.size() + 1 );

    emplace_back( extEdges.begin()->first );
    for ( edge_set_t::iterator e = extEdges.begin(); e != extEdges.end() && e->first == back();
          e                      = extEdges.lower_bound( edge_t( back(), std::numeric_limits<idx_t>::min() ) ) ) {
        emplace_back( e->second );
        extEdges.erase( *e );
    }
    ATLAS_ASSERT( front() == back() );

    // exhaust remaining edges, should represent additional cycles, if any
    while ( !extEdges.empty() ) {
        operator+=( Polygon( extEdges ) );
    }
}

Polygon::operator bool() const {
    return !Polygon::empty();
}

Polygon& Polygon::operator+=( const Polygon& other ) {
    if ( empty() ) {
        return operator=( other );
    }

    // polygon can have multiple cycles, but must be connected graphs
    // Note: a 'cycle' is handled by repeating the indices, excluding (repeated)
    // last index
    ATLAS_ASSERT( other.front() == other.back() );
    const difference_type N = difference_type( other.size() ) - 1;

    container_t cycle( 2 * size_t( N ) );
    std::copy( other.begin(), other.begin() + N, cycle.begin() );
    std::copy( other.begin(), other.begin() + N, cycle.begin() + N );

    for ( const_iterator c = cycle.begin(); c != cycle.begin() + N; ++c ) {
        iterator here = std::find( begin(), end(), *c );
        if ( here != end() ) {
            insert( here, c, c + N );
            return *this;
        }
    }

    throw_AssertionFailed( "Polygon: could not merge polygons, they are not connected", Here() );
}

void Polygon::print( std::ostream& s ) const {
    char z = '{';
    for ( auto n : static_cast<const container_t&>( *this ) ) {
        s << z << n;
        z = ',';
    }
    s << '}';
}

//------------------------------------------------------------------------------------------------------

PolygonCoordinates::PolygonCoordinates( const Polygon& poly, const atlas::Field& lonlat, bool removeAlignedPoints ) {
    ATLAS_ASSERT( poly.size() > 2 );
    ATLAS_ASSERT( poly.front() == poly.back() );

    // Point coordinates
    // - use a bounding box to quickly discard points,
    // - except when that is above/below bounding box but poles should be included

    coordinates_.clear();
    coordinates_.reserve( poly.size() );

    auto ll         = array::make_view<double, 2>( lonlat );
    coordinatesMin_ = PointLonLat( ll( poly[0], LON ), ll( poly[0], LAT ) );
    coordinatesMax_ = coordinatesMin_;

    size_t nb_removed_points_due_to_alignment = 0;

    for ( size_t i = 0; i < poly.size(); ++i ) {
        PointLonLat A( ll( poly[i], LON ), ll( poly[i], LAT ) );
        coordinatesMin_ = PointLonLat::componentsMin( coordinatesMin_, A );
        coordinatesMax_ = PointLonLat::componentsMax( coordinatesMax_, A );

        // if new point is aligned with existing edge (cross product ~= 0) make the
        // edge longer
        if ( ( coordinates_.size() >= 2 ) && removeAlignedPoints ) {
            const PointLonLat& B = coordinates_.back();
            const PointLonLat& C = coordinates_[coordinates_.size() - 2];
            if ( eckit::types::is_approximately_equal( 0., cross_product_analog( A, B, C ) ) ) {
                coordinates_.back() = A;
                ++nb_removed_points_due_to_alignment;
                continue;
            }
        }

        coordinates_.emplace_back( A );
    }

    ATLAS_ASSERT( coordinates_.size() == poly.size() - nb_removed_points_due_to_alignment );
}

PolygonCoordinates::PolygonCoordinates( const std::vector<PointLonLat>& points ) : coordinates_( points ) {
    ATLAS_ASSERT( coordinates_.size() > 2 );
    ATLAS_ASSERT( eckit::geometry::points_equal( coordinates_.front(), coordinates_.back() ) );

    coordinatesMin_ = coordinates_.front();
    coordinatesMax_ = coordinatesMin_;
    for ( const PointLonLat& P : coordinates_ ) {
        coordinatesMin_ = PointLonLat::componentsMin( coordinatesMin_, P );
        coordinatesMax_ = PointLonLat::componentsMax( coordinatesMax_, P );
    }
}

namespace {
PointLonLat compute_centroid(const std::vector<PointLonLat>& points) {

    ATLAS_ASSERT( eckit::geometry::points_equal( points.front(), points.back() ) );

    PointLonLat centroid = {0, 0};
    double signed_area = 0.;
    double a = 0.;  // Partial signed area

    for (size_t i=0; i<points.size()-1; ++i) {
        const PointLonLat& p0 = points[i];
        const PointLonLat& p1 = points[i+1];
        a = p0[0]*p1[1] - p1[0]*p0[1];
        signed_area += a;
        centroid[0] += (p0[0] + p1[0])*a;
        centroid[1] += (p0[1] + p1[1])*a;
    }
    signed_area *= 0.5;
    centroid[0] /= (6.*signed_area);
    centroid[1] /= (6.*signed_area);

    return centroid;
}

double compute_inner_radius_squared(const std::vector<PointLonLat>& points, const PointLonLat& centroid) {
  auto distance2 = [](const PointLonLat& a, const PointLonLat& b) {
    double dx = (a[0]-b[0]);
    double dy = (a[1]-b[1]);
    return dx*dx + dy*dy;
  };
  auto dot = [](const PointLonLat& a, const PointLonLat& b) {
    return a[0]*b[0]+a[1]*b[1];
  };
  double R2 = std::numeric_limits<double>::max();
  PointLonLat projection;
  for( size_t i=0; i<points.size()-1; ++i ) {
    double d2 = distance2(points[i],points[i+1]);
    double t = std::max(0.,std::min(1.,dot(centroid-points[i],points[i+1]-points[i])/d2));
    projection[0] = points[i][0] + t * (points[i+1][0]-points[i][0]);
    projection[1] = points[i][1] + t * (points[i+1][1]-points[i][1]);
    R2 = std::min( R2, distance2(projection,centroid) );
    //Log::info() << "Segment " << points[i] << " - " << points[i+1] << " :  projection = " << projection << "   \t distance = "  << std::sqrt(distance2(projection,centroid) ) << std::endl;
  }
  return R2;
}
} // namespace

PolygonCoordinates::PolygonCoordinates( const std::vector<PointLonLat>& points,  bool removeAlignedPoints ) {
    coordinates_.clear();
    coordinates_.reserve( points.size() );

    coordinatesMin_ = PointLonLat( std::numeric_limits<double>::max(), std::numeric_limits<double>::max() );
    coordinatesMax_ = PointLonLat( std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest() );

    size_t nb_removed_points_due_to_alignment = 0;

    for ( size_t i = 0; i < points.size(); ++i ) {
        const PointLonLat& A = points[i];
        coordinatesMin_ = PointLonLat::componentsMin( coordinatesMin_, A );
        coordinatesMax_ = PointLonLat::componentsMax( coordinatesMax_, A );

        // if new point is aligned with existing edge (cross product ~= 0) make the
        // edge longer
        if ( ( coordinates_.size() >= 2 ) && removeAlignedPoints ) {
            const PointLonLat& B = coordinates_.back();
            const PointLonLat& C = coordinates_[coordinates_.size() - 2];
            if ( eckit::types::is_approximately_equal( 0., cross_product_analog( A, B, C ) ) ) {
                coordinates_.back() = A;
                ++nb_removed_points_due_to_alignment;
                continue;
            }
        }
        coordinates_.emplace_back( A );
    }

    ATLAS_ASSERT( coordinates_.size() > 2 );
    ATLAS_ASSERT( eckit::geometry::points_equal( coordinates_.front(), coordinates_.back() ) );

    centroid_ = compute_centroid( coordinates_ );
    inner_radius_squared_ = compute_inner_radius_squared( coordinates_, centroid_ );
}

PolygonCoordinates::~PolygonCoordinates() = default;

const PointLonLat& PolygonCoordinates::coordinatesMax() const {
    return coordinatesMax_;
}

const PointLonLat& PolygonCoordinates::coordinatesMin() const {
    return coordinatesMin_;
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas