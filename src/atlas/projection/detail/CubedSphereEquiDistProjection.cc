/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereEquiDistProjection.h"

#include <cmath>
#include <iostream>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

CubedSphereEquiDistProjection::CubedSphereEquiDistProjection( const eckit::Parametrisation& params )
    : CubedSphereProjectionBase(params) {
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistProjection::lonlat2xy( double crd[] ) const {

    Log::info() << "lonlat2xy start : lonlat = " << crd[LON] << " " << crd[LAT] << std::endl;

    idx_t t;
    double ab[2]; // alpha-beta coordinate
    double xyz[3]; // on Cartesian grid

    CubedSphereProjectionBase::lonlat2xypre(crd, t, xyz);

    //now should be tile 0 - now calculate (alpha, beta) in radians.
    // should be between - pi/4 and pi/4
    ab[0] =   M_PI_4 * xyz[YY] / xyz[XX] ;
    ab[1] = - M_PI_4 * xyz[ZZ] / xyz[XX];

    Log::debug() << "lonlat2xy xyz ab : "
       << xyz[0] << " " << xyz[1]  << " " << xyz[2] << " "
       << ab[0] << " " << ab[1] << std::endl;

    CubedSphereProjectionBase::alphabetatt2xy(t, ab, crd);

    Log::info() << "lonlat2xy end : xy = " << crd[LON] << " " << crd[LAT] << std::endl;

}


// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistProjection::xy2lonlat( double crd[] ) const {

    const double rsq3 = 1.0/sqrt(3.0);
    double xyz[3];
    double ab[2]; // alpha-beta coordinate
    idx_t t;  // tile index

    // calculate xy (in degrees) to alpha beta (in radians) and t - tile index.
    CubedSphereProjectionBase::xy2alphabetat(crd, t, ab);

    Log::info() << "xy2lonlat:: crd t ab  : "  << crd[LON] << " " << crd[1] << " " << t << " " << ab[0] << " " << ab[1] << std::endl;

    xyz[0] = -rsq3;
    xyz[1] = -rsq3 * ab[0] / M_PI_4;
    xyz[2] = -rsq3 * ab[1] / M_PI_4;

    CubedSphereProjectionBase::xy2lonlatpost(xyz, t, crd);

    Log::info() << "end of equidistant xy2lonlat lonlat = " << crd[LON] << " " << crd[LAT] << std::endl;
}


// -------------------------------------------------------------------------------------------------

ProjectionImpl::Jacobian CubedSphereEquiDistProjection::jacobian(const PointLonLat& ) const {
    ATLAS_NOTIMPLEMENTED;
}

// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------

CubedSphereEquiDistProjection::Spec CubedSphereEquiDistProjection::spec() const {
    // Fill projection specification
    Spec proj;
    proj.set( "type", static_type() );
    return proj;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistProjection::hash( eckit::Hash& h ) const {
    // Add to hash
    h.add( static_type() );
    CubedSphereProjectionBase::hash(h);
}

// -------------------------------------------------------------------------------------------------

namespace {
static ProjectionBuilder<CubedSphereEquiDistProjection>
register_1( CubedSphereEquiDistProjection::static_type() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas