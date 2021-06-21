/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <array>
#include <string>
#include <iostream>

#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
class Hash;
}  // namespace eckit

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
class PointXY;
class PointLonLat;

namespace util {
class Config;
}  // namespace util

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace cubedspheretiles {
class CubedSphereTiles;
class FV3CubedSphereTiles;
class LFRicCubedSphereTiles;
}  // namespace cubespheretiles
#endif

//---------------------------------------------------------------------------------------------------------------------

class CubedSphereTiles : DOXYGEN_HIDE(
    public util::ObjectHandle<atlas::cubedspheretiles::CubedSphereTiles> ) {
public:
    using Spec = util::Config;

public:
    using Handle::Handle;
    CubedSphereTiles() = default;
    CubedSphereTiles( const eckit::Parametrisation& );

    /// Type of the cubed-sphere tiles:
    std::string type() const;

    void tile0Rotate( double xyz[] ) const;

    void tile1Rotate( double xyz[] ) const;

    void tile2Rotate( double xyz[] ) const;

    void tile3Rotate( double xyz[] ) const;

    void tile4Rotate( double xyz[] ) const;

    void tile5Rotate( double xyz[] ) const;

    void tile0RotateInverse( double xyz[] ) const;

    void tile1RotateInverse( double xyz[] ) const;

    void tile2RotateInverse( double xyz[] ) const;

    void tile3RotateInverse( double xyz[] ) const;

    void tile4RotateInverse( double xyz[] ) const;

    void tile5RotateInverse( double xyz[] ) const;

    idx_t tileFromXY( const double xy[] ) const;

    idx_t tileFromLonLat( const double lonlat[] ) const;

    void enforceXYdomain( double xy[] ) const;

private:
    /// Output to stream
    void print( std::ostream& ) const;

    friend std::ostream& operator<<( std::ostream& s, const CubedSphereTiles& cst );

};


class FV3CubedSphereTiles : public CubedSphereTiles {
public:
    using CubedSphereTiles::CubedSphereTiles;
    FV3CubedSphereTiles() : CubedSphereTiles() {
      std::cout << "grid/Tiles.h FV3CubedSphereTiles"  << std::endl;
    }

    FV3CubedSphereTiles( const CubedSphereTiles& );

private:
    const ::atlas::cubedspheretiles::FV3CubedSphereTiles* cubedspheretiles_;
};

class LFRicCubedSphereTiles : public CubedSphereTiles {
public:
    using CubedSphereTiles::CubedSphereTiles;
    LFRicCubedSphereTiles() : CubedSphereTiles() {
      std::cout << "grid/Tiles.h LFRicCubedSphereTiles"  << std::endl;
    }

    LFRicCubedSphereTiles( const CubedSphereTiles& );

private:
    const ::atlas::cubedspheretiles::LFRicCubedSphereTiles* cubedspheretiles_;
};


}  // namespace atlas