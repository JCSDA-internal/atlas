/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/domain.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Point.h"



namespace atlas {
namespace projection {
namespace detail {


template <typename Rotation>
class StretchLAM final : public ProjectionImpl {
public:
    using Spec = ProjectionImpl::Spec;
    atlas::PointXY point_reg;

    // constructor uses parametrisation and point to stretch
    StretchLAM(const eckit::Parametrisation& );
    // projection name
    static std::string static_type() { return Rotation::typePrefix() + "stretch"; }
    std::string type() const override { return static_type(); }

    // projection and inverse projection, I will use the POINTXY in the function
    void xy2lonlat( double crd[] ) const override;
    void lonlat2xy( double crd[] ) const override;


    // specification for stretching
    Spec spec() const override;

    std::string units() const override { return "meters"; }

    void hash( eckit::Hash& ) const override;

    // NOT for the moment
    Jacobian jacobian( const PointLonLat& ) const override;

    bool strictlyRegional() const override { return true; }  // Stretch projection cannot be used for global grids
    RectangularLonLatDomain lonlatBoundingBox( const Domain& domain ) const override {
        return ProjectionImpl::lonlatBoundingBox( domain );
    }

    //void stretch_LAM( eckit::Parametrisation&, PointXY );
    //atlas::PointXY stretch_LAM() const;
    void checkvalue(const double&, const double&) const;
    double general_stretch (double&, const bool&, int&, const int&) const;
    //atlas::PointXY stretch_LAM_gen(atlas::PointXY&);
    //void stretch_LAM_gen(double crd[]);



protected:
    //Normalise normalise_;
    //util::NormaliseLongitude normalise_stretch_;

    // INPUT used in stretch
    //util::Config Domainconfig;
    //util::Config reg_grid_config;
    //util::Config stretch_config;


    double delta_low_; //resolution of the host model
    double delta_high_; //resolution of the regional model (regular grid)
    double var_ratio_; //power used for the stretching
    double x_reg_start_; //xstart of the internal regional grid
    double y_reg_start_; //ystart of the internal regional grid
    double x_reg_end_; //xend of the internal regional grid
    double y_reg_end_; ////yend of the internal regional grid
    double startx_; //domain startx
    double endx_; //domain endx
    double starty_; //domain starty
    double endy_; //domain endy


    void setup( const eckit::Parametrisation& p );

private:
    //PointXY pointxy_;
    Rotation rotation_;
};

using StretchProjection        = StretchLAM<NotRotated>;
using RotatedStretchProjection = StretchLAM<Rotated>;

}  // namespace detail
}  // namespace projection
}  // namespace atlas
