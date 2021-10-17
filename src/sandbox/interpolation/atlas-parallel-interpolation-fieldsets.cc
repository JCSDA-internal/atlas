/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <string>

#include "eckit/filesystem/LocalPathName.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "eckit/config/Resource.h"
#include "eckit/config/Configuration.h"

using namespace atlas;

auto vortex_rollup = []( double lon, double lat, double t ) {
    // lon and lat in radians!

    // Formula found in "A Lagrangian Particle Method with Remeshing for Tracer Transport on the Sphere"
    // by Peter Bosler, James Kent, Robert Krasny, CHristiane Jablonowski, JCP 2015

    auto sqr           = []( const double x ) { return x * x; };
    auto sech          = []( const double x ) { return 1. / std::cosh( x ); };
    const double T     = 1.;
    const double Omega = 2. * M_PI / T;
    t *= T;
    const double lambda_prime = std::atan2( -std::cos( lon - Omega * t ), std::tan( lat ) );
    const double rho          = 3. * std::sqrt( 1. - sqr( std::cos( lat ) ) * sqr( std::sin( lon - Omega * t ) ) );
    double omega              = 0.;
    double a                  = util::Earth::radius();
    if ( rho != 0. ) {
        omega = 0.5 * 3 * std::sqrt( 3 ) * a * Omega * sqr( sech( rho ) ) * std::tanh( rho ) / rho;
    }
    double q = 1. - std::tanh( 0.2 * rho * std::sin( lambda_prime - omega / a * t ) );
    return q;
};


class FieldConfiguration  {
public:

    FieldConfiguration(eckit::Configuration & conf) ;

    const std::string & getSourceGridName(){ return source_grid_name_;};
    const std::string & getSourceParitionerName(){ return source_partitioner_name_;};
    const std::size_t & getSourceLevels(){ return source_levels_;};
    const std::string & getTargetGridName(){ return target_grid_name_;};
    const std::string & getTargetPartitionerName(){ return target_partitioner_name_;};

private:
    std::string source_grid_name_;
    std::string source_partitioner_name_;
    std::size_t source_levels_;
    std::string target_grid_name_;
    std::string target_partitioner_name_;
};


FieldConfiguration::FieldConfiguration(eckit::Configuration & conf) :
    source_grid_name_(conf.getString("source grid name")),
    source_partitioner_name_(conf.getString("source partitioner name")),
    source_levels_(conf.getUnsigned("source levels")),
    target_grid_name_(conf.getString("target grid name")),
    target_partitioner_name_(conf.getString("target partitioner name")) {};


class FieldSetConfiguration  {
public:
    FieldSetConfiguration(eckit::Configuration & conf);

    const bool & getIncludeAdjoint(){ return include_adjoint_;};
    const bool & getIncludeRedistribution(){ return include_redistribution_;};
    const std::string & getInterpolationMethod(){ return interpolation_method_;};
    std::size_t getNoFields(){ return field_configurations_.size(); };
    eckit::LocalConfiguration & getFieldConfiguration(const std::size_t & i){return field_configurations_[i];};

private:
    bool include_adjoint_;
    bool include_redistribution_;
    std::string interpolation_method_;
    std::vector<eckit::LocalConfiguration> field_configurations_;
};

FieldSetConfiguration::FieldSetConfiguration(eckit::Configuration & conf) :
    include_adjoint_(conf.getBool("include adjoint")),
    include_redistribution_(conf.getBool("include redistribution")),
    interpolation_method_(conf.getString("interpolation method")),
    field_configurations_(conf.getSubConfigurations("field configurations")) {};

class AtlasParallelInterpolationFieldSet : public AtlasTool {
    int execute( const AtlasTool::Args& args ) override;

  //  int numberOfPositionalArguments() override { return -1; }
  //  int minimumPositionalArguments() override { return 0; }
public:
    AtlasParallelInterpolationFieldSet( int argc, char* argv[] ) : AtlasTool( argc, argv ) {
        add_option( new SimpleOption<std::string>( "yaml filename", "file name of yaml" ) );
    }
};

int AtlasParallelInterpolationFieldSet::execute(const AtlasTool::Args &args) {
    ATLAS_TRACE( "AtlasParallelInterpolationFieldSet::execute" );
    eckit::LocalPathName path(args.getString("yaml filename", "interp.yaml" ));
    atlas::util::Config config(path);

    FieldSetConfiguration fs(config);

    std::cout << "Config setting" << fs.getIncludeAdjoint() << " "
              << fs.getIncludeRedistribution() << " "
              << fs.getInterpolationMethod() << std::endl;

    for (std::size_t i = 0; i < fs.getNoFields(); ++i ) {
        auto f = FieldConfiguration(fs.getFieldConfiguration(i));
        std::cout << i << " " << f.getSourceGridName() << " " << f.getTargetGridName() << " "
                  << f.getSourceParitionerName() << " " << f.getTargetPartitionerName() << " "
                  << f.getSourceLevels() << std::endl;
     }

    // get main program to read the configuration and print the values

    // get main program to allocate input and output FieldSets.

    // populate the input FieldSet with variants of the vortex roll up.

    // get main progran to interface to wrapper.

    // gmsh output from interpolation.


    // apply adjoint interpolation

    // test?
    return 1;
}

int main( int argc, char* argv[] ) {
    AtlasParallelInterpolationFieldSet tool( argc, argv );
    return tool.start();
}
