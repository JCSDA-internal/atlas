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
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"

#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/grid/detail/partitioner/CheckerboardPartitioner.h"

#include "atlas/grid.h"

#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"

#include "atlas/interpolation.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "eckit/config/Resource.h"
#include "eckit/config/Configuration.h"

#include "InterpRedistrWrapper.h"

using namespace atlas;


double testFunction( double lon, double lat ) {
    return std::sin( 3 * lon * M_PI / 180 ) * std::sin( 2 * lat * M_PI / 180 );
}


class FieldConfiguration  {
public:

    FieldConfiguration(eckit::Configuration & conf) ;

    const std::string & getSourceGridName(){ return source_grid_name_;};
    const std::string & getSourcePartitionerName(){ return source_partitioner_name_;};
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
public:
    AtlasParallelInterpolationFieldSet( int argc, char* argv[] ) : AtlasTool( argc, argv ) {
        add_option( new SimpleOption<std::string>( "yaml filename", "file name of yaml" ) );
    }
};

int AtlasParallelInterpolationFieldSet::execute(const AtlasTool::Args &args) {
    ATLAS_TRACE( "AtlasParallelInterpolationFieldSet::execute" );

    // get main program to read the configuration and print the values
    eckit::LocalPathName path(args.getString("yaml filename", "interp.yaml" ));
    atlas::util::Config config(path);

    FieldSetConfiguration fs(config);

    std::cout << "Config setting" << fs.getIncludeAdjoint() << " "
              << fs.getIncludeRedistribution() << " "
              << fs.getInterpolationMethod() << std::endl;

    for (std::size_t i = 0; i < fs.getNoFields(); ++i ) {
        auto f = FieldConfiguration(fs.getFieldConfiguration(i));
        std::cout << i << " " << f.getSourceGridName() << " " << f.getTargetGridName() << " "
                  << f.getSourcePartitionerName() << " " << f.getTargetPartitionerName() << " "
                  << f.getSourceLevels() << std::endl;
     }

    using atlas::grid::detail::partitioner::TransPartitioner;
    using atlas::grid::detail::partitioner::CheckerboardPartitioner;

    // get main program to allocate input and output FieldSets.
    std::vector<atlas::functionspace::StructuredColumns> srcFS;
    std::vector<atlas::functionspace::StructuredColumns> tgtFS;
    std::vector<atlas::StructuredGrid> srcG;
    std::vector<atlas::StructuredGrid> tgtG;
    std::vector<TransPartitioner> srcP;
    std::vector<CheckerboardPartitioner> tgtP;
    atlas::FieldSet srcFieldSet;
    atlas::FieldSet tgtFieldSet;

    for (std::size_t i= 0; i < fs.getNoFields(); ++i) {
        auto f = FieldConfiguration(fs.getFieldConfiguration(i));

        // create grid
        srcG.emplace_back(f.getSourceGridName());
        tgtG.emplace_back(f.getTargetGridName());

        // create partitioner
      //  srcP.emplace_back(CheckerboardPartitioner());
      //  tgtP.emplace_back(TransPartitioner());

        // extract SLat S L Slon G F O
        if ( (f.getSourceGridName().compare(0, 1, "S") == 0)|| (f.getSourceGridName().compare(0, 1, "L") == 0) ||
             (f.getSourceGridName().compare(0, 1, "G") == 0) || (f.getSourceGridName().compare(0, 1, "F") == 0)  ||
             (f.getSourceGridName().compare(0, 1, "O") == 0) ) {

           const auto funConfig =  util::Config( "halo", 1 ) |
                                   util::Config( "levels", f.getSourceLevels());

           srcFS.push_back(atlas::functionspace::StructuredColumns(srcG[i], new TransPartitioner(), funConfig));


           srcFieldSet.add(srcFS[i].createField<double>(
                 util::Config( "name",  "field " + std::to_string(i) + " with levels " + std::to_string(f.getSourceLevels()) ) ) );

        }

        if ( (f.getTargetGridName().compare(0, 1, "S") == 0) || (f.getTargetGridName().compare(0, 1, "L") == 0) ||
             (f.getTargetGridName().compare(0, 1, "G") == 0) || (f.getTargetGridName().compare(0, 1, "F") == 0) ||
             (f.getTargetGridName().compare(0, 1, "O") == 0) ) {

            const auto funConfig =  util::Config( "halo", 1 ) |
                                    util::Config( "levels", f.getSourceLevels());

            tgtFS.push_back(atlas::functionspace::StructuredColumns(tgtG[i], new CheckerboardPartitioner(), funConfig));

            tgtFieldSet.add( tgtFS[i].createField<double>(
                  util::Config( "name",  "field " + std::to_string(i) + " with levels " + std::to_string(f.getSourceLevels()) ) ) );
        }

        if ( (f.getSourceGridName().compare(0, 2, "CS") == 0) ||  (f.getTargetGridName().compare(0, 2, "CS") == 0) ){
            // create CubedSphere mesh

            // Set mesh config.
            const auto meshConfig = util::Config( "partitioner", f.getSourcePartitionerName() ) |
                                    util::Config( "halo", 1 ) |
                                    util::Config( "levels", f.getSourceLevels());

            // Set mesh generator.
            const auto meshGen = MeshGenerator( "cubedsphere", meshConfig );

            if ( f.getSourceGridName().compare(0, 2, "CS") ) {

                const auto mesh = meshGen.generate( srcG[i] );

                srcFieldSet.add( atlas::functionspace::CubedSphereCellColumns(mesh).createField<double>(
                  util::Config( "name", "field " + std::to_string(i) + " with levels " + std::to_string(f.getSourceLevels()) ) ) );
            }

            if ( f.getTargetGridName().compare(0, 2, "CS") ) {

                const auto mesh = meshGen.generate( tgtG[i] );

                tgtFieldSet.add( atlas::functionspace::CubedSphereCellColumns(mesh).createField<double>(
                  util::Config( "name", "field " + std::to_string(i) + " with levels " + std::to_string(f.getSourceLevels()) ) ) );
            }
        }
    }

    // populate the source fields and take halo exchange;
    idx_t i(0);
    for (auto & Field : srcFieldSet) {

        auto f = FieldConfiguration(fs.getFieldConfiguration(i));

        auto testView1 = array::make_view<double, 2>( Field );

        // Willem it would be nice to be able to collapse this.
        // Ideally some like
        //  auto funcS = ( f.getSourceGridName().compare(0, 2, "CS") ?
        //           functionspace::CubedSphereCellColumns(Field.functionspace()) :
        //           functionspace::StructuredColumns(Field.functionspace()) )
        //
        //
        // Maybe a template solution?


        if ( (f.getSourceGridName().compare(0, 2, "CS") == 0)) {

           auto funcS = functionspace::CubedSphereCellColumns(Field.functionspace());

           auto lonLatView = array::make_view<double, 2>(funcS.lonlat() );

           funcS.parallel_for([&]( idx_t index, idx_t k )
             { testView1( index, k ) = testFunction(lonLatView(index, LON),
                                                    lonLatView(index, LAT)); } );

           funcS.haloExchange( Field );


        } else {

            auto lonLatView = array::make_view<double, 2>(srcFS[i].lonlat() );

            srcFS[i].parallel_for([&]( idx_t index, idx_t k )
              { testView1( index, k ) = testFunction(lonLatView(index, LON),
                                                     lonLatView(index, LAT)); } );

            srcFS[i].haloExchange( Field );

        }

        ++i;
   }


  {
      std::cout << " MATCHNIG MESH" << std::endl;

      auto f = FieldConfiguration(fs.getFieldConfiguration(0));

      auto funConfig =  util::Config( "halo", 1 ) |
                        util::Config( "levels", f.getSourceLevels());

      atlas::StructuredGrid sg(f.getSourceGridName());
      atlas::StructuredGrid tg(f.getTargetGridName());

      atlas::functionspace::StructuredColumns tgtFS_1(tg, new CheckerboardPartitioner(), funConfig);

      atlas::functionspace::StructuredColumns matchFS(sg,
                  atlas::grid::MatchingPartitioner(tgtFS_1),
                  atlas::option::levels(f.getSourceLevels()) | atlas::option::halo(1));


       std::cout << " MATCHING MESH END" << std::endl;
   }

   std::cout << " MATCHNIG MESH" << std::endl;

   // get main program to interface to wrapper.
   // note that (unlike here) the instantiation of the wrapper object may use
   // FieldSets that are out of scope by the time they are used to do interpolation

   eckit::LocalConfiguration conf;
   conf.set("interpolation method", "linear");

   std::cout << "create interpolation obj" << std::endl;

   InterpRedistr interp(conf, srcFieldSet, tgtFieldSet);

   std::cout << "execute interpolation obj" << std::endl;

   interp.execute(srcFieldSet, tgtFieldSet);

   std::cout << " done  "  << std::endl;

   // dump out gmsh
   return 1;
}

int main( int argc, char* argv[] ) {
   AtlasParallelInterpolationFieldSet tool( argc, argv );
   return tool.start();
}
