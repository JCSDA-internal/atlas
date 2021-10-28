/*
 * (C) Britsh Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <functional>
#include <iterator>
#include <memory>
#include <vector>
#include <tuple>

#include "eckit/exception/Exceptions.h"

#include "atlas/array/MakeView.h"

#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"

#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"


#include "atlas/redistribution/Redistribution.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "atlas/interpolation.h"

#include "InterpRedistrWrapper2.h"

using atlas::grid::detail::partitioner::TransPartitioner;
using atlas::grid::detail::partitioner::MatchingMeshPartitioner;

namespace {

// -----------------------------------------------------------------------------
// atlas wrapper method
// output functionspace keys are product of the grid name and the name of vertical levels
// there is no metadata for the partitioner
// For now do not include partitionerName - we maybe want to include it in the future as part of the string key.
const std::vector<std::tuple<std::string, std::string, std::size_t> > createFunctionSpaceKeys(const atlas::FieldSet & fieldset) {

    std::vector<std::tuple<std::string, std::string, std::size_t> > key;
    for (auto & field : fieldset) {
        key.push_back(std::make_tuple( atlas::functionspace::StructuredColumns(
                                      field.functionspace()).grid().name(),
                                      field.functionspace().distribution(),
                                      field.levels() ) );

    }

    for (auto & k : key) {
        std::cout << atlas::mpi::rank() << " " << std::get<0>(k) << " " << std::get<1>(k)  << " " << std::get<2>(k) << std::endl;
    }

    return key;
}

// -----------------------------------------------------------------------------
const std::set<std::tuple<std::string, std::string, std::size_t>>
createDifferingFunctionSpaceKeys(const std::vector<std::tuple<std::string, std::string, std::size_t>> & FSkeys ) {

   std::set<std::tuple<std::string, std::string, std::size_t>> differingFSKeys;
   std::for_each(FSkeys.begin(), FSkeys.end(),
                 [&](const std::tuple<std::string, std::string, std::size_t> elem){differingFSKeys.insert(elem);});
   return differingFSKeys;
}

// -----------------------------------------------------------------------------
std::map<std::tuple<std::string, std::string, std::size_t>,
atlas::functionspace::StructuredColumns>
    createFS(const std::vector<std::tuple<std::string, std::string, std::size_t>> & FSkeys,
             const std::set<std::tuple<std::string, std::string, std::size_t>> & differingFSkeys,
             const atlas::FieldSet & fieldset) {

    std::map<std::tuple<std::string, std::string, std::size_t>, atlas::functionspace::StructuredColumns>
        functionSpaces;

    std::set<std::tuple<std::string, std::string, std::size_t>> tmpKeys(differingFSkeys);

    atlas::idx_t i(0);

    atlas::functionspace::StructuredColumns fs;
    for (auto & f : fieldset) {
        // get name of functionspace associated with f
        if (tmpKeys.erase(FSkeys[i])) {

           atlas::functionspace::StructuredColumns fs(f.functionspace());

           std::cout << "output grids " << fs.grid().name() << std::endl;

           functionSpaces[FSkeys[i]] = fs;
        }
        ++i;
    }


    return functionSpaces;
}

// -----------------------------------------------------------------------------------
//matching mesh functionspaces
std::map<std::tuple<std::string, std::string, std::size_t, std::string>, atlas::functionspace::StructuredColumns>
createMatchingMeshFunctionSpaces(
      const std::vector<std::tuple<std::string, std::string, std::size_t>> inputKeys,
      const std::vector<std::tuple<std::string, std::string, std::size_t>> outputKeys,
      const std::map<std::tuple<std::string, std::string, std::size_t>,
          atlas::functionspace::StructuredColumns> & inputFS,
      const std::map<std::tuple<std::string, std::string, std::size_t>,
          atlas::functionspace::StructuredColumns> & outputFS) {

    std::map<std::tuple<std::string, std::string, std::size_t, std::string>, atlas::functionspace::StructuredColumns>
        functionSpaces;

    // we expect that the number of Fields in the input FieldSet and the output FieldSet are the same.
    ASSERT(inputKeys.size() == outputKeys.size());

    // calculate and set the number of matching mesh FS keys
    std::set<std::tuple<std::string, std::string, std::size_t, std::string>> matchingMeshKeys;
    for (std::size_t i = 0; i < inputKeys.size(); ++i) {
       matchingMeshKeys.insert(std::make_tuple( std::get<0>(inputKeys[i]), std::get<1>(inputKeys[i]),
                                                std::get<2>(inputKeys[i]), std::get<0>(outputKeys[i])));
    }


    // Note I am recreating the output grid here.
    // This might be more costly than it should
    for (auto & k : matchingMeshKeys) {
        auto inputKey = std::make_tuple(std::get<0>(k), std::get<1>(k), std::get<2>(k));
        auto p = atlas::grid::MatchingPartitioner(inputFS.at(inputKey));
        auto g = atlas::Grid(std::get<3>(k));

        functionSpaces[k] = atlas::functionspace::StructuredColumns
                         ( g, p,
                          atlas::option::levels(std::get<2>(k)) | atlas::option::halo(1) );
    }

    return functionSpaces;
}


// -----------------------------------------------------------------------------
//atlas createAtlasInterpolations
std::map<std::tuple<std::string, std::string, std::size_t, std::string>, atlas::Interpolation>
createAtlasInterpolations(const std::map< std::tuple<std::string, std::string, std::size_t>,
                          atlas::functionspace::StructuredColumns> & inputFS,
                          const std::map< std::tuple<std::string, std::string, std::size_t, std::string>,
                          atlas::functionspace::StructuredColumns> & matchingFS) {

    std::cout << " atlas wrapper construct" << std::endl;

    // maybe we should expose this configuration
    atlas::util::Config interp;
    interp.set( "type", "structured-linear2D" );

    std::map<std::tuple<std::string, std::string, std::size_t, std::string>, atlas::Interpolation> interps;

    for (auto & m : matchingFS) {
        auto matchingKey = m.first;
        auto inputKey = std::make_tuple(std::get<0>(matchingKey), std::get<1>(matchingKey), std::get<2>(matchingKey));

        interps[matchingKey] =
          atlas::Interpolation(
                    interp |
                    atlas::util::Config( "adjoint", true ),
                    inputFS.at(inputKey), matchingFS.at(matchingKey));

    }

    return interps;
};

//-----------------------------------------------------------------------------
//atlas createAtlasRedistributions
std::map<std::tuple<std::string, std::string, std::size_t, std::string, std::string>, atlas::Redistribution>
createAtlasRedistributions(const std::vector<std::tuple<std::string, std::string, std::size_t>> inputKeys,
                           const std::vector<std::tuple<std::string, std::string, std::size_t>> outputKeys,
                           const std::map< std::tuple<std::string, std::string, std::size_t, std::string>,
                           atlas::functionspace::StructuredColumns> & matchingFS,
                           const std::map< std::tuple<std::string, std::string, std::size_t>,
                           atlas::functionspace::StructuredColumns> & outputFS) {

    std::map< std::tuple<std::string, std::string, std::size_t, std::string, std::string>, atlas::Redistribution>
            redistributions;

    //create a set of unique redistribution keys
    std::set<std::tuple<std::string, std::string, std::size_t, std::string, std::string>> redistributionKeys;
    for (std::size_t i = 0; i < inputKeys.size(); ++i) {
       redistributionKeys.insert(std::make_tuple( std::get<0>(inputKeys[i]), std::get<1>(inputKeys[i]),
                                                  std::get<2>(inputKeys[i]), std::get<0>(outputKeys[i]),
                                                  std::get<1>(outputKeys[i])));
    }

    for (auto & k : redistributionKeys) {
        auto matchingKey = std::make_tuple( std::get<0>(k), std::get<1>(k),
                                            std::get<2>(k), std::get<3>(k) );
        auto outputKey = std::make_tuple( std::get<3>(k), std::get<4>(k), std::get<2>(k) );

        redistributions[k] = atlas::Redistribution( matchingFS.at(matchingKey), outputFS.at(outputKey) );

        auto inverseKey = std::make_tuple( std::get<0>(k) + "_inverse", std::get<1>(k),
                                           std::get<2>(k), std::get<3>(k), std::get<4>(k) );

        redistributions[inverseKey] = atlas::Redistribution( outputFS.at(outputKey), matchingFS.at(matchingKey) );

    }

    return redistributions;
};


} // anonymous namespace
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// ATLAS INTERPOLATION WRAPPER
//-------------------------------------------------------------------------------------------------
InterpRedistr2::InterpRedistr2( const eckit::Configuration & conf,
                              const atlas::FieldSet & srcFieldset,
                              const atlas::FieldSet & tarFieldset) :
     inputFSKeys_(createFunctionSpaceKeys(srcFieldset)),
     outputFSKeys_(createFunctionSpaceKeys(tarFieldset)),
     differingInputFSKeys_(createDifferingFunctionSpaceKeys(inputFSKeys_)),
     differingOutputFSKeys_(createDifferingFunctionSpaceKeys(outputFSKeys_)),
     inputFS_(createFS(inputFSKeys_, differingInputFSKeys_, srcFieldset)),
     outputFS_(createFS(outputFSKeys_, differingOutputFSKeys_, tarFieldset)),
     matchingFS_(createMatchingMeshFunctionSpaces(inputFSKeys_, outputFSKeys_, inputFS_, outputFS_)),
     interps_(createAtlasInterpolations(inputFS_, matchingFS_)),
     redistr_(createAtlasRedistributions(inputFSKeys_, outputFSKeys_, matchingFS_,outputFS_))
{
};

//-------------------------------------------------------------------------------------------------
void InterpRedistr2::execute( const atlas::Field & srcField,
                              atlas::Field & tgtField ) const {

    // create a field of matchingmesh functionspace
    auto srcFS =
        atlas::functionspace::StructuredColumns(srcField.functionspace());
    auto tgtFS =
        atlas::functionspace::StructuredColumns(tgtField.functionspace());

    auto src_v = atlas::array::make_view<double, 2>( srcField );

    auto srcTmp = srcFS.createField<double>(atlas::option::name(srcField.name()) | atlas::option::levels(srcField.levels()));
    auto srcTmp_v = atlas::array::make_view<double, 2>( srcTmp );

    for ( atlas::idx_t t = 0; t < srcField.shape( 0 ); ++t ) {
        for ( atlas::idx_t k = 0; k < srcField.shape( 1 ); ++k ) {
            srcTmp_v( t, k ) = src_v( t, k );
        }
    }
    srcTmp.haloExchange();

    std::tuple<std::string, std::string, std::size_t, std::string>
       matchingKey = std::make_tuple( srcFS.grid().name(), srcFS.distribution(), srcFS.levels(),
                                      tgtFS.grid().name() );

    std::tuple<std::string, std::string, std::size_t, std::string, std::string>
       redistributionKey = std::make_tuple( srcFS.grid().name(), srcFS.distribution(), srcFS.levels(),
                                            tgtFS.grid().name(), tgtFS.distribution() );

    auto tmpField = (*(matchingFS_.at(matchingKey))).createField<double>(atlas::option::name(tgtField.name()));

    auto tmp_v = atlas::array::make_view<double, 2>( tmpField );
    tmp_v.assign(0.0);

    (*(interps_.at(matchingKey))).execute(srcTmp, tmpField);

    (*(redistr_.at(redistributionKey))).execute(tmpField, tgtField);

    tgtField.haloExchange();

}

void InterpRedistr2::executeAdjoint( atlas::Field & srcField,
                                     const atlas::Field & tgtField ) const {

    // create a field that is a copy of target
    // take halo adjoint of temp field
    atlas::Field tmpField = tgtField.functionspace().createField<double>(atlas::option::name(tgtField.name()));
    auto tmp_v = atlas::array::make_view<double, 2>( tmpField );
    auto tgt_v = atlas::array::make_view<double, 2>( tgtField );

    for ( atlas::idx_t t = 0; t < tgtField.shape( 0 ); ++t ) {
        for ( atlas::idx_t k = 0; k < tgtField.shape( 1 ); ++k ) {
            tmp_v( t, k ) = tgt_v( t, k );
        }
    }
    tmpField.adjointHaloExchange();

    auto srcFS =
        atlas::functionspace::StructuredColumns(srcField.functionspace());
    auto tgtFS =
        atlas::functionspace::StructuredColumns(tgtField.functionspace());

    std::tuple<std::string, std::string, std::size_t, std::string>
       matchingKey = std::make_tuple( srcFS.grid().name(), srcFS.distribution(), srcFS.levels(),
                                      tgtFS.grid().name() );

    std::tuple<std::string, std::string, std::size_t, std::string, std::string>
       redistributionKey = std::make_tuple( srcFS.grid().name() + "_inverse", srcFS.distribution(), srcFS.levels(),
                                            tgtFS.grid().name(), tgtFS.distribution() );

    auto tmp2Field = (*(matchingFS_.at(matchingKey))).createField<double>(atlas::option::name(tgtField.name()));
    auto tmp2_v = atlas::array::make_view<double, 2>( tmp2Field );
    tmp2_v.assign(0.0);

    // redistribution

    (*(redistr_.at(redistributionKey))).execute(tmpField, tmp2Field);

    (*(interps_.at(matchingKey))).execute_adjoint(srcField, tmp2Field);

}

void InterpRedistr2::execute( const atlas::FieldSet & srcFieldSet,
                                  atlas::FieldSet & targetFieldSet ) const {

    for (auto srcField : srcFieldSet) {
        std::string nameStr = srcField.name();
        execute(srcField, targetFieldSet[nameStr]);
    }

}

void InterpRedistr2::executeAdjoint( atlas::FieldSet & srcFieldSet,
                                         const atlas::FieldSet & targetFieldSet ) const {

    for (auto srcField : srcFieldSet) {
        std::string nameStr = srcField.name();
        executeAdjoint(srcField, targetFieldSet[nameStr]);
    }

}

// ------------------------------------------------------------------------------------------------

