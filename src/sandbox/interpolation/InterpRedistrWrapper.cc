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

#include "eckit/exception/Exceptions.h"

#include "atlas/array/MakeView.h"

#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"

#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"


#include "atlas/redistribution/Redistribution.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "atlas/interpolation.h"

#include "InterpRedistrWrapper.h"

using atlas::grid::detail::partitioner::TransPartitioner;


namespace {


// -----------------------------------------------------------------------------
atlas::StructuredGrid createGaussGrid(std::shared_ptr<const atlas::FieldSet> & fieldset) {
    // Atlas uses "F" (followed by a number # denoting wavenumber) to denote a classic Gaussian grid
    // with # wavenumbers.
    // The maximum wavenumber is defined as the number of points around a latitude circle
    // divided by 4. ie.(frequencies that are not affected by aliasing)
    //
    // We are assuming here that there is a single input grid that is used for all interpolations
    // and that only the number of levels can change.
    //
    std::vector<std::string> fieldnames = fieldset->field_names();

    auto field_fs = atlas::functionspace::StructuredColumns((*fieldset)[fieldnames[0]].functionspace());
    int longitudes = field_fs.grid().nx(0);

    std::string grid_uid( "F" + std::to_string( longitudes/4 ) );
    atlas::StructuredGrid gaussGrid( grid_uid );
    return gaussGrid;
}


// -----------------------------------------------------------------------------
const std::vector<std::string> prependFieldNamesWithGauss(
        std::vector<std::string> & fieldNames)  {

    std::vector<std::string> gaussNames;
    // create field names prepended with "gauss_"
    std::string gaussPrepend("gauss_");
    for (std::string & name : fieldNames) {
        gaussNames.push_back(gaussPrepend + name);
    }

    return gaussNames;
}


// -----------------------------------------------------------------------------
// atlas wrapper method
// output functionspace keys are product of the grid name and the name of vertical levels
// there is no metadata for the partitioner
// For now do not include partitionerName - we maybe want to include it in the future as part of the string key.
const std::vector<std::pair<std::string, std::size_t> > createOutputFunctionSpaceKeys(const atlas::FieldSet & fieldset) {

    std::vector<std::pair<std::string, std::size_t> > key;
    for (auto field : fieldset) {
        key.push_back(std::make_pair( atlas::functionspace::StructuredColumns(
                                      field.functionspace()).grid().name(),
                                      field.levels() ) );
    }

    return key;
}

// -----------------------------------------------------------------------------
const std::set<std::pair<std::string, std::size_t>>
createDifferingOutputFunctionSpaceKeys(const std::vector<std::pair<std::string, std::size_t>> & outputFSkeys  ) {

   std::set<std::pair<std::string, std::size_t>> differingOutputFSKeys;
   std::for_each(outputFSkeys.begin(), outputFSkeys.end(),
                 [&](const std::pair<std::string, std::size_t> elem){differingOutputFSKeys.insert(elem);});
   return differingOutputFSKeys;
}

// -----------------------------------------------------------------------------
std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::StructuredGrid>>
manageGridSubset(const std::vector<std::pair<std::string, std::size_t>> & outputFSkeys,
                 const std::set<std::pair<std::string, std::size_t>> & differingOutputFSkeys,
                 const atlas::FieldSet & fieldset) {

    std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::StructuredGrid>> keysGrids;

    std::set<std::pair<std::string, std::size_t>> tmpKeys(differingOutputFSkeys);

    atlas::idx_t i(0);

    atlas::StructuredGrid grid;
    for (auto f : fieldset) {
        if (tmpKeys.erase(outputFSkeys[i])) {

            grid = atlas::StructuredGrid(atlas::functionspace::StructuredColumns(f.functionspace()).grid());

            keysGrids[outputFSkeys[i]] =
                 std::unique_ptr<const atlas::StructuredGrid>( new const atlas::StructuredGrid(
                    grid) );

        }
        ++i;
    }

    return keysGrids;
}


// -----------------------------------------------------------------------------
std::map<std::pair<std::string, std::size_t>,
std::unique_ptr<const atlas::functionspace::StructuredColumns>>
    createInputFS( const atlas::StructuredGrid & gaussGrid,
    const std::set<std::pair<std::string, std::size_t>> & differingOutputFSKeys) {

    std::map<std::pair<std::string, std::size_t>,
            std::unique_ptr<const atlas::functionspace::StructuredColumns> >
        functionSpaces;

    // extract a set of levels from the output keys
    std::set<std::size_t> uniqueLevels;
    for ( auto & p : differingOutputFSKeys) {
        uniqueLevels.insert(p.second);
    }

    std::map< std::pair<std::string, std::size_t>,
        std::unique_ptr<const atlas::StructuredGrid> > keyInputGrids;

    for (auto l : uniqueLevels) {
       std::pair<std::string, std::size_t> key = std::make_pair(gaussGrid.name(), l);

       keyInputGrids[key] =
           std::unique_ptr<atlas::StructuredGrid>(new atlas::StructuredGrid(gaussGrid) );
    }


    // Note that below is very similar to Matching Mesh ...
    // ... except for the Partitioner.
    for (auto & s : keyInputGrids) {

        auto key = s.first;

        atlas::functionspace::StructuredColumns inputFS(
                    *(s.second),
                    atlas::grid::Partitioner(new TransPartitioner() ),
                    atlas::option::levels(key.second) | atlas::option::halo(1) );

        functionSpaces[key] =
             std::unique_ptr<const atlas::functionspace::StructuredColumns>(
                    new const atlas::functionspace::StructuredColumns(inputFS));
    }

    return functionSpaces;

}

// -----------------------------------------------------------------------------
std::map<std::pair<std::string, std::size_t>,
std::unique_ptr<const atlas::functionspace::StructuredColumns>>
    createOutputFS(const std::vector<std::pair<std::string, std::size_t>> & outputFSkeys,
                   const std::set<std::pair<std::string, std::size_t>> & differingOutputFSkeys,
                   const atlas::FieldSet & fieldset) {

    std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::functionspace::StructuredColumns>>
        functionSpaces;

    std::set<std::pair<std::string, std::size_t>> tmpKeys(differingOutputFSkeys);

    atlas::idx_t i(0);

    atlas::functionspace::StructuredColumns fs;
    for (auto & f : fieldset) {
        // get name of functionspace associated with f
        if (tmpKeys.erase(outputFSkeys[i])) {

            fs = atlas::functionspace::StructuredColumns(f.functionspace());

            functionSpaces[outputFSkeys[i]] =
               std::unique_ptr<const atlas::functionspace::StructuredColumns>( new
                  const atlas::functionspace::StructuredColumns(fs) );
        }
        ++i;
    }


    return functionSpaces;
}

// -----------------------------------------------------------------------------------
//matching mesh functionspaces
std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::functionspace::StructuredColumns>>
createMatchingMeshFunctionSpaces(
        const std::map<std::pair<std::string, std::size_t>,
        std::unique_ptr<const atlas::functionspace::StructuredColumns>> & inputFS,
        const std::map<std::pair<std::string, std::size_t>,
        std::unique_ptr<const atlas::StructuredGrid>> & keyOutputGrids) {

    std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::functionspace::StructuredColumns>>
        functionSpaces;

    // assuming that inputFS has the same grid throughout.
    // but has differing sets of levels.
    // so we create an error trap for this and extract the grid name.
    std::set<std::string> inputGridNames;
    for (auto & f : inputFS) {
        auto key = f.first;
        inputGridNames.insert(key.first);
    }
    if (inputGridNames.size() != 1 ) {
        std::string err_message("input functionspaces with different horizontal grid are not currently supported ");
        throw eckit::UnexpectedState(err_message);
    }

    std::pair<std::string, std::size_t> inputKey;
    inputKey.first = *(inputGridNames.begin());
    for (auto & s : keyOutputGrids) {
        auto key = s.first;

        // create input key
        inputKey.second = key.second;

        atlas::functionspace::StructuredColumns
                outputFS( *(s.second),
                          atlas::grid::MatchingPartitioner( *(inputFS.at(inputKey)) ),
                          atlas::option::levels(key.second) | atlas::option::halo(1) );

        functionSpaces[key] =
            std::unique_ptr<const atlas::functionspace::StructuredColumns>(
                new const atlas::functionspace::StructuredColumns(outputFS) );

    }

    return functionSpaces;
}


// -----------------------------------------------------------------------------
//atlas createAtlasInterpolations
std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::Interpolation>>
createAtlasInterpolations(const std::map< std::pair<std::string, std::size_t>,
                          std::unique_ptr<const atlas::functionspace::StructuredColumns>> & inputFS,
                          const std::map< std::pair<std::string, std::size_t>,
                          std::unique_ptr<const atlas::functionspace::StructuredColumns>> & matchingFS,
                          const std::map< std::pair<std::string, std::size_t>,
                          std::unique_ptr<const atlas::StructuredGrid>> & keysGrids) {

    std::cout << " atlas wrapper construct" << std::endl;

    // maybe we should expose this configuration
    atlas::util::Config interp;
    interp.set( "type", "structured-linear2D" );

    std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::Interpolation>> interps;

    auto inputKeyFirst = (*inputFS.begin()).first;

    std::pair<std::string, std::size_t> inputKey;

    // here inputKey.first refers to the grid name
    // we are assuming that for the input FS the grid name is unchanged
    // That is why the line below works.
    inputKey.first = inputKeyFirst.first;
    for (auto & s : keysGrids) {
        auto key = s.first;

        //here inputKey.second refers to the model levels
        inputKey.second = key.second;

       // std::unique_ptr<MyClass> test(new MyClass(data));

        interps[key] =
          std::unique_ptr<const atlas::Interpolation>( new const atlas::Interpolation(
                    interp |
                    atlas::util::Config( "adjoint", true ),
                    *(inputFS.at(inputKey)), *(matchingFS.at(key))));
    }
    return interps;
};

//-----------------------------------------------------------------------------
//atlas createAtlasRedistributions
std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::Redistribution>>
createAtlasRedistributions(const std::map< std::pair<std::string, std::size_t>,
                           std::unique_ptr<const atlas::functionspace::StructuredColumns>> & matchingFS,
                           const std::map< std::pair<std::string, std::size_t>,
                           std::unique_ptr<const atlas::functionspace::StructuredColumns>> & outputFS) {

    std::map< std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::Redistribution>>
            redistributions;

    for (auto & f : matchingFS) {
        std::pair<std::string, std::size_t> key = f.first;

        redistributions[key] =

            std::unique_ptr<const atlas::Redistribution>( new const atlas::Redistribution(
                *(f.second), *(outputFS.at(key)) ) );

        std::pair<std::string, std::size_t> key2 = make_pair(key.first + "_inverse", key.second);

        redistributions[key2] =
            std::unique_ptr<const atlas::Redistribution>( new const atlas::Redistribution(
                *(outputFS.at(key)), *(f.second) ) );
    }

    return redistributions;
};


} // anonymous namespace
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// ATLAS INTERPOLATION WRAPPER
//-------------------------------------------------------------------------------------------------
InterpRedistr::InterpRedistr( const eckit::Configuration & conf,
                                        std::shared_ptr<const atlas::FieldSet> & fieldset) :
     fieldsetNames_((*fieldset).field_names()),
     outputFSKeys_(createOutputFunctionSpaceKeys(*fieldset)),
     differingOutputFSKeys_(createDifferingOutputFunctionSpaceKeys(outputFSKeys_)),
     keyOutputGrids_(manageGridSubset(outputFSKeys_, differingOutputFSKeys_, *fieldset)),
     gaussNames_(prependFieldNamesWithGauss(fieldsetNames_)),
     gaussGrid_(createGaussGrid(fieldset)),
     inputFS_(createInputFS(gaussGrid_, differingOutputFSKeys_)),
     outputFS_(createOutputFS(outputFSKeys_, differingOutputFSKeys_, *fieldset)),
     matchingFS_(createMatchingMeshFunctionSpaces(inputFS_, keyOutputGrids_)),
     interps_(createAtlasInterpolations(inputFS_, matchingFS_, keyOutputGrids_)),
     redistr_(createAtlasRedistributions(matchingFS_,outputFS_))
{
};

//-------------------------------------------------------------------------------------------------
void InterpRedistr::execute( const atlas::Field & srcField,
                                  atlas::Field & targetField ) const {

    // create a field of matchingmesh functionspace
    auto srcFS =
        atlas::functionspace::StructuredColumns(srcField.functionspace());
    auto targetFS =
        atlas::functionspace::StructuredColumns(targetField.functionspace());

    auto src_v = atlas::array::make_view<double, 2>( srcField );

    auto srcTmp = srcFS.createField<double>(atlas::option::name(srcField.name()) | atlas::option::levels(srcField.levels()));
    auto srcTmp_v = atlas::array::make_view<double, 2>( srcTmp );

    for ( atlas::idx_t t = 0; t < srcField.shape( 0 ); ++t ) {
        for ( atlas::idx_t k = 0; k < srcField.shape( 1 ); ++k ) {
            srcTmp_v( t, k ) = src_v( t, k );
        }
    }
    srcTmp.haloExchange();

    std::pair<std::string, std::size_t> key = std::make_pair( targetFS.grid().name(), targetFS.levels() );

    auto tmpField = (*(matchingFS_.at(key))).createField<double>(atlas::option::name(targetField.name()));

    auto tmp_v = atlas::array::make_view<double, 2>( tmpField );
    tmp_v.assign(0.0);

    (*(interps_.at(key))).execute(srcTmp, tmpField);

    (*(redistr_.at(key))).execute(tmpField, targetField);

    targetField.haloExchange();

}

void InterpRedistr::executeAdjoint( atlas::Field & srcField,
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

    auto targetFS =
        atlas::functionspace::StructuredColumns(tgtField.functionspace());
    std::pair<std::string, std::size_t> key = std::make_pair( targetFS.grid().name(), targetFS.levels() );
    std::pair<std::string, std::size_t> key2 = std::make_pair( key.first + "_inverse", key.second );

    auto tmp2Field = (*(matchingFS_.at(key))).createField<double>(atlas::option::name(tgtField.name()));
    auto tmp2_v = atlas::array::make_view<double, 2>( tmp2Field );
    tmp2_v.assign(0.0);

    // redistribution

    (*(redistr_.at(key2))).execute(tmpField, tmp2Field);

    (*(interps_.at(key))).execute_adjoint(srcField, tmp2Field);

}

void InterpRedistr::execute( const atlas::FieldSet & srcFieldSet,
                                  atlas::FieldSet & targetFieldSet ) const {

    for (auto srcField : srcFieldSet) {
        std::string nameStr = srcField.name();
        execute(srcField, targetFieldSet[nameStr.erase(0,6)]);
    }

}

void InterpRedistr::executeAdjoint( atlas::FieldSet & srcFieldSet,
                                         const atlas::FieldSet & targetFieldSet ) const {

    for (auto srcField : srcFieldSet) {
        std::string nameStr = srcField.name();
        executeAdjoint(srcField, targetFieldSet[nameStr.erase(0,6)]);
    }

}

// ------------------------------------------------------------------------------------------------

