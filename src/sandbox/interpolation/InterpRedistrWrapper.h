/*
 *
 * (C) British Crown Copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_UMJEDI_UTILITIES_ATLASINTERPWRAPPER_H
#define SRC_UMJEDI_UTILITIES_ATLASINTERPWRAPPER_H


#include <string>
#include <map>
#include <memory>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/interpolation.h"
#include "atlas/redistribution/Redistribution.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "um-jedi/Geometry/Geometry.h"
#include "um-jedi/Utilities/BespokeInterpWrapper.h"

namespace atlas {
class Field;
class FieldSet;
}

namespace eckit {
  class Configuration;
}


namespace unifiedmodel {

// I will need to remove the dependency util::Printable and Object Counter
class AtlasInterpWrapper :  public util::Printable,
        private util::ObjectCounter<AtlasInterpWrapper> {

public:
    static const std::string classname() {return "unifiedmodel::AtlasInterpWrapper";}

    AtlasInterpWrapper( const eckit::Configuration & conf,
                        std::shared_ptr<const atlas::FieldSet> & fieldset );

    AtlasInterpWrapper(){};

    void execute( const atlas::Field & srcField,
                        atlas::Field & targetField ) const ;

    void executeAdjoint( atlas::Field & srcField,
                         const atlas::Field & targetField ) const ;

    void execute( const atlas::FieldSet & srcFieldSet,
                        atlas::FieldSet & targetFieldSet ) const ;

    void executeAdjoint( atlas::FieldSet & srcFieldSet,
                         const atlas::FieldSet & targetFieldSet ) const ;


private:
    void print(std::ostream &) const {};
    std::vector<std::string> fieldsetNames_;
    std::vector<std::pair<std::string, std::size_t>> outputFSKeys_;
    std::set<std::pair<std::string, std::size_t>> differingOutputFSKeys_;
    std::map< std::pair< std::string, std::size_t>,
              std::unique_ptr<const atlas::StructuredGrid> > keyOutputGrids_;
    const std::vector<std::string>  gaussNames_;
    atlas::StructuredGrid gaussGrid_;

    std::map<std::pair<std::string, std::size_t>,
             std::unique_ptr<const atlas::functionspace::StructuredColumns>> inputFS_;
    std::map<std::pair<std::string, std::size_t>,
             std::unique_ptr<const atlas::functionspace::StructuredColumns>> outputFS_;
    std::map<std::pair<std::string, std::size_t>,
             std::unique_ptr<const atlas::functionspace::StructuredColumns>> matchingFS_;

    std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::Interpolation>>
    interps_;

    std::map<std::pair<std::string, std::size_t>, std::unique_ptr<const atlas::Redistribution>>
    redistr_;

};

} // namespace unifiedmodel

#endif // SRC_UMJEDI_UTILITIES_ATLASINTERPWRAPPER_H
