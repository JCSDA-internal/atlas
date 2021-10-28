/*
 *
 * (C) Crown Copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <map>
#include <memory>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "atlas/redistribution/Redistribution.h"

namespace atlas {
class Field;
class FieldSet;
}

namespace eckit {
  class Configuration;
}

// This is based on InterpRedistr, but instead uses the inputFieldSet to
// work out the appropriate keys (based on grid, partitioner_name, level)
// and is no longer hardwired to a Gaussian grid functionspace based on trans partitioner.

class InterpRedistr2  {

public:
    static const std::string classname() {return "unifiedmodel::AtlasInterpWrapper";}

   InterpRedistr2( const eckit::Configuration & conf,
                   const atlas::FieldSet & srcFieldset,
                   const atlas::FieldSet & tarFieldset );

   InterpRedistr2(){};

   void execute( const atlas::Field & srcField,
                       atlas::Field & targetField ) const ;

   void executeAdjoint( atlas::Field & srcField,
                         const atlas::Field & targetField ) const ;

   void execute( const atlas::FieldSet & srcFieldSet,
                        atlas::FieldSet & targetFieldSet ) const ;

   void executeAdjoint( atlas::FieldSet & srcFieldSet,
                         const atlas::FieldSet & targetFieldSet ) const ;


private:
    std::vector<std::tuple<std::string, std::string, std::size_t>> inputFSKeys_;
    std::vector<std::tuple<std::string, std::string, std::size_t>> outputFSKeys_;
    std::set<std::tuple<std::string, std::string, std::size_t>> differingInputFSKeys_;
    std::set<std::tuple<std::string, std::string, std::size_t>> differingOutputFSKeys_;

    std::map<std::tuple<std::string, std::string, std::size_t>,
             atlas::functionspace::StructuredColumns> inputFS_;
    std::map<std::tuple<std::string, std::string, std::size_t>,
             atlas::functionspace::StructuredColumns> outputFS_;
    std::map<std::tuple<std::string, std::string, std::size_t, std::string>,
             atlas::functionspace::StructuredColumns> matchingFS_;

    std::map<std::tuple<std::string, std::string, std::size_t, std::string>, atlas::Interpolation>
    interps_;

    std::map<std::tuple<std::string, std::string, std::size_t, std::string, std::string>, atlas::Redistribution>
    redistr_;

};

