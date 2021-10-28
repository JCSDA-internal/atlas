/*
 *
 * (C) Crown Copyright, Met Office.
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

// This is based on InterpRedistr2, but instead allows for
//  for the target to be not fixed to StructuredColumns
//  but instead be also CubedSphere.

class InterpRedistr3  {

public:
    static const std::string classname() {return "unifiedmodel::AtlasInterpWrapper";}

   InterpRedistr3( const eckit::Configuration & conf,
                   const atlas::FieldSet & srcFieldset,
                   const atlas::FieldSet & tarFieldset );

   InterpRedistr3(){};

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
             atlas::FunctionSpace> inputFS_;
    std::map<std::tuple<std::string, std::string, std::size_t>,
             atlas::FunctionSpace> outputFS_;
    std::map<std::tuple<std::string, std::string, std::size_t, std::string>,
             atlas::FunctionSpace> matchingFS_;

    std::map<std::tuple<std::string, std::string, std::size_t, std::string>, atlas::Interpolation>
    interps_;

    std::map<std::tuple<std::string, std::string, std::size_t, std::string, std::string>, atlas::Redistribution>
    redistr_;

};

