// (C) Copyright 1996-2014 ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation nor
// does it submit to any jurisdiction.


#include "atlas/Gmsh.hpp"
#include "atlas/BuildEdges.hpp"
#include "atlas/BuildDualMesh.hpp"
#include "atlas/BuildPeriodicBoundaries.hpp"
#include "atlas/Partitioner.hpp"
#include "atlas/MPL.hpp"

#define DATADIR ATLAS_DATADIR

using namespace atlas;
int main(int argc, char *argv[])
{
  MPL::init();

  Mesh& mesh = Gmsh::read( std::string(DATADIR) + "/meshes/T47.msh");

  build_periodic_boundaries(mesh);
  build_edges(mesh);
  build_dual_mesh(mesh);

  Gmsh::write(mesh,"bla.msh");
  
  MPL::finalize();
  return 0;
}
