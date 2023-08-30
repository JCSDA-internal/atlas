/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

namespace option {
    struct mpi_comm : public util::Config {
        mpi_comm(std::string_view name) {
            set("mpi_comm",std::string(name));
        }
    };
}

int color() {
    static int c = mpi::comm("world").rank()%2;
    return c;
}

Grid grid() {
    static Grid g (color() == 0 ? "O32" : "N32" );
    return g;
}

struct Fixture {
    Fixture() {
        mpi::comm().split(color(),"split");
    }
    ~Fixture() {
        if (eckit::mpi::hasComm("split")) {
            eckit::mpi::deleteComm("split");
        }
    }
};

CASE("StructuredMeshGenerator") {
    Fixture fixture;

    StructuredMeshGenerator meshgen{option::mpi_comm("split")};
    Mesh mesh = meshgen.generate(grid());
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
}

CASE("Mesh constructor") {
    Fixture fixture;

    auto mesh = Mesh(grid(), option::mpi_comm("split"));
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
