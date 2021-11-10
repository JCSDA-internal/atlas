/**
* (C) Crown copyright 2021, Met Office
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/


#include <algorithm>
#include <vector>

#include "atlas/array.h"
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/projection.h"
#include "atlas/projection/detail/StretchProjection.h"
#include "tests/AtlasTestEnvironment.h"
#include <iomanip> ///< for setprecision

using namespace atlas::util;
using namespace atlas::grid;

namespace {

// use vectors of double for testing
const double lon_LAM_str[] = {
                                346.9838795150959,  347.9838795150959,
                                348.9838795150959,  349.8677242998692,  350.6659835640297,  351.3869448275862,
                                352.0380931034483,  352.6892413793103,  353.3403896551724,  353.9915379310345,
                                354.6426862068965,  355.2938344827586,  355.9449827586207,  356.5961310344828,
                                357.2472793103448,  357.8984275862069,  358.5495758620689,  359.200724137931,
                                359.8518724137931,  360.5030206896552,  361.1541689655172,  361.8053172413793,
                                362.4564655172414,  363.1076137931034,  363.7587620689655,  364.4099103448276,
                                365.0610586206897,  365.7122068965517,  366.3633551724138,  367.0843164359702,
                                367.8825757001305,  368.7664204849037,
                                369.7664204849037,  370.7664204849037
    };

const double lat_LAM_str[] = {  -9.411818911205892,  -8.411818911205892,
                                -7.411818911205892,  -6.52797499345702,  -5.729716448839834,  -5.008755777777777,
                                -4.357608037037037,  -3.706460296296296,  -3.055312555555555,  -2.404164814814814,
                                -1.753017074074074,  -1.101869333333333,  -0.4507215925925925,  0.2004261481481482,
                                0.8515738888888897,  1.50272162962963,  2.153869370370371,  2.805017111111111,
                                3.456164851851852,  4.107312592592594,  4.758460333333334,  5.409608074074075,
                                6.060755814814815,  6.711903555555557,  7.363051296296296, 8.014199037037038,
                                8.66534677777778,  9.386306244003876,  10.18456345460838,  11.06840589531815,
                                12.06840589531815, 13.06840589531815
    };
/*
const double lon_LAM_reg[] = { 348.13120344827576, 348.78235172413787,
                              349.4335,  350.084648275862,  350.7357965517241, 351.3869448275862,
                              352.0380931034483,  352.6892413793103,  353.3403896551724,  353.9915379310345,
                              354.6426862068965,  355.2938344827586,  355.9449827586207,  356.5961310344828,
                              357.2472793103448,  357.8984275862069,  358.5495758620689,  359.200724137931,
                              359.8518724137931,  360.5030206896552,  361.1541689655172,  361.8053172413793,
                              362.4564655172414,  363.1076137931034,  363.7587620689655,  364.4099103448276,
                              365.0610586206897,  365.7122068965517,  366.3633551724138,  367.0145034482759,
                              367.6656517241379,  368.3168, 368.9679482758621, 369.6190965517242
                              };
*/

const double lon_LAM_reg[] = {348.13120344827576, 348.78235172413787,
                              349.4335,   350.084648275862094,   350.735796551724206,   351.386944827586319,   352.038093103448432,
                              352.689241379310545,   353.340389655172658,   353.991537931034770,   354.642686206896883,   355.293834482758996,
                              355.944982758621109,   356.596131034483221,   357.247279310345334,   357.898427586207447,   358.549575862069560,
                              359.200724137931672,   359.851872413793785,   360.503020689655898,   361.154168965518011,   361.805317241380124,
                              362.456465517242236,   363.107613793104349,   363.758762068966462,   364.409910344828575,   365.061058620690687,
                              365.712206896552800,   366.363355172414913,   367.014503448277026,   367.665651724139138,   368.3168,
                              368.9679482758621, 369.6190965517242};

/*
const double lat_LAM_reg[] = { -8.264494481481481, -7.613346740740741,
                              -6.962199,  -6.311051259259259,  -5.659903518518519,  -5.008755777777777,
                              -4.357608037037037,  -3.706460296296296,  -3.055312555555555,  -2.404164814814814,
                              -1.753017074074074,  -1.101869333333333,  -0.4507215925925925,  0.2004261481481482,
                              0.8515738888888897,  1.50272162962963,  2.153869370370371,  2.805017111111111,
                              3.456164851851852,  4.107312592592594,  4.758460333333334,  5.409608074074075,
                              6.060755814814815,  6.711903555555557,  7.363051296296296,  8.014199037037038,
                              8.66534677777778,  9.316494518518519,  9.967642259259259,  10.61879,
                              11.26993774074074, 11.921085481481482
};
*/


const double lat_LAM_reg[] = { -8.264495551724226, -7.613347275862113,
                              -6.962199,   -6.311050724137887,   -5.659902448275774,   -5.008754172413662,   -4.357605896551549,
                              -3.706457620689436,   -3.055309344827323,   -2.404161068965211,   -1.753012793103098,   -1.101864517240985,
                              -0.450716241378872,   0.200432034483240,   0.851580310345353,   1.502728586207466,   2.153876862069579,
                              2.805025137931692,   3.456173413793804,   4.107321689655917,   4.758469965518030,   5.409618241380143,
                              6.060766517242255,   6.711914793104368,   7.363063068966481,   8.014211344828594,   8.665359620690706,
                              9.316507896552819,   9.967656172414932,   10.618804448277045,
                              11.269952724139157, 11.92110100000127
                              };

};




namespace atlas {
namespace test {
CASE( "LAMstretch" ) {

    auto proj = Projection( "stretch", Config( "delta_low", 1. ) | Config( "delta_hi", 0.6511482758621128 ) |
                                                 Config( "var_ratio", 1.13 ) | Config( "x_reg_start", 351.386944827586319 ) |
                                                 Config( "y_reg_start", -5.008754172413662 ) | Config( "x_reg_end", 366.363355172414913 ) |
                                                 Config( "y_reg_end", 8.665359620690706 ) | Config( "startx", 348.13120344827576) |
                                                 Config( "starty",  -8.264495551724226) | Config( "endx", 369.6190965517242) | Config( "endy", 11.92110100000127) |
                                                 Config( "north_pole", {0.0, 90.0} ) | Config( "rim_widthx", 4. ) | Config( "rim_widthy", 4. ) );





    atlas::util::Config XSpaceConfig;
    XSpaceConfig.set( "type", "linear" );
    XSpaceConfig.set( "N", 34 );
    XSpaceConfig.set("start", 348.13120344827576 );
    XSpaceConfig.set("end", 369.6190965517242 );

    atlas::grid::detail::grid::Structured::XSpace XS(XSpaceConfig);

    atlas::util::Config YSpaceConfig;
    YSpaceConfig.set( "type", "linear" );
    YSpaceConfig.set( "N", 32 );
    YSpaceConfig.set("start",  -8.264494481481481 );
    YSpaceConfig.set("end", 11.921085481481482 );

    atlas::grid::detail::grid::Structured::YSpace YS(YSpaceConfig);

    // definition of stretched grid
    auto grid_st = StructuredGrid(XS, YS, proj );

    // create regular grid
    atlas::util::Config reg_grid_config;
    reg_grid_config.set( "type", "structured" );
    reg_grid_config.set ("xspace", []() {
      atlas::util::Config config;
      config.set( "type", "linear" );
      config.set( "N", 34 );
      config.set( "start", 348.13120344827576 );
      config.set( "end", 369.6190965517242 );
      return config;
    }() );

    reg_grid_config.set ("yspace", []() {
      atlas::util::Config config;
      config.set( "type", "linear" );
      config.set( "N", 32 );
      config.set( "start", -8.264494481481481 );
      config.set( "end", 11.921085481481482);
      return config;
    }() );

    atlas::StructuredGrid reg_grid(reg_grid_config);

    atlas::functionspace::StructuredColumns nodes_reg(
           reg_grid, atlas::grid::Partitioner( "equal_regions" ), reg_grid_config );

    int sizei = nodes_reg.i_end(0) - nodes_reg.i_begin(0);
    int sizej = nodes_reg.j_end() - nodes_reg.j_begin();

    std::vector<double> lon_p_arr_st(sizei);
    std::vector<double> lon_p_arr(sizei);
    std::vector<double> lat_p_arr_st(sizej);
    std::vector<double> lat_p_arr(sizej);

    // check over regular grid points stretched using new atlas object and check using look-up table
    for (atlas::idx_t j = nodes_reg.j_begin(); j < nodes_reg.j_end(); ++j) {
        for (atlas::idx_t i = nodes_reg.i_begin(j); i < nodes_reg.i_end(j); ++i) {
            auto ll1lon = lon_LAM_str[i];
            auto ll1lat = lat_LAM_str[j];
            auto ll2 = grid_st.lonlat( i, j );
            auto ll2lon = ll2.lon();
            auto ll2lat = ll2.lat();
            std::cout << " i lon: " << i << " j lat: " << j << std::endl;
            std::cout << " i lon ar: " << ll1lon << " j lat ar: " << ll1lat << std::endl;
            std::cout << " i lon comp: " << ll2lon << " j lat comp: " << ll2lat << std::endl;
            std::cout <<  std::setprecision(18) << " dlon : " << ll1lon - ll2lon << std::endl;
            std::cout <<  std::setprecision(18) << " dlat : " << ll1lat - ll2lat << std::endl;
            std::cout << " "  << std::endl;
            EXPECT_APPROX_EQ( ll1lon , ll2lon, 1.e-10 );
            EXPECT_APPROX_EQ( ll1lat, ll2lat, 1.e-10 );
        }
    }


    auto proj_reg = Projection( "stretch", Config( "delta_low", 1. ) | Config( "delta_hi", 0.6511482758621128 ) |
                                                 Config( "var_ratio", 1.0 ) | Config( "x_reg_start", 351.386944827586319 ) |
                                                 Config( "y_reg_start", -5.008754172413662 ) | Config( "x_reg_end", 366.363355172414913 ) |
                                                 Config( "y_reg_end", 8.665359620690706 ) | Config( "startx", 348.13120344827576) |
                                                 Config( "starty",  -8.264495551724226) | Config( "endx", 369.6190965517242) | Config( "endy", 11.92110100000127) |
                                                 Config( "north_pole", {0.0, 90.0} ) | Config( "rim_widthx", 4. ) | Config( "rim_widthy", 4. ) );



    // definition of stretched grid
    auto grid_reg = StructuredGrid(XS, YS, proj_reg );
    // Check if bounding box is correct
        {
            RectangularLonLatDomain bb{grid_reg.lonlatBoundingBox()};
            EXPECT( RectangularLonLatDomain( bb ) );
            const double tolerance = 1.e-6;
            EXPECT_APPROX_EQ( bb.west(), 348.13120344827576, tolerance );
            EXPECT_APPROX_EQ( bb.east(), 369.6190965517242, tolerance );
            EXPECT_APPROX_EQ( bb.south(), -8.264494481481481, tolerance );
            EXPECT_APPROX_EQ( bb.north(), 11.921085481481482, tolerance );
            for ( PointLonLat p : grid_reg.lonlat() ) {
                EXPECT( bb.contains( p ) );
            }
          }




    // check over regular grid points stretched using new atlas object and check using look-up table
    for (atlas::idx_t j = nodes_reg.j_begin(); j < nodes_reg.j_end(); ++j) {
        for (atlas::idx_t i = nodes_reg.i_begin(j); i < nodes_reg.i_end(j); ++i) {
            auto ll1lon = lon_LAM_reg[i];
            auto ll1lat = lat_LAM_reg[j];
            auto ll2 = grid_reg.lonlat( i, j );
            auto ll2lon = ll2.lon();
            auto ll2lat = ll2.lat();
            EXPECT_APPROX_EQ( ll1lon , ll2lon, 1.e-10 );
            EXPECT_APPROX_EQ( ll1lat, ll2lat, 1.e-10 );
        }
    }

}




}  // namespace test
}  // namespace atlas


int main( int argc, char* argv[] ) {
    return atlas::test::run( argc, argv );
}
