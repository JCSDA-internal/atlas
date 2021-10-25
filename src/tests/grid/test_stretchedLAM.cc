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

using namespace atlas::util;
using namespace atlas::grid;

namespace {
// use a vector of double instead
const double lon_LAM_str[] = {
                                348.9838795150959,  349.8677242998692,  350.6659835640297,  351.3869448275862,
                                352.0380931034483,  352.6892413793103,  353.3403896551724,  353.9915379310345,
                                354.6426862068965,  355.2938344827586,  355.9449827586207,  356.5961310344828,
                                357.2472793103448,  357.8984275862069,  358.5495758620689,  359.200724137931,
                                359.8518724137931,  360.5030206896552,  361.1541689655172,  361.8053172413793,
                                362.4564655172414,  363.1076137931034,  363.7587620689655,  364.4099103448276,
                                365.0610586206897,  365.7122068965517,  366.3633551724138,  367.0843164359702,
                                367.8825757001305,  368.7664204849037
    };

const double lat_LAM_str[] = {
                                -7.411818911205892,  -6.52797499345702,  -5.729716448839834,  -5.008755777777777,
                                -4.357608037037037,  -3.706460296296296,  -3.055312555555555,  -2.404164814814814,
                                -1.753017074074074,  -1.101869333333333,  -0.4507215925925925,  0.2004261481481482,
                                0.8515738888888897,  1.50272162962963,  2.153869370370371,  2.805017111111111,
                                3.456164851851852,  4.107312592592594,  4.758460333333334,  5.409608074074075,
                                6.060755814814815,  6.711903555555557,  7.363051296296296, 8.014199037037038,
                                8.66534677777778,  9.386306244003876,  10.18456345460838,  11.06840589531815
    };

const double lon_LAM_reg[] = {349.4335,  350.084648275862,  350.7357965517241, 351.3869448275862,
                              352.0380931034483,  352.6892413793103,  353.3403896551724,  353.9915379310345,
                              354.6426862068965,  355.2938344827586,  355.9449827586207,  356.5961310344828,
                              357.2472793103448,  357.8984275862069,  358.5495758620689,  359.200724137931,
                              359.8518724137931,  360.5030206896552,  361.1541689655172,  361.8053172413793,
                              362.4564655172414,  363.1076137931034,  363.7587620689655,  364.4099103448276,
                              365.0610586206897,  365.7122068965517,  366.3633551724138,  367.0145034482759,
                              367.6656517241379,  368.3168};

const double lat_LAM_reg[] = {-6.962199,  -6.311051259259259,  -5.659903518518519,  -5.008755777777777,
                              -4.357608037037037,  -3.706460296296296,  -3.055312555555555,  -2.404164814814814,
                              -1.753017074074074,  -1.101869333333333,  -0.4507215925925925,  0.2004261481481482,
                              0.8515738888888897,  1.50272162962963,  2.153869370370371,  2.805017111111111,
                              3.456164851851852,  4.107312592592594,  4.758460333333334,  5.409608074074075,
                              6.060755814814815,  6.711903555555557,  7.363051296296296,  8.014199037037038,
                              8.66534677777778,  9.316494518518519,  9.967642259259259,  10.61879
};

};




namespace atlas {
namespace test {
CASE( "LAMstretch" ) {

    /*
    double delta_low_; //resolution of the host model
    double delta_high_; //resolution of the regional model (regular grid)
    double var_ratio_; //power used for the stretching
    double x_reg_start_; //xstart of the internal regional grid
    double y_reg_start_; //ystart of the internal regional grid
    double x_reg_end_; //xend of the internal regional grid
    double y_reg_end_; ////yend of the internal regional grid
    double startx_; //domain startx
    double endx_; //domain endx
    double starty_; //domain starty
    double endy_; //domain endy
    */



    auto proj = Projection( "stretch", Config( "delta_low", 0.6545448275862076 ) | Config( "delta_hi", 0.6511482758621128 ) |
                                                 Config( "var_ratio", 1.13 ) | Config( "x_reg_start", 351.3869448275862 ) |
                                                 Config( "y_reg_start", -5.008755777777777 ) | Config( "x_reg_end", 366.3633551724138 ) |
                                                 Config( "y_reg_end", 8.66534677777778 ) | Config( "startx", 349.4335) |
                                                 Config( "starty",  -6.962199) | Config( "endx", 368.3168) | Config( "endy", 10.61879) |
                                                 Config( "north_pole", {0.0, 90.0} ) );




    std::cout<< "projection name " << proj.type() <<std::endl;

    atlas::util::Config XSpaceConfig;
    XSpaceConfig.set( "type", "linear" );
    XSpaceConfig.set( "N", 30 );
    XSpaceConfig.set("start", 349.4335 );
    XSpaceConfig.set("end", 368.3168 );

    atlas::grid::detail::grid::Structured::XSpace XS(XSpaceConfig);

    atlas::util::Config YSpaceConfig;
    YSpaceConfig.set( "type", "linear" );
    YSpaceConfig.set( "N", 28 );
    YSpaceConfig.set("start", -6.962199 );
    YSpaceConfig.set("end", 10.61879 );

    atlas::grid::detail::grid::Structured::YSpace YS(YSpaceConfig);

    /*
    atlas::util::Config DomainConfig;
    DomainConfig.set("type", "rectangular"); // can be "global"(default) "zonal_band"
                                             // "empty"
    DomainConfig.set("xmin", 349.4335);
    DomainConfig.set("xmax", 368.3168);
    DomainConfig.set("ymin", -6.962199);
    DomainConfig.set("ymax", 10.61879);
    DomainConfig.set("units", "degrees"); // "degrees" is the default.
    */

    //definition of stretched grid
    //auto grid = ReducedGaussianGrid( nx, proj );
    auto grid_st = StructuredGrid(XS, YS, proj );

    // create regular grid
    atlas::util::Config reg_grid_config;
    reg_grid_config.set( "type", "structured" );
    reg_grid_config.set ("xspace", []() {
      atlas::util::Config config;
      config.set( "type", "linear" );
      config.set( "N", 30 );
      config.set( "start", 349.4335 );
      config.set( "end", 368.3168 );
      return config;
    }() );

    reg_grid_config.set ("yspace", []() {
      atlas::util::Config config;
      config.set( "type", "linear" );
      config.set( "N", 28 );
      config.set( "start", -6.962199 );
      config.set( "end", 10.61879);
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




    //check over regular grid points stretched using new atlas object and check using look-up table
    for (atlas::idx_t j = nodes_reg.j_begin(); j < nodes_reg.j_end(); ++j) {
        for (atlas::idx_t i = nodes_reg.i_begin(j); i < nodes_reg.i_end(j); ++i) {
            auto ll1lon = lon_LAM_str[i];
            auto ll1lat = lat_LAM_str[j];
            auto ll2 = grid_st.lonlat( i, j );
            auto ll2lon = ll2.lon();
            auto ll2lat = ll2.lat();
            //std::cout << " i lon: " << i << " j lat: " << j << std::endl;
            //std::cout << "lon from array ll1: " << ll1lon << " lon computed ll2: " << ll2lon << std::endl;
            EXPECT_APPROX_EQ( ll1lon , ll2lon, 1.e-10 );
            EXPECT_APPROX_EQ( ll1lat, ll2lat, 1.e-10 );
        }
    }



    auto proj_reg = Projection( "stretch", Config( "delta_low", 0.6545448275862076 ) | Config( "delta_hi", 0.6511482758621128 ) |
                                                 Config( "var_ratio", 1.0 ) | Config( "x_reg_start", 351.3869448275862 ) |
                                                 Config( "y_reg_start", -5.008755777777777 ) | Config( "x_reg_end", 366.3633551724138 ) |
                                                 Config( "y_reg_end", 8.66534677777778 ) | Config( "startx", 349.4335) |
                                                 Config( "starty",  -6.962199) | Config( "endx", 368.3168) | Config( "endy", 10.61879) |
                                                 Config( "north_pole", {0.0, 90.0} ) );
    //definition of stretched grid
    auto grid_reg = StructuredGrid(XS, YS, proj_reg );


    std::cout<< "projection name " << proj.type() <<std::endl;

    //definition of stretched grid
    //auto grid_reg = StructuredGrid(proj_grid_config );


    //check over regular grid points stretched using new atlas object and check using look-up table
    for (atlas::idx_t j = nodes_reg.j_begin(); j < nodes_reg.j_end(); ++j) {
        for (atlas::idx_t i = nodes_reg.i_begin(j); i < nodes_reg.i_end(j); ++i) {
            auto ll1lon = lon_LAM_reg[i];
            auto ll1lat = lat_LAM_reg[j];
            auto ll2 = grid_reg.lonlat( i, j );
            auto ll2lon = ll2.lon();
            auto ll2lat = ll2.lat();
            //std::cout << " i lon: " << i << " j lat: " << j << std::endl;
            //std::cout << "lon from array ll1: " << ll1lon << " lon computed ll2: " << ll2lon << std::endl;
            EXPECT_APPROX_EQ( ll1lon , ll2lon, 1.e-10 );
            EXPECT_APPROX_EQ( ll1lat, ll2lat, 1.e-10 );
        }
    }

    /*
    //atlas::util::Config proj_grid_config;
    //proj_grid_config.set( "type", "stretch" );
    proj_reg( "type", "stretch" );
    proj_reg.set( "delta_low",  0.6545448275862076);
    proj_reg.set( "delta_hi", 0.6545448275862076 );
    proj_reg.set( "var_ratio", 1.0 );
    proj_reg.set( "x_reg_start", 351.3869448275862 );
    proj_reg.set( "y_reg_start", -5.008755777777777 );
    proj_reg.set( "x_reg_end", 366.3633551724138 );
    proj_reg.set( "y_reg_end", 8.66534677777778 );
    proj_reg.set( "startx", 349.4335 );
    proj_reg.set( "starty", -6.962199 );
    proj_reg.set( "endx", 368.3168 );
    proj_reg.set( "endy", 10.61879 );
    proj_reg.set( "north_pole", {0.0, 90.0} );
    */

}




}  // namespace test
}  // namespace atlas


int main( int argc, char* argv[] ) {
    return atlas::test::run( argc, argv );
}
