/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <vector>

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/KDTree.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::util;

namespace atlas {
namespace test {

//------------------------------------------------------------------------------------------------
namespace {  // helpers

template <typename T>
std::string to_string(const std::vector<T>& vector) {
    std::stringstream s;
    s << vector;
    return s.str();
}
template <typename T>
std::string to_string(const T& v) {
    return std::to_string(v);
}

class PayloadGenerator {
public:
    class iterator {
        friend class PayloadGenerator;

    public:
        long int operator*() const { return i_; }
        const iterator& operator++() {
            ++i_;
            return *this;
        }
        /*
        // unused
        iterator operator++( int ) {
            iterator copy( *this );
            ++i_;
            return copy;
        }*/

        bool operator==(const iterator& other) const { return i_ == other.i_; }
        bool operator!=(const iterator& other) const { return i_ != other.i_; }

    protected:
        iterator(long int start): i_(start) {}

    private:
        unsigned long i_;
    };

    iterator begin() const { return begin_; }
    iterator end() const { return end_; }
    PayloadGenerator(long int end): begin_(0), end_(end) {}

    template <typename Container, typename Value = typename Container::value_type>
    void fill(Container& container) {
        std::iota(container.begin(), container.end(), Value(*begin_));
    }

    template <typename Value, size_t size>
    std::array<Value, size> make_array() {
        std::array<Value, size> array;
        fill(array);
        return array;
    }

private:
    iterator begin_;
    iterator end_;
};

static double radius() {
    return util::Earth::radius();
}

static Geometry& geometry() {
    static Geometry _geometry(radius());
    return _geometry;
}

static PointXYZ make_xyz(const PointLonLat& lonlat) {
    return geometry().xyz(lonlat);
}

static std::array<double, 7>& test_lon() {
    static std::array<double, 7> lon = {0., 30., 60., 90., 120., 150., 180.};
    return lon;
}

static std::array<double, 7>& test_lat() {
    static std::array<double, 7> lat = {90., 60., 30., 0., -30., -60., -90.};
    return lat;
}

static std::array<PointLonLat, 7>& test_lonlat() {
    static std::array<PointLonLat, 7> lonlat{PointLonLat{0., 90.},   PointLonLat{30., 60.},   PointLonLat{60., 30.},
                                             PointLonLat{90., 0.},   PointLonLat{120., -30.}, PointLonLat{150., -60.},
                                             PointLonLat{180., -90.}};
    return lonlat;
}

static std::array<PointXYZ, 7>& test_xyz() {
    static std::array<PointXYZ, 7> xyz{make_xyz(PointLonLat{0., 90.}),    make_xyz(PointLonLat{30., 60.}),
                                       make_xyz(PointLonLat{60., 30.}),   make_xyz(PointLonLat{90., 0.}),
                                       make_xyz(PointLonLat{120., -30.}), make_xyz(PointLonLat{150., -60.}),
                                       make_xyz(PointLonLat{180., -90.})};
    return xyz;
}

static std::array<idx_t, 7>& test_payloads() {
    static auto payloads = PayloadGenerator(7).make_array<idx_t, 7>();
    EXPECT_EQ(std::distance(payloads.begin(), payloads.end()), 7);
    return payloads;
}

void validate(const KDTree<idx_t>& tree) {
    EXPECT_NO_THROW(tree.closestPoint(PointLonLat{180., 45.}));
    // Search 4 nearest neighbours (k=4), sorted by shortest distance
    auto neighbours        = tree.closestPoints(PointLonLat{89.9, 44.9}, 4);
    auto expected_payloads = std::vector<idx_t>{2, 1, 3, 0};
    EXPECT_EQ(neighbours.payloads(), expected_payloads);
}

static const IndexKDTree& nanoflann_search() {
    static IndexKDTree kdtree = []() {
        util::IndexKDTree kdtree{new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry())};
        auto grid = Grid{"O32"};
        kdtree.build(grid.lonlat(), PayloadGenerator(grid.size()));
        return kdtree;
    }();
    return kdtree;
}

template <typename TreeImpl>
auto time_kdtree_build_and_query(
    const std::vector<PointLonLat>& index_points,
    const std::vector<idx_t>& payloads,
    const std::vector<PointLonLat>& query_points,
    int num_closet_points
) {
    util::IndexKDTree tree_{new TreeImpl(geometry())};

    auto build_t0 = std::chrono::steady_clock::now();
    tree_.build(index_points, payloads);
    auto build_t1 = std::chrono::steady_clock::now();

    auto query_t0 = std::chrono::steady_clock::now();
    for (const auto& p : query_points) {
        tree_.closestPoints(p, num_closet_points);
    }
    auto query_t1 = std::chrono::steady_clock::now();

    return std::make_pair(
        std::chrono::duration_cast<std::chrono::milliseconds>(build_t1 - build_t0).count(),
        std::chrono::duration_cast<std::chrono::milliseconds>(query_t1 - query_t0).count());
}

template <typename TreeImpl1, typename TreeImpl2>
std::pair<int, int> compare_kdtree_query_results(
    const std::vector<PointLonLat>& index_points,
    const std::vector<idx_t>& payloads,
    const std::vector<PointLonLat>& query_points,
    int num_closet_points,
    bool print_differences = false
) {

    util::IndexKDTree ek{new TreeImpl1(geometry())};
    util::IndexKDTree nf{new TreeImpl2(geometry())};
    ek.build(index_points, payloads);
    nf.build(index_points, payloads);
    int count_differences = 0;
    int count_errors = 0;
    for (const auto& p : query_points) {
        
        auto ek_point = ek.closestPoint(p);
        auto ek_payload = ek_point.payload();
        auto ek_distance = ek_point.distance();
        
        auto nf_point = nf.closestPoint(p);
        auto nf_payload = nf_point.payload();
        auto nf_distance = nf_point.distance();
        
        if (ek_payload != nf_payload) {

            if (ek_distance > nf_distance) {
                if (print_differences) {
                    std::cout << "Warning: different distances causes a difference in results ek_payload=" << ek_payload << " ek_distance=" << ek_distance
                              << " nf_payload=" << nf_payload << " nf_distance=" << nf_distance << std::endl;
                }
                ++count_differences;
                continue;
            }

            auto nf_points = nf.closestPoints(p, num_closet_points);
            
            std::optional<size_t> maybe_payload_idx = std::nullopt;
            for (size_t i = 0; i < nf_points.size(); ++i) {
                if (nf_points[i].payload() == ek_payload) {
                    maybe_payload_idx = std::make_optional<size_t>(i);
                    break;
                }
            }

            if (maybe_payload_idx.has_value() && nf_points[maybe_payload_idx.value()].distance() == ek_distance) {
                if (print_differences) {
                    std::cout << "Warning: tie-breaker issue (ek_payload found in closest points search with same distance):"
                                << " ek_payload=" << ek_payload << " ek_distance=" << ek_distance
                                << " nf_payload=" << nf_payload << " nf_distance=" << nf_distance << std::endl;
                }
                ++count_differences;
            } else {
                if (print_differences) {
                    std::cout << "Error: closest points search didn't have target payload,"
                              << " ek_payload=" << ek_payload << " ek_distance=" << ek_distance
                              << ", in the results." << std::endl;
                }
                ++count_errors;
            }
        }
    }

    return std::make_pair(count_errors, count_differences);
}

}  // namespace
//------------------------------------------------------------------------------------------------

static bool ECKIT_515_implemented = ATLAS_ECKIT_VERSION_AT_LEAST(1, 13, 2);
// --> implements eckit::KDTree::size() and eckit::KDTree::empty()

CASE("compare eckit and nanoflann kdtrees - cs to random") {
    const auto grid = Grid("CS-LFR-C-300");
    const auto config = util::Config("partitioner", "equal_regions") | util::Config("halo", 3);
    const auto mesh = atlas::MeshGenerator("cubedsphere", config).generate(grid);
    const auto csGrid = CubedSphereGrid(mesh.grid());

    // Get views to cell data.
    const auto lonlatView   = atlas::array::make_view<double, 2>(mesh.cells().field("lonlat"));
    const auto haloView = atlas::array::make_view<int, 1>(mesh.cells().halo());

    // make points and payloads vectors.
    auto points   = std::vector<PointLonLat>{};
    auto payloads = std::vector<idx_t>{};

    // Iterate over cells.
    auto halo = config.getInt("halo", 0);
    for (idx_t i = 0; i < mesh.cells().size(); ++i) {
        if (haloView(i) <= halo) {
            points.emplace_back(lonlatView(i, LON), lonlatView(i, LAT));
            payloads.emplace_back(i);
        }
    }

    // Create a random lon-lat list of the same size as points
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> lon_dist(-180.0, 180.0);
    std::uniform_real_distribution<> lat_dist(-90.0, 90.0);
    std::vector<PointLonLat> random_points;
    random_points.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        random_points.emplace_back(lon_dist(gen), lat_dist(gen));
    }

    const int num_closet_points = 3;

    // time nanoflann kd-tree.
    const auto [nanoflann_build_duration, nanoflann_query_duration] = 
        time_kdtree_build_and_query<util::detail::KDTree_nanoflann<idx_t, Point3>>(
            points, payloads, random_points, num_closet_points
        );

    // time eckit kd-tree (KDTreeMemory).
    const auto [eckit_build_duration, eckit_query_duration] =
        time_kdtree_build_and_query<util::detail::KDTreeMemory<idx_t, Point3>>(
            points, payloads, random_points, num_closet_points
        );

    std::cout << "KDTree Grid: " << grid.name() << " size=" << grid.size() << std::endl;
    std::cout << "Query Grid: random size=" << random_points.size() << std::endl;
    std::cout << "Eckit KDTree:\n"
              << "- build time: " << eckit_build_duration << " ms\n"
              << "- query time: " << eckit_query_duration << " ms\n";
    std::cout << "Nanoflann KDTree:\n"
              << "- build time = " << nanoflann_build_duration
              << " ms (speed-up = " << (double(eckit_build_duration) / double(nanoflann_build_duration)) << ")\n"
              << "- query time = " << nanoflann_query_duration
              << " ms (speed-up = " << (double(eckit_query_duration) / double(nanoflann_query_duration)) << ")\n"
              << "- total time = " << (nanoflann_build_duration + nanoflann_query_duration)
              << " ms (speed-up = " << (double(eckit_build_duration + eckit_query_duration) / double(nanoflann_build_duration + nanoflann_query_duration)) << ")\n";

    auto [errors, differences] = 
        compare_kdtree_query_results<util::detail::KDTreeMemory<idx_t, Point3>, util::detail::KDTree_nanoflann<idx_t, Point3>>(
            points, payloads, random_points, num_closet_points
        );
    
    std::cout << "Total errors: " << errors << " out of " << random_points.size() << " queries." << std::endl;
    std::cout << "Total differences: " << differences << " out of " << random_points.size() << " queries." << std::endl;
    EXPECT(errors == 0);
}

CASE("compare eckit and nanoflann kdtrees - cs to gauss") {
    const auto grid = Grid("CS-LFR-C-300");
    const auto config = util::Config("partitioner", "equal_regions") | util::Config("halo", 3);
    const auto mesh = atlas::MeshGenerator("cubedsphere", config).generate(grid);
    const auto csGrid = CubedSphereGrid(mesh.grid());

    // Get views to cell data.
    const auto lonlatView   = atlas::array::make_view<double, 2>(mesh.cells().field("lonlat"));
    const auto haloView = atlas::array::make_view<int, 1>(mesh.cells().halo());

    // make points and payloads vectors.
    auto points   = std::vector<PointLonLat>{};
    auto payloads = std::vector<idx_t>{};

    // Iterate over cells.
    auto halo = config.getInt("halo", 0);
    for (idx_t i = 0; i < mesh.cells().size(); ++i) {
        if (haloView(i) <= halo) {
            points.emplace_back(lonlatView(i, LON), lonlatView(i, LAT));
            payloads.emplace_back(i);
        }
    }

    const auto gaussGrid = Grid("F320");
    const auto gaussLonLat = gaussGrid.lonlat();
    std::vector<PointLonLat> gauss_points;
    gauss_points.reserve(gaussGrid.size());
    for (const auto& p : gaussLonLat) {
        gauss_points.emplace_back(p);
    }

    int num_closet_points = 5;

    // time nanoflann kd-tree.
    const auto [nanoflann_build_duration, nanoflann_query_duration] = 
        time_kdtree_build_and_query<util::detail::KDTree_nanoflann<idx_t, Point3>>(
            points, payloads, gauss_points, num_closet_points
        );

    // time eckit kd-tree.
    const auto [eckit_build_duration, eckit_query_duration] =
        time_kdtree_build_and_query<util::detail::KDTreeMemory<idx_t, Point3>>(
            points, payloads, gauss_points, num_closet_points
        );

    std::cout << "KDTree Grid: " << grid.name() << " size=" << grid.size() << std::endl;
    std::cout << "Query Grid: " << gaussGrid.name() << " size=" << gaussGrid.size() << std::endl;
    std::cout << "Eckit KDTree:\n"
              << "- build time: " << eckit_build_duration << " ms\n"
              << "- query time: " << eckit_query_duration << " ms\n";
    std::cout << "Nanoflann KDTree:\n"
              << "- build time = " << nanoflann_build_duration
              << " ms (speed-up = " << (double(eckit_build_duration) / double(nanoflann_build_duration)) << ")\n"
              << "- query time = " << nanoflann_query_duration
              << " ms (speed-up = " << (double(eckit_query_duration) / double(nanoflann_query_duration)) << ")\n"
              << "- total time = " << (nanoflann_build_duration + nanoflann_query_duration)
              << " ms (speed-up = " << (double(eckit_build_duration + eckit_query_duration) / double(nanoflann_build_duration + nanoflann_query_duration)) << ")\n";

    auto [errors, differences] =
        compare_kdtree_query_results<util::detail::KDTreeMemory<idx_t, Point3>, util::detail::KDTree_nanoflann<idx_t, Point3>>(
            points, payloads, gauss_points, num_closet_points
        );
    
    std::cout << "Total errors: " << errors << " out of " << gauss_points.size() << " queries." << std::endl;
    std::cout << "Total differences: " << differences << " out of " << gauss_points.size() << " queries." << std::endl;
    
    EXPECT(errors == 0);
}

CASE("test kdtree") {
    auto grid = Grid{"O32"};

    util::IndexKDTree search{new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry())};
    EXPECT(search.empty());

    search.reserve(grid.size());
    idx_t n{0};
    for (auto& point : grid.lonlat()) {
        search.insert(point, n++);
    }
    search.build();
    EXPECT_EQ(search.size(), grid.size());
    EXPECT_NO_THROW(search.closestPoint(PointLonLat{180., 45.}));
    auto neighbours          = search.closestPoints(PointLonLat{180., 45.}, 4).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761};
    EXPECT_EQ(neighbours, expected_neighbours);
}

CASE("test assertion") {
    auto grid = Grid{"O32"};

    IndexKDTree search(new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry()));
    search.reserve(grid.size());
    idx_t n{0};
    for (auto& point : grid.lonlat()) {
        search.insert(point, n++);
    }
    // Forgot to call search.build() --> assertion thrown when trying to access
    EXPECT_THROWS_AS(search.closestPoint(PointLonLat{180., 45.}), eckit::AssertionFailed);
}

CASE("test no assertion") {
    // Like case "test assertion", but without reserving size
    auto grid = Grid{"O32"};

    IndexKDTree search(new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry()));
    // No search.reserve() --> build() will not be necessary.
    EXPECT(search.empty());
    idx_t n{0};
    for (auto& point : grid.lonlat()) {
        search.insert(point, n++);
        if (ECKIT_515_implemented) {
            EXPECT_EQ(search.size(), n);
        }
    }
    EXPECT_EQ(search.size(), grid.size());
    // search.build() Not required
    EXPECT_NO_THROW(search.closestPoint(PointLonLat{180., 45.}));
}

CASE("test kdtree building with separate lon and lat and payload arrays") {
    IndexKDTree search(new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry()));
    search.build(test_lon(), test_lat(), test_payloads());
    validate(search);
}

CASE("test kdtree building with separate lon and lat and raw payload iterators") {
    IndexKDTree search(new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry()));
    auto lon       = test_lon();
    auto lat       = test_lat();
    auto payloads_ = test_payloads();
    search.build(lon.begin(), lon.end(), lat.begin(), lat.end(), payloads_.begin(), payloads_.end());
    validate(search);
}

CASE("test kdtree building with separate PointLonLat and payload containers") {
    IndexKDTree search(new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry()));
    search.build(test_lonlat(), test_payloads());
    validate(search);
}

CASE("test kdtree building with separate PointXYZ and payload containers") {
    IndexKDTree search(new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry()));
    search.build(test_xyz(), test_payloads());
    validate(search);
}

CASE("test assignment") {
    IndexKDTree search;
    search = IndexKDTree(new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry()));
    search.build(test_lonlat(), test_payloads());
    validate(search);
}

CASE("test closestPoint") {
    auto neighbour          = nanoflann_search().closestPoint(PointLonLat{180., 45.}).payload();
    auto expected_neighbour = 760;
    EXPECT_EQ(neighbour, expected_neighbour);
}

CASE("test closestPoints") {
    auto neighbours          = nanoflann_search().closestPoints(PointLonLat{180., 45.}, 4).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761};
    EXPECT_EQ(neighbours, expected_neighbours);
}

CASE("test closestPointsWithinRadius") {
    double km                = 1000. * radius() / util::Earth::radius();
    auto neighbours          = nanoflann_search().closestPointsWithinRadius(PointLonLat{180., 45.}, 500 * km).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761, 841, 843, 682};
    EXPECT_EQ(neighbours, expected_neighbours);
}

CASE("test IndexKDTree 2D vs 3D") {
    IndexKDTree2D search2d(new util::detail::KDTree_nanoflann<idx_t, Point2>(geometry()));
    IndexKDTree3D search3d(new util::detail::KDTree_nanoflann<idx_t, Point3>(geometry()));
    search2d.build(test_lonlat(), test_payloads());
    search3d.build(test_lonlat(), test_payloads());
    auto payloads2d = search2d.closestPoints(PointLonLat{89.9, 44.9}, 4).payloads();
    auto payloads3d = search3d.closestPoints(PointLonLat{89.9, 44.9}, 4).payloads();
    EXPECT_EQ(payloads2d, (std::vector<idx_t>{2, 3, 1, 4}));
    EXPECT_EQ(payloads3d, (std::vector<idx_t>{2, 1, 3, 0}));
    // Note that the expected values are different whether 2D search or 3D search is used
}


// CASE("test kdtree with configured geometry") {
//     auto grid = Grid{"O32"};

//     IndexKDTree search_unit(util::Config("geometry","UnitSphere"));
//     IndexKDTree search_earth(util::Config("geometry","Earth"));

//     search_unit.build(grid.lonlat(),PayloadGenerator(grid.size()));
//     search_earth.build(grid.lonlat(),PayloadGenerator(grid.size()));

//     double km_unit                = 1000. / util::Earth::radius();
//     double km_earth               = 1000.;
//     auto neighbours_unit          = search_unit. closestPointsWithinRadius(PointLonLat{180., 45.}, 500 * km_unit ).payloads();
//     auto neighbours_earth         = search_earth.closestPointsWithinRadius(PointLonLat{180., 45.}, 500 * km_earth).payloads();

//     auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761, 841, 843, 682};

//     EXPECT_EQ(neighbours_unit,  expected_neighbours);
//     EXPECT_EQ(neighbours_earth, expected_neighbours);
// }

//------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
