/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "atlas/grid/Partitioner.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/option/Options.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {
namespace {

struct FieldStats {
    double mean{};
    double stddev{};
    double min{};
    double max{};
};

util::Config get_test_config() {
    const auto config_path =
        eckit::Resource<std::string>("--config", "testinput/interpolation_approx_inverse_local_pseudoinverse.yaml");
    return util::Config{eckit::PathName{config_path}};
}

Field get_vortex_field(const FunctionSpace& function_space, const std::string& name) {
    const auto lonlat = array::make_view<double, 2>(function_space.lonlat());

    auto field      = function_space.createField<double>(option::name{name});
    auto field_view = array::make_view<double, 1>(field);
    for (idx_t idx = 0; idx < field.shape(0); ++idx) {
        field_view(idx) = util::function::vortex_rollup(lonlat(idx, atlas::LON), lonlat(idx, atlas::LAT), 1.);
    }

    return field;
}

template <typename InverseInterpolationType>
std::pair<FieldSet, FieldSet> get_fields_and_residuals(const Interpolation& forward_interpolation,
                                                       const InverseInterpolationType& inverse_interpolation) {
    const auto x_field = get_vortex_field(inverse_interpolation.target(), "x_field");
    const auto y_field = get_vortex_field(inverse_interpolation.source(), "y_field");

    const auto abs_difference = [](const Field& lhs, const Field& rhs, const std::string& name) {
        auto diff_field      = lhs.functionspace().createField<double>(option::name(name));
        auto diff_field_view = array::make_view<double, 1>(diff_field);
        const auto lhs_view  = array::make_view<double, 1>(lhs);
        const auto rhs_view  = array::make_view<double, 1>(rhs);

        for (idx_t idx = 0; idx < diff_field.shape(0); ++idx) {
            diff_field_view(idx) = lhs_view(idx) - rhs_view(idx);
        }

        return diff_field;
    };

    const auto do_interpolation = [](const auto& interpolation, const Field& input_field,
                                     const std::string& output_field_name) {
        auto output_field = interpolation.target().template createField<double>(option::name(output_field_name));
        interpolation.execute(input_field, output_field);
        output_field.haloExchange();
        return output_field;
    };

    const auto do_interpolation_adjoint =
        [](const auto& interpolation, const Field& input_field, const std::string& output_field_name) {
            auto output_field = interpolation.source().template createField<double>(option::name(output_field_name));
            array::make_view<double, 1>(output_field).assign(0.);
            interpolation.execute_adjoint(output_field, input_field);
            interpolation.source().haloExchange(output_field);
            return output_field;
        };

    const auto w_x_field       = do_interpolation(forward_interpolation, x_field, "w_x_field");
    const auto w_x_residual= abs_difference(w_x_field, y_field, "w_x_residual");

    const auto wplus_y_field       = do_interpolation(inverse_interpolation, y_field, "wplus_y_field");
    const auto wplus_y_residual = abs_difference(wplus_y_field, x_field, "wplus_y_residual");

    const auto wplus_w_x_field          = do_interpolation(inverse_interpolation, w_x_field, "wplus_w_x_field");
    const auto wplus_w_x_residual = abs_difference(wplus_w_x_field, x_field, "wplus_w_x_residual");

    const auto w_wplus_y_field          = do_interpolation(forward_interpolation, wplus_y_field, "w_wplus_y_field");
    const auto w_wplus_y_residual = abs_difference(w_wplus_y_field, y_field, "w_wplus_y_residual");

    const auto wT_y_field = do_interpolation_adjoint(forward_interpolation, y_field, "wT_y_field");
    const auto wT_wplus_y_residual =
        do_interpolation_adjoint(forward_interpolation, w_wplus_y_residual, "wT_wplus_y_residual");


    auto x_space_fields = FieldSet{};
    x_space_fields.add(x_field);
    x_space_fields.add(wplus_y_field);
    x_space_fields.add(wplus_y_residual);
    x_space_fields.add(wplus_w_x_field);
    x_space_fields.add(wplus_w_x_residual);
    x_space_fields.add(wT_y_field);
    x_space_fields.add(wT_wplus_y_residual);

    auto y_space_fields = FieldSet{};
    y_space_fields.add(y_field);
    y_space_fields.add(w_x_field);
    y_space_fields.add(w_x_residual);
    y_space_fields.add(w_wplus_y_field);
    y_space_fields.add(w_wplus_y_residual);

    return {x_space_fields, y_space_fields};
}

void make_gmsh_output(const std::string& file_name_prefix, const FieldSet& fields) {
    const auto node_columns = functionspace::NodeColumns(fields[0].functionspace());
    ATLAS_ASSERT(node_columns.valid(), "Expected NodeColumns functionspace for Gmsh output");

    const auto gmsh_config =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto gmsh = output::Gmsh(file_name_prefix, gmsh_config);
    gmsh.write(node_columns.mesh());
    gmsh.write(fields, fields[0].functionspace());
}

FieldStats compute_field_stats(const Field& field) {
    const auto field_view = array::make_view<double, 1>(field);
    const auto ghost_view = array::make_view<int, 1>(field.functionspace().ghost());

    double local_sum         = 0.;
    double local_sum_squares = 0.;
    double local_min         = std::numeric_limits<double>::max();
    double local_max         = std::numeric_limits<double>::lowest();
    int local_count          = 0;

    for (idx_t idx = 0; idx < field.shape(0); ++idx) {
        if (ghost_view(idx)) {
            continue;
        }
        const auto value = field_view(idx);
        local_sum += value;
        local_sum_squares += value * value;
        local_min = std::min(local_min, value);
        local_max = std::max(local_max, value);
        ++local_count;
    }

    auto& comm = eckit::mpi::comm();
    comm.allReduceInPlace(local_sum, eckit::mpi::Operation::SUM);
    comm.allReduceInPlace(local_sum_squares, eckit::mpi::Operation::SUM);
    comm.allReduceInPlace(local_min, eckit::mpi::Operation::MIN);
    comm.allReduceInPlace(local_max, eckit::mpi::Operation::MAX);
    comm.allReduceInPlace(local_count, eckit::mpi::Operation::SUM);

    FieldStats stats{};
    stats.mean          = local_sum / local_count;
    const auto variance = local_sum_squares / local_count - stats.mean * stats.mean;
    stats.stddev        = std::sqrt(variance);
    stats.min           = local_min;
    stats.max           = local_max;
    return stats;
}

void log_field_stats(const std::string& label, const Field& field) {
    const auto stats = compute_field_stats(field);
    if (eckit::mpi::comm().rank() == 0) {
        Log::info() << "  " << label << ": mean=" << stats.mean << ", stddev=" << stats.stddev << ", min=" << stats.min
                    << ", max=" << stats.max << std::endl;
    }
}

std::pair<FunctionSpace, FunctionSpace> make_functionspaces(const util::Config& run_config) {
    const util::Config source_grid_config        = run_config.getSubConfiguration("source grid");
    const util::Config target_grid_config        = run_config.getSubConfiguration("target grid");
    const util::Config source_partitioner_config = run_config.getSubConfiguration("source partitioner");
    const util::Config target_partitioner_config = run_config.getSubConfiguration("target partitioner");
    const util::Config source_mesh_config        = run_config.getSubConfiguration("source mesh generator");
    const util::Config target_mesh_config        = run_config.getSubConfiguration("target mesh generator");
    const auto source_halo_size                  = static_cast<size_t>(run_config.getInt("source halo"));
    const auto target_halo_size                  = static_cast<size_t>(run_config.getInt("target halo"));

    const auto target_grid           = Grid(target_grid_config);
    const auto target_partitioner    = grid::Partitioner{target_partitioner_config};
    const auto target_mesh_generator = MeshGenerator{target_mesh_config | option::halo{target_halo_size}};
    const auto target_mesh           = target_mesh_generator.generate(target_grid, target_partitioner);
    const auto target_function_space = functionspace::NodeColumns{target_mesh, option::halo{target_halo_size}};

    const auto source_grid        = Grid(source_grid_config);
    const auto source_partitioner = [&]() -> grid::Partitioner {
        bool matching_mesh = source_partitioner_config.getBool("matching mesh", false);
        if (matching_mesh) {
            return grid::MatchingPartitioner{target_mesh, source_partitioner_config};
        }
        return grid::Partitioner{source_partitioner_config};
    }();
    const auto source_mesh_generator = MeshGenerator{source_mesh_config | option::halo{source_halo_size}};
    const auto source_mesh           = source_mesh_generator.generate(source_grid, source_partitioner);
    const auto source_function_space = functionspace::NodeColumns{source_mesh, option::halo{source_halo_size}};

    return {source_function_space, target_function_space};
}

void run_approx_inverse_case(const util::Config& global_config, const util::Config& run_config) {
    const auto test_name            = run_config.getString("test name");
    const auto interpolation_config = run_config.getSubConfiguration("interpolation");
    const auto inverse_interpolation_config =
        global_config.getSubConfiguration("inverse interpolation").set("scheme", interpolation_config);

    const auto write_gmsh = global_config.getBool("write gmsh");

    const auto [source_function_space, target_function_space] = make_functionspaces(run_config);

    const auto forward_interpolation =
        Interpolation{interpolation_config, target_function_space, source_function_space};

    const auto approx_inverse_scheme = util::Config{"scheme", interpolation_config} | inverse_interpolation_config;
    const auto approx_inverse_interpolation =
        Interpolation{approx_inverse_scheme, source_function_space, target_function_space};

    auto [target_field_set, source_field_set] =
        get_fields_and_residuals(forward_interpolation, approx_inverse_interpolation);

    if (eckit::mpi::comm().rank() == 0) {
        Log::info() << "Approx-inverse run: method=" << inverse_interpolation_config << ", test=" << test_name
                    << std::endl;
    }

    Log::info() << "Target functionspace statistics:" << std::endl;
    for (const auto& field : target_field_set) {
        log_field_stats(field.name(), field);
    }
    Log::info() << std::endl;

    Log::info() << "Source functionspace statistics:" << std::endl;
    for (const auto& field : source_field_set) {
        log_field_stats(field.name(), field);
    }
    Log::info() << std::endl;

    if (!write_gmsh) {
        return;
    }

    const auto output_prefix = test_name + "_approx_inverse_" + inverse_interpolation_config.getString("type");
    make_gmsh_output(output_prefix + "_source.msh", source_field_set);
    make_gmsh_output(output_prefix + "_target.msh", target_field_set);
}

}  // namespace

CASE("Run approx-inverse scenarios from YAML") {
    const auto global_config = get_test_config();
    std::vector<util::Config> runs;
    ATLAS_ASSERT(global_config.get("runs", runs), "Config is missing top-level 'runs' list");

    for (const auto& run_config : runs) {
        run_approx_inverse_case(global_config, run_config);
    }
}

CASE("plot inverse interpolation kernel") {
    const auto global_config        = get_test_config();
    const auto run_config           = global_config.getSubConfiguration("visualise inverse interpolation kernel");
    const auto interpolation_config = run_config.getSubConfiguration("interpolation");
    const auto inverse_interpolation_config =
        global_config.getSubConfiguration("inverse interpolation").set("scheme", interpolation_config);

    const auto [source_function_space, target_function_space] = make_functionspaces(run_config);

    const auto interp_inv = Interpolation{inverse_interpolation_config, source_function_space, target_function_space};

    const auto interp_inv_matrix_storage = interpolation::MatrixCache(interp_inv).matrix();
    const auto interp_inv_matrix         = linalg::make_non_owning_eckit_sparse_matrix(interp_inv_matrix_storage);

    const auto target_ghost_view = array::make_view<int, 1>(target_function_space.ghost());
    const auto target_ridx_view  = array::make_indexview<idx_t, 1>(target_function_space.remote_index());
    const auto target_part_view  = array::make_view<int, 1>(target_function_space.partition());

    std::vector<int> target_remote_indices;
    ATLAS_ASSERT(run_config.get("plot rows", target_remote_indices), "Config is missing 'plot rows' list");
    const auto target_partition = run_config.getInt("plot rank");

    const auto kernel_field_name = inverse_interpolation_config.getString("type") + "_kernel_field";
    auto kernel_field            = source_function_space.createField<double>(option::name(kernel_field_name));
    auto kernel_field_view       = array::make_view<double, 1>(kernel_field);
    kernel_field_view.assign(0.);

    for (auto matrix_itr = interp_inv_matrix.begin(); matrix_itr != interp_inv_matrix.end(); ++matrix_itr) {
        const auto row_idx = matrix_itr.row();
        const auto col_idx = matrix_itr.col();

        if (target_ghost_view(row_idx)) {
            continue;
        }

        for (const auto target_remote_index : target_remote_indices) {
            if (target_ridx_view(row_idx) == target_remote_index && target_part_view(row_idx) == target_partition) {
                kernel_field_view(col_idx) += *matrix_itr;
            }
        }
    }

    // Note: adjointHaloExchange followed by haloExchange is for visualisation purposes only!
    source_function_space.adjointHaloExchange(kernel_field);
    source_function_space.haloExchange(kernel_field);

    auto kernel_field_set = FieldSet{};
    kernel_field_set.add(kernel_field);
    make_gmsh_output(kernel_field_name + ".msh", kernel_field_set);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    // pause for ten seconds
    // std::this_thread::sleep_for(std::chrono::seconds(10));
    return atlas::test::run(argc, argv);
}
