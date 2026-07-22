/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <string>
#include <utility>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/interpolation.h"
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


/// function to generate a synthetic field (vortex field)
Field getVortexField(const FunctionSpace& functionSpace) {
    const auto lonlat = array::make_view<double, 2>(functionSpace.lonlat());

    auto vField     = functionSpace.createField<double>(option::name{"vortex field"});
    auto vFieldView = array::make_view<double, 1>(vField);
    for (idx_t idx = 0; idx < vField.shape(0); ++idx) {
        vFieldView(idx) = util::function::vortex_rollup(lonlat(idx, atlas::LON), lonlat(idx, atlas::LAT), 1.);
    }

    return vField;
}

/// Generate and the difference between a target field and the analytic solution
Field getErrorField(const Field& targetField) {
    auto errorField               = targetField.functionspace().createField<double>(option::name("error field"));
    const auto referenceField     = getVortexField(targetField.functionspace());
    auto errorFieldView           = array::make_view<double, 1>(errorField);
    const auto referenceFieldView = array::make_view<double, 1>(referenceField);
    const auto targetFieldView    = array::make_view<double, 1>(targetField);
    for (idx_t idx = 0; idx < errorField.shape(0); ++idx) {
        errorFieldView(idx) = std::abs(targetFieldView(idx) - referenceFieldView(idx));
    }
    return errorField;
}

/// Generate the residual between low_res -> high_res -> low_res transformation
std::tuple<Field, Field, Field> getResidualField(const Interpolation& interpScheme, const Interpolation& binningScheme) {
    const auto lowResField = getVortexField(interpScheme.source());

    auto highResField = interpScheme.target().createField<double>(option::name("high res field"));
    interpScheme.execute(lowResField, highResField);
    highResField.haloExchange();

    auto residualField = lowResField.functionspace().createField<double>(option::name("residual field"));
    binningScheme.execute(highResField, residualField);
    residualField.haloExchange();

    auto residualFieldView     = array::make_view<double, 1>(residualField);
    const auto lowResFieldView = array::make_view<double, 1>(lowResField);
    for (idx_t idx = 0; idx < residualField.shape(0); ++idx) {
        residualFieldView(idx) = std::abs(residualFieldView(idx) - lowResFieldView(idx));
    }

    return {lowResField, highResField, residualField};
}

/// Generate the residual between high_res -> low_res -> high_res transformation
std::tuple<Field, Field, Field> getResidualFieldHighToLowToHigh(const Interpolation& interpScheme, const Interpolation& binningScheme) {
    const auto highResField = getVortexField(interpScheme.target());

    auto lowResField = interpScheme.source().createField<double>(option::name("low res field"));
    binningScheme.execute(highResField, lowResField);
    lowResField.haloExchange();

    auto residualField = highResField.functionspace().createField<double>(option::name("residual field"));
    interpScheme.execute(lowResField, residualField);
    residualField.haloExchange();

    auto residualFieldView      = array::make_view<double, 1>(residualField);
    const auto highResFieldView = array::make_view<double, 1>(highResField);
    for (idx_t idx = 0; idx < residualField.shape(0); ++idx) {
        residualFieldView(idx) = std::abs(residualFieldView(idx) - highResFieldView(idx));
    }

    return {lowResField, highResField, residualField};
}

/// function to write a field set in a Gmsh file
void makeGmshOutput(const std::string& fileNamePrefix, const FieldSet& fields) {
    const auto& fs = fields[0].functionspace();

    const auto mesh = [&]() {
        if (const auto nc = functionspace::NodeColumns(fs)) {
            return nc.mesh();
        }
        if (const auto sc = functionspace::StructuredColumns(fs)) {
            return MeshGenerator{"structured", option::halo(sc.halo())}.generate(sc.grid());
        }
        throw_Exception("Unsupported function space type for Gmsh output: " + fs.type());
    }();

    const auto gmshConfig =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto gmsh = output::Gmsh(fileNamePrefix, gmshConfig);
    gmsh.write(mesh);
    gmsh.write(fields, fs);
}


/// function to carry out a dot product
double dotProd(const Field& fieldA, const Field& fieldB) {
    double dprod{};

    const auto field01_view = array::make_view<double, 1>(fieldA);
    const auto field02_view = array::make_view<double, 1>(fieldB);

    for (idx_t i = 0; i < field01_view.shape(0); ++i) {
        dprod += field01_view(i) * field02_view(i);
    }
    eckit::mpi::comm().allReduceInPlace(dprod, eckit::mpi::Operation::SUM);

    return dprod;
}


void regriddingTest(const util::Config& config) {
    const auto sourceGridName = config.getString("source_grid");
    const auto targetGridName = config.getString("target_grid");


    const auto sourceGrid = Grid(sourceGridName);
    const auto targetGrid = Grid(targetGridName);

    const auto [sourceFunctionSpace, targetFunctionSpace] = [&]() -> std::pair<FunctionSpace, FunctionSpace> {
        const auto fSpaceName = config.getString("functionspace");
        const size_t haloSize = config.getInt("halo");
        if (fSpaceName == "NodeColumns") {
            const auto meshGenerator =
                MeshGenerator{config.getString("mesh_generator"),
                              option::halo{haloSize} | util::Config("partitioner", "equal_regions")};

            const auto sourceMesh = meshGenerator.generate(sourceGrid);
            const auto targetMesh = meshGenerator.generate(targetGrid);

            return {functionspace::NodeColumns{sourceMesh, option::halo{haloSize}},
                    functionspace::NodeColumns{targetMesh, option::halo{haloSize}}};
        }
        if (fSpaceName == "StructuredColumns") {
            return {functionspace::StructuredColumns{sourceGrid, option::halo{haloSize}},
                    functionspace::StructuredColumns{targetGrid, option::halo{haloSize}}};
        }
        throw_Exception("Unsupported functionspace type: " + fSpaceName);
    }();

    const auto sourceField = getVortexField(sourceFunctionSpace);
    auto targetField       = targetFunctionSpace.createField<double>(option::name("field_target"));

    const auto interpScheme = config.getSubConfiguration("scheme");
    const auto interp       = Interpolation{interpScheme, targetFunctionSpace, sourceFunctionSpace};

    const auto binningScheme = option::type{"local-pseudoinverse"} | util::Config{"scheme", interpScheme};
    const auto binning       = Interpolation{binningScheme, sourceFunctionSpace, targetFunctionSpace};

    sourceField.haloExchange();
    binning.execute(sourceField, targetField);
    targetField.haloExchange();

    auto targetFieldSet = FieldSet{};
    targetFieldSet.add(targetField);
    targetFieldSet.add(getErrorField(targetField));

    auto [lowResField, interpedField, lowResResidualField] = getResidualField(interp, binning);
    auto [highResField, binnedField, highResResidualField] = getResidualFieldHighToLowToHigh(interp, binning);

    targetFieldSet.add(lowResField);
    targetFieldSet.add(lowResResidualField);

    auto sourceFieldSet = FieldSet{};
    sourceFieldSet.add(sourceField);
    sourceFieldSet.add(interpedField);
    sourceFieldSet.add(highResResidualField);

    const auto binningType = binningScheme.getString("scheme.type");
    makeGmshOutput(sourceGridName + "_to_" + targetGridName + "_" + binningType + "_source.msh", sourceFieldSet);
    makeGmshOutput(sourceGridName + "_to_" + targetGridName + "_" + binningType + "_target.msh", targetFieldSet);
}

CASE("Regridding from high to low resolution: cubed sphere, bilinear") {
    const auto config = util::Config{"source_grid", "CS-LFR-140"} | util::Config{"target_grid", "CS-LFR-28"} |
                        util::Config{"functionspace", "NodeColumns"} | util::Config{"halo", 5} |
                        util::Config{"mesh_generator", "cubedsphere_dual"} |
                        util::Config{"scheme", option::type{"spherical-mean-value"}};

    regriddingTest(config);
}

CASE("Regridding from high to low resolution: cubed sphere, bilinear") {
    const auto config = util::Config{"source_grid", "CS-LFR-112"} | util::Config{"target_grid", "CS-LFR-28"} |
                        util::Config{"functionspace", "NodeColumns"} | util::Config{"halo", 3} |
                        util::Config{"mesh_generator", "cubedsphere_dual"} |
                        util::Config{"scheme", option::type{"spherical-mean-value"}};

    regriddingTest(config);
}

CASE("Regridding from high to low resolution: cubed sphere, bilinear") {
    const auto config = util::Config{"source_grid", "CS-LFR-84"} | util::Config{"target_grid", "CS-LFR-28"} |
                        util::Config{"functionspace", "NodeColumns"} | util::Config{"halo", 2} |
                        util::Config{"mesh_generator", "cubedsphere_dual"} |
                        util::Config{"scheme", option::type{"spherical-mean-value"}};

    regriddingTest(config);
}

CASE("Regridding from high to low resolution: cubed sphere, bilinear") {
    const auto config = util::Config{"source_grid", "CS-LFR-56"} | util::Config{"target_grid", "CS-LFR-28"} |
                        util::Config{"functionspace", "NodeColumns"} | util::Config{"halo", 1} |
                        util::Config{"mesh_generator", "cubedsphere_dual"} |
                        util::Config{"scheme", option::type{"spherical-mean-value"}};

    regriddingTest(config);
}


CASE("Regridding from high to low resolution: gaussian, spherical-mean-value") {
    const auto config = util::Config{"source_grid", "O96"} | util::Config{"target_grid", "O24"} |
                        util::Config{"functionspace", "NodeColumns"} | util::Config{"halo", 8} |
                        util::Config{"mesh_generator", "structured"} |
                        util::Config{"scheme", option::type{"spherical-mean-value"}};

    regriddingTest(config);
}

CASE("plot binning kernel") {
    const auto sourceGrid = Grid{"CS-LFR-35"};
    const auto targetGrid = Grid{"CS-LFR-5"};

    const auto sourceHaloOption = option::halo(3);
    const auto targetHaloOption = option::halo(0);

    const auto sourceMesh = MeshGenerator("cubedsphere_dual", sourceHaloOption).generate(sourceGrid);
    const auto targetMesh = MeshGenerator("cubedsphere_dual", targetHaloOption).generate(targetGrid);

    const auto sourceFunctionSpace = functionspace::NodeColumns{sourceMesh, sourceHaloOption};
    const auto targetFunctionSpace = functionspace::NodeColumns{targetMesh, targetHaloOption};

    const auto binningScheme =
        option::type{"local-pseudoinverse"} | util::Config{"scheme", option::type{"spherical-mean-value"}};
    const auto binning = Interpolation{binningScheme, sourceFunctionSpace, targetFunctionSpace};

    const auto binningMatrixStorage = interpolation::MatrixCache(binning).matrix();
    const auto binningMatrix        = linalg::make_non_owning_eckit_sparse_matrix(binningMatrixStorage);

    const auto targetGhostView = array::make_view<int, 1>(targetFunctionSpace.ghost());
    const auto targetRidxView  = array::make_indexview<idx_t, 1>(targetFunctionSpace.remote_index());
    const auto targetPartView  = array::make_view<int, 1>(targetFunctionSpace.partition());

    const auto targetRemoteIndices = std::vector{0, 2, 12};
    const auto targetPartition     = 0;

    auto kernelField     = sourceFunctionSpace.createField<double>(option::name("kernel_field"));
    auto kernelFieldView = array::make_view<double, 1>(kernelField);
    kernelFieldView.assign(0.);

    for (auto matrixItr = binningMatrix.begin(); matrixItr != binningMatrix.end(); ++matrixItr) {
        const auto rowIdx = matrixItr.row();
        const auto colIdx = matrixItr.col();

        if (targetGhostView(rowIdx)) {
            continue;
        }

        for (const auto targetRemoteIndex : targetRemoteIndices) {
            if (targetRidxView(rowIdx) == targetRemoteIndex && targetPartView(rowIdx) == targetPartition) {
                kernelFieldView(colIdx) += *matrixItr;
            }
        }
    }
    // Force halo weights into owned region of field.
    sourceFunctionSpace.adjointHaloExchange(kernelField);
    sourceFunctionSpace.haloExchange(kernelField);

    makeGmshOutput("binning_kernel.msh", kernelField);
}

// /// test to carry out the 'dot-product' test for the rigridding from
// /// 'high' to 'low' resolution, for a given type of grid (CS-LFR)
// ///
// CASE("dot-product test for the rigridding from high to low resolution; grid type: Cubed Sphere") {
//     // source grid (high res.)
//     const auto sourceGrid          = Grid("CS-LFR-100");
//     const auto sourceMesh          = MeshGenerator("cubedsphere_dual").generate(sourceGrid);
//     const auto sourceFunctionSpace = functionspace::NodeColumns(sourceMesh);

//     // target grid (low res.)
//     const auto targetGrid          = Grid("CS-LFR-50");
//     const auto targetMesh          = MeshGenerator("cubedsphere_dual").generate(targetGrid);
//     const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh);

//     // source field
//     const auto sourceField = getVortexField(sourceFunctionSpace);

//     auto sourceFieldSet = FieldSet{};
//     sourceFieldSet.add(sourceField);

//     // target field
//     auto targetField = targetFunctionSpace.createField<double>(option::name("field_01_t"));

//     auto targetFieldSet = FieldSet{};
//     targetFieldSet.add(targetField);

//     const auto scheme = util::Config("type", "local-pseudoinverse") | util::Config("scheme", option::type("spherical-mean-value")) |
//                         util::Config("adjoint", true);

//     Interpolation binning(scheme, sourceFunctionSpace, targetFunctionSpace);

//     // performing the regridding from high to low resolution
//     binning.execute(sourceFieldSet, targetFieldSet);


//     targetFieldSet["field_01_t"].haloExchange();

//     // target field (adjoint)
//     auto targetFieldStar = targetFunctionSpace.createField<double>(option::name("field_01_ad_t"));
//     array::make_view<double, 1>(targetFieldStar).assign(array::make_view<double, 1>(targetField));
//     targetFieldStar.adjointHaloExchange();

//     auto targetFieldSetStar = FieldSet{};
//     targetFieldSetStar.add(targetFieldStar);

//     // source field (adjoint)
//     auto sourceFieldStar = sourceFunctionSpace.createField<double>(option::name("field_01_ad_s"));
//     array::make_view<double, 1>(sourceFieldStar).assign(0.);

//     auto sourceFieldSetStar = FieldSet{};
//     sourceFieldSetStar.add(sourceFieldStar);

//     // performing adjoint operation
//     binning.execute_adjoint(sourceFieldSetStar, targetFieldSetStar);


//     const auto targetDotTarget     = dotProd(targetField, targetField);
//     const auto sourceDotSourceStar = dotProd(sourceField, sourceFieldStar);

//     double scaled_diff = std::abs(targetDotTarget - sourceDotSourceStar) / std::abs(targetDotTarget);

//     // carrrying out a dot-product test ...
//     Log::info() << "\n- dot-product test:\n"
//                 << "(Ax) . (Ax) = " << targetDotTarget << "; "
//                 << "x . (A^t A x) = " << sourceDotSourceStar << "; "
//                 << "scaled difference = " << scaled_diff << "\n"
//                 << std::endl;

//     EXPECT(scaled_diff < 1e-12);
// }

}  // namespace test
}  // namespace atlas


//--

int main(int argc, char** argv) {
    // std::this_thread::sleep_for(std::chrono::seconds(10));
    return atlas::test::run(argc, argv);
}
