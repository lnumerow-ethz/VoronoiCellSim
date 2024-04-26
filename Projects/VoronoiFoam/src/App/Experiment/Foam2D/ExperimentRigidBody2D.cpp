#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentRigidBody2D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RigidBody2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergy2D.h"

#include <iostream>

void ExperimentRigidBody2D::setup(FoamSubApp* foam_sub_app) {
    foam_sub_app->tessellation_selector.first = 0;  /// Start with Voronoi diagram

    foam_sub_app->scenario_selector.first = 3;  /// Rigid Body
    auto scenario = std::dynamic_pointer_cast<RigidBody2D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);
    scenario->num_sites = 2000;

    foam_sub_app->energy_selector.first = 0;
    auto energy = std::dynamic_pointer_cast<CellEnergy2D>(
        foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
    energy->weights[AREA] = 0.05;
    energy->weights[PERIMETER_TARGET] = 0;
    energy->weights[PERIMETER_MINIMIZATION] = 0;
    energy->weights[CENTROID] = 1.0;
    energy->weights[SECOND_MOMENT] = 100.0;

    std::fill(std::begin(foam_sub_app->site_param_free), std::end(foam_sub_app->site_param_free), false);
    foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = true;

    foam_sub_app->use_dynamics = false;  // Start with no dynamics to find initial equilibrium.
    foam_sub_app->dynamics_data.dt = 0.01;
    foam_sub_app->dynamics_data.tolerance = 1e-10;
    foam_sub_app->dynamics_mass_pos = 0.0003;
    foam_sub_app->dynamics_mass_boundary = 0.0003;
    foam_sub_app->dynamics_viscosity_pos = 0.0003;
    foam_sub_app->dynamics_viscosity_boundary = 0.003;

    foam_sub_app->applyAllSettings();

    foam_sub_app->model_definition.boundary_free_param_indices =
        VectorXI::Zero(0);  // Fix rigid body to find initial equilibrium.

    iter_start_time = std::chrono::high_resolution_clock::now();
}

void ExperimentRigidBody2D::loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) {
    if (optimization_status == Optimization::SUCCESS) {
        num_frame_iters++;
    }
    if (optimization_status == Optimization::CONVERGED) {
        if (!pre_optimized_0) {
            // foam_sub_app->tessellation_selector.first = 1;  /// Switch to Power diagram.
            foam_sub_app->use_dynamics = true;
            foam_sub_app->optimize = true;  // Turn optimization back on.

            foam_sub_app->model_definition.boundary_free_param_indices = Vector3I(0, 1, 2);
            foam_sub_app->clearState();

            pre_optimized_0 = true;
        }

        auto iter_end_time = std::chrono::high_resolution_clock::now();
        auto iter_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(iter_end_time - iter_start_time);

        num_frame_iters++;

        foam_sub_app->write_png = true;
        foam_sub_app->png_file_name = "coarsening_frame_" + std::to_string(frame) + ".png";
        foam_sub_app->write_rendering_file = true;
        foam_sub_app->render_file_name = "coarsening_frame_render_" + std::to_string(frame) + ".txt";

        std::string frame_info_file_name = "rigid_body_frame_info_" + std::to_string(frame) + ".txt";
        std::ofstream info_file(frame_info_file_name);

        info_file << frame << " " << iter_time_ms.count() << " " << num_frame_iters << "\n";
        info_file.close();

        foam_sub_app->write_png = true;
        foam_sub_app->png_file_name = "rigid_body_frame_" + std::to_string(frame) + ".png";
        foam_sub_app->write_rendering_file = true;
        foam_sub_app->render_file_name = "render_rigid_body_frame_" + std::to_string(frame) + ".txt";

        iter_start_time = iter_end_time;
        num_frame_iters = 0;
        frame++;
    }
}
