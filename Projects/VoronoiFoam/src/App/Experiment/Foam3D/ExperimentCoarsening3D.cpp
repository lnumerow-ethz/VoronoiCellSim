#include "Projects/VoronoiFoam/include/App/Experiment/Foam3D/ExperimentCoarsening3D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario3D/RandomSitesInBox3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy3D/CellEnergyCoarsening3D.h"

#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellVolume3D.h"
#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"

#define TBB_SUPPRESS_DEPRECATED_MESSAGES (true)

#include <tbb/tbb.h>

#include <iostream>
#include <fstream>

void ExperimentCoarsening3D::setup(FoamSubApp* foam_sub_app) {
    foam_sub_app->optimizer = Optimization::LBFGS;
    foam_sub_app->tessellation_selector.first = 1;  /// Power diagram

    foam_sub_app->scenario_selector.first = 0;
    auto scenario = std::dynamic_pointer_cast<RandomSitesInBox3D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);
    scenario->num_sites = num_cells;
    scenario->dimensions_independent = true;
    std::fill(std::begin(scenario->dimensions_free), std::end(scenario->dimensions_free), false);
    std::fill(std::begin(scenario->dimension_defaults), std::end(scenario->dimension_defaults), 1.0);
    scenario->dimension_defaults[0] = 0.1;

    foam_sub_app->energy_selector.first = 2;  /// Coarsening energy
    auto energy = std::dynamic_pointer_cast<CellEnergyCoarsening3D>(
        foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
    energy->weights[VOLUME] = 0.002;
    energy->weights[SURFACE_TARGET] = 0.0;
    energy->weights[SURFACE_MINIMIZATION] = 0.01;
    energy->weights[CENTROID] = 0.1;
    energy->weights[SECOND_MOMENT] = 0.0;

    std::fill(std::begin(foam_sub_app->site_param_free), std::end(foam_sub_app->site_param_free), false);
    foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = false;
    foam_sub_app->site_param_free[SITE_PARAM_SIZE_TARGET] = false;  // Start with Voronoi to find initial equilibrium.

    foam_sub_app->use_dynamics = false;  // Start with no dynamics to find initial equilibrium.
    foam_sub_app->dynamics_data.dt = 1.0 / pow(num_cells, 1.0);
    foam_sub_app->dynamics_mass_pos = 0;
    std::fill(foam_sub_app->dynamics_mass_param.begin(), foam_sub_app->dynamics_mass_param.end(), 0);
    foam_sub_app->dynamics_mass_param[SITE_PARAM_SIZE_TARGET] = 0.0;
    foam_sub_app->dynamics_viscosity_pos = 0.0001;
    std::fill(foam_sub_app->dynamics_viscosity_param.begin(), foam_sub_app->dynamics_viscosity_param.end(), 0.001);
    foam_sub_app->dynamics_viscosity_param[SITE_PARAM_SIZE_TARGET] = 1.0;

    foam_sub_app->applyAllSettings();

    iter_start_time = std::chrono::high_resolution_clock::now();
}

void ExperimentCoarsening3D::loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) {
    if (optimization_status == Optimization::SUCCESS) {
        F avg_size = 0.8 / num_cells;

        int remaining_cells = 0;
        for (Site& site : foam_sub_app->degrees_of_freedom.sites) {
            if (site.param(SITE_PARAM_SIZE_TARGET) < 1e-3 * avg_size || site.is_removed) {
                site.is_removed = true;
            } else {
                remaining_cells++;
            }
        }
        num_cells = remaining_cells;

        num_frame_iters++;
    }
    if (optimization_status == Optimization::CONVERGED) {
        saveFrame(foam_sub_app);

        foam_sub_app->dynamics_data.dt = 1.0 / pow(num_cells, 1.0);

        if (!foam_sub_app->use_dynamics) {
            foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = true;
            foam_sub_app->site_param_free[SITE_PARAM_SIZE_TARGET] = true;
            foam_sub_app->use_dynamics = true;
            foam_sub_app->optimize = true;
            foam_sub_app->optimizer = Optimization::NEWTON;
            foam_sub_app->clearState();
        }

        if (frame > 500) {
            foam_sub_app->optimize = false;
        }
    }
    if (optimization_status == Optimization::FAILURE) {
        foam_sub_app->optimize = true;
    }
}

void ExperimentCoarsening3D::saveFrame(FoamSubApp* foam_sub_app) {
    auto iter_end_time = std::chrono::high_resolution_clock::now();
    auto iter_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(iter_end_time - iter_start_time);

    num_frame_iters++;

    foam_sub_app->write_png = true;
    foam_sub_app->png_file_name = "coarsening_frame_" + std::to_string(frame) + ".png";
    foam_sub_app->write_rendering_file = true;
    foam_sub_app->render_file_name = "coarsening_frame_render_" + std::to_string(frame) + ".txt";

    std::string frame_info_file_name = "coarsening_frame_info_" + std::to_string(frame) + ".txt";
    std::ofstream info_file(frame_info_file_name);

    VectorXF dof;
    ModelHelper::getDOFVector(foam_sub_app->model_definition, foam_sub_app->degrees_of_freedom, dof);

    info_file << frame << " " << iter_time_ms.count() << " " << num_frame_iters << "\n";
    for (int i = 0; i < dof.rows(); i++) {
        info_file << dof(i) << " ";
    }
    info_file << "\n";
    info_file.close();

    iter_start_time = iter_end_time;
    num_frame_iters = 0;
    frame++;
}
