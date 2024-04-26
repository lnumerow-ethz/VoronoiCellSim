#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentDCMComparison2D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInBox2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergy2D.h"

#include <iostream>
#include <fstream>

void ExperimentDCMComparison2D::setup(FoamSubApp* foam_sub_app) {
    foam_sub_app->tessellation_selector.first = 0;  /// Start with Voronoi diagram

    foam_sub_app->scenario_selector.first = 0;  /// Sites in Box
    auto scenario = std::dynamic_pointer_cast<RandomSitesInBox2D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);
    scenario->num_sites = 30;

    foam_sub_app->energy_selector.first = 0;
    auto energy = std::dynamic_pointer_cast<CellEnergy2D>(
        foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
    energy->weights[AREA] = 5.0;
    energy->weights[PERIMETER_TARGET] = 0.01;
    energy->weights[PERIMETER_MINIMIZATION] = 0;
    energy->weights[CENTROID] = 0;
    energy->weights[SECOND_MOMENT] = 0;

    std::fill(std::begin(foam_sub_app->site_param_free), std::end(foam_sub_app->site_param_free), false);
    foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = true;

    foam_sub_app->use_dynamics = false;
    foam_sub_app->applyAllSettings();
}

void ExperimentDCMComparison2D::loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) {
    if (optimization_status == Optimization::SUCCESS && stage == 1) {
        num_frame_iters++;
    }
    if (optimization_status == Optimization::CONVERGED) {
        auto energy = std::dynamic_pointer_cast<CellEnergy2D>(
            foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
        switch (stage) {
            case 0:
                foam_sub_app->write_rendering_file = true;
                foam_sub_app->render_file_name = "cells.txt";

                foam_sub_app->tessellation_selector.first = 1;  // Switch to power diagram
                energy->weights[PERIMETER_TARGET] = 0.012;
                energy->weights[CENTROID] = 0.001;     // Add small centroid regularizer
                energy->weights[SECOND_MOMENT] = 1.5;  // Add small centroid regularizer

                foam_sub_app->clearState();
                foam_sub_app->optimize = true;  // Turn optimization back on.

                stage = 1;
                break;
            case 1:
                foam_sub_app->write_png = true;
                foam_sub_app->png_file_name = "dcm_comparison_frame_" + std::to_string(frame) + ".png";
                foam_sub_app->write_rendering_file = true;
                foam_sub_app->render_file_name = "render_dcm_comparison_frame_" + std::to_string(frame) + ".txt";

                {
                    auto iter_end_time = std::chrono::high_resolution_clock::now();
                    auto iter_time_ms =
                        std::chrono::duration_cast<std::chrono::milliseconds>(iter_end_time - iter_start_time);
                    num_frame_iters++;

                    std::string frame_info_file_name = "dcm_comparison_frame_info_" + std::to_string(frame) + ".txt";
                    std::ofstream info_file(frame_info_file_name);

                    info_file << frame << " " << iter_time_ms.count() << " " << num_frame_iters << "\n";
                    info_file.close();

                    iter_start_time = iter_end_time;
                    num_frame_iters = 0;
                }

                if (frame < 100) {
                    frame++;
                    foam_sub_app->degrees_of_freedom.boundary_param =
                        Vector2F(1.0 + 0.005 * frame, 1.0 / (1.0 + 0.005 * frame));

                    foam_sub_app->optimize = true;  // Turn optimization back on.
                }
                break;
            default:
                assert(0);
        }
    }
}
