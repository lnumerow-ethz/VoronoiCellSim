#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentRuntimeScaling2D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInBox2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergy2D.h"

#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"

#define TBB_SUPPRESS_DEPRECATED_MESSAGES (true)

#include <tbb/tbb.h>

#include <iostream>
#include <fstream>

void ExperimentRuntimeScaling2D::setup(FoamSubApp* foam_sub_app) {
    foam_sub_app->optimizer = Optimization::NEWTON;
    foam_sub_app->tessellation_selector.first = 1;  /// Power diagram

    foam_sub_app->scenario_selector.first = 0;
    auto scenario = std::dynamic_pointer_cast<RandomSitesInBox2D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);
    scenario->dimensions_independent = false;
    std::fill(std::begin(scenario->dimensions_free), std::end(scenario->dimensions_free), false);

    foam_sub_app->energy_selector.first = 0;  /// Normal energy
    auto energy = std::dynamic_pointer_cast<CellEnergy2D>(
        foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
    energy->weights[AREA] = 1.0;
    energy->weights[PERIMETER_TARGET] = 0.0;
    energy->weights[PERIMETER_MINIMIZATION] = 0.0;
    energy->weights[CENTROID] = 0.1;
    energy->weights[SECOND_MOMENT] = 10.0;

    std::fill(std::begin(foam_sub_app->site_param_free), std::end(foam_sub_app->site_param_free), false);
    foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = true;

    foam_sub_app->use_dynamics = false;

    scenario->num_sites = std::round(num_sites_float);
    foam_sub_app->applyAllSettings();

    start_time = std::chrono::high_resolution_clock::now();
}

void ExperimentRuntimeScaling2D::loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    auto scenario = std::dynamic_pointer_cast<RandomSitesInBox2D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);

    std::cout << "Runtime for " << scenario->num_sites << " cells: " << runtime.count() << " microseconds."
              << std::endl;

    num_sites_float *= pow(10.0, 1.0 / 10.0);
    scenario->num_sites = std::round(num_sites_float);
    foam_sub_app->applyAllSettings();

    if (num_sites_float > 100001) {
        foam_sub_app->optimize = false;
    }
    start_time = std::chrono::high_resolution_clock::now();
}
