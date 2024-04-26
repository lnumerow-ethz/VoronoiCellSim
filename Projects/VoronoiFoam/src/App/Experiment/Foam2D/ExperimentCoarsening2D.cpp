#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentCoarsening2D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInBox2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergyCoarsening2D.h"

#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellArea2D.h"

#define TBB_SUPPRESS_DEPRECATED_MESSAGES (true)

#include <tbb/tbb.h>

#include <iostream>

void ExperimentCoarsening2D::setup(FoamSubApp* foam_sub_app) {
    foam_sub_app->tessellation_selector.first = 1;  /// Power diagram

    foam_sub_app->scenario_selector.first = 0;
    auto scenario = std::dynamic_pointer_cast<RandomSitesInBox2D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);
    scenario->num_sites = num_cells;
    scenario->dimensions_independent = false;
    std::fill(std::begin(scenario->dimensions_free), std::end(scenario->dimensions_free), false);
    std::fill(std::begin(scenario->dimension_defaults), std::end(scenario->dimension_defaults), 1.0);

    foam_sub_app->energy_selector.first = 2;  /// Coarsening energy
    auto energy = std::dynamic_pointer_cast<CellEnergyCoarsening2D>(
        foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
    energy->weights[AREA] = 0.001;
    energy->weights[PERIMETER_TARGET] = 0.0;
    energy->weights[PERIMETER_MINIMIZATION] = 0.001;
    energy->weights[CENTROID] = 0.1;
    energy->weights[SECOND_MOMENT] = 0.0;

    std::fill(std::begin(foam_sub_app->site_param_free), std::end(foam_sub_app->site_param_free), false);
    foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = true;
    foam_sub_app->site_param_free[SITE_PARAM_SIZE_TARGET] =
        false;  // Start with no dynamics to find initial equilibrium.

    foam_sub_app->use_dynamics = false;  // Start with no dynamics to find initial equilibrium.
    foam_sub_app->dynamics_data.dt = 3.0 / pow(num_cells, 1.0);
    foam_sub_app->dynamics_mass_pos = 0;
    std::fill(foam_sub_app->dynamics_mass_param.begin(), foam_sub_app->dynamics_mass_param.end(), 0);
    foam_sub_app->dynamics_mass_param[SITE_PARAM_SIZE_TARGET] = 0.00;
    foam_sub_app->dynamics_viscosity_pos = 0;
    std::fill(foam_sub_app->dynamics_viscosity_param.begin(), foam_sub_app->dynamics_viscosity_param.end(), 0);
    foam_sub_app->dynamics_viscosity_param[SITE_PARAM_SIZE_TARGET] = 0.2;

    foam_sub_app->applyAllSettings();
}

void ExperimentCoarsening2D::loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) {
    if (optimization_status == Optimization::SUCCESS) {
        F avg_size = 4.0 / num_cells;

        // bool success;
        // Model model(foam_sub_app->model_definition, foam_sub_app->degrees_of_freedom, 0, success);
        // assert(success);
        // PerCellArea2D cell_area_function;
        // tbb::concurrent_vector<int> remaining_cells;
        // tbb::parallel_for_each(model.cells.begin(), model.cells.end(), [&](const TessellationCell& cell) {
        //     PerCellValue cell_area(model, cell, 0);
        //     cell_area_function.getValue(model, cell_area);

        //     if (cell_area.value < 1e-3 * avg_size) {
        //         foam_sub_app->degrees_of_freedom.sites[cell.site_index].is_removed = true;
        //     } else {
        //         remaining_cells.push_back(cell.site_index);
        //     }
        // });
        // num_cells = remaining_cells.size();

        int remaining_cells = 0;
        for (Site& site : foam_sub_app->degrees_of_freedom.sites) {
            if (site.param(SITE_PARAM_SIZE_TARGET) < 1e-3 * avg_size) {
                site.is_removed = true;
            } else {
                remaining_cells++;
            }
        }
        num_cells = remaining_cells;
    }
    if (optimization_status == Optimization::CONVERGED) {
        foam_sub_app->dynamics_data.dt = 3.0 / pow(num_cells, 1.0);

        if (!foam_sub_app->use_dynamics) {
            foam_sub_app->site_param_free[SITE_PARAM_SIZE_TARGET] = true;
            foam_sub_app->use_dynamics = true;
            foam_sub_app->optimize = true;
            foam_sub_app->clearState();
        }
    }
    if (optimization_status == Optimization::FAILURE) {
        foam_sub_app->optimize = true;
    }
}
