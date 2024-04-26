#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentConvergenceTest2D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/ConvergenceTest2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergy2D.h"

#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"

#include <iostream>
#include <iomanip>

ExperimentConvergenceTest2D::ExperimentConvergenceTest2D() {
    F area = 0.00;
    F perimeter = 0.0008;
    F centroid = 0.02;
    F initial_dt = 0.007;
    std::vector<int> div = {1, 2, 3, 4, 6, 8, 12, 16, 24, 32};

    // initial_dt = 0.008;
    for (int a : div) {
        run_params.emplace_back(0, FIRST_ORDER, 0, perimeter, 0, 0, initial_dt / a, 2 * a);
    }
    for (int a : div) {
        run_params.emplace_back(0, FIRST_ORDER, area, 0, centroid, 0, initial_dt / a, 2 * a);
    }
    // initial_dt = 0.08;
    for (int a : div) {
        run_params.emplace_back(0, SECOND_ORDER, 0, perimeter, 0, 0, initial_dt / a, 2 * a);
    }
    for (int a : div) {
        run_params.emplace_back(0, SECOND_ORDER, area, 0, centroid, 0, initial_dt / a, 2 * a);
    }

    output_string << std::setprecision(15);
}

void ExperimentConvergenceTest2D::setup(FoamSubApp* foam_sub_app) {
    const ConvergenceTest2DParams& params = run_params[run_index];

    foam_sub_app->tessellation_selector.first = params.tessellation;

    foam_sub_app->scenario_selector.first = 2;
    auto scenario = std::dynamic_pointer_cast<ConvergenceTest2D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);

    foam_sub_app->energy_selector.first = 0;
    auto energy = std::dynamic_pointer_cast<CellEnergy2D>(
        foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
    energy->weights[AREA] = params.area_weight;
    energy->weights[PERIMETER_TARGET] = params.perimeter_weight;
    energy->weights[PERIMETER_MINIMIZATION] = 0.0;
    energy->weights[CENTROID] = params.centroid_weight;
    energy->weights[SECOND_MOMENT] = params.second_moment_weight;

    std::fill(std::begin(foam_sub_app->site_param_free), std::end(foam_sub_app->site_param_free), false);
    foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = true;

    F mass = 0.01;
    F visc = 0.01;

    foam_sub_app->use_dynamics = true;
    foam_sub_app->dynamics_data.discretization = params.discretization;
    foam_sub_app->dynamics_data.tolerance = 1e-14;
    foam_sub_app->dynamics_data.dt = params.time_step;
    foam_sub_app->dynamics_mass_pos = mass;
    std::fill(foam_sub_app->dynamics_mass_param.begin(), foam_sub_app->dynamics_mass_param.end(), mass);
    foam_sub_app->dynamics_mass_param[SITE_PARAM_SIZE_TARGET] = mass;
    foam_sub_app->dynamics_viscosity_pos = visc;
    std::fill(foam_sub_app->dynamics_viscosity_param.begin(), foam_sub_app->dynamics_viscosity_param.end(), visc);
    foam_sub_app->dynamics_viscosity_param[SITE_PARAM_SIZE_TARGET] = visc;

    foam_sub_app->applyAllSettings();
}

void ExperimentConvergenceTest2D::loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) {
    if (optimization_status == Optimization::CONVERGED) {
        current_time_step++;
        if (current_time_step == run_params[run_index].num_time_steps) {
            VectorXF dof;
            ModelHelper::getDOFVector(foam_sub_app->model_definition, foam_sub_app->degrees_of_freedom, dof);

            for (int i = 0; i < dof.rows(); i++) {
                output_string << dof(i) << " ";
            }
            output_string << std::endl;

            if (run_index == (int)run_params.size() - 1) {
                std::cout << "All runs finished." << std::endl;
                foam_sub_app->optimize = false;

                std::ofstream output_file;
                output_file.open("convergence_test.csv");
                output_file << output_string.rdbuf();
                output_file.close();
            } else {
                std::cout << "Run " << run_index << " finished." << std::endl;
                run_index++;
                current_time_step = 0;

                setup(foam_sub_app);
                foam_sub_app->applyAllSettings();
            }
        }
    }
}
