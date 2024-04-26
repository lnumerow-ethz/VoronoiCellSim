#include "Projects/VoronoiFoam/include/App/Experiment/Foam3D/ExperimentCleavage3D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario3D/RandomSitesInMembrane3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy3D/CellEnergyCleavage3D.h"

#include <fstream>

void ExperimentCleavage3D::setup(FoamSubApp* foam_sub_app) {
    foam_sub_app->tessellation_selector.first = 1;  /// Power diagram

    foam_sub_app->scenario_selector.first = 1;  /// Membrane
    auto scenario = std::dynamic_pointer_cast<RandomSitesInMembrane3D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);
    scenario->num_sites = 1;
    scenario->subdivision_level = 4;
    scenario->spring_constant = 0.3;

    foam_sub_app->energy_selector.first = 1;  /// Nondimensionalized energy
    auto energy = std::dynamic_pointer_cast<CellEnergyCleavage3D>(
        foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
    energy->weights[VOLUME] = 50.0;
    energy->weights[CENTROID] = 1.0;
    energy->weights[SECOND_MOMENT] = 5.0;
    energy->weights[GRAVITY] = 0.0;

    std::fill(std::begin(foam_sub_app->site_param_free), std::end(foam_sub_app->site_param_free), false);
    foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = true;

    foam_sub_app->use_dynamics = false;  // Start with no dynamics to find initial equilibrium.
    foam_sub_app->dynamics_data.dt = 0.01;
    foam_sub_app->dynamics_data.tolerance = 1e-10;
    foam_sub_app->dynamics_mass_pos = 0;
    std::fill(foam_sub_app->dynamics_mass_param.begin(), foam_sub_app->dynamics_mass_param.end(), 0);
    foam_sub_app->dynamics_mass_param[SITE_PARAM_SIZE_TARGET] = 0.0;
    foam_sub_app->dynamics_viscosity_pos = 5.0;
    std::fill(foam_sub_app->dynamics_viscosity_param.begin(), foam_sub_app->dynamics_viscosity_param.end(), 5.0);
    foam_sub_app->dynamics_viscosity_boundary = 0.5;

    foam_sub_app->applyAllSettings();

    iter_start_time = std::chrono::high_resolution_clock::now();
}

void ExperimentCleavage3D::loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) {
    if (optimization_status == Optimization::SUCCESS) {
        num_frame_iters++;
    }
    if (optimization_status == Optimization::CONVERGED) {
        foam_sub_app->use_dynamics = true;  // Enable dynamics after convergence to static equilibrium.
        foam_sub_app->optimize = true;      // Turn optimization back on.

        int curr_n_cells = foam_sub_app->degrees_of_freedom.sites.size();
        for (int i = 0; i < curr_n_cells; i++) {
            if (growing_cells.find(i) == growing_cells.end()) {
                std::uniform_real_distribution<F> dist(0, 1);
                F rando = dist(rng);
                F size_target = foam_sub_app->degrees_of_freedom.sites[i].param(SITE_PARAM_SIZE_TARGET);
                // if (rando < 0.003 * size_target || (curr_n_cells == 1 && frame == 1)) {
                if (frame % 60 == 1) {
                    VectorXF rando_displacement(3);

                    rando_displacement(0) = (dist(rng) - 0.5);
                    rando_displacement(1) = (dist(rng) - 0.5);
                    rando_displacement(2) = (dist(rng) - 0.5);
                    rando_displacement = rando_displacement.normalized() * 0.1 * pow(size_target, 1.0 / 3.0);

                    if (frame == 1) {
                        rando_displacement = Vector3F(1, 0, 0) * 0.1 * pow(size_target, 1.0 / 3.0);
                    } else if (frame == 61) {
                        rando_displacement =
                            Vector3F(0.1 * dist(rng), 1, 0.1 * dist(rng)) * 0.1 * pow(size_target, 1.0 / 3.0);
                    } else if (frame == 121) {
                        rando_displacement =
                            Vector3F(0.1 * dist(rng), 0.1 * dist(rng), 1) * 0.1 * pow(size_target, 1.0 / 3.0);
                    }

                    foam_sub_app->degrees_of_freedom.sites.emplace_back();

                    Site& new_site = foam_sub_app->degrees_of_freedom.sites.back();

                    VectorXF old_site_pos = foam_sub_app->degrees_of_freedom.sites[i].pos;
                    new_site.pos = old_site_pos + rando_displacement;
                    foam_sub_app->degrees_of_freedom.sites[i].pos = old_site_pos - rando_displacement;

                    new_site.param(SITE_PARAM_POWER_WEIGHT) =
                        foam_sub_app->degrees_of_freedom.sites[i].param(SITE_PARAM_POWER_WEIGHT);

                    new_site.param(SITE_PARAM_SIZE_TARGET) = 0.5 * size_target;
                    foam_sub_app->degrees_of_freedom.sites[i].param(SITE_PARAM_SIZE_TARGET) = 0.5 * size_target;

                    int new_site_index = foam_sub_app->degrees_of_freedom.sites.size() - 1;
                    growing_cells[i] = 0.5 * size_target;
                    growing_cells[new_site_index] = 0.5 * size_target;
                }
            } else {
                Site& site = foam_sub_app->degrees_of_freedom.sites[i];
                F growth_target = growing_cells[i];

                if (site.param(SITE_PARAM_SIZE_TARGET) < growth_target) {
                    site.param(SITE_PARAM_SIZE_TARGET) =
                        std::min(site.param(SITE_PARAM_SIZE_TARGET) + 0.07 * growth_target, growth_target);
                } else {
                    growing_cells.erase(i);
                }
            }
        }
        foam_sub_app->clearState();
        foam_sub_app->write_png = true;
        foam_sub_app->png_file_name = "cleavage_membrane_frame_" + std::to_string(frame) + ".png";
        foam_sub_app->write_rendering_file = true;
        foam_sub_app->render_file_name = "render_cleavage_membrane_frame_" + std::to_string(frame) + ".txt";

        auto iter_end_time = std::chrono::high_resolution_clock::now();
        auto iter_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(iter_end_time - iter_start_time);

        num_frame_iters++;

        std::string frame_info_file_name = "cleavage_membrane_frame_info_" + std::to_string(frame) + ".txt";
        std::ofstream info_file(frame_info_file_name);

        info_file << frame << " " << iter_time_ms.count() << " " << num_frame_iters << "\n";
        info_file.close();

        frame++;
    }
}
