#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"
#include "Projects/VoronoiFoam/include/Model/Energy/PotentialEnergy.h"
#include "Projects/VoronoiFoam/include/Model/Energy/DynamicsObjective.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryGenerator.h"
#include "Projects/VoronoiFoam/include/Model/Energy/PerCellFunction.h"
#include "Projects/VoronoiFoam/include/App/Rendering/Rendering.h"

#include "CRLHelper/Optimization.h"

#include <iostream>

void FoamSubApp::resetOptimization() {
    int dims_space = model_definition.boundary_generator->getDims();
    int dims_site_params_free = model_definition.site_free_param_indices.rows();
    int dims_site_dof = dims_space + dims_site_params_free;
    int num_dof =
        (int)(degrees_of_freedom.sites.size() * dims_site_dof + model_definition.boundary_free_param_indices.rows());

    /// Reset BGFS Hessian to initial guess.
    // bfgs_inverse_hessian = MatrixXF::Identity(num_dof, num_dof);
    lbfgs_mat.reset(num_dof, lbfgs_m);
}

void FoamSubApp::clearState() {
    /// Set tessellation.
    model_definition.tessellation_generator = tessellation_selector.second[tessellation_selector.first];

    /// Set cell energy function.
    model_definition.cell_energy_function = energy_selector.second[energy_selector.first].second;

    /// Set free site parameters.
    ModelHelper::setFreeSiteParam(model_definition, site_param_free);

    resetOptimization();

    /// Reset Dynamics velocities to zero (previous states = current state).
    ModelHelper::getDOFVector(model_definition, degrees_of_freedom, dynamics_data.prev_dof_1);
    dynamics_data.prev_dof_2 = dynamics_data.prev_dof_1;
    dynamics_data.prev_dof_3 = dynamics_data.prev_dof_1;
    setDynamicsMatrices();
}

void FoamSubApp::applyAllSettings() {
    scenario_selector.second[scenario_selector.first]->assignScenario(model_definition, degrees_of_freedom);
    clearState();
}

void FoamSubApp::setDynamicsMatrices() {
    int dims_space = model_definition.boundary_generator->getDims();
    int dims_site_params_free = model_definition.site_free_param_indices.rows();
    int dims_c = dims_space + dims_site_params_free;
    int n_sites = degrees_of_freedom.sites.size();
    int nc = degrees_of_freedom.sites.size() * dims_c;
    int np = model_definition.boundary_free_param_indices.rows();
    int num_dof = nc + np;

    dynamics_data.mass_matrix_diagonal = VectorXF::Zero(num_dof);
    dynamics_data.viscosity_matrix_diagonal = VectorXF::Zero(num_dof);
    for (int i = 0; i < n_sites; i++) {
        dynamics_data.mass_matrix_diagonal.segment(i * dims_c, dims_space).setConstant(dynamics_mass_pos);
        for (int j = 0; j < dims_site_params_free; j++) {
            dynamics_data.mass_matrix_diagonal(i * dims_c + dims_space + j) =
                dynamics_mass_param[model_definition.site_free_param_indices(j)];
        }

        dynamics_data.viscosity_matrix_diagonal.segment(i * dims_c, dims_space).setConstant(dynamics_viscosity_pos);
        for (int j = 0; j < dims_site_params_free; j++) {
            dynamics_data.viscosity_matrix_diagonal(i * dims_c + dims_space + j) =
                dynamics_viscosity_param[model_definition.site_free_param_indices(j)];
        }
    }
    dynamics_data.mass_matrix_diagonal.tail(np).setConstant(dynamics_mass_boundary);
    dynamics_data.viscosity_matrix_diagonal.tail(np).setConstant(dynamics_viscosity_boundary);
}

void FoamSubApp::initializeSubApp() {
    applyAllSettings();
}

void FoamSubApp::optimizationConverged(const VectorXF &dof_vector) {
    std::cout << "Optimization converged." << std::endl;
    if (use_dynamics) {
        DynamicsObjective::nextDynamicsStep(dynamics_data, dof_vector);
        resetOptimization();
    } else {
        optimize = false;
    }
}

Optimization::OptimizationStatus FoamSubApp::energyMinimizationStep() {
    auto objective_function = [&](const VectorXF &y, F &energy) {
        return DynamicsObjective::computePotentialEnergyOrDynamicsObjective(model_definition, degrees_of_freedom,
                                                                            dynamics_data, use_dynamics, y, energy);
    };

    auto gradient_function = [&](const VectorXF &y, F &energy, VectorXF &gradient) {
        return DynamicsObjective::computePotentialEnergyOrDynamicsObjective(
            model_definition, degrees_of_freedom, dynamics_data, use_dynamics, y, energy, gradient);
    };

    auto hessian_function = [&](const VectorXF &y, F &energy, VectorXF &gradient, HessianF &hessian) {
        return DynamicsObjective::computePotentialEnergyOrDynamicsObjective(
            model_definition, degrees_of_freedom, dynamics_data, use_dynamics, y, energy, gradient, hessian);
    };

    VectorXF dof_vector;
    ModelHelper::getDOFVector(model_definition, degrees_of_freedom, dof_vector);
    Optimization::OptimizationStatus status;
    F tolerance = use_dynamics ? dynamics_data.tolerance : 1e-14;
    switch (optimizer) {
        case Optimization::GRADIENT_DESCENT:
            status = Optimization::stepGradientDescent(dof_vector, objective_function, gradient_function, tolerance);
            break;
        case Optimization::NEWTON:
            status = Optimization::stepNewton(dof_vector, objective_function, hessian_function, false, tolerance);
            break;
        case Optimization::NEWTON_WOODBURY:
            status = Optimization::stepNewton(dof_vector, objective_function, hessian_function, true, tolerance);
            break;
        // case Optimization::BFGS:
        //     status = Optimization::stepBFGS(dof_vector, bfgs_inverse_hessian, objective_function, gradient_function,
        //                                     tolerance);
        //     break;
        case Optimization::LBFGS:
            status = Optimization::stepLBFGS(dof_vector, lbfgs_mat, objective_function, gradient_function, tolerance);
            break;
        default:
            assert(0);
    }
    ModelHelper::setDOFFromVector(model_definition, degrees_of_freedom, dof_vector);

    switch (status) {
        case Optimization::FAILURE:
            std::cout << "Optimization step failed." << std::endl;
            optimize = false;
            break;
        case Optimization::CONVERGED:
            optimizationConverged(dof_vector);
            break;
        case Optimization::SUCCESS:
            std::cout << "Optimization step success." << std::endl;
            break;
    }

    return status;
}

void FoamSubApp::checkEnergyGradients(int order, F epsilon, bool print_all) {
    auto objective_function = [&](const VectorXF &y, F &energy) {
        return DynamicsObjective::computePotentialEnergyOrDynamicsObjective(model_definition, degrees_of_freedom,
                                                                            dynamics_data, use_dynamics, y, energy);
    };

    auto gradient_function = [&](const VectorXF &y, VectorXF &gradient) {
        F energy;
        return DynamicsObjective::computePotentialEnergyOrDynamicsObjective(
            model_definition, degrees_of_freedom, dynamics_data, use_dynamics, y, energy, gradient);
    };

    auto hessian_function = [&](const VectorXF &y, MatrixXF &hessian_dense) {
        F energy;
        VectorXF gradient;
        HessianF hessian;
        bool success = DynamicsObjective::computePotentialEnergyOrDynamicsObjective(
            model_definition, degrees_of_freedom, dynamics_data, use_dynamics, y, energy, gradient, hessian);
        if (success) {
            hessian_dense = hessian.evalDense();
        }
        return success;
    };

    VectorXF dof_vector;
    ModelHelper::getDOFVector(model_definition, degrees_of_freedom, dof_vector);

    int print_level = print_all ? 2 : 1;
    if (order == 1) {
        VectorXF gradient_error;
        Optimization::checkGradient(dof_vector, gradient_error, objective_function, gradient_function, epsilon,
                                    print_level);
    } else if (order == 2) {
        MatrixXF hessian_error;
        Optimization::checkHessian(dof_vector, hessian_error, gradient_function, hessian_function, epsilon,
                                   print_level);
    }
}

void FoamSubApp::mainLoop() {
    if (optimize) {
        Optimization::OptimizationStatus status = energyMinimizationStep();
        experiment_selector.second[experiment_selector.first]->loop(this, status);
        generateRenderOutput();
    }
}

void FoamSubApp::generateRenderOutput() {
    if (write_rendering_file) {
        bool model_generation_success;
        Model model(model_definition, degrees_of_freedom, 0, model_generation_success);
        if (!model_generation_success) {
            std::cout << "Broken model encountered in FoamSubApp3D::getViewerData." << std::endl;
            assert(0 && "Broken model encountered in FoamSubApp3D::getViewerData.");
        }

        Rendering::writeFile(model, render_file_name);
        write_rendering_file = false;
    }
}

void FoamSubApp::runExperiment(int experiment) {
    experiment_selector.first = experiment;
    experiment_selector.second[experiment_selector.first]->setup(this);
    optimize = true;
}

void FoamSubApp::makeConfigWindow() {
    if (ImGui::Combo("Optimizer", reinterpret_cast<int *>(&optimizer), Optimization::optimizer_names)) {
        clearState();
    }
    if (ImGui::Checkbox("Optimize", &optimize)) {
        clearState();
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Experiments", ImGuiTreeNodeFlags_None)) {
        std::vector<std::string> experiment_names;
        for (const auto &experiment : experiment_selector.second) {
            experiment_names.push_back(experiment->getName());
        }
        if (ImGui::Combo("Experiment##Selector", &experiment_selector.first, experiment_names)) {
            experiment_selector.second[experiment_selector.first]->setup(this);
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Scenario", ImGuiTreeNodeFlags_DefaultOpen)) {
        std::vector<std::string> scenario_names;
        for (const auto &scenario : scenario_selector.second) {
            scenario_names.push_back(scenario->getName());
        }
        if (ImGui::Combo("Scenario##Selector", &scenario_selector.first, scenario_names)) {
            applyAllSettings();
        }
        scenario_selector.second[scenario_selector.first]->makeConfigMenu();
        if (ImGui::Button("Regenerate")) {
            applyAllSettings();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Tessellation / DOF", ImGuiTreeNodeFlags_DefaultOpen)) {
        std::vector<std::string> tessellation_names;
        for (const auto &tessellation : tessellation_selector.second) {
            tessellation_names.push_back(tessellation->getName());
        }
        if (ImGui::Combo("Tessellation##Selector", &tessellation_selector.first, tessellation_names)) {
            clearState();
        }
        bool site_free_param_changed = false;
        if (model_definition.tessellation_generator->hasSiteParam(SITE_PARAM_POWER_WEIGHT))
            site_free_param_changed |=
                ImGui::Checkbox("Power weight##free", &(site_param_free[SITE_PARAM_POWER_WEIGHT]));
        site_free_param_changed |= ImGui::Checkbox("Size target##free", &(site_param_free[SITE_PARAM_SIZE_TARGET]));
        site_free_param_changed |=
            ImGui::Checkbox("Surface target##free", &(site_param_free[SITE_PARAM_SURFACE_TARGET]));
        if (site_free_param_changed) {
            clearState();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Energy", ImGuiTreeNodeFlags_DefaultOpen)) {
        std::vector<std::string> energy_names;
        for (const auto &energy : energy_selector.second) {
            energy_names.push_back(energy.first);
        }
        if (ImGui::Combo("Energy##Selector", &energy_selector.first, energy_names)) {
            applyAllSettings();
        }
        model_definition.cell_energy_function->makeConfigMenu();
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Dynamics", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Checkbox("Use Dynamics", &use_dynamics)) {
            clearState();
        }
        ImGui::InputDouble("Time Step##Dynamics", &dynamics_data.dt, DYNAMICS_STEP_DT, 0, "%.4f");

        bool dynamic_matrices_changed = false;
        ImGui::Text("Mass");
        dynamic_matrices_changed |=
            ImGui::InputDouble("Spatial##Mass", &dynamics_mass_pos, DYNAMICS_STEP_MASS, 0, "%.4f");
        if (model_definition.tessellation_generator->hasSiteParam(SITE_PARAM_POWER_WEIGHT) &&
            site_param_free[SITE_PARAM_POWER_WEIGHT]) {
            dynamic_matrices_changed |= ImGui::InputDouble(
                "Power Weight##Mass", &dynamics_mass_param[SITE_PARAM_POWER_WEIGHT], DYNAMICS_STEP_MASS, 0, "%.4f");
        }
        if (site_param_free[SITE_PARAM_SIZE_TARGET]) {
            dynamic_matrices_changed |= ImGui::InputDouble(
                "Size Target##Mass", &dynamics_mass_param[SITE_PARAM_SIZE_TARGET], DYNAMICS_STEP_MASS, 0, "%.4f");
        }
        if (site_param_free[SITE_PARAM_SURFACE_TARGET]) {
            dynamic_matrices_changed |= ImGui::InputDouble(
                "Surface Target##Mass", &dynamics_mass_param[SITE_PARAM_SURFACE_TARGET], DYNAMICS_STEP_MASS, 0, "%.4f");
        }
        if (model_definition.boundary_free_param_indices.rows() > 0) {
            dynamic_matrices_changed |=
                ImGui::InputDouble("Boundary##Mass", &dynamics_mass_boundary, DYNAMICS_STEP_MASS, 0, "%.4f");
        }

        ImGui::Text("Viscosity");
        dynamic_matrices_changed |=
            ImGui::InputDouble("Spatial##Viscosity", &dynamics_viscosity_pos, DYNAMICS_STEP_VISCOSITY, 0, "%.4f");
        if (model_definition.tessellation_generator->hasSiteParam(SITE_PARAM_POWER_WEIGHT) &&
            site_param_free[SITE_PARAM_POWER_WEIGHT]) {
            dynamic_matrices_changed |=
                ImGui::InputDouble("Power Weight##Viscosity", &dynamics_viscosity_param[SITE_PARAM_POWER_WEIGHT],
                                   DYNAMICS_STEP_VISCOSITY, 0, "%.4f");
        }
        if (site_param_free[SITE_PARAM_SIZE_TARGET]) {
            dynamic_matrices_changed |=
                ImGui::InputDouble("Size Target##Viscosity", &dynamics_viscosity_param[SITE_PARAM_SIZE_TARGET],
                                   DYNAMICS_STEP_VISCOSITY, 0, "%.4f");
        }
        if (site_param_free[SITE_PARAM_SURFACE_TARGET]) {
            dynamic_matrices_changed |=
                ImGui::InputDouble("Surface Target##Viscosity", &dynamics_viscosity_param[SITE_PARAM_SURFACE_TARGET],
                                   DYNAMICS_STEP_VISCOSITY, 0, "%.4f");
        }
        if (model_definition.boundary_free_param_indices.rows() > 0) {
            dynamic_matrices_changed |= ImGui::InputDouble("Boundary##Viscosity", &dynamics_viscosity_boundary,
                                                           DYNAMICS_STEP_VISCOSITY, 0, "%.4f");
        }

        if (dynamic_matrices_changed) {
            clearState();
        }
    }
}

void FoamSubApp::makeAnalysisWindow() {
    if (ImGui::Button("Check Gradient")) {
        checkEnergyGradients(1, pow(10.0, check_gradient_epsilon_exponent), check_gradient_print_all);
    }
    ImGui::SameLine();
    ImGui::Text("Epsilon: 10^");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(70);
    ImGui::InputInt("##check_gradient_epsilon", &check_gradient_epsilon_exponent);
    ImGui::SameLine();
    ImGui::Checkbox("Print All##0", &check_gradient_print_all);
    if (ImGui::Button("Check Hessian")) {
        checkEnergyGradients(2, pow(10.0, check_hessian_epsilon_exponent), check_hessian_print_all);
    }
    ImGui::SameLine();
    ImGui::Text("Epsilon: 10^");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(70);
    ImGui::InputInt("##check_hessian_epsilon", &check_hessian_epsilon_exponent);
    ImGui::SameLine();
    ImGui::Checkbox("Print All##1", &check_hessian_print_all);
}

bool FoamSubApp::callbackKeyPressed(const CRLControlState &control_state, int key) {
    switch (key) {
        case GLFW_KEY_SPACE:
            optimize = !optimize;
            clearState();
            break;
        default:
            break;
    }
    return false;
}
